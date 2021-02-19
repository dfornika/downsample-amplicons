#!/usr/bin/env python

# Adapted from artic `align_trim.py`, Written by Nick Loman
# https://github.com/artic-network/fieldbioinformatics/blob/master/artic/align_trim.py

import argparse
import json
import pandas as pd
import pysam
import re
import sys

from copy import copy
from collections import defaultdict
from operator import itemgetter

def getPrimerDirection(primerID):
    """Infer the primer direction based on it's ID containing LEFT/RIGHT
    Taken from: https://github.com/artic-network/fieldbioinformatics/blob/master/artic/vcftagprimersites.py
    Parameters
    ----------
    primerID : string
        The primer ID from the 4th field of the primer scheme
    """
    if 'LEFT' in primerID:
        return '+'
    elif 'RIGHT':
        return '-'
    else:
        print("LEFT/RIGHT must be specified in Primer ID", file=sys.stderr)
        raise SystemExit(1)


def read_bed_file(fn):
    """Parses a bed file and collapses alts into canonical primer sites
    Taken from: https://github.com/artic-network/fieldbioinformatics/blob/master/artic/vcftagprimersites.py
    Parameters
    ----------
    fn : str
        The bedfile to parse
    Returns
    -------
    list
        A list of dictionaries, where each dictionary contains a row of the parsed bedfile.
        The available dictionary keys are - chrom, start, end, primer_id, pool_name, direction
        eg:
        [
          {
              "chrom": "MN908947.3",
              "start": 30,
              "end": 54,
              "primer_id": "nCoV-2019_1_LEFT",
              "pool_name": "1",
              "direction": "+"
          },
          {
              "chrom": "MN908947.3",
              "start": 1183,
              "end": 1205,
              "primer_id": "nCoV-2019_1_RIGHT",
              "pool_name": "1",
              "direction": "-"
          },
          {
              "chrom": "MN908947.3",
              "start": 1100,
              "end": 1128,
              "primer_id": "nCoV-2019_2_LEFT",
              "pool_name": "2",
              "direction": "+"
          },
          ...
        ]
    """

    # read the primer scheme into a pandas dataframe and run type, length and null checks
    primers = pd.read_csv(fn, sep='\t', header=None,
                          names=['chrom', 'start', 'end',
                                 'primer_id', 'pool_name'],
                          dtype={'chrom': str, 'start': int, 'end': int,
                                 'primer_id': str, 'pool_name': str},
                          usecols=(0, 1, 2, 3, 4),
                          skiprows=0)
    if len(primers.index) < 1:
        print("primer scheme file is empty", file=sys.stderr)
        raise SystemExit(1)
    if primers.isnull().sum().sum():
        print("malformed primer scheme file", file=sys.stderr)
        raise SystemExit(1)

    # compute the direction
    primers['direction'] = primers.apply(
        lambda row: getPrimerDirection(row.primer_id), axis=1)

    # separate alt primers into a new dataframe
    altFilter = primers['primer_id'].str.contains('_alt')
    alts = pd.DataFrame(
        columns=('chrom', 'start', 'end', 'primer_id', 'pool_name', 'direction'))
    alts = pd.concat([alts, primers[altFilter]])
    primers = primers.drop(primers[altFilter].index.values)

    # convert the primers dataframe to dictionary, indexed by primer_id
    #  - verify_integrity is used to prevent duplicate primer_ids being processed
    bedFile = primers.set_index('primer_id', drop=False,
                                verify_integrity=True).T.to_dict()

    # if there were no alts, return the bedfile as a list of dicts
    if len(alts.index) == 0:
        return list(bedFile.values())

    # merge alts
    for _, row in alts.iterrows():
        primerID = row['primer_id'].split('_alt')[0]

        # check the bedFile if another version of this primer exists
        if primerID not in bedFile:

            # add to the bed file and continue
            bedFile[primerID] = row.to_dict()
            continue

        # otherwise, we've got a primer ID we've already seen so merge the alt
        mergedSite = merge_sites(bedFile[primerID], row)

        # update the bedFile
        bedFile[primerID] = mergedSite

    # return the bedFile as a list
    return [value for value in bedFile.values()]


def midpoint(start, end):
    mid = round(start + ((end - start) / 2))
    return mid


def get_amplicon_midpoints(parsed_bed):
    """
    returns: {
               "nCoV-2019_1": {
                   "start": 30,
                   "midpoints": [
                       618
                   ],
                   "end": 1205
               },
               "nCoV-2019_2": {
                   "start": 1100,
                   "midpoints": [
                       1683
                   ],
                   "end": 2266
               },
               ...
             }
    """
    amplicon_midpoints = {}
    for primer in parsed_bed:
        amplicon_id = re.match("(.+)_LEFT", primer['primer_id'])
        if amplicon_id:
            amplicon_midpoints[amplicon_id.group(1)] = {}

    for amplicon_id in amplicon_midpoints.keys():
        left_primer_id = amplicon_id + "_LEFT"
        right_primer_id = amplicon_id + "_RIGHT"
        left_primer = list(filter(lambda p: p['primer_id'] == left_primer_id, parsed_bed))[0]
        right_primer = list(filter(lambda p: p['primer_id'] == right_primer_id, parsed_bed))[0]
        amplicon_start = left_primer['start']
        amplicon_end = right_primer['end']
        amplicon_midpoint = midpoint(amplicon_start, amplicon_end)
        amplicon_midpoints[amplicon_id]['start'] = amplicon_start
        amplicon_midpoints[amplicon_id]['midpoints'] = [amplicon_midpoint]
        amplicon_midpoints[amplicon_id]['end'] = amplicon_end
    
    return amplicon_midpoints

def subdivide_amplicon_midpoints(amplicon_midpoints):
    """
    """
    new_amplicon_midpoints = {}
    for amplicon_id in amplicon_midpoints:
        amplicon_start = amplicon_midpoints[amplicon_id]['start']
        amplicon_end = amplicon_midpoints[amplicon_id]['end']

        original_points = [amplicon_start] + amplicon_midpoints[amplicon_id]['midpoints'] + [amplicon_end]

        #first_midpoint = amplicon_midpoints[amplicon_id]['midpoints'][0]
#
#        if len(amplicon_midpoints[amplicon_id]['midpoints']) > 2:
#            internal_midpoints = amplicon_midpoints[amplicon_id]['midpoints'][1:-1]
#        else:
#            internal_midpoints = amplicon_midpoints[amplicon_id]['midpoints']
#
#        last_midpoint = amplicon_midpoints[amplicon_id]['midpoints'][-1]
#
#        new_first_midpoint = midpoint(amplicon_start, first_midpoint)
#
        new_midpoints = []
        for idx in range(len(original_points)):
            try:
                new_midpoint = midpoint(original_points[idx], original_points[idx + 1])
                new_midpoints.append(new_midpoint)
            except IndexError as e:
                pass

        new_midpoints = original_points[1:-1] + new_midpoints
        new_midpoints.sort()
        
        new_amplicon_midpoints[amplicon_id] = { 'start': amplicon_start,
                                                'midpoints': new_midpoints,
                                                'end': amplicon_end }
    return new_amplicon_midpoints


def find_primer(bed, pos, direction):
    """Given a reference position and a direction of travel, walk out and find the nearest primer site.
    Parameters
    ----------
    bed : list
        A list of dictionaries, where each dictionary contains a row of bedfile data
    pos : int
        The position in the reference sequence to start from
    direction : string
        The direction to search along the reference sequence
    Returns
    -------
    tuple
        The offset, distance and bed entry for the closest primer to the query position. eg:
        [8, -8, {"chrom": "MN908947.3", "start": 30, "end": 54, "primer_id": "nCoV-2019_1_LEFT", "pool_name": "nCoV-2019_1", "direction": "+"}]
        [1010, -1010, {"chrom": "MN908947.3", "start": 28677, "end": 28699, "primer_id": "nCoV-2019_29_LEFT", "pool_name": "nCoV-2019_1", "direction": "+"}]
    """
    

    if direction == '+':
        closest = min([(abs(p['start'] - pos), p['start'] - pos, p)
                       for p in bed if p['direction'] == direction], key=itemgetter(0))
    else:
        closest = min([(abs(p['end'] - pos), p['end'] - pos, p)
                       for p in bed if p['direction'] == direction], key=itemgetter(0))
    return closest


def main(args):
    """
    """
    # prepare the report outfile
    if args.report:
        reportfh = open(args.report, "w")
        print("QueryName\tReferenceStart\tReferenceEnd\tPrimerPair\tPrimer1\tPrimer1Start\tPrimer2\tPrimer2Start\tIsSecondary\tIsSupplementary\tStart\tEnd\tCorrectlyPaired", file=reportfh)

    # set up a counter to track amplicon abundance
    counter = defaultdict(int)

    # open the primer scheme and get the pools
    bed = read_bed_file(args.bedfile)
    # print(json.dumps(bed))

    amplicon_midpoints = get_amplicon_midpoints(bed)

    for _ in range(args.num_amplicon_subdivisions - 1):
        amplicon_midpoints = subdivide_amplicon_midpoints(amplicon_midpoints)

    checkpoints = []
    for amplicon_id in amplicon_midpoints.keys():
        checkpoints.append(amplicon_midpoints[amplicon_id]['start'])
        checkpoints.extend(amplicon_midpoints[amplicon_id]['midpoints'])
        checkpoints.append(amplicon_midpoints[amplicon_id]['end'])
    checkpoints.sort()

    required_depth_achieved = [False] * len(checkpoints)
        
    depths = [0] * 30000

    # open the input SAM file and process read groups
    infile = pysam.AlignmentFile("-", "rb")
    bam_header = infile.header.copy().to_dict()

    # prepare the alignment outfile
    outfile = pysam.AlignmentFile("-", "wh", header=bam_header)


    # iterate over the alignment segments in the input SAM file
    for segment in infile:
        checkpoints_under_segment = list(filter(lambda cp: segment.reference_start <= cp and segment.reference_end >= cp, checkpoints))
        checkpoints_achieved_required_depth = [True if depths[checkpoint] >= args.depth else False for checkpoint in checkpoints_under_segment]
        
        if not all(checkpoints_achieved_required_depth):
            outfile.write(segment)
            for checkpoint in checkpoints_under_segment:
                depths[checkpoint] += 1

    checkpoint_depths = [depths[checkpoint] for checkpoint in checkpoints]
    # print(json.dumps(checkpoint_depths))

    # close up the file handles
    infile.close()
    outfile.close()
    if args.report:
        reportfh.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Downsample alignments from an amplicon scheme.')
    parser.add_argument('bedfile', help='BED file containing the amplicon scheme')
    parser.add_argument('--depth', type=int, help='Subsample to n coverage')
    parser.add_argument('--no-read-groups', dest='no_read_groups',
                        action='store_true', help='Do not divide reads into groups in SAM output')
    parser.add_argument('--report', type=str, help='Output report to file')
    parser.add_argument('--num-amplicon-subdivisions', type=int, default=3, help='How many times to divide amplicons to detect coverage')
    parser.add_argument('--remove-incorrect-pairs', action='store_true')
    parser.add_argument('--verbose', action='store_true', help='Debug mode')
    args = parser.parse_args()
    main(args)
