#!/usr/bin/env python

# Adapted from artic `align_trim.py`, Written by Nick Loman
# https://github.com/artic-network/fieldbioinformatics/blob/master/artic/align_trim.py

import argparse
import json
import pandas as pd
import pysam
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
        The available dictionary keys are - Primer_ID, direction, start, end
    """

    # read the primer scheme into a pandas dataframe and run type, length and null checks
    primers = pd.read_csv(fn, sep='\t', header=None,
                          names=['chrom', 'start', 'end',
                                 'Primer_ID', 'PoolName'],
                          dtype={'chrom': str, 'start': int, 'end': int,
                                 'Primer_ID': str, 'PoolName': str},
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
        lambda row: getPrimerDirection(row.Primer_ID), axis=1)

    # separate alt primers into a new dataframe
    altFilter = primers['Primer_ID'].str.contains('_alt')
    alts = pd.DataFrame(
        columns=('chrom', 'start', 'end', 'Primer_ID', 'PoolName', 'direction'))
    alts = pd.concat([alts, primers[altFilter]])
    primers = primers.drop(primers[altFilter].index.values)

    # convert the primers dataframe to dictionary, indexed by Primer_ID
    #  - verify_integrity is used to prevent duplicate Primer_IDs being processed
    bedFile = primers.set_index('Primer_ID', drop=False,
                                verify_integrity=True).T.to_dict()

    # if there were no alts, return the bedfile as a list of dicts
    if len(alts.index) == 0:
        return list(bedFile.values())

    # merge alts
    for _, row in alts.iterrows():
        primerID = row['Primer_ID'].split('_alt')[0]

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
        [8, -8, {"chrom": "MN908947.3", "start": 30, "end": 54, "Primer_ID": "nCoV-2019_1_LEFT", "PoolName": "nCoV-2019_1", "direction": "+"}]
        [1010, -1010, {"chrom": "MN908947.3", "start": 28677, "end": 28699, "Primer_ID": "nCoV-2019_29_LEFT", "PoolName": "nCoV-2019_1", "direction": "+"}]
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
    
    pools = set([row['PoolName'] for row in bed])
    pools.add('unmatched')

    # open the input SAM file and process read groups
    infile = pysam.AlignmentFile("-", "rb")
    bam_header = infile.header.copy().to_dict()

    if not args.no_read_groups:
        bam_header['RG'] = []
        for pool in pools:
            read_group = {}
            read_group['ID'] = pool
            bam_header['RG'].append(read_group)

    # prepare the alignment outfile
    outfile = pysam.AlignmentFile("-", "wh", header=bam_header)

    # iterate over the alignment segments in the input SAM file
    for segment in infile:

        # locate the nearest primers to this alignment segment
        p1 = find_primer(bed, segment.reference_start, '+')
        if p1[0] > 1000:
            print(json.dumps(p1))
            exit(0)
        p2 = find_primer(bed, segment.reference_end, '-')

        # check if primers are correctly paired and then assign read group
        # NOTE: removed this as a function as only called once
        # TODO: will try improving this / moving it to the primer scheme processing code
        correctly_paired = p1[2]['Primer_ID'].replace('_LEFT', '') == p2[2]['Primer_ID'].replace('_RIGHT', '')
        # if not args.no_read_groups:
        #     if correctly_paired:
        #         segment.set_tag('RG', p1[2]['PoolName'])
        #     else:
        #         segment.set_tag('RG', 'unmatched')
        # if args.remove_incorrect_pairs and not correctly_paired:
        #     print("%s skipped as not correctly paired" %
        #           (segment.query_name), file=sys.stderr)
        #     continue

        # update the report with this alignment segment + primer details
        report = "%s\t%s\t%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d" % (segment.query_name, segment.reference_start, segment.reference_end, p1[2]['Primer_ID'], p2[2]['Primer_ID'], p1[2]['Primer_ID'], abs(
            p1[1]), p2[2]['Primer_ID'], abs(p2[1]), segment.is_secondary, segment.is_supplementary, p1[2]['start'], p2[2]['end'], correctly_paired)
        if args.report:
            print(report, file=reportfh)
        if args.verbose:
            print(report, file=sys.stderr)


        # normalise if requested
        if args.normalise:
            pair = "%s-%s-%d" % (p1[2]['Primer_ID'],
                                 p2[2]['Primer_ID'], segment.is_reverse)
            counter[pair] += 1
            if counter[pair] > args.normalise:
                print("%s dropped as abundance theshold reached" %
                      (segment.query_name), file=sys.stderr)
                continue

        # current alignment segment has passed filters, send it to the outfile
        outfile.write(segment)

    # close up the file handles
    infile.close()
    outfile.close()
    if args.report:
        reportfh.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Downsample alignments from an amplicon scheme.')
    parser.add_argument('bedfile', help='BED file containing the amplicon scheme')
    parser.add_argument('--normalise', type=int, help='Subsample to n coverage per strand')
    parser.add_argument('--no-read-groups', dest='no_read_groups',
                        action='store_true', help='Do not divide reads into groups in SAM output')
    parser.add_argument('--report', type=str, help='Output report to file')
    parser.add_argument('--remove-incorrect-pairs', action='store_true')
    parser.add_argument('--verbose', action='store_true', help='Debug mode')
    args = parser.parse_args()
    main(args)
