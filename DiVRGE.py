"""
Written by: Kate Johnson
Email: kej310@nyu.edu

# 11.21.2022:
Removed the two 'N' section

# 01.31.2023: Adjusted to output empty file if no deletions were found

# 02.21.2023: Included flags for 0 read groupings (cant count)
run:
python3 DiVRGE.v3.py   \
    --strain ${STRAIN}   \
    --ref ${REF}   \
    --file "/home/kate/Lab/DiVRGE/testing/N_files/New_W17_1415_d02_rep2.H9N2.split.txt"   \
    --idx "/home/kate/Lab/DiVRGE/testing/DVG/H9N2.FeaturesOutput.csv"   \
    --save_dir "/home/kate/Lab/DiVRGE/testing/DVG/"   \
    --align_length ${ALIGN_LEN}   \
    --total_mismatch ${MIS}   \
    --group_bw 5 \
    --nbams 1 \
    --njob 2

"""

import argparse
import pandas as pd
import gc
import os
from joblib import Parallel, delayed
import time, math
from DiVRGE_functions import (PrepDVGFrame, PrepOneDeletions, PrepTwoDeletions,
                           PrepFreq, GroupOneDeletionDVGs,
                           GroupTwoDeletionDVGs, MergeReadsGroups,
                           CountGroupedDVGs,
                           ReduceDF, CleanFeatureCounts)

import warnings

start = time.time()

parser = argparse.ArgumentParser()

parser.add_argument('--strain', '-s', required=True, help='Indicate strain')

parser.add_argument('--ref', '-r', required=True,
                    help='Directory path + ref name')

parser.add_argument('--file', '-f', required=True, help='Split CIGAR txt file')

parser.add_argument('--idx', '-i', required=True, help='Features file')

parser.add_argument('--save_dir', '-d', default='.',
                    help='Save directory (default: ./) ')

parser.add_argument('--align_length', '-l', type=int, default=28,
                    help='Give the length that each portion of the read must \
                    be (default 28 bp)')

parser.add_argument('--gap_size', '-g', type=int, default=5,
                    help='Give the minimum size of the gap (default 5)')

parser.add_argument('--total_mismatch', '-m', type=int, default=2,
                    help='Only 2 INDELS allowed in alignment portion')

parser.add_argument('--group_bw', '-w', type=int, default=5,
                    help='Provide the bandwidth for the grouping script (default: 5)')

parser.add_argument('--njob', '-j', type=int, default=2,
                    help='Number of threads')

parser.add_argument('--nbams', '-b', type=int, default=0,
                    help='Provide number of samples being used in idx file')


args = parser.parse_args()

if __name__ == '__main__':
    print("")
    print("Loading in parameters and files")
    bandwidth = args.group_bw
    infile = args.file  # read in split read file
    idxstats = args.idx  # read in index stat file
    njob = args.njob
    
    try:
        fdf = pd.read_csv(infile, sep=',', keep_default_na=False)  # read split read file
    except pd.errors.EmptyDataError:
        fdf = pd.DataFrame()

    if fdf.shape[0] > 0 :
        prefix = list(set(list(fdf['name'])))[0] if args.nbams == 1 else args.strain
        print("INPUT FILE: {0}".format(infile))
        print("INPUT FEATURES: {0}".format(idxstats))
        print("OUTPUT FILENAME: {0}".format(prefix))
        print('STRAIN: {0}'.format(args.strain))
        print("MATCH LENGTH: {0}".format(args.align_length))
        print("GROUPING SIZE: {0}".format(args.group_bw))
        print("MIN. DELETION SIZE: {0}".format(args.gap_size))
        print("")
        print("")

        # file names for files to be made, candidate 2 gaps output but not
        Reads1N = "{0}/{1}.CandidateDVG_OneGap.{6}.N{2}.Mis{3}.M{4}.G{5}.csv".format(
            args.save_dir, prefix, args.gap_size,
            args.total_mismatch, args.align_length, args.group_bw,
            args.strain)

        Reads2N = "{0}/{1}.CandidateDVG_TwoGap.{6}.N{2}.Mis{3}.M{4}.G{5}.csv".format(
            args.save_dir, prefix, args.gap_size,
            args.total_mismatch, args.align_length, args.group_bw,
            args.strain)

        Group1N = "{0}/{1}.DVG.grouped.OneGap.{6}.N{2}.Mis{3}.M{4}.G{5}.csv".format(
                args.save_dir, prefix, args.gap_size,
                args.total_mismatch, args.align_length, args.group_bw,
                args.strain)

        Final1N = "{0}/{1}.DVG.FINAL.OneGap.{6}.N{2}.Mis{3}.M{4}.G{5}.csv".format(
                args.save_dir, prefix, args.gap_size,
                args.total_mismatch, args.align_length, args.group_bw,
                args.strain)

        # working on data:
        print("Cleaning the feature counts file")
        idf = CleanFeatureCounts(pd.read_csv(idxstats, sep=',', keep_default_na=False), args.strain)  # clean feat cts

        print("Extracting CIGAR string information for filtering steps")
        df = PrepDVGFrame(fdf, idf, args.align_length)  # prep for filt

        print("Filtering for single deletions that pass set thresholds")
        one_gap = df[(df.number_N == 1) & (df.tup == bool('True')) & (
            df.total_indel <= args.total_mismatch)]  # temp df of group data

        print("Filtering for double-deletions that pass set thresholds")
        two_gap = df[(df.number_N == 2) & (df.tup == bool('True')) & (
            df.total_indel <= args.total_mismatch)]  # temp df of group data

        # remove from memory
        del df
        del fdf
        del idf
        gc.collect()

        # SINGLE DELETIONS
        print("Prepping single-deletions")
        one_count = PrepOneDeletions(one_gap, Reads1N, args.strain,
                                    args.gap_size, args.total_mismatch,
                                    args.align_length)

        # DOUBLE DELETIONS
        print("Prepping double-deletions")
        two_count = PrepTwoDeletions(two_gap, Reads2N,
                                    args.strain, args.gap_size,
                                    args.total_mismatch, args.align_length)

        del one_gap
        del two_gap
        gc.collect()

        if one_count > 0:
            del one_count
            gc.collect()

            # groups files which will be output from grouping
            df_merge = PrepFreq(Reads1N, 1)

            print("")
            print("Grouping single-deletions")
            GroupOneDeletionDVGs(df_merge, bandwidth, Group1N, njob)

            # Merging read information with NewGap information
            d1g_m = MergeReadsGroups(Reads1N, Group1N, 1)

            if d1g_m is not None: 
                # Count the New Gap information for each sample and segment
                d1c = CountGroupedDVGs(d1g_m, 1)

                print("Final 1N output: {0}".format(Final1N))
                d1c = ReduceDF(d1c, 1)
                d1c.to_csv(Final1N, index=False)

                del d1g_m
                del d1c
                gc.collect()
            else:
                print("No reads passed input requirements for {0}".format(prefix))
        else:
            print("NOTE: No reads with a single-deletion passed for {0}".format(prefix))

    elif fdf.shape[0] == 0:
        print('NOTE: {0} was empty. Skipping.'.format(infile))
        #continue # will skip the rest of the block and move to next file

print("finished running")

end = time.time()

print('{:.4f} s'.format(end-start))
