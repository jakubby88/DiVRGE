"""
Written by: Kate Johnson

Functions for finding DVGs in NGS data

Updated 02.21.2023: Fail checks
Updated: 11.16.2022
Updated the naming
Removed 2N section

"""

import pandas as pd
from scipy.stats.distributions import binom  # scikit-learn
import re
import numpy as np
import os
from sklearn.cluster import MeanShift # scikit-learn
from joblib import Parallel, delayed
import time, math
import warnings


#ORF_flag, ORF_check, CDS ESTIMATED SIZE
def frame_check(del_estimated_size):
    """
    INPUT: the estimated frag. size
           will not work without CDS
    OUTPUT: if orf is maintained (T=true or F=false)
    """
    if del_estimated_size % 3 == 0:  # if the remainder is zero-ORF maintained
        frame_flag = 'T'
    elif del_estimated_size % 3 != 0:  # elif the remainder is not zero
        frame_flag = "F"
    return frame_flag


def EstSize(segment_length, gap_size):
    """
    INPUT: chrom/seg size, size of deletion
    OUTPUT: estimated frag. size remaining after deletion
    """
    CDS_size = int(segment_length) - int(gap_size)
    frame_flag = frame_check(int(gap_size))
    return frame_flag, CDS_size


def MNM(start, M1, N, M2):
    """
    INPUT: CIGAR string, Left-most mapping position
    OUTPUT: The deletion location information
    """
    gap_start = start + M1
    gap_end = gap_start + (N-1)
    gap = '{0}_{1}'.format(gap_start, gap_end)
    return gap, gap_start, gap_end


def MNMNM(start, M1, N1, M2, N2, M3):
    """
    INPUT: CIGAR string, left-most mapping position
    OUTPUT: where the deletions start and end
    """
    gap_start1 = start + M1
    gap_end1 = gap_start1 + (N1 - 1)
    after_gap1 = gap_start1 + N1
    gap1 = '{0}_{1}'.format(gap_start1, gap_end1)
    gap_start2 = after_gap1 + M2
    gap_end2 = gap_start2 + (N2 - 1)
    gap2 = '{0}_{1}'.format(gap_start2, gap_end2)
    return gap1, gap2, gap_start1, gap_end1, gap_start2, gap_end2


def printer1(outfile, name, segment, readname, gap, gap_start, gap_end,
             estimated_size, length, segment_length, frame_flag, totalreads,
             readflags, strain, cigar):
    """
    INPUT: info on reads that span deletion
    OUTPUT: print out info to file
    """
    print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13}'.format(
        name, segment, readname, gap, gap_start, gap_end,
        estimated_size, length, segment_length, frame_flag, totalreads,
        readflags, strain, cigar), end="\n", file=outfile)


def printer2(outfile, name, segment, readname, gap1, gap2, gap_start1,
             gap_end1, gap_start2, gap_end2, gap_size1, gap_size2,
             align_start, cigar, totalreads, readflags, strain,
             segment_length):
    """
    INPUT: info for two gaps
    OUTPUT: print info for two gaps to a specific file
    """
    print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16}'.format(
            name, segment, readname, gap1, gap2, gap_start1, gap_end1,
            gap_start2, gap_end2, gap_size1, gap_size2, align_start, cigar,
            totalreads, readflags, strain, segment_length), end="\n",
          file=outfile)

def PrepDVGFrame(split_df, features_df, min_align_length):
    """
    INPUT: Split read file (split_df) and features count file (features_df)
    OUTPUT: DF with length/cigar checks to be filtered
    """

    split_df = split_df.astype({"name": str})  # updated 08.04.2022
    features_df = features_df.astype({"name": str})  # updated 08.04.2022
    # merge the two files to combine information (will use this to work on)
    df = pd.merge(split_df, features_df, how='left',
                  left_on=['segment', 'name'], right_on=['segment', 'name'])
    cigs = list(df['cigar'])

    break_cigar = [re.findall(r'(\d+)(\w)', x) for x in cigs]
    # Set empty lists - will be added to df
    # At end, lists and df should be same len. Can check with print(len(list))
    total_indel = []  # list of total indels
    Where_N = []  # Where deletion located in cigar
    Cigar_len = []  # Length of cigar
    tup = []  # List containing T/F info on Length pass check
    no_soft_list = []

    for cig_tuple in break_cigar:
        indel = [int(x[0]) for x in cig_tuple if x[1] == 'I' or x[1] == 'D']
        total_indel_sum = sum(indel)  # sum total deletions and insertions
        total_indel.append(total_indel_sum)
        # doesn't include soft clips or insertions in alignment)
        no_soft = [i for i in cig_tuple if i[1] != 'S' and i[1] != 'I']
        no_soft_list.append(no_soft)
        # check to see if alignment length passes input parameter: T or F
        align_pass = all(int(i[0]) >= min_align_length for i in no_soft if i[1] == 'M')
        tup.append(align_pass)  # used as a filter for alignment L
        cigar_len = len(no_soft)
        Cigar_len.append(cigar_len)
        # determine number of 'N'
        where_n = [x for x, y in enumerate(no_soft) if y[1] == 'N']
        Where_N.append(len(where_n))  # add to list should be as long as df

    df['number_N'] = Where_N  # add lists to DF
    df['Cigar_len'] = Cigar_len
    df['total_indel'] = total_indel
    df['no_soft'] = no_soft_list
    df['tup'] = tup

    return(df)


def PrepOneDeletions(one_gap, Reads1N, strain, gap_size, total_mismatch, align_length):
    """
    INPUT: the estimated size using the CDS coordinates
    OUTPUT: Whether the estimated size is divisible by 3 and
    therefore maintaining open read frame
    """
    OneGapFile = open(Reads1N, 'w')

    HEADER1 = 'name,segment,readname,gap,gap_start,gap_end,gap_size,estimated_length,segment_size,frame_flag,totalreads,readflags,strain,cigar'

    print(HEADER1, end="\n", file=OneGapFile)

    # Iterate through rows to generate gap info for one and two gapped reads
    one_count = 0

    for index, row in one_gap.iterrows():
        seg_length = row['SegmentLength']
        start = row['left_pos']
        no_soft = row['no_soft']
        Cigar_len = row['Cigar_len']
        # determine the number of 'N'
        Where_N = [x for x, y in enumerate(no_soft) if y[1] == 'N']
        N_idx = Where_N[0]  # find the position where the gap is
        N = int(no_soft[N_idx][0])  # the length of the gap
        # incase there are INDELS- identify everything before and after 'N'
        Before_idx = range(0, N_idx)  # find the idx pos of everything before N
        After_idx = range(N_idx+1, Cigar_len)

        five_length = []  # how many idx pos are there before the 'N'
        for x in Before_idx:
            five_length.append(int(no_soft[x][0]))
        five_length = (sum(five_length))  # sum up all of the pos M,D,I, bef N

        three_length = []  # how many idx positions are there after 'N'
        for x in After_idx:
            three_length.append(int(no_soft[x][0]))
        three_length = sum(three_length)  # sum up all of the M/D/I aft. "N"

        M1 = five_length  # new 'match' before N
        M2 = three_length  # new "match" length after 'N'

        # calculate the following useing the function MNM for the pattern
        gap, gap_start, gap_end = MNM(start, M1, N, M2)

        if N > gap_size:
            one_count += 1
            frame_flag, est_length = EstSize(seg_length, N)
            printer1(OneGapFile, row['name'], row['segment'], row['readname'],
                     gap, gap_start, gap_end, N, est_length,
                     row['SegmentLength'], frame_flag, row['totalreads'],
                     row['flags'], strain, row['cigar'])

    OneGapFile.close()

    return(one_count)


def PrepTwoDeletions(two_gap, Reads2N, strain, gap_size, total_mismatch,
                     align_length):
    """
    INPUT: the estimated size using the CDS coordinates
    OUTPUT: Whether the estimated size is divisible by 3 and
            therefore maintaining open read frame
    """
    TwoGapFile = open(Reads2N, 'w')

    HEADER2 = 'name,segment,readname,gap1,gap2,gap_start1,gap_end1,gap_start2,gap_end2,gap_size1,gap_size2,align_start,cigar,totalreads,readflags,strain,segment_size'

    print(HEADER2, end="\n", file=TwoGapFile)

    two_count = 0

    for index, row in two_gap.iterrows():
        seg_length = row['SegmentLength']
        start = row['left_pos']
        no_soft = row['no_soft']
        Cigar_len = row['Cigar_len']
        # determine the number of 'N'
        Where_N = [x for x, y in enumerate(no_soft) if y[1] == 'N']
        N_idx1 = Where_N[0]  # find the first gap in the pattern
        N1 = int(no_soft[N_idx1][0])  # N1 is the length of first gap
        N_idx2 = Where_N[1]  # find the second gap in the pattern
        N2 = int(no_soft[N_idx2][0])  # length of second gap
        match1 = range(0, N_idx1)  # match positions before N1
        match2 = range(N_idx1 + 1, N_idx2)  # match positions between N1 and
        match3 = range(N_idx2 + 1, Cigar_len)

        match_1 = []
        for x in match1:
            match_1.append(int(no_soft[x][0]))
        match_1 = (sum(match_1))

        match_2 = []
        for x in match2:
            match_2.append(int(no_soft[x][0]))
        match_2 = sum(match_2)

        match_3 = []
        for x in match3:
            match_3.append(int(no_soft[x][0]))
        match_3 = sum(match_3)

        gap1, gap2, gap_start1, gap_end1, gap_start2, gap_end2 = MNMNM(start,
                                                                       match_1,
                                                                       N1,
                                                                       match_2,
                                                                       N2,
                                                                       match_3)

        if N1 > gap_size and N2 > gap_size:
            two_count += 1

            printer2(TwoGapFile, row['name'], row['segment'],
                     row['readname'], gap1, gap2,
                     gap_start1, gap_end1, gap_start2, gap_end2, N1, N2,
                     start, row['cigar'], row['totalreads'],
                     row['flags'], strain, seg_length)

    TwoGapFile.close()

    return(two_count)


def PrepFreq(filename, GapNumber):
    """
    INPUT: The candidate gap dataframe file
    ({0}/{1}.CandidateDI_OneGap.N{2}.Mis{3}.csv)

    OUTPUT: Outputs dataframe to be written as a csv file

    """
    if GapNumber == 1:

        groups = ['name', 'segment', 'gap_start', 'gap_end', 'gap', 'gap_size',
                  'estimated_length', 'frame_flag']

        df = pd.read_csv(filename, sep=',', keep_default_na=False)  # read file
        # count number of reads that match particular dvg
        dvg_freq = df.groupby(groups).size().reset_index(name="freq")
        # drop any duplicates to not double count
        dvg_freq = dvg_freq.drop_duplicates()
        return dvg_freq  # return dataframe to be used in grouping script

    if GapNumber == 2:
        # columns used to group and merge with
        groups = ['name','segment', 'gap_start1', 'gap_end1', 'gap1',
                  'gap_start2', 'gap_end2', 'gap2', 'gap_size1', 'gap_size2']

        # read file
        df = pd.read_csv(filename, sep=',', keep_default_na=False)
        # count the number of reads that match particular dvg across samples
        dvg_freq = df.groupby(groups).size().reset_index(name="freq")
        # drop read information to count number of samples with given dvg
        #df2 = df.drop(['readname', 'totalreads', 'readflags'], axis=1)
        # drop any duplicates
        dvg_freq = dvg_freq.drop_duplicates()
        # count number of samples with DVG types
        #samp_freq = df2.groupby(groups).size().reset_index(
        #        name='number_of_samples')
        # merge the two dataframes to be used in grouping script
        #df_merge = pd.merge(dvg_freq, samp_freq, how='left',
        #                    left_on=groups, right_on=groups)

        return dvg_freq

def GroupingOne(s, n, df, bandwidth=5):
    """
    INPUT: Freq. information for each sample
    OUTPUT: Grouping information for each sample/segment
    """

    warnings.filterwarnings("ignore", category=UserWarning)

    print("Grouping {0} 1N DVGs for {1}".format(s, n))
    dfs = df[(df.segment == s) & (df.name == n)].copy()  # subset dataframe for given segment

    if not dfs.empty:
        s = list(dfs['gap_start'])  # list of starts
        e = list(dfs['gap_end'])  # list of ends
        # zip and make into an array
        DVGs = np.array(list(zip(s, e)))  # zip will work left to right
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(DVGs)
        labels = ms.labels_
        labels_unique = list(np.unique(labels))
        n_clusters_ = len(labels_unique)

        #print("Number of estimated clusters : {0}".format(n_clusters_))

        grouping_info = list(zip(list(dfs['segment']), s, e, list(dfs['gap']),
            list(dfs['gap_size']),
            list(dfs['estimated_length']),
            list(dfs['freq']),
            list(labels)))

        g_df = pd.DataFrame(grouping_info,
            columns=['segment', 'gap_start', 'gap_end',
                    'gap', 'gap_size', 'estimated_length',
                    'freq', 'DVG_group'])

        # now we need to iterate through all of the groups and add in info
        for group in labels_unique:
            t = g_df[g_df.DVG_group == group]  # temp df of group data
            g_df.loc[g_df.DVG_group == group, 'Deletion'] = list(t[t.freq == t.freq.max()]['gap'])[0]
            g_df.loc[g_df.DVG_group == group, 'DeletionStart'] = list(t[t.freq == t.freq.max()]['gap_start'])[0]
            g_df.loc[g_df.DVG_group == group, 'DeletionEnd'] = list(t[t.freq == t.freq.max()]['gap_end'])[0]
            g_df.loc[g_df.DVG_group == group, 'GroupBoundaries'] = '{0}-{1}_{2}-{3}'.format(t['gap_start'].min(), t['gap_start'].max(), t['gap_end'].min(), t['gap_end'].max())
            g_df.loc[g_df.DVG_group == group, 'DeletionSize'] = list(t[t.freq == t.freq.max()]['gap_size'])[0]
            g_df.loc[g_df.DVG_group == group, 'EstimatedFragLength'] = list(t[t.freq == t.freq.max()]['estimated_length'])[0]

        return(g_df)
    else:
        print('no deletions {0} {1}'.format(s, n))

# MAKING INTO IT'S OWN FUNCTION! added 09.16.2021
#def GroupOneDeletionDVGs(infile, bandwidth, outfilename, njob):
def GroupOneDeletionDVGs(indataframe, bandwidth, outfilename, njob):
    """
    ADD IN PARAMETER FOR NUMBER OF THREADS/JOBS
    INPUT: freq. dataframe for all samples
    OUTPUT: group dataframe for all samples
    """
    #df = pd.read_csv(infile, sep=',', keep_default_na=False)
    df = indataframe
    SEGMENTS = list(set(list(df['segment'])))
    NAMES = list(set(list(df['name']))) # ADDED
    masterDF = pd.concat(Parallel(n_jobs=njob)(delayed(GroupingOne)(s, n, df, bandwidth) for s in SEGMENTS for n in NAMES))
    print("Grouping 1N finished")
    masterDF.to_csv(outfilename, index=False)



def GroupingTwo(s, n, df, bandwidth=5):
    print("Grouping {0} 2N DVGs for {1}".format(s, n))
    warnings.filterwarnings("ignore", category=UserWarning)
    dfs = df[(df.segment == s) & (df.name == n)].copy()  # subset dataframe for given segment

    if not dfs.empty:
        s = list(dfs['gap_start1'])  # list of starts
        e = list(dfs['gap_end1'])  # list of ends
        s2 = list(dfs['gap_start2'])
        e2 = list(dfs['gap_end2'])
        DVGs = np.array(list(zip(s, e, s2, e2)))
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(DVGs)
        labels = ms.labels_
        labels_unique = list(np.unique(labels))
        n_clusters_ = len(labels_unique)
        grouping_info = list(zip(list(dfs['segment']), s, e,
                                 list(dfs['gap1']), s2, e2,
                                 list(dfs['gap2']),
                                 list(dfs['freq']),
                                 list(labels)))
        g_df = pd.DataFrame(grouping_info, columns=['segment', 'gap_start1',
                                                    'gap_end1', 'gap1',
                                                    'gap_start2', 'gap_end2',
                                                    'gap2', 'freq',
                                                    'DVG_group'])

        # now we need to iterate through all of the groups and add in info
        for group in labels_unique:
            t = g_df[g_df.DVG_group == group]  # temp df of group data
            g_df.loc[g_df.DVG_group == group, 'Deletion1'] = list(t[t.freq == t.freq.max()]['gap1'])[0]
            g_df.loc[g_df.DVG_group == group, 'DeletionStart1'] = list(t[t.freq == t.freq.max()]['gap_start1'])[0]
            g_df.loc[g_df.DVG_group == group, 'DeletionEnd1'] = list(t[t.freq == t.freq.max()]['gap_end1'])[0]
            g_df.loc[g_df.DVG_group == group, 'Deletion2'] = list(t[t.freq == t.freq.max()]['gap2'])[0]
            g_df.loc[g_df.DVG_group == group, 'DeletionStart2'] = list(t[t.freq == t.freq.max()]['gap_start2'])[0]
            g_df.loc[g_df.DVG_group == group, 'DeletionEnd2'] = list(t[t.freq == t.freq.max()]['gap_end2'])[0]
            g_df.loc[g_df.DVG_group == group, 'group1_boundaries'] = '{0}-{1}_{2}-{3}'.format(t['gap_start1'].min(), t['gap_start1'].max(), t['gap_end1'].min(), t['gap_end1'].max())
            g_df.loc[g_df.DVG_group == group, 'group2_boundaries'] = '{0}-{1}_{2}-{3}'.format(t['gap_start2'].min(), t['gap_start2'].max(), t['gap_end2'].min(), t['gap_end2'].max())

        return(g_df)
    else:
        print('no deletions {0} {1}'.format(s, n))


def GroupTwoDeletionDVGs(infile, bandwidth, outfilename, njob):
    """
    INPUT: two deletion freq. information
    OUTPUT: grouped two deletion samples GROUPING ACROSS SAMPLES, NOT INDIV.
    """
    df = infile
    SEGMENTS = list(set(list(df['segment'])))
    NAMES = list(set(list(df['name']))) # ADDED
    masterDF = pd.concat(Parallel(n_jobs=njob)(delayed(GroupingTwo)(s, n, df, bandwidth) for s in SEGMENTS for n in NAMES))
    print("Grouping 2N finished")
    masterDF.to_csv(outfilename, index=False)


def MergeReadsGroups(ReadsFile, GroupedFile, GapNumber):
    """
    INPUT:Read file info, group file, and number of deletions
    OUTPUT: df with read and group info, used to calc: rpkm, percentage, etc.
    """

    if GapNumber == 1:
        # read the candidate reads file
        df = pd.read_csv(ReadsFile, sep=',', keep_default_na=False)

        # read the grouped/new gap files from grouping script
        d1g = pd.read_csv(GroupedFile, sep=',', keep_default_na=False)

        r, c = d1g.shape

        if r > 1:
            # merge the two dataframes
            d1g_m = pd.merge(df, d1g, how='left', on=['segment', 'gap',
                                                      'gap_start', 'gap_end',
                                                      'gap_size',
                                                      'estimated_length'])

            # drop any duplicates that were generated by merge
            d1g_m = d1g_m.drop_duplicates()

            # rename the freq column to keep for later use
            d1g_m = d1g_m.rename(columns={"freq": "FreqAcrossSamples"})

            # return to be used for rpkm, percentages, etc.
            return d1g_m

        else:
            print('No gap 1 files to merge')

    if GapNumber == 2:
        # read the candidate two gap file
        df2 = pd.read_csv(ReadsFile, sep=',', keep_default_na=False)

        # read the grouped New gap file made using grouping script
        d2g = pd.read_csv(GroupedFile, sep=',', keep_default_na=False)

        r, c = d2g.shape

        if r > 1:
            # merge the two dataframes
            d2g_m = pd.merge(df2, d2g,
                             how='left',
                             on=['segment', 'gap1', 'gap_start1',
                                 'gap_end1', 'gap2', 'gap_start2', 'gap_end2'])

            # drop any duplicates potentially generated from merging
            d2g_m = d2g_m.drop_duplicates()

            # rename freq column to be used later on and not get confused
            d2g_m = d2g_m.rename(columns={"freq": "FreqAcrossSamples"})

            # return to be used for counting, rpkm, percentage, etc.
            return d2g_m

        else:
            print("No files for two gap")


def CountGroupedDVGs(df, GapNumber):
    """
    INPUT: Take in grouping/read merged dataframe from MergedReadGroups,
    count the new gap information for each sample

    OUTPUT: new dataframe that has the freq info for total number of gaps and
    Count information for each individual dvg type within the samples.
    This will then be used as input into the binomial check.
    Counts calculated with "New" gaps generated by grouping script
    """
    if GapNumber == 1 and df is not None:
        # count num reads / dvg
        dvg_freq = df.groupby(['name', 'segment',
                               'DeletionStart', 'DeletionEnd',
                               'Deletion']).size().reset_index(name="deletion_count")

        # count num reads /sample segment
        total_dvg = df.groupby(['name', 'segment']).size().reset_index(name="SegTotalDVG")

        # merge the dvg freq and total dvg counts into df
        dvg_c = pd.merge(dvg_freq, total_dvg, on=['name', 'segment'])

        # merge the dvg freq, total dvg counts with the full df with read info
        dvg_f = pd.merge(df, dvg_c, on=['name', 'segment', 'DeletionStart', 'DeletionEnd', 'Deletion'])

        # drop any duplicates generated during merging
        dvg_f.drop_duplicates()

        # return for rpkm, percentage calcs
        return dvg_f
    else: 
        print("No single gaps to count")

    if GapNumber == 2 and df is not None:
        # count num reads / dvg
        dvg_freq = df.groupby(['name', 'segment',
                               'DeletionStart1', 'DeletionEnd1',
                               'Deletion1','DeletionStart2',
                               'DeletionEnd2',
                               'Deletion2']).size().reset_index(name="deletion_count")

        # count num reads/sample segment
        total_dvg = df.groupby(['name', 'segment']).size().reset_index(name="SegTotalDVG")

        # merge the dvg counts and total dvg counts
        dvg_c = pd.merge(dvg_freq, total_dvg, on=['name', 'segment'])

        # merge count information with all data which includes read info
        dvg_f = pd.merge(df, dvg_c, on=['name', 'segment',
                                        'DeletionStart1', 'DeletionEnd1', 'Deletion1',
                                        "DeletionStart2", "DeletionEnd2", "Deletion2"])

        # drop any duplicates potentially generated during merge
        dvg_f.drop_duplicates()

        # return dataframe to be used for rpkm, binocheck, percentage, etc.
        return dvg_f
    else:
        print("No 2 deletions to count")




def ReduceDF(df, GapNumber):
    """
    INPUT: Read dataframe
    OUTPUT: A reduced dataframe, with info that we care about for figures
    """
    if GapNumber == 1:
        dropcols = ['readname', 'gap', 'gap_start', 'gap_end',
                    'gap_size', 'estimated_length', 'frame_flag',
                    'readflags', 'FreqAcrossSamples', 'cigar']

        df = df.drop(dropcols, axis=1)

        df = df.drop_duplicates()

        return df

    if GapNumber == 2:
        dropcols = ['readname', 'align_start', 'cigar', 'readflags',
                    'gap_start1', 'gap_end1',
                    'gap1', 'gap_start2', 'gap_end2', 'gap2', 'gap_size1',
                    'gap_size2', 'FreqAcrossSamples']

        df = df.drop(dropcols, axis=1)

        df = df.drop_duplicates()

        return df


def CleanFeatureCounts(df, strain):
    idcols = list(df.columns[0:6])
    sampcols = list(df.columns[6:])

    fdf = pd.melt(df,
                  id_vars=idcols, value_vars=sampcols,
                  var_name="bam_name", value_name='totalreads')

    fdf['bam_name'] = fdf['bam_name'].astype(str)

    # updating to work with multiple differe naming types - not just our pipeline
    fdf['name'] = fdf["bam_name"].apply(lambda x: os.path.splitext(os.path.basename(x))[0].split('.')[0])

    # clean up and change column names
    fdf = fdf[['Chr', 'Length', 'name', 'totalreads']]

    # rename columns to fit in script
    fdf = fdf.rename(columns={'Chr': 'segment',
                              'Length': 'SegmentLength',
                              'name': 'name',
                              'totalreads': 'totalreads'})

    return fdf
