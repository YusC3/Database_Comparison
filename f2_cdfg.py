"""
This is code for the research project conducted at the University of Washinton Bothell.
This code is meant to r
Authors: YeeMey Seah, Yusuf Corr."""

import pandas as pd
import numpy as np
import seaborn as sns
from argparse import ArgumentParser
import matplotlib.pyplot as plt

def make_user_interface():
    #defining help string
    help_string = \
    """
    DISCLAIMER: This code is meant to deal with outputs from Dr.
    Yee Mey Seah’s code found here:
    https://github.com/yeemey/dvh_mms2/blob/master/notebooks/data_exploration/mutation_coverage.ipynb
    IN PARTICULAR, the researcher should run this line:
    SAMPLE_LINEmutN, SAMPLE_LINEevidenceN = subset_gd_to_df('/path/to/annotated.gd', 'D2-0', 'D2', '0', cov=True),
    for each generational time stamp in a particular sample,
    using N + 1 to create different titles for each df creation.
    Creation of csv file for all generational time stamps for a particular sample line:
    For each SAMPLE_LINEmutN, use pd.concat([SAMPLE_LINEmut1, SAMPLE_LINEmut2, SAMPLE_LINEmutN], ignore_index= True)
    to get a combined df of all ACCEPTED mutations. Then use: Accepted_mut.to_csv to create a CSV file
    USE SAME METHOD FOR SAMPLE_LINEevidenceN, to create Evidence_df.csv

    USE THESE FINAL CSV FILES FOR THIS PROGRAM.

    This code has the sole purpose of filtering through a compiled data frame of
    ALL accepted mutation across all sequenced generational time points,
    find the mutations that do not have “accepted data” across all generational time points
    (and this can be specified by number of matches that SHOULD exist in the data frame),
    and finally, using another data frame (contain ALL the evidence for every generational time stamp,
    both accepted and rejected evidence for mutations)
    and create a new data frame with the rejected evidence that matches the data missing for the mutations that are not covered throughout
    each time stamp.
    """
    #making parser option
    parser = ArgumentParser(description = help_string)
    parser.add_argument('accepted_df', help=\
    "input path to csv for a sample line of ACCEPTED mutations, ALL GEN")
    parser.add_argument('evidence_df', help=\
    "input path to csv for evidence for a specific sample line, ALL GEN")
    parser.add_argument('-i', '--int', help=\
    "this is the number of expected lines for each mutation.", type = int)
    parser.add_argument('-q', '--que1', help=\
    "if 'yes' or 'y', will create an output file.", default= "none", type = str)
    parser.add_argument('-p', '--pos', help=\
    "position number of a specific mutation in a specific part of genome.", default= 0, type = int)
    parser.add_argument('-g', '--g_id', help=\
    "genome id for part of genome the mutation resides.", default= "none", type = str)
    return parser

def main():
    """Run the analysis"""
    parser = make_user_interface()
    #Interpret user interface arguments
    args = parser.parse_args()
    answer1 = args.que1
    number = args.int
    df1 = args.accepted_df
    df2 = args.evidence_df
    position = args.pos
    genomeid = args.g_id

    #open user files
    accepted_mut = pd.read_csv(df1)
    evid_df = pd.read_csv(df2)

    def missing_column(a_df, number):
        """
        This is the part of the code that will build a list of rows(by a unique key value made up of two column vlaues)
        that will be used to search through an evidence data frame to filter out said rows into a new dataframe
        """
        filter_list =[]
        filter_list2 = []
        for row in a_df.index:
            u = a_df.iloc[row]['genome_id']
            t = a_df.iloc[row]['position']
            df_same = a_df[(a_df['position'] == t) & (a_df['genome_id'] == u)]
            if int(len(df_same.index)) != number: #have user input number
                filter_list.append([t])
                filter_list2.append([u])
        s1 = pd.Series(v[0] for v in filter_list)
        s2 = pd.Series(v[0] for v in filter_list2)
        filtered_column = pd.concat([s1, s2], axis = 1, ignore_index = True)
        filtered_column.rename(index=str, columns = {0: 'position', 1: 'genome_id'}, inplace=True)

        return filtered_column

    def missing_data_df(filter_column, df_evidence):
        """
        This is the part that will take the filtered data and find said rows using two key values as markers.
        It will then return a new df with those specific rows only.
        f_evi is equal to the evidence dataframe PRE-FILTERED TO FOCUS ON REJECTED EVIDENCE only.
        """
        f_evi = df_evidence[df_evidence['reject'].isna() == False]
        final_df = pd.merge(filter_column, f_evi, how = 'inner', on = ['genome_id', 'position'])

        return final_df

    def create_output(answer, missingdf):
        """
        This creates a CSV file that contains all of the evidence that matches the blanks found
        in the accepted data base across given generational time stamps.
        """
        if answer == 'yes' or answer == 'y':
            missingdf.to_csv('Missing_Mutations')
        elif answer != 'yes' or answer != 'y':
            print('No output')

    missing = missing_column(accepted_mut, number)
    missing_evid = missing_data_df(missing, evid_df)
    output = create_output(answer1, missing_evid)


if __name__== "__main__":
    main()
