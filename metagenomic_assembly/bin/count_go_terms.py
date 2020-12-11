#!/bin/env python

#$ -cwd
#$ -V
#$ -jc long
#$ -mods l_hard mem_free 16G
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID

"""
count_go_terms.py: Counts occurrences of GO terms in each sample
based on cDNA alignment read counts

Run `count_go_terms.py --help` for usage information
"""

import os
import pandas as pd
import argparse
import sys

def process_row(row):

    """
    process_row: Counts occurrences of each GO term weighted by number of mapped reads per sample 
        in a row from a dataframe

    Required params:
        row (series): Row from df obtained using df.apply()

    """
    go_terms=row['GOs'].split(',')
    for term in go_terms:
        for sample in samples:
            count=row[sample]
            if sample in term_counts:
                if term in term_counts[sample]:
                    term_counts[sample][term]+=count
                else:
                    term_counts[sample][term]=count
            else:
                term_counts[sample]={term:count}

parser = argparse.ArgumentParser(
	description="Create counts of GO terms per-sample from per-read summary count tables"
)

parser.add_argument('--dir', help='Path to directory containing full_cdna_count_summary.txt file', required=True)

args = parser.parse_args()

global samples
global term_counts
samples=('2000', '2001', '2002', '2006', '2007', '2009', '2011', '2012', '2013', '2023', '2024', '2025')

term_counts={}
merged=pd.read_csv('{}/full_cdna_count_summary.txt'.format(args.dir),sep='\t')
merged.apply(lambda x: process_row(x),axis=1)
terms_df=pd.DataFrame.from_dict(term_counts)
terms_df.to_csv('{}/go_count_table.txt'.format(args.dir),sep='\t',index=True)

