#!/bin/env python

#$ -cwd
#$ -pe smp 64
#$ -jc long
#$ -mods l_hard mfree 200G
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -j y

"""
Determines counts of GO slim terms from cDNA mapping results

Note this is parallised using pandardallel, which has considerable 
memory demands (~200 Gb). Reducing the number of parallel threads 
will reduce memory usage at the expense of runtime.

Run `count_go_slim_terms.py --help` for usage 
"""

import pandas as pd
import numpy as np
import os
import sys
import pprint
import argparse
from pandarallel import pandarallel

def process_row(row):

    """
    Processes a row of results. Each row contains a queryname column (cDNA id), a column of comma-separated mapped GO terms ('GOs'), followed
    by a column for each sample containing the number of reads mapped to the cDNA in that sample.

    Required parameters:
        row (pd.Series): As described above

    Returns:
        df (pd.dataFrame): Dataframe of slim counts per sample for cDNA
    """
    term_counts={}
    go_terms=row['GOs'].split(',')
    slims=[]
    for term in go_terms:
        if term in slim_mapping['GO'].unique():
            slim=slim_mapping.loc[slim_mapping['GO']==term,'Slim'].values[0]
            slims.append(slim)
    # need to make the list unique to avoid overcounting...only want 1 representation of a slim per gene
    slims=list(set(slims))    
    for slim in slims:
        for sample in samples:
            count=row[sample]
            if sample in term_counts:
                if slim in term_counts[sample]:
                    term_counts[sample][slim]+=count
                else:
                    term_counts[sample][slim]=count
            else:
                term_counts[sample]={slim:count}
    df=pd.DataFrame.from_dict(term_counts)
    return(df)

def read_slim_mapping(goterm_gaf,slim_gaf):

    """
    Reads input gaf file used in slim mapping, and resulting output file and 
    merges these to create a mapping table

    Required parameters: 
        goterm_gaf(str): path to gaf containing all GO terms provided to owltools for slim mapping
        slim_gaf(str): path to output gaf from owltools containing slim terms

    Returns:
        slim_mapping(pd.DataFrame): Slim mapping table
    """

    orig_terms=pd.read_csv('{}'.format(goterm_gaf),sep="\t",header=None,skiprows=1)
    slim_terms=pd.read_csv('{}'.format(slim_gaf),sep="\t",header=None,skiprows=5)
    orig_terms=orig_terms.iloc[:,[1,4]]
    slim_terms=slim_terms.iloc[:,[1,4]]
    slim_mapping=pd.merge(orig_terms,slim_terms,how='right',on=1)
    slim_mapping=slim_mapping.iloc[:,[1,2]]
    slim_mapping.columns=['GO','Slim']

    return(slim_mapping)

def main():
    parser = argparse.ArgumentParser(
            description="Create counts of GO slim terms per-sample from per-read summary count tables"
    )

    parser.add_argument('--goterm_gaf', help='Path to .gaf file containing all GO terms provided to slim mapping', required=True)
    parser.add_argument('--slim_gaf', help='Path to .gaf file created during slim mapping', required=True)
    parser.add_argument('--counts', help='Path to file containing cDNA counts', required=True)
    parser.add_argument('--out', help='Path to output file to write', required=True)

    args = parser.parse_args()

    global slim_mapping, samples
    slim_mapping=read_slim_mapping(args.goterm_gaf, args.slim_gaf)

    pandarallel.initialize(nb_workers=64)
    merged=pd.read_csv('{}'.format(args.counts),sep='\t')
    samples=('2000', '2001', '2002', '2006', '2007', '2009', '2011', '2012', '2013', '2023', '2024', '2025')
    res=merged.parallel_apply(process_row,axis=1)
    terms_df=None
    for df in res:
        df.reset_index(level=0,inplace=True)
        if terms_df is None:
            terms_df=df
        else:
            terms_df=pd.concat([terms_df,df],ignore_index=True)

    summed=terms_df.groupby('index').sum()
    summed.index.rename('GO Term',inplace=True)
    summed.to_csv('{}'.format(args.out),sep='\t',index=True)

if __name__=='__main__':
    main()
