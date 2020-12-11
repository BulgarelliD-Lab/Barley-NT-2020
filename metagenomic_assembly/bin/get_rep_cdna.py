#!/bin/env python

#$ -cwd
#$ -V
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -mods l_hard m_mem_free 48G
#$ -adds l_hard ramdisk 10G

from Bio import SeqIO
from shutil import copyfile
import os
import argparse

parser = argparse.ArgumentParser(
	description="Obtain cDNA sequences representative of MMseqs2 clusters based on IDs of matching protein sequences"
)

parser.add_argument('--cdna', help='Fasta file of cDNA sequences from protein predictions',
	required=True)
parser.add_argument('--repprot', help='Fasta file of cluster representative proteins', 
	required=True)
parser.add_argument('--repcdna', help='Output file of cluster representative cDNAs to create', 
	required=True)

args = parser.parse_args()

RAMDISK=os.environ['TMPDIR1']
copyfile(args.cdna, "{}/cDNA.fa".format(RAMDISK))
cdna_dict=SeqIO.index('{}/cDNA.fa'.format(RAMDISK),'fasta')

with open(args.repcdna, "w") as handle:
	for record in SeqIO.parse(args.repprot, 'fasta'):
		cdna = cdna_dict[record.id]
		SeqIO.write([cdna], handle, "fasta")

