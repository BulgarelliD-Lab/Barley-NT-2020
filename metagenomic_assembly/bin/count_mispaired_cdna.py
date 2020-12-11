#!/bin/env python

#$ -cwd
#$ -jc inf
#$ -pe smp 24
#$ -mods l_hard mfree 200G
#$ -adds l_hard ramdisk 4G
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID.$TASK_ID

"""

count_mispaired_cdna.py: Generates counts of reads mapped to different reference sequences

Each read is assigned a count of 0.5 to ensure correct weighting, whereas paired reads are
and singletons are assigned a count of 1 since they represent a full sequence fragment.

Read are only considered applicable if they share n GO terms in column, defined by the value 
passed with the --common_terms argument.

BAM files are split into subsections and subsequently processed in parallel.

For usage information, run `count_mispaired_cdna.py --help`

"""

import pandas as pd
from multiprocessing import Pool,Lock
import pysam
import io
import os
import sys
import argparse
import glob
from shutil import copyfile
import pprint

lock=Lock()

def split_bam(infile):

	"""
	split_bam: Subdivides a bam file into chunks of 200000 alignments
	
	Required parameters:
		infile (str):	Name of BAM file

	Returns:
		chunk (int):	Count of number of chunks written
	"""

	tmpdir=os.environ['TMPDIR1']

	chunk=1
	read_count=0

	basename=os.path.basename(infile)
	base=os.path.splitext(basename)[0]
	outfile='{}.{}.bam'.format(base,chunk)

	bam=pysam.Samfile("{}/{}".format(tmpdir,basename),"rb")
	header=bam.header

	outbam=pysam.Samfile("{}/{}".format(tmpdir,outfile),'wb', header=header)

	print('Writing chunk {}'.format(chunk))
	for reads in bam.fetch():
		outbam.write(reads)
		read_count+=1
		if(read_count==200000):
			print('last chunk: {} reads\n'.format(read_count))
			chunk+=1
			print('Writing chunk {}'.format(chunk))
			read_count=0
			outbam.close()
			pysam.index("{}/{}".format(tmpdir,outfile))

			outfile='{}.{}.bam'.format(base,chunk)
			outbam=pysam.Samfile("{}/{}".format(tmpdir,outfile),'wb', header=header)
	outbam.close()
	print('last chunk: {} reads\n'.format(read_count))
	pysam.index("{}/{}".format(tmpdir,outfile))
	os.remove('{}/{}'.format(tmpdir,basename))
	return(chunk)

def count_chunk(chunk):

	"""
	count_chunk: Determines count of reads mapped to separate reference sequences
	within BAM chunk.

	Required parameters:
		chunk (int): Number of chunk to operate on

	Returns:
		counts_df (pd.DataFrame): dataframe of read counts
	"""

	tmpdir=os.environ['TMPDIR1']
	base=os.path.splitext(bam_file)[0]
	infile='{}.{}.bam'.format(base,chunk)		

	with lock:
		sys.stdout.write("Processing {}/{}\n".format(tmpdir,infile))

	bam=pysam.AlignmentFile('{}/{}'.format(tmpdir,infile), 'rb')

	# Certain high-level GO terms are excluded from comparison since they are so overly abundant
	# that they would lead to false associations of read pairs
	excluded_terms=('GO:0008150','GO:0007582', 'GO:0044699', 'GO:0000004','GO:0009987',
					'GO:0050875', 'GO:0044763', 'GO:0008151','GO:0008152','GO:0044710',
					'GO:0044236','GO:0044237')
	read_counts={}
	for read in bam.fetch():

		# counter for GO terms which are present on the references of both reads of the pair
		go_match=0 

		# just look at first read in pair to avoid double counting
		if read.is_read1:
			if read.reference_name in go_map['#query_name'].values and \
		   		read.next_reference_name in go_map['#query_name'].values:

				read1_terms=go_map['GOs'].loc[go_map['#query_name']==read.reference_name].to_string(index=False)
				read2_terms=go_map['GOs'].loc[go_map['#query_name']==read.next_reference_name].to_string(index=False)
				read1_terms=read1_terms.split(',')
				read2_terms=read2_terms.split(',')
            
				for term in read1_terms:
					if term not in excluded_terms:
						term=term.strip()
						if term in read2_terms:
							go_match+=1

				# Accept read-pairing if we have more terms than value defined by --common_terms argument
				if go_match>count_threshold:
					# count each read as 0.5 since we normally count read pairs
					if read.reference_name in read_counts:
						read_counts[read.reference_name]+=0.5
					else:
						read_counts[read.reference_name]=0.5

					if read.next_reference_name in read_counts:
						read_counts[read.next_reference_name]+=0.5
					else:
						read_counts[read.next_reference_name]=0.5

	counts_df=pd.DataFrame.from_dict(read_counts,orient='index')
	with lock:
		sys.stdout.write("{} returning {} counts\n".format(infile,len(counts_df.index)))

	return(counts_df)

def main():

	sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)

	parser = argparse.ArgumentParser(
		description="Count mispaired reads alignments where references share GO terms"
	)

	parser.add_argument('--eggnog', help='Output file from eggnog-mapper',
        required=True)
	parser.add_argument('--bam_dir', help='Path to directory of cdna alignment bam files. Files should be name *.diff_cdna.bam',
        required=True)
	parser.add_argument('--out_dir', help='Path to output directory', required=True)
	parser.add_argument('--common_terms', help='Number of common terms required for a match to be accepted', 
		required=False, default=0)

	args = parser.parse_args()

	global bam_file
	global go_map
	global count_threshold
	count_threshold=int(args.common_terms)

	go_map=pd.read_csv(args.eggnog, sep="\t",header=1,skiprows=2)
	go_map=go_map[['#query_name','GOs']]
	go_map=go_map[go_map['GOs'].notna()]

	task_id=int(os.environ['SGE_TASK_ID'])
	if not task_id or task_id=='undefined':
		print('This script must be submitted as an array job')
		sys.exit(1)

	if not os.path.isdir(args.out_dir):
		os.mkdir(args.out_dir)

	bam_files=glob.glob('{}/*diff_cdna.bam'.format(args.bam_dir))
	bam_path=bam_files[task_id-1]
	bam_file=os.path.basename(bam_path)
	sample=os.path.basename(bam_path).split('.')[0]

	tmp_file='{}/{}'.format(os.environ['TMPDIR1'],bam_file)
	index_path='{}.bai'.format(bam_path)
	tmp_index='{}.bai'.format(tmp_file)

	copyfile(bam_path,tmp_file)
	copyfile(index_path,tmp_index)

	chunks=split_bam(bam_path)

	pool=Pool(24)
	results=pool.map(count_chunk,[i for i in range(1,chunks+1)])
	pool.close()
	pool.join()

	sample_df=results.pop(0)
	for df in results:
		if (len(df.index)>0):
			sample_df.merge(df,how='outer',left_index=True,right_index=True)
	sample_df.columns=[sample]
	sample_df.to_csv('{}/{}.invalid_pairs.read_count_summary.txt'.format(args.out_dir,sample),sep='\t',index=True)

if __name__=='__main__':
	main()

