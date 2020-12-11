#!/bin/env python

"""
Parses obo format ontologies into a tab-delimited text file

Run `obo_to_tab.py` --help for usage
"""

import argparse
import re

def parse_file(infile:str, outfile:str):

	"""
	Parses obo file specified as infile, writing outputs to outfile

	Required parameters:
		infile (str): path to obo input file
		outfile (str): path to output file (txt)
	"""

	term=None
	namespace=None
	name=None

	term_re=re.compile('^id: (GO:\d+)')
	namespace_re=re.compile('^namespace: (\S+)')
	name_re=re.compile('^name: (.+)')

	file = open(infile,'r')
	outfile = open(outfile,'w')

	for line in file:
		if '[Term]' in line:
			if term is not None:
				outfile.write("{}\t{}\t{}\n".format(term, namespace, name))
				term=None
				namespace=None
				name=None

		m=term_re.match(line)
		if m:
			term=m.group(1)

		m=namespace_re.match(line)
		if m:
			namespace=m.group(1)

		m=name_re.match(line)
		if m:
			name=m.group(1)

def main():

	parser = argparse.ArgumentParser(
		description="Create counts of GO terms per-sample from per-read summary count tables"
	)

	parser.add_argument('--infile', help='Path to obo format ontology file', required=True)
	parser.add_argument('--outfile', help='Path to output file', required=True)

	args = parser.parse_args()

	parse_file(args.infile, args.outfile)

if __name__=='__main__':
	main()