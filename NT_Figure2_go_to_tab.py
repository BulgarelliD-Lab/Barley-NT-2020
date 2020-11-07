#!/bin/env python

import re

term=None
namespace=None
name=None

term_re=re.compile('^id: (GO:\d+)')
namespace_re=re.compile('^namespace: (\S+)')
name_re=re.compile('^name: (.+)')

file = open('go_mapping/go-basic.obo','r')
outfile = open('go_mapping/go_table.txt','w')

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
