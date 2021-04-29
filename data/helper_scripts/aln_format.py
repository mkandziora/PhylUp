# -*- coding: utf-8 -*-
"""
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

change format of alignment to fasta
"""
import sys
import os
import dendropy
import argparse

input_aln = sys.argv[1]
input_frmt = sys.argv[2]


# create variables that can be entered in the command line
parser = argparse.ArgumentParser(usage='''
aln_format.py [<args>]
The aln_format.py arguments are:
  -aln   STR     path to alignment (required)
  -f     STR     format of alignment (required)
''')
parser.add_argument('-aln', type = str, required = True, help = 'path to alignment ')
parser.add_argument('-f', type = str, required = True, help = 'format of alignment')


args = parser.parse_args()

input_aln = os.path.abspath(args.aln)

aln = dendropy.DnaCharacterMatrix.get(path=input_aln, schema=args.f)
for tax in aln.taxon_namespace:
    tax.label = tax.label.replace(" ", "_")

aln.write(path='{}.fasta'.format(input_aln), schema="fasta", preserve_underscores=True)
