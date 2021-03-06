"""
Martha Kandziora
martha.kandziora@mailbox.org

Split concatenated alignment into loci - to be run where the file is located

"""

import sys
import os
import dendropy
import argparse

# create variables that can be entered in the command line
parser = argparse.ArgumentParser(description='Split concatenated alignment into loci - to be run where the file is located')

parser.add_argument('-aln', required=True, help=' path to alignment')
parser.add_argument('-f', required=True, help='format alignment')
parser.add_argument('-start', required=True, help='start of locus')
parser.add_argument('-stop', required=True, help='end of locus')
parser.add_argument('-name', required=True, help='name of locus')

args = parser.parse_args()

input_aln = os.path.abspath(args.aln)

prune = []
aln = dendropy.DnaCharacterMatrix.get(path=input_aln, schema=args.f)
for tax in aln.taxon_namespace:
    aln[tax] = aln[tax][args.start:args.stop]
    tax.label = tax.label.replace(" ", "_")

    seq = aln[tax].symbols_as_string().replace('?', '')
    if len(seq) == 0:
        prune.append(tax)
aln.discard_sequences(prune)     
aln.write(path='{}.fasta'.format(args.name), schema="fasta", preserve_underscores=True)

