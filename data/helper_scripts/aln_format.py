"""
Martha Kandziora
martha.kandziora@mailbox.org

change format of alignment to fasta
"""
import os
import dendropy
import argparse


# create variables that can be entered in the command line
parser = argparse.ArgumentParser(description='change format of an alignment to fasta')

parser.add_argument('-aln', type=str, required=True, help='path to alignment ')
parser.add_argument('-f', type=str, required=True, help='format of alignment')


args = parser.parse_args()

input_aln = os.path.abspath(args.aln)

aln = dendropy.DnaCharacterMatrix.get(path=input_aln, schema=args.f)
for tax in aln.taxon_namespace:
    tax.label = tax.label.replace(" ", "_")

aln.write(path='{}.fasta'.format(input_aln), schema="fasta", preserve_underscores=True)
