"""
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

make a spn name file from alignment - test if all names are accepted using check_spn.py
"""
from ncbiTAXONparser import ncbi_data_parser
import argparse
import csv
import os
from Bio import SeqIO


# create variables that can be entered in the command line
parser = argparse.ArgumentParser(usage='''
check_spn.py [<args>]
The check_spn.py arguments is:
  -aln   STR     path to file with alignment (required)
''')

parser.add_argument('-aln', required = True, help = 'path to file with alignment')


args = parser.parse_args()

ncbi_parser = ncbi_data_parser.Parser(names_file='./data/names.dmp', nodes_file="./data/nodes.dmp")

output_spn = os.path.abspath('spn_list.csv')

print('Species name file is called spn_list.csv')
with open(output_spn, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    records = list(SeqIO.parse(args.aln, "fasta"))
    for record in records:
        print(record.name)
        spn = record.name.replace(':', '_').replace('.','').replace('-','_').replace(' ','_')
        writer.writerow(["{}, {}".format(record.name, spn)])
