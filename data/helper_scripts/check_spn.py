"""
Martha Kandziora
martha.kandziora@mailbox.org

check if names provided are accepted by ncbi - prints names that are not
"""

from ncbiTAXONparser import ncbi_data_parser
import pandas as pd
import argparse

# create variables that can be entered in the command line
parser = argparse.ArgumentParser(description='check if names in files are accepted by ncbi - prints names that are not')

parser.add_argument('-txn_list', required=True, help='path to file with sample name and ncbi name')
args = parser.parse_args()

tax_col = ["seq_name",
           "ncbi_name"]
taxonname_list = pd.read_csv(args.txn_list, header=None, index_col=False, names=tax_col)

ncbi_parser = ncbi_data_parser.Parser(names_file='./data/names.dmp',
                                      nodes_file='./data/nodes.dmp', interactive=False)

for name in taxonname_list['ncbi_name']:
    tax_id = ncbi_parser.get_id_from_name(name)
    valid = ncbi_parser.taxid_is_valid(tax_id)

    if valid == False:
        print(name)