"""
Martha Kandziora
martha.kandziora@mailbox.org

Get ncbi taxon name from species name

"""
from ncbiTAXONparser import ncbi_data_parser
import argparse

# create variables that can be entered in the command line
parser = argparse.ArgumentParser(description='Get ncbi taxon name from species name')

parser.add_argument('-spn', required=True, help='name of species')

args = parser.parse_args()
ncbi_parser = ncbi_data_parser.Parser(names_file="./data/nodes.dmp",
                                      nodes_file="./data/names.dmp")

taxid = ncbi_parser.get_id_from_name(args.spn)
print(taxid)
