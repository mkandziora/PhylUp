import sys
from ncbiTAXONparser import ncbi_data_parser

spn = sys.argv[1]
nodes = "./data/nodes.dmp"
names = "./data/names.dmp"


ncbi_parser = ncbi_data_parser.Parser(names_file=names,
                                      nodes_file=nodes)

taxid = ncbi_parser.get_id_from_name(spn)
return taxid
