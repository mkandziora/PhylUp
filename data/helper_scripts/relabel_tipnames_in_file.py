"""
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

Relabel tip names in files using the updated.table information
"""

import pandas as pd
import os
import argparse
from PhylUp.phylogenetic_helpers import replace_uid_with_name


# create variables that can be entered in the command line
parser = argparse.ArgumentParser(usage='''
check_spn.py [<args>]
The check_spn.py arguments is:
  -wd   STR     path to working directory (required)
  -fn   STR     file name that is to be relabeled (required)
  -f    STR     type of file: aln (alignment) or tree (required)
''')

parser.add_argument('-wd', required = True, help = 'path to working directory ')
parser.add_argument('-fn', required = True, help = 'file name that is to be relabeled')
parser.add_argument('-f', required = True, help = 'type of file: aln (alignment) or tree (required)')


args = parser.parse_args()

assert args.f in ['tree', 'aln']

fullpath = os.path.join(args.wd, args.fn)
table = pd.read_csv(os.path.join(args.wd, 'table.updated'))
replace_uid_with_name(fullpath, table, args.f)
