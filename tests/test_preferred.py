from PhylUp import wrapper
from distutils.dir_util import copy_tree
from dendropy import Tree, DnaCharacterMatrix
import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser

import os
from distutils.dir_util import copy_tree

from PhylUp import phyl_up, config, blast



def test_preferred():
    copy_tree('tests/output/test_runs/', 'tests/output/test_run_preferred/test_its')
    copy_tree('./data/blast_for_tests', 'tests/output/test_run_preferred/test_its/blast')
    copy_tree('./data/blast_for_test_ets', 'tests/output/test_run_preferred/test_ets/blast')
    copy_tree('./data/tmp_for_test_ETS', 'tests/output/test_run_preferred/test_ets/tmp')

    workdir = "tests/output/test_run_preferred/test_its"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast_test_preferred.config"

    ignore_acc1 = []  # sequence accession numbers that shall not be added

    # define the different loci to be updated and matched
    dict1 = {'workdir': workdir,
             'idtospn': id_to_spn,
             'seqaln': seqaln,
             'mattype': mattype,
             'trfn': trfn,
             'schema_trf': schema_trf,
             'ignore_acc_list': ignore_acc1
             }

    dict2 = {'workdir': "tests/output/test_run_preferred/test_ets/",
             'idtospn': "data/tiny_test_ETS/nicespl.csv",
             'seqaln': "data/tiny_test_ETS/test_ets.fasta",
             'mattype': "fasta",
             'trfn': 'data/tiny_test_ETS/test_ets.tre',
             'schema_trf': 'newick',
             'ignore_acc_list': ignore_acc1
             }

    # define the configuration files to use, min one has to be provided
    confs = [configfi
             ]

    data = {
        'its': dict1,
        'ets': dict2,
    }

    overlap_folder = "tests/output/test_run_preferred/"

    ############################################
    # end defines the number of blast rounds. is limited to one if you want to match taxa across loci
    wrapper.run_multiple(data, confs, end=0, overlap_folder=overlap_folder)

    print('do the actual test')
    new_seqs_its = DnaCharacterMatrix.get(path="{}/new_seqs.fasta".format(workdir), schema='fasta')
    new_seqs_ets = DnaCharacterMatrix.get(path="tests/output/test_run_preferred/test_ets/new_seqs.fasta",
                                          schema='fasta')

    print(len(new_seqs_its))
    taxid_its = []
    taxid_ets = []
    for tax, seq in new_seqs_its.items():
        txid = blast.get_taxid_from_acc(tax.label, '/media/blubb/schmuh/local_blast_db_new/nt', workdir)
        taxid_its.append(txid[0])

    for tax, seq in new_seqs_ets.items():
        txid = blast.get_taxid_from_acc(tax.label, '/media/blubb/schmuh/local_blast_db_new/nt',
                                        'tests/output/test_run_preferred/test_ets')
        taxid_ets.append(txid[0])
    print(taxid_its)

    with open("data/tiny_test_example/preferred_taxa_test.txt", 'r') as content:
        preferred_taxa_ids = content.read().splitlines()
        if preferred_taxa_ids[1].isnumeric():
            for i in range(0, len(preferred_taxa_ids)):

                if '.' in preferred_taxa_ids[i]:
                    int(preferred_taxa_ids[i].split('.')[0])
                else:
                    preferred_taxa_ids[i] = int(preferred_taxa_ids[i])
        else:
            ncbi_parser = ncbi_data_parser.Parser(names_file="./data/names.dmp",
                                                  nodes_file="./data/nodes.dmp")

            preferred_taxa_ids = list(map(ncbi_parser.get_id_from_name, preferred_taxa_ids))
    print(preferred_taxa_ids)

    for i in taxid_its:
        assert i in preferred_taxa_ids, i
    for i in taxid_ets:
        assert i in preferred_taxa_ids, i

