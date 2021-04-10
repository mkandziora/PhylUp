from PhylUp import wrapper, phylogenetic_helpers
from distutils.dir_util import copy_tree
from dendropy import Tree, DnaCharacterMatrix
import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser

import os
from distutils.dir_util import copy_tree

from PhylUp import phyl_up, config, blast



def test_prefer():
    workdir = "tests/output/test_run_preferred"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast_test.config"

    tmp_folder = os.path.join(workdir, 'tmp')
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    # call(['cp', '-a', 'data/tmp_for_test/', tmp_folder])
    copy_tree('data/tmp_for_test/', tmp_folder)

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    conf.preferred_taxa = True
    conf.allow_parent = False

    conf.downtorank = "species"
    conf.preferred_taxa_fn = "data/tiny_test_example/preferred_taxa_test.txt"
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)
    new_seqs = test.extend()
    aln = phylogenetic_helpers.read_in_aln(test.aln_fn, test.aln_schema)
    taxids_original = set(test.table.ncbi_txid)
    new_seqs = new_seqs[~new_seqs['accession'].isin(test.table['accession'])]  # ~ is the pd not in/!
    new_seqs = test.basic_filters(aln, test.mrca, new_seqs)
    # next filter need infos in table
    new_seqs = test.add_new_seqs(new_seqs)

    assert len(new_seqs.index) > 0, new_seqs.index
    # print(new_seqs)
    for i in range(0, len(new_seqs.ncbi_txid)):
        try:
            # print(i)
            new_seqs.ncbi_txid[i] = int(new_seqs.ncbi_txid[i])
        except:
            print('wrong index')

    # print(new_seqs.columns)
    taxids_new_seqs = set(new_seqs.ncbi_txid.tolist())

    found = False
    for item in [462523, 1268580, 1268581, 1268589]:
        if item in taxids_new_seqs:
            found = True
    assert found == True

    f = phyl_up.FilterPreferredTaxa(test.config, test.table, test.status)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs
    assert len(new_seqs.index) > 0, new_seqs.index

    taxids_after = set(new_seqs.ncbi_txid.tolist())

    preferred_taxa = f.get_preferred_ids()
    assert len(taxids_new_seqs) >= len(taxids_after), (len(taxids_new_seqs), len(taxids_after))
    assert len(taxids_after) <= len(preferred_taxa), (len(taxids_after), len(taxids_after))
    print(preferred_taxa)
    for i in taxids_after:
        assert i in preferred_taxa, (i, preferred_taxa)

def fixtest_preferred():
    copy_tree('tests/output/test_runs/', 'tests/output/test_run_preferred/test_its')
    copy_tree('./data/blast_for_tests', 'tests/output/test_run_preferred/test_its/blast')
    copy_tree('./data/blast_for_test_ets', 'tests/output/test_run_preferred/test_ets/blast')
    copy_tree('./data/tmp_for_test_ETS', 'tests/output/test_run_preferred/test_ets/tmp')
    copy_tree('./data/tiny_test_ETS/', 'tests/output/test_run_preferred/test_ets')


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
    conf = config.ConfigObj(confs[0], workdir, interactive=False)

    data = {
        'its': dict1,
        'ets': dict2,
    }

    overlap_folder = "tests/output/test_run_preferred/"

    ############################################
    #todo fix run multiple - need to copy from somewhere
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
        txid = blast.get_taxid_from_acc(tax.label, conf.blastdb_path, workdir)
        taxid_its.append(txid[0])

    for tax, seq in new_seqs_ets.items():
        txid = blast.get_taxid_from_acc(tax.label, conf.blastdb_path,
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

