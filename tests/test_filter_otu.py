import os
import pandas as pd
import pytest

from distutils.dir_util import copy_tree

from PhylUp import phyl_up, config, phylogenetic_helpers
from copy import deepcopy


@pytest.fixture(autouse=True)
def configure():
    workdir = "tests/output/test_runs_filter_otu"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"

    configfi = "data/localblast_test.config"
    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    conf.preferred_taxa = True
    conf.allow_parent = False
    conf.preferred_taxa_fn = "data/tiny_test_example/preferred_taxa_test.txt"
    pytest.conf = conf

    tmp_folder = os.path.join(workdir, 'tmp')
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    # call(['cp', '-a', 'data/tmp_for_test/', tmp_folder])
    copy_tree('data/tmp_for_test/', tmp_folder)

    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, pytest.conf)
    print(test)
    pytest.test = test

    new_seqs = test.extend()
    assert len(new_seqs) > 0, len(new_seqs)

    aln = phylogenetic_helpers.read_in_aln(test.aln_fn, test.aln_schema)

    new_seqs = new_seqs[~new_seqs['accession'].isin(test.table['accession'])]  # ~ is the pd not in/!
    assert len(new_seqs) > 0, len(new_seqs)

    #new_seqs = test.basic_filters(aln, test.mrca, new_seqs)
    orig_len = new_seqs
    columns = ['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date",
               'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title']
    all_del = pd.DataFrame(columns=columns)

    remove_basics = [phyl_up.FilterUniqueAcc(test.config, test.table),
                     # FilterBLASTThreshold(self.config),
                     phyl_up.FilterLength(test.config, aln),
                     phyl_up.FilterMRCA(test.config, test.mrca)
                     ]
    for f in remove_basics:
        internal_new_seq = new_seqs
        f.filter(new_seqs)
        new_seqs = f.upd_new_seqs
        del_seq = f.del_table
        # all_del = all_del.append(del_seq, ignore_index=True)
        all_del = pd.concat([all_del, del_seq], ignore_index=True, sort=True)
        phyl_up.check_filter_numbers(del_seq, new_seqs, internal_new_seq)
        assert len(new_seqs) > 0, len(new_seqs)

        if len(new_seqs) == 0:
            return new_seqs  # stop filtering if new seqs is 0
    phyl_up.check_filter_numbers(all_del, new_seqs, orig_len)
    assert len(new_seqs) > 0, len(new_seqs)

    # next filter need infos in table
    new_seqs = test.add_new_seqs(new_seqs)
    pytest.new_seqs = new_seqs

    assert len(pytest.new_seqs) > 0, len(pytest.new_seqs)


def test_filter_otu_no_rank():
    test = pytest.test
    print(test.table)
    new_seqs = pytest.new_seqs
    assert len(new_seqs) > 0, len(new_seqs)

    conf = pytest.conf


    before = len(new_seqs)
    print(before)
    new_seqs_org = new_seqs
    before_table = deepcopy(test.table)

    f = phyl_up.FilterNumberOtu(conf, test.table, test.status)
    #assert len(test.table[test.table['status'] == test.status]) == 0, test.table[test.table['status'] == test.status]
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    del_tab = len(f.del_table)
    after = len(new_seqs)


    assert before > after
    assert del_tab > 0

    # with downtorank
    new_seqs = new_seqs_org
    before_dtr = len(new_seqs)

    f = phyl_up.FilterNumberOtu(conf, before_table, test.status)
    f.filter(new_seqs_org, 'genus')
    new_seqs = f.upd_new_seqs

    del_tab_dtr = len(f.del_table)
    after_dtr = len(new_seqs)

    print(before, after, del_tab)
    print(before_dtr, after_dtr, del_tab_dtr)

    print(new_seqs[['ncbi_txn', 'accession', 'ncbi_txid']])
    assert before_dtr > after_dtr
    assert del_tab_dtr > 0

    assert after > after_dtr
    assert del_tab < del_tab_dtr


def test_get_preferred_taxa():
    test = pytest.test
    new_seqs = pytest.new_seqs

    taxids_before = test.table.ncbi_txid.tolist()


    f = phyl_up.FilterNumberOtu(test.config, test.table, test.status)
    # assert len(test.table[test.table['status'] == test.status]) == 0, test.table[test.table['status'] == test.status]
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    taxids_after = new_seqs.ncbi_txid.tolist()

    assert taxids_before not in taxids_after

def test_prefer_different_OTU():
    test = pytest.test
    new_seqs = pytest.new_seqs

    taxids_before = test.table.ncbi_txid.tolist()

    f = phyl_up.FilterNumberOtu(test.config, test.table, test.status)
    # assert len(test.table[test.table['status'] == test.status]) == 0, test.table[test.table['status'] == test.status]
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    taxids_after = new_seqs.ncbi_txid.tolist()

    assert taxids_before not in taxids_after

