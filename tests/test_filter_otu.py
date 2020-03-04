import os
from PhylUp import phyl_up, config, phylogenetic_helpers
from copy import deepcopy


def test_filter_otu_no_rank():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")

    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)

    new_seqs = test.extend()
    aln = phylogenetic_helpers.read_in_aln(test.aln_fn, test.aln_schema)

    new_seqs = new_seqs[~new_seqs['accession'].isin(test.table['accession'])]  # ~ is the pd not in/!
    new_seqs = test.basic_filters(aln, test.mrca, new_seqs)
    # next filter need infos in table
    new_seqs = test.add_new_seqs(new_seqs)

    before = len(new_seqs)
    new_seqs_org = new_seqs
    before_table = deepcopy(test.table)

    f = phyl_up.FilterNumberOtu(test.config, test.table, test.status)
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

    f = phyl_up.FilterNumberOtu(test.config, before_table, test.status)
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

