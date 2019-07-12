from pandas_filter import pandas_numpy_try1, config, aln_updater

def test_filter_otu_no_rank():
    workdir = "test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

    new_seqs = None
    new_seqs = test.extend(new_seqs)

    aln = test.read_in_aln()
    f = pandas_numpy_try1.FilterUniqueAcc(test.config)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    before = len(new_seqs)
    new_seqs_org = new_seqs
    f = pandas_numpy_try1.FilterNumberOtu(test.config, test.table, aln)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    del_tab = len(f.del_table)
    after = len(new_seqs)
    print(before, after, del_tab)
    assert before > after
    assert del_tab > 0

    # with downtorank
    new_seqs = new_seqs_org
    before_dtr = len(new_seqs)
    aln = test.read_in_aln()

    f = pandas_numpy_try1.FilterNumberOtu(test.config, test.table, aln)
    f.filter(new_seqs, 'genus')
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
