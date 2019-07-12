from pandas_filter import pandas_numpy_try1, config

def test_filter_thresh():
    print('test_filter_unique')
    workdir = "tests/test_runs"
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

    before = len(new_seqs)

    f = pandas_numpy_try1.FilterBLASTThreshold(test.config)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    del_tab = len(f.del_table)
    after = len(new_seqs)
    print(before, after, del_tab)
    assert before >= after
    assert del_tab >= 0
    assert after + del_tab == before

    # note: test never finds seq below threshold...