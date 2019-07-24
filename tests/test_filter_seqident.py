from pandas_filter import pandas_numpy_try1, config, aln_updater

def test_filter_seqident():
    print('test_filter_seqident')
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

    present = test.table['status'] != 'deleted'
    old_status = test.table['status'] < test.status
    comp_table = test.table[present & old_status]

    seqs_to_comp = comp_table[['ncbi_txid', 'tip_name', 'sseq']]


    f = pandas_numpy_try1.FilterSeqIdent(test.config, seqs_to_comp, test.table, test.status)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    del_tab = len(f.del_table)
    after = len(new_seqs)
    print(before, after, del_tab)
    assert before >= after
    assert del_tab >= 0
    assert after + del_tab == before
