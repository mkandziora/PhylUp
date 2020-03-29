import os
from distutils.dir_util import copy_tree

from PhylUp import phyl_up, config

def test_filter_thresh():
    print('test_filter_unique')
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast_test.config"

    tmp_folder = os.path.join(workdir, 'tmp')
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    #call(['cp', '-a', 'data/tmp_for_test/', tmp_folder])
    copy_tree('data/tmp_for_test/', tmp_folder)

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)

    new_seqs = test.extend()

    before = len(new_seqs)

    f = phyl_up.FilterBLASTThreshold(test.config)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    del_tab = len(f.del_table)
    after = len(new_seqs)
    print(before, after, del_tab)
    assert before >= after
    assert del_tab >= 0
    assert after + del_tab == before

    # note: test never finds seq below threshold...
