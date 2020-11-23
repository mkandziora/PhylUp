import os
from distutils.dir_util import copy_tree

from PhylUp import phyl_up, config

def test_add_unpublished():
    print('test unpublished seq')
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast_test.config"

    if not os.path.exists(workdir):
        os.mkdir(workdir)
    tmp_folder = os.path.join(workdir, 'tmp')
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    #call(['cp', '-a', 'data/tmp_for_test/', tmp_folder])
    copy_tree('data/tmp_for_test/', tmp_folder)


    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.unpublished = True
    conf.unpubl_data = 'data/unpublished_seqs/'
    conf.unpubl_names = 'data/unpublished_names.csv'
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")

    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)

    new_seqs = test.extend_with_unpublished()
    print(new_seqs)
    assert len(new_seqs) > 0
