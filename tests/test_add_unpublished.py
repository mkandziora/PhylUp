import os
from PhylUp import phyl_up, config

def test_add_unpublished():
    print('test unpublished seq')
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.unpublished = True
    conf.unpubl_data = 'data/unpublished_seqs/'
    conf.unpubl_names = 'data/unpublished_names.csv'

    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)


    new_seqs = test.extend_with_unpublished()
    print(new_seqs)
    assert len(new_seqs) > 0