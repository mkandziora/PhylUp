import os
from pandas_filter import pandas_numpy_try1, config

print('test_filter_unique')
workdir = "tests/test_runs"
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

test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

new_seqs = test.run()




def test_add_unpublished():
    print('test_filter_unique')
    workdir = "tests/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.unpublished = True
    conf.unpubl_data = 'data/unpublished_seqs/'
    conf.unpubl_names = 'data/unpublished_seqs/unpublished_names.csv'

    test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)
