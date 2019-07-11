
import pandas as pd

from pandas_filter import pandas_numpy_try1, config, aln_updater

def test_blacklist():
    workdir = "tests/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    blacklist = ['JX895419.1']

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

    new_seqs = None
    new_seqs = test.extend(new_seqs)
    len_no_bl = len(new_seqs)

    test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794, blacklist=blacklist)
    new_seqs = None
    new_seqs = test.extend(new_seqs)
    len_bl = len(new_seqs)

    assert len_no_bl > len_bl






