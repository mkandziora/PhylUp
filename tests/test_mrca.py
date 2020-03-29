# package import
import os
from distutils.dir_util import copy_tree

from PhylUp import phyl_up, config



def test_mrca():
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
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)
    print(test.mrca)

    assert test.mrca == {18794}

# #################################
def test_mrca_list():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast_test.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    conf.mrca_input = "18794, 422320, 422327, 422329, 422331"
    # ingroup_mrca = [senecio, culcitium, hasteola, iocenes, lasiocephalus]
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)

    assert test.mrca == {18794, 422320, 422327, 422329, 422331}  # Senecionineae

def test_no_mrca():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast_test.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.mrca_input = None
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)
    print(test.mrca)

    assert test.mrca == {18794}



