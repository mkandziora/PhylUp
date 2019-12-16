# package import
from PhylUp import phyl_up, config



def test_mrca():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)
    print(test.mrca)

    test.mrca == 18794

# #################################
def test_mrca_list():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    ingroup_mrca = [18794, 422320, 422327, 422329, 422331]
    # ingroup_mrca = [senecio, culcitium, hasteola, iocenes, lasiocephalus]
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=ingroup_mrca)

    test.mrca == 795077  # Senecionineae

def test_no_mrca():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=None)
    print(test.mrca)

    assert test.mrca == 18794

def test_mrca_outgroup():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    ingroup_mrca = [18794, 422320, 422327, 422329, 422331, 84584]
    # ingroup_mrca = [senecio, culcitium, hasteola, iocenes, lasiocephalus, ABROTANELLA]
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=ingroup_mrca)

    test.mrca != 795077  # Senecioninea
    test.mrca == 102812  # Senecioneae



