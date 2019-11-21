from PhylUp import phyl_up, config, phylogen_updater

def test_blacklist():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    blacklist = ['JX895419.1']

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

    new_seqs = test.extend()
    len_no_bl = len(new_seqs)

    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794, blacklist=blacklist)
    new_seqs = test.extend()
    len_bl = len(new_seqs)

    assert len_no_bl > len_bl






