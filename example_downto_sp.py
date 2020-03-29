from PhylUp import phyl_up, config

workdir = "tests/output/downtorank_genus"
trfn = "data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = r"data/tiny_test_example/test_nicespl.csv"
seqaln = "data/tiny_test_example/test.fas"
mattype = "fasta"

configfi = "data/localblast_genusexample.config"

conf = config.ConfigObj(configfi, workdir, interactive=False)

test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)
test.run()


