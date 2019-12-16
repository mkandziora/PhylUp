from PhylUp import phyl_up, config

workdir = "tests/output/downtorank_sp"
trfn = "data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = r"data/tiny_test_example/test_nicespl.csv"
seqaln = "data/tiny_test_example/test.fas"
mattype = "fasta"

configfi = "data/localblast.config"

conf = config.ConfigObj(configfi, workdir, interactive=False)
conf.downtorank = 'genus'

test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=102812)
test.run()


