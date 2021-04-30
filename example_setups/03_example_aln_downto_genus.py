from PhylUp import phyl_up, config

workdir = "./example/downtorank_genus"
trfn = "data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = r"data/tiny_test_example/test_nicespl.csv"
seqaln = "data/tiny_test_example/test.fas"
mattype = "fasta"

# in this configuration file the threshold per OTU is based on the genus rank
# sequences will be added until the defined number of samples per genus are reached
configfi = "example_setups/localblast_rankgenus.config"

conf = config.ConfigObj(configfi, workdir, interactive=False)

test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)
test.run(status_end=2)


