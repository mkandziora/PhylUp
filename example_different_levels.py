from PhylUp import phyl_up, config


# round 1 Senecio
workdir = "tests/output/test_different_level"
trfn = "data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = "data/tiny_test_example/test_nicespl.csv"
seqaln = "data/tiny_test_example/test.fas"
mattype = "fasta"
configfi = "data/localblast_Senecio.config"
blacklist = ['JX895419.1']

conf = config.ConfigObj(configfi, workdir, interactive=False)
test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, blacklist=blacklist)
test.run()



# round 2 Senecioneae

seqaln = "tests/output/different_level/updt_aln.fasta"
configfi = "data/localblast_Senecioneae.config"
blacklist = ['JX895419.1']

conf = config.ConfigObj(configfi, workdir, interactive=False)
test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, blacklist=blacklist)
test.run()



# round 3 Asteroideae

seqaln = "tests/output/different_level/updt_aln.fasta"
configfi = "data/localblast_Asteroideae.config"
blacklist = ['JX895419.1']

conf = config.ConfigObj(configfi, workdir, interactive=False)
test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, blacklist=blacklist)
test.run()


