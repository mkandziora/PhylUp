from PhylUp import phyl_up, config

workdir = "tests/output/test_runs"  # working directory
trfn = "data/tiny_test_example/test.tre"  # phylogeny
schema_trf = "newick"  # format of phylogeny
id_to_spn = "data/tiny_test_example/test_nicespl.csv"  # tab-delimited file where tip names correspond to ncbi names
seqaln = "data/tiny_test_example/test.fas"  # alignment
mattype = "fasta"  # format of alignment
configfi = "data/localblast.config"  # configuration file

blacklist = ['JX895419.1']  # sequence accession numbers that shall not be added

# update the data:
conf = config.ConfigObj(configfi, workdir, interactive=False)

test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, blacklist=blacklist)
test.run(status_end=1)

