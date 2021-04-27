from PhylUp import phyl_up, config

workdir = "./example/run_singleseq"  # working directory
id_to_spn = "data/single_seq_input/single_seq_spn.csv"  # tab-delimited file where tip names correspond to ncbi names
seqaln = "data/single_seq_input/single_seq.fas"  # alignment
mattype = "fasta"  # format of alignment
configfi = "example_setups/localblast.config"  # configuration file

ignore_acc = ['JX895419.1']  # sequence accession numbers that shall not be added


# update the data:
conf = config.ConfigObj(configfi, workdir, interactive=False)
example = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, None, None, conf, ignore_acc_list=ignore_acc)
example.run(status_end=1) # number of blast rounds: here a single blast run
