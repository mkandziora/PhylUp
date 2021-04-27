from PhylUp import wrapper

workdir = "./example/different_ranks"
trfn = "data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = "data/tiny_test_example/test_nicespl.csv"
seqaln = "data/tiny_test_example/test.fas"
mattype = "fasta"
ignore_acc = ['JX895419.1']

dict1 = {'workdir': workdir,
         'idtospn': id_to_spn,
         'seqaln': seqaln,
         'mattype': mattype,  # format of alignment,
         'trfn': trfn,  # phylogeny
         'schema_trf': schema_trf,  # format of phylogeny,
         'ignore_acc_list': ignore_acc
         }

data = {
        'ITS':  dict1
        }

# define the configuration files to use, minimum one has to be provided -
# if there are different sampling strategies across different taxonomic ranks,
# there need to be multiple config files- one per each rank
confs = ["example_setups/localblast_Senecio.config",
         "example_setups/localblast_Senecioneae.config",
         "example_setups/localblast_Asteroideae.config"
         ]

# end defines the number of blast rounds. is limited to one if you want to match taxa across loci
wrapper.run_multiple(data, confs, end=None, overlap_folder=None)
