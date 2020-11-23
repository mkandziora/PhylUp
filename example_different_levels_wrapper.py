from PhylUp import wrapper

workdir = "tests/output/test_different_level"
trfn = "data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = "data/tiny_test_example/test_nicespl.csv"
seqaln = "data/tiny_test_example/test.fas"
mattype = "fasta"
configfi = "data/localblast_Senecio.config"
ignore_acc = ['JX895419.1']

dict1 = {'workdir': workdir,
          'idtospn': id_to_spn,
          'seqaln': seqaln,
          'mattype': mattype ,  # format of alignment,
          'trfn': trfn, # phylogeny
          'schema_trf': schema_trf,  # format of phylogeny,
          'ignore_acc_list': ignore_acc
          }

data = {
        'ITS':  dict1
        }

# define the configuration files to use, min one has to be provided
confs = ["data/localblast_Senecio.config",
         "data/localblast_Senecioneae.config",
         "data/localblast_Asteroideae.config"
         ]

# end defines the number of blast rounds. is limited to one if you want to match taxa across loci
wrapper.run_multiple(data, confs, end=None, overlap_folder=None)
