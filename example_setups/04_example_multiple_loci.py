from PhylUp import wrapper

# define the different loci to be updated and matched

ignore_acc = []  # sequence accession numbers that shall not be added
workdir = "./example/multiple_loci/ITS"
id_to_spn = "data/tiny_test_example/test_nicespl.csv"
seqaln = "data/tiny_test_example/test.fas"
mattype = "fasta"

dict1 = {'workdir': workdir,
         'idtospn': id_to_spn,
         'seqaln': seqaln,
         'mattype': mattype,  # format of alignment,
         'trfn': None,  # phylogeny
         'schema_trf': None,  # format of phylogeny,
         'ignore_acc_list': ignore_acc
         }

workdir = "./example/multiple_loci/ETS"
id_to_spn = "data/tiny_test_ETS/nicespl.csv"
seqaln = "data/tiny_test_ETS/test_ets.fasta"

dict2 =  {'workdir': workdir,
         'idtospn': id_to_spn,
         'seqaln': seqaln,
         'mattype': mattype,  # format of alignment,
         'trfn': None,  # phylogeny
         'schema_trf': None,  # format of phylogeny,
         'ignore_acc_list': ignore_acc
         }


# define the configuration files to use, min one has to be provided
confs = ["example_setups/localblast_preferred.config"]


data = {
        'ITS':dict1,
        'ETS': dict2,
        }

# folder defines where to store information that are needed to match taxa across loci
overlap_folder = "./example/multiple_loci/"


############################################
# end defines the number of blast rounds. Currently, it is automatically limited to one if you want to match taxa across loci
wrapper.run_multiple(data, confs, end=None, overlap_folder=overlap_folder, first_locus=True)
