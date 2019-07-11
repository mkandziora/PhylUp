# package import
from pandas_filter import pandas_numpy_try1, config



def test_mrca():
    workdir = "tests/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)
    print(test.mrca)

    test.mrca == 18794

# #################################
def test_mrca_list():
    workdir = "test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    ingroup_mrca = [18794, 422320, 422327, 422329, 422331]
    # ingroup_mrca = [senecio, culcitium, hasteola, iocenes, lasiocephalus]
    test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=ingroup_mrca)

    test.mrca == 795077  # Senecionineae

def test_no_mrca():
    workdir = "test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)

    test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=None)
    print(test.mrca)

    assert test.mrca == 18794

def test_mrca_outgroup():
    workdir = "test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)

    ingroup_mrca = [18794, 422320, 422327, 422329, 422331, 84584]
    # ingroup_mrca = [senecio, culcitium, hasteola, iocenes, lasiocephalus, ABROTANELLA]
    test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=ingroup_mrca)

    test.mrca != 795077  # Senecioninea
    test.mrca == 102812  # Senecioneae





#
#
# # todo: reimplement
# # needed for palms idea
# def test_mrca_modify_later():
#     seqaln = "tests/data/tiny_test_example/test.fas"
#     mattype = "fasta"
#     trfn = "tests/data/tiny_test_example/test.tre"
#     schema_trf = "newick"
#     id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
#     workdir = "tests/output/test_mrcalist_modify"
#     configfi = "tests/data/test.config"
#     otu_jsonfi = "{}/otu_dict.json".format(workdir)
#
#     ingroup_mrca = 723076
#
#     # setup the run
#     if not os.path.exists("{}".format(workdir)):
#         os.makedirs("{}".format(workdir))
#
#     conf = ConfigObj(configfi)
#     ids = IdDicts(conf, workdir=workdir)
#
#     # print(ids.mrca_ott, ids.mrca_ncbi)
#
#     data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
#     filteredScrape = FilterBlast(data_obj, ids, ingroup_mrca)
#
#
#     ingroup_mrca = 4728090
#
#     mrca_1 = filteredScrape.mrca_ncbi_list
#     filteredScrape.data.ott_mrca = int(ingroup_mrca)
#     mrca_id = filteredScrape.ids.ott_to_ncbi[filteredScrape.data.ott_mrca]
#     print(mrca_id)
#     filteredScrape.mrca_ncbi_list = set()
#     filteredScrape.mrca_ncbi_list.add(mrca_id)
#     mrca_2 = filteredScrape.mrca_ncbi_list
#     assert mrca_1 != mrca_2
#
#     if not os.path.exists("{}/tmp".format(filteredScrape.workdir)):
#         os.mkdir("{}/tmp".format(filteredScrape.workdir))
#     src = "/home/blubb/Documents/gitdata/physcraper/tests/data/precooked/full_seqs"
#     src_files = os.listdir(src)
#     for file_name in src_files:
#         # print(src_files)
#         full_file_name = os.path.join(src, file_name)
#         if (os.path.isfile(full_file_name)):
#             # shutil.copy(full_file_name, dest)
#             shutil.copy(full_file_name, "{}/tmp".format(filteredScrape.workdir))
#
#     blast_dir = "tests/data/precooked/fixed/tte_blast_files"
#
#     filteredScrape.read_blast_wrapper(blast_dir=blast_dir)
#     filteredScrape.remove_identical_seqs()

