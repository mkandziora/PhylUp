from pandas_filter import pandas_numpy_try1, config

workdir = "test_runs"
trfn = "data/tiny_test_example/test.tre"
schema_trf = "newick"
id_to_spn = "data/tiny_test_example/test_nicespl.csv"


seqaln = "data/tiny_test_example/test.fas"
mattype = "fasta"

configfi = "data/localblast.config"

print("test")
print("data/tiny_test_example/test_nicespl.csv")



conf = config.ConfigObj(configfi, workdir, interactive=False)
test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)
print()

# test.build_table(fn)
test.run()
tab = test.table
print(test.table[['tip_name', 'status']])

print(len(test.table[test.table['status'] != 'original']))
