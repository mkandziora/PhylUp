
from pandas_filter import pandas_numpy_try1, config, aln_updater

def test_reconcile():
    workdir = "tests/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    seqalnmiss = "data/tiny_test_example/test_missingseq.fas"
    treefilemiss = "data/tiny_test_example/test_missingtip.tre"


    conf = config.ConfigObj(configfi, workdir, interactive=False)
    print('###################### RUN 1 ######################')
    test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)


    aln = test.read_in_aln()
    tre = test.read_in_tree(aln.taxon_namespace)

    cleaner = aln_updater.InputCleaner(tre, aln, seqaln, test.table, test.config, test.mrca)

    len_aln = (len(cleaner.aln.taxon_namespace))
    # ################################
    print('###################### RUN 2 ######################')

    test = pandas_numpy_try1.Update_data(id_to_spn, seqalnmiss, mattype, trfn, schema_trf, conf, mrca=18794)

    aln = test.read_in_aln()
    tre = test.read_in_tree(aln.taxon_namespace)

    cleaner = aln_updater.InputCleaner(tre, aln, seqalnmiss, test.table, test.config, test.mrca)

    len_aln_missaln = (len(cleaner.aln))
    #len_tre_missaln = (len(cleaner.tre))
    tre_asstring = cleaner.tre.as_string(schema='newick')
    print('###################### RUN 3 ######################')


    test = pandas_numpy_try1.Update_data(id_to_spn, seqaln, mattype, treefilemiss, schema_trf, conf, mrca=18794)

    aln = test.read_in_aln()
    tre = test.read_in_tree(aln.taxon_namespace)

    cleaner = aln_updater.InputCleaner(tre, aln, seqaln, test.table, test.config, test.mrca)

    len_aln_misstre = (len(cleaner.aln))
    #len_tre_misstre = (len(cleaner.tre))
    tre_asstring_misstre = cleaner.tre.as_string(schema='newick')



    print(len_aln, len_aln_missaln, len_aln_misstre)
    assert len_aln != len_aln_missaln
    assert len_aln == len_aln_misstre
    assert tre_asstring != tre_asstring_misstre
