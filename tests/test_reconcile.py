
from PhylUp import phyl_up, config, phylogen_updater

def test_reconcile():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    seqalnmiss = "data/tiny_test_example/test_missingseq.fas"
    treefilemiss = "data/tiny_test_example/test_missingtip.tre"


    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    print('###################### RUN 1 ######################')
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

    cleaner = phylogen_updater.InputCleaner(trfn, seqaln, test.table, test.config, test.mrca)

    len_tre = len(cleaner.tre.taxon_namespace)
    len_aln = (len(cleaner.aln.taxon_namespace))
    tre_asstring = cleaner.tre.as_string(schema='newick')

    # ################################
    print('###################### RUN 2 ######################')

    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqalnmiss, mattype, trfn, schema_trf, conf, mrca=18794)

    cleaner = phylogen_updater.InputCleaner(trfn, seqalnmiss, test.table, test.config, test.mrca)

    len_aln_missaln = (len(cleaner.aln))

    print('###################### RUN 3 ######################')


    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, treefilemiss, schema_trf, conf, mrca=18794)
    test.config.update_tree == True

    cleaner = phylogen_updater.InputCleaner(treefilemiss, seqaln, test.table, test.config, test.mrca)

    len_aln_misstre = (len(cleaner.aln))
    tre_asstring_misstre = cleaner.tre.as_string(schema='newick')



    print(len_aln, len_aln_missaln, len_aln_misstre)
    assert len_aln != len_aln_missaln
    assert len_aln != len_aln_misstre
    assert len_aln == len_tre

    assert tre_asstring != tre_asstring_misstre
