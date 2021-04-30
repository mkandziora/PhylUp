import os
from distutils.dir_util import copy_tree
from dendropy import Tree, DnaCharacterMatrix

from PhylUp import phyl_up, config, phylogen_updater

def test_reconcile():
    workdir = "tests/output/test_runs_reconcile"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast_test.config"

    seqalnmiss = "data/tiny_test_example/test_missingseq.fas"
    treefilemiss = "data/tiny_test_example/test_missingtip.tre"

    if not os.path.exists(workdir):
        os.mkdir(workdir)
    tmp_folder = os.path.join(workdir, 'tmp')
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    #call(['cp', '-a', 'data/tmp_for_test/', tmp_folder])
    copy_tree('data/tmp_for_test/', tmp_folder)


    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    print('###################### RUN 1 ######################')
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)

    cleaner = phylogen_updater.InputCleaner(trfn, schema_trf, seqaln, mattype, test.table, test.config) # test.mrca
    len_tre = len(cleaner.tre.taxon_namespace)
    len_aln = (len(cleaner.aln.taxon_namespace))
    tre_asstring = cleaner.tre.as_string(schema='newick')

    # ################################
    print('###################### RUN 2 ######################')

    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqalnmiss, mattype, trfn, schema_trf, conf)
    cleaner = phylogen_updater.InputCleaner(trfn, schema_trf, seqalnmiss, mattype, test.table, test.config) # test.mrca

    len_aln_missaln = (len(cleaner.aln))

    print('###################### RUN 3 ######################')


    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, treefilemiss, schema_trf, conf)
    test.config.update_tree == True
    cleaner = phylogen_updater.InputCleaner(treefilemiss, schema_trf, seqaln, mattype, test.table, test.config) # test.mrca

    len_aln_misstre = (len(cleaner.aln))
    tre_asstring_misstre = cleaner.tre.as_string(schema='newick')



    print(len_aln, len_aln_missaln, len_aln_misstre)
    assert len_aln != len_aln_missaln
    assert len_aln != len_aln_misstre
    assert len_aln == len_tre

    assert tre_asstring != tre_asstring_misstre


    print('###################### RUN 4 ######################')

    seqaln = "data/tiny_test_example/test_emptyseq.fas"
    trfn = "data/tiny_test_example/test.tre"

    aln = DnaCharacterMatrix.get(path=os.path.abspath(seqaln),
                                 schema='fasta')
    tre = Tree.get(path=trfn, schema="newick", preserve_underscores=True)
    len_aln_before = (len(aln))
    print(len(tre.taxon_namespace))
    len_tre_before = (len(tre.taxon_namespace))

    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)
    test.config.update_tree == True
    cleaner = phylogen_updater.InputCleaner(treefilemiss, schema_trf, seqaln, mattype, test.table,
                                            test.config)  # test.mrca

    len_aln_after = (len(cleaner.aln))
    len_tre_after = (len(cleaner.tre.taxon_namespace))

    print(len_aln, len_aln_missaln, len_aln_misstre)
    assert len_aln_after < len_aln_before
    assert len_tre_after < len_tre_before


    print('###################### RUN 5 ######################')

    seqaln = "data/tiny_test_example/test_emptyseq.fas"
    trfn = "data/tiny_test_example/test.tre"

    aln = DnaCharacterMatrix.get(path=os.path.abspath(seqaln),
                                 schema='fasta')
    tre = Tree.get(path=trfn, schema="newick", preserve_underscores=True)

    conf.backbone = True
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)
    test.config.update_tree == True
    cleaner = phylogen_updater.InputCleaner(treefilemiss, schema_trf, seqaln, mattype, test.table,
                                            test.config)  # test.mrca

    assert os.path.exists(os.path.join(workdir, 'backbone.tre'))


