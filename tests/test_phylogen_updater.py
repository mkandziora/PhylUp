import os
from distutils.dir_util import copy_tree
import pandas as pd
from PhylUp import phyl_up, config, phylogenetic_helpers, phylogen_updater
from copy import deepcopy

def test_remove_short_fromaln():
    workdir = "tests/output/test_runs_length"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test_extrashortseq.fas"
    mattype = "fasta"
    configfi = "data/localblast_test.config"

    if not os.path.exists(workdir):
        os.mkdir(workdir)
    tmp_folder = os.path.join(workdir, 'tmp')
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    # call(['cp', '-a', 'data/tmp_for_test/', tmp_folder])
    copy_tree('data/tmp_for_test/', tmp_folder)

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    conf.minlen = 0.9

    aln0 = phylogenetic_helpers.read_in_aln(seqaln, mattype).taxon_namespace
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)


    alnbefore = deepcopy(test.aln.taxon_namespace)

    aln_upd = phylogen_updater.AlnUpdater(conf, test.aln, test.table, 2, None)

    #aln_upd.delete_short_seqs()
    print(aln0)
    print(alnbefore)
    print(aln_upd.aln.taxon_namespace)
    assert alnbefore == aln0, (alnbefore, aln0)
    assert alnbefore != aln_upd.aln.taxon_namespace, (alnbefore,  aln_upd.aln.taxon_namespace)





