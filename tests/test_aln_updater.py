# package import
import os, shutil
from distutils.dir_util import copy_tree
import pandas as pd
from PhylUp import phyl_up, config, blast, phylogen_updater
import sys
import numpy
import datetime
import pandas as pd
import numpy as np
import os
from dendropy import Tree, DnaCharacterMatrix

def test_alnupdater():
    workdir = "tests/output/test_alnupdt"  # working directory
    trfn = "data/tiny_test_example/test.tre"  # phylogeny
    schema_trf = "newick"  # format of phylogeny
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"  # tab-delimited file where tip names correspond to ncbi names
    seqaln = "data/tiny_test_example/test.fas"  # alignment
    mattype = "fasta"  # format of alignment
    configfi = "data/localblast_test.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    conf.identical_seqs = False

    if not os.path.exists(workdir):
        os.rename(workdir)
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    tmp_folder = os.path.join(workdir, 'tmp')
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    # call(['cp', '-a', 'data/tmp_for_test/', tmp_folder])
    copy_tree('data/tmp_for_test/', tmp_folder)
    # shutil.copyfile('data/tiny_test_example/updt_aln.fasta', os.path.join(workdir, 'updt_aln.fasta'))
    # shutil.copyfile('data/tiny_test_example/updt_tre.tre', os.path.join(workdir, 'updt_tre.tre'))
    shutil.copyfile('data/tiny_test_example/updt_aln.fasta', os.path.join(workdir, 'updt_aln.fasta'))
    shutil.copyfile('data/tiny_test_example/updt_tre.tre', os.path.join(workdir, 'updt_tre.tre'))

    ###################################
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, ignore_acc_list=None)
    # test.run()
    # tre = Tree.get(path=os.path.abspath(os.path.join(workdir, 'updt_tre.tre')), schema='newick',
    #               taxon_namespace=test.aln.taxon_namespace, preserve_underscores=True)


    aln = DnaCharacterMatrix.get(path=os.path.abspath(os.path.join(workdir, 'updt_aln.fasta')),
                                 schema='fasta')

    print('prepare input')
    test.call_input_cleaner()
    # assert len(table) > 1, (len(table), table)  # not the case if single seq is used as input

    new_seqs = test.extend()  # todo rename to find new seqs

    new_seqs = test.call_filter(new_seqs, test.aln)

    print(len(new_seqs.index))
    assert len(new_seqs.index) >= 171, len(new_seqs.index)
    aln_before = aln
    aln_upd = phylogen_updater.AlnUpdater(test.config, test.aln, test.status, test.table, test.status, test.tre)
    assert aln_before != aln_upd.aln

    aln_upd.write_labelled('updt_aln.fasta')