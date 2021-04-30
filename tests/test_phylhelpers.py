# package import
import os, shutil
from distutils.dir_util import copy_tree
import pandas as pd
from PhylUp import phyl_up, config, blast, phylogenetic_helpers, phylogen_updater
import sys
import numpy
import datetime
import pandas as pd
import numpy as np
import os
from dendropy import Tree, DnaCharacterMatrix


def test_truncate_papara_aln():
    workdir = "tests/output/test_phylhelper"  # working directory
    seqaln = "data/tiny_test_example/test.fas"  # alignment
    mattype = "fasta"  # format of alignment
    trfn = "data/tiny_test_example/test.tre"  # phylogeny
    schema_trf = "newick"  # format of phylogeny
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"  # tab-delimited file where tip names correspond to ncbi names
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
    # shutil.copyfile('data/tiny_test_example/updt_aln.fasta', os.path.join(workdir, 'updt_aln.fasta'))
    # shutil.copyfile('data/tiny_test_example/updt_tre.tre', os.path.join(workdir, 'updt_tre.tre'))
    tmp_folder = os.path.join(workdir, 'tmp')
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    copy_tree('data/tmp_for_test/', tmp_folder)


    ###################

    aln = DnaCharacterMatrix.get(path=seqaln,
                                 schema='fasta')
    lenaln = len(aln[0])

    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, ignore_acc_list=None)
    print('prepare input')
    test.call_input_cleaner()
    print('extend')
    new_seqs = test.extend()  # todo rename to find new seqs
    new_seqs = test.call_filter(new_seqs, test.aln)
    tre = Tree.get(path=os.path.abspath(os.path.join(workdir, 'updt_tre.tre')), schema='newick',
                   taxon_namespace=test.aln.taxon_namespace, preserve_underscores=True)

    aln = DnaCharacterMatrix.get(path=os.path.abspath(os.path.join(workdir, 'updt_aln.fasta')),
                                 schema='fasta')

    aln_upd = phylogen_updater.AlnUpdater(test.config, aln, test.table, test.status, tre)
    aln_upd.write_papara_queryseqs()


    phylogenetic_helpers.write_papara_alnfile(aln, workdir)
    phylogenetic_helpers.write_papara_trefile(tre, workdir)

    cwd = os.getcwd()

    os.chdir(workdir)
    phylogenetic_helpers.run_papara()
    #path = os.path.join(os.getcwd(), "papara_alignment.phylip".format())
    #assert os.path.exists(path), "{} does not exists".format(path)
    os.chdir(cwd)

    aln = aln_upd.trim(os.path.join(workdir, 'papara_alignment.phylip'), 'phylip')

    phylogenetic_helpers.truncate_papara_aln(aln)

    aln = DnaCharacterMatrix.get(path=os.path.abspath(os.path.join(workdir, 'updt_aln.fasta')),
                                 schema='fasta')
    lenaln_after = len(aln[0])


    assert lenaln > lenaln_after


def xtest_run_modeltest():
    workdir = "tests/output/test_phylhelper"  # working directory
    phylogenetic_helpers.run_modeltest('updt_aln.fasta', workdir, 'AIC')
    #phylogenetic_helpers.find_best_model(seqaln, 'AICc', None)

# def test_best_model():
#     workdir = "tests/output/test_phylhelper"  # working directory
#     seqaln = "data/tiny_test_example/test.fas"  # alignment
#     mattype = "fasta"  # format of alignment
#     trfn = "data/tiny_test_example/test.tre"  # phylogeny
#     schema_trf = "newick"  # format of phylogeny
#     id_to_spn = "data/tiny_test_example/test_nicespl.csv"  # tab-delimited file where tip names correspond to ncbi names
#     configfi = "data/localblast_test.config"
#
#     phylogenetic_helpers.find_best_model(seqaln, 'AICc', None)

def test_read_tree():
    trfn = "data/tiny_test_example/test.tre"  # phylogeny
    phylogenetic_helpers.read_in_tree(trfn)


def test_writemafft():
    workdir = "tests/output/test_phylhelper"  # working directory
    seqaln = "data/tiny_test_example/test.fas"  # alignment
    mattype = "fasta"  # format of alignment
    aln = DnaCharacterMatrix.get(path=os.path.abspath(seqaln),
                                 schema='fasta')
    phylogenetic_helpers.make_mafft_aln(aln, workdir)

    assert os.path.exists(os.path.join(workdir, "mafft.fasta"))


def test_rewrite():
    workdir = "tests/output/test_phylhelper"  # working directory
    seqaln = "data/tiny_test_example/test.fas"  # alignment
    mattype = "fasta"  # format of alignment
    configfi = "data/localblast_test.config"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    trfn = "data/tiny_test_example/test.tre"  # phylogeny
    schema_trf = "newick"  # format of phylogeny


    conf = config.ConfigObj(configfi, workdir, interactive=False)
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)

    shutil.copyfile(seqaln, os.path.join(workdir, 'updt_aln.fasta'))
    shutil.copyfile(trfn, os.path.join(workdir, 'updt_tre.tre'))

    phylogenetic_helpers.replace_uid_with_name( os.path.join(workdir, 'updt_aln.fasta'), test.table, 'aln')
    phylogenetic_helpers.replace_uid_with_name(os.path.join(workdir, 'updt_tre.tre'), test.table, 'tree')

    assert os.path.exists(os.path.join(workdir, 'updt_aln.fasta_relabel'))
    assert os.path.exists(os.path.join(workdir, 'updt_tre.tre_relabel'))
