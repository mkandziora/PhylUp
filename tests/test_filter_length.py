import os
from distutils.dir_util import copy_tree

import pandas as pd
from PhylUp import phyl_up, config, phylogenetic_helpers


def test_filterlen():
    workdir = "tests/output/test_runs_length"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test_extralongseq.fas"
    mattype = "fasta"
    configfi = "data/localblast_test.config"

    if not os.path.exists(workdir):
        os.mkdir(workdir)
    tmp_folder = os.path.join(workdir, 'tmp')
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    #call(['cp', '-a', 'data/tmp_for_test/', tmp_folder])
    copy_tree('data/tmp_for_test/', tmp_folder)


    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)

    aln = phylogenetic_helpers.read_in_aln(test.aln_fn, test.aln_schema)

    extra_long = {'ncbi_txn': ['extra_long'],
                  'ncbi_txid': [123],
                  'org_sp_name': ['extra_long'],
                  'tip_name': ['extra_long'],
                  'status': ['test'], 'status_note': ['note'], "date": ['2001'],
                  'index': [5], 'accession': ['532342.1'],
                  'pident': [99.533], 'evalue': [0], 'bitscore': [1175],
                  'sseq': [
                      'TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTTCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCTTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGC'],
                  'title': ['no title']
                  }

    new_seqs = pd.DataFrame(extra_long)
    before = len(new_seqs)

    f = phyl_up.FilterLength(test.config, aln)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    del_tab = len(f.del_table)
    after = len(new_seqs)

    assert before > after
    assert del_tab > 0

    extra_short = {'ncbi_txn': ['extra_long'],
                  'ncbi_txid': [123],
                  'org_sp_name': ['extra_long'],
                  'tip_name': ['extra_long'],
                  'status': ['test'], 'status_note': ['note'], "date": ['2001'],
                  'index': [5], 'accession': ['532342.1'],
                  'pident': [99.533], 'evalue': [0], 'bitscore': [1175],
                  'sseq': [
                      'TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGA'],
                  'title': ['no title']
                  }

    new_seqs = pd.DataFrame(extra_short)
    before = len(new_seqs)

    f = phyl_up.FilterLength(test.config, aln)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    del_tab = len(f.del_table)
    after = len(new_seqs)

    assert before > after
    assert del_tab > 0
    assert before == after + del_tab
