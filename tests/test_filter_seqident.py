from PhylUp import phyl_up, config, phylogen_updater
import pandas as pd
from copy import deepcopy
#
# print('test_filter_seqident')
# workdir = "tests/output/test_runs"
# trfn = "data/tiny_test_example/test.tre"
# schema_trf = "newick"
# id_to_spn = "data/tiny_test_example/test_nicespl.csv"
# seqaln = "data/tiny_test_example/test.fas"
# mattype = "fasta"
# configfi = "data/localblast.config"
#
# conf = config.ConfigObj(configfi, workdir, interactive=False)
# conf.threshold = 2
# test = pandas_numpy_try1.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)
#
# new_seqs = test.extend()
#
# new_existing_seq = """TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAG
# TATCGGGCTCTTGTTCGATTTCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCC
# TTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCC
# AAGGAAATATAAACTTAAGAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGA
# TTGCATTGAAACTTGCTTCTTTATAATTCATAAACGACTCTCGGCAACGGATATCTCG
# GCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCC
# GTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGT
# CTGCCTGGGCGTCACATATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTT
# TGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGCCTGAATTTTGA
# GTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCT
# TATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTG
# TTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAG"""
#
# new_seqs = new_seqs[~new_seqs['accession'].isin(test.table['accession'])]  # ~ is the pd not in/!
#
# f = pandas_numpy_try1.FilterUniqueAcc(test.config, test.table)
# f.filter(new_seqs)
# new_seqs = f.upd_new_seqs
# new_seqs.loc[1, 'sseq'] = new_existing_seq
# assert new_seqs.loc[1, 'sseq'] == new_existing_seq
# new_seqs_before = new_seqs
#
# len_table_before = len(test.table)
# new_seqs = test.add_new_seqs(new_seqs)
# assert new_seqs.index.tolist() != new_seqs_before.index.tolist()
# # assert new_existing_seq in [new_seqs['sseq']]
# assert new_seqs['sseq'].str.contains(new_existing_seq).any()
# print(test.table.index)
# len_table_after = len(test.table)
# assert len_table_before < len_table_after, (len_table_before, len_table_after)
# before = len(new_seqs)
#
# assert test.table['sseq'].hasnans == False, test.table['sseq']
#
# print(new_seqs.index)
# f = pandas_numpy_try1.FilterSeqIdent(test.config, test.table, test.status)
# f.filter(new_seqs)
# new_seqs = f.upd_new_seqs
#
# print(f.del_table)
# del_tab = len(f.del_table)
# after = len(new_seqs)
# print(before, after, del_tab)
# assert before >= after
# assert del_tab >= 0
# assert after + del_tab == before, (after, del_tab, before)
#
# seq_table = test.table['sseq']
# found = seq_table[seq_table.str.contains(new_existing_seq)]
# print(found)
# print(found.count())
# assert found.count() == 1, found.count()


# df[~df["col"].str.contains(word)]
def test_contain():

    present = {'ncbi_txn':['1', '2'],
        'sseq':['TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTTCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCTTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGC',
                'GCAGAACGACCTGTGAACATGTAACAACAATCGGGTGTTCTAAGTATCGGGCTCTTGTCCGATTCCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGACCCTTAGGTTTTCTAAGGACGTCGCGTCGACACAACAACAAACCCCCGGCACGGAATGTGCCAAGGAAATATAAACTTAAGAAGGGCTTGTTCCATGCTTTGCCGTTTTCGCGGTGATTGYGTTGAATCTTGCTTCTTTATAAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGCCGAGGGCACGTCTGCCTGGGCGTCACACATCGCGWCGCCCCCATCACACCCCTTGACGGGGATGTATGAATGGGGAGCGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGATGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTCGTGTGTTCCTAGGAGTAAGGAAGATCTCTT']}
    df = pd.DataFrame(present)
    longer = 'TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTTCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCTTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCCGCGCGCGC'
    shorter = 'TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTTCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCTTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGAC'
    same = 'TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTTCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCTTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAGTATTGGGCTCTTGTYCGATTYCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAATAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGATTGCATTGAAACTTGCTTCTTTATAATCATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGC'
    assert df['sseq'].str.contains(longer).any() == False
    assert df['sseq'].str.contains(shorter).any() == True
    assert df['sseq'].str.contains(same).any() == True
    txid_compare = 1
    same_old = df[(df.sseq.str.contains(same)) & (df['ncbi_txn'] == txid_compare)]
    print(len(same_old.index))
    print(same)
    found = 0
    if len(same_old.index) > 0:
        print("identical old")
        found =1
    assert found == 0


def test_filter_compare():
    print('test_filter_seqident')
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

    new_seqs = test.extend()

    new_existing_seq = """TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAACAATTGGGTGTTCTAAGTATCGGGCTCTTGTCCGATTCCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAAGAAGGGCTTGTTCCATGCATTGCCGTTCGTGGTGACTGCATTGAAACTTGCTTCTCTATAATTAATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAACCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGAAGATTGGTCTCCTGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGC"""

    new_seqs = new_seqs[~new_seqs['accession'].isin(test.table['accession'])]  # ~ is the pd not in/!

    f = phyl_up.FilterUniqueAcc(test.config, test.table)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs
    new_seqs.loc[1, 'sseq'] = new_existing_seq
    assert new_seqs.loc[1, 'sseq'] == new_existing_seq
    assert new_seqs['sseq'].str.contains(new_existing_seq).any()
    assert test.table['sseq'].str.contains(new_existing_seq).any()

    new_seqs_before = new_seqs

    len_table_before = len(test.table)
    new_seqs = test.add_new_seqs(new_seqs)
    assert new_seqs.index.tolist() != new_seqs_before.index.tolist()
    # assert new_existing_seq in [new_seqs['sseq']]
    len_table_after = len(test.table)
    assert len_table_before < len_table_after, (len_table_before, len_table_after)

    assert test.table['sseq'].hasnans == False, test.table['sseq']

    new_seqs = test.compare_filter(new_seqs)

    assert new_seqs['sseq'].str.contains(new_existing_seq).any() == False


    all_avail_data = test.table[test.table['status'].between(0, test.status, inclusive=True)]

    seq_table = all_avail_data['sseq']
    found = seq_table[seq_table.str.contains(new_existing_seq)]
    print(found)
    print(found.count())
    assert found.count() == 1, found.count()



def test_filter_compare_shorter():
    print('test_filter_seqident')
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

    new_seqs = test.extend()

    new_existing_seq = """AGAACGACCTGTGAACATGTAAAACAATTGGGTGTTCTAAGTATCGGGCTCTTGTCCGATTCCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAAGAAGGGCTTGTTCCATGCATTGCCGTTCGTGGTGACTGCATTGAAACTTGCTTCTCTATAATTAATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAACCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGAAGATTGGTCTCCTGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGC"""

    new_seqs = new_seqs[~new_seqs['accession'].isin(test.table['accession'])]  # ~ is the pd not in/!

    f = phyl_up.FilterUniqueAcc(test.config, test.table)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs
    new_seqs.loc[1, 'sseq'] = new_existing_seq
    assert new_seqs.loc[1, 'sseq'] == new_existing_seq
    assert new_seqs['sseq'].str.contains(new_existing_seq).any()
    assert test.table['sseq'].str.contains(new_existing_seq).any()

    new_seqs_before = new_seqs

    len_table_before = len(test.table)
    new_seqs = test.add_new_seqs(new_seqs)
    assert new_seqs.index.tolist() != new_seqs_before.index.tolist()
    # assert new_existing_seq in [new_seqs['sseq']]
    len_table_after = len(test.table)
    assert len_table_before < len_table_after, (len_table_before, len_table_after)

    assert test.table['sseq'].hasnans == False, test.table['sseq']

    new_seqs = test.compare_filter(new_seqs)

    assert new_seqs['sseq'].str.contains(new_existing_seq).any() == False


    all_avail_data = test.table[test.table['status'].between(0, test.status, inclusive=True)]

    seq_table = all_avail_data['sseq']
    found = seq_table[seq_table.str.contains(new_existing_seq)]
    print(found)
    print(found.count())
    assert found.count() == 1, found.count()




def test_filter_seqident_newexist():
    print('test_filter_seqident')
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

    new_seqs = test.extend()

    new_existing_seq = """TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAG
    TATCGGGCTCTTGTTCGATTTCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCC
    TTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCC
    AAGGAAATATAAACTTAAGAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGA
    TTGCATTGAAACTTGCTTCTTTATAATTCATAAACGACTCTCGGCAACGGATATCTCG
    GCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCC
    GTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGT
    CTGCCTGGGCGTCACATATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTT
    TGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGCCTGAATTTTGA
    GTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCT
    TATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTG
    TTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAG"""

    new_seqs = new_seqs[~new_seqs['accession'].isin(test.table['accession'])]  # ~ is the pd not in/!

    f = phyl_up.FilterUniqueAcc(test.config, test.table)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs
    new_seqs.loc[1, 'sseq'] = new_existing_seq
    assert new_seqs.loc[1, 'sseq'] == new_existing_seq
    new_seqs_before = new_seqs

    len_table_before = len(test.table)
    new_seqs = test.add_new_seqs(new_seqs)
    assert new_seqs.index.tolist() != new_seqs_before.index.tolist()
    # assert new_existing_seq in [new_seqs['sseq']]
    assert new_seqs['sseq'].str.contains(new_existing_seq).any()
    print(test.table.index)
    len_table_after = len(test.table)
    assert len_table_before < len_table_after, (len_table_before, len_table_after)
    before = len(new_seqs)

    assert test.table['sseq'].hasnans == False, test.table['sseq']

    print(new_seqs.index)
    f = phyl_up.FilterSeqIdent(test.config, test.table, test.status)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    print(f.del_table)
    del_tab = len(f.del_table)
    after = len(new_seqs)
    print(before, after, del_tab)
    assert before >= after
    assert del_tab >= 0
    assert after + del_tab == before, (after, del_tab, before)

    all_avail_data = test.table[test.table['status'].between(0, test.status, inclusive=True)]

    seq_table = all_avail_data['sseq']
    found = seq_table[seq_table.str.contains(new_existing_seq)]
    print(found)
    print(found.count())
    assert found.count() == 1, found.count()



def test_oldseq():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

    new_seqs = test.extend()

    old_seq = """TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAAACAATTGGGTGTTCTAAG
       TATCGGGCTCTTGTTCGATTTCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCC
       TTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCC
       AAGGAAATATAAACTTAAGAAGGGCTTGTTCCATGCATTGCCGTTCGCGGTGA
       TTGCATTGAAACTTGCTTCTTTATAATTCATAAACGACTCTCGGCAACGGATATCTCG
       GCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCC
       GTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGT
       CTGCCTGGGCGTCACATATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTT
       TGAATGGGGACGGAGATTGGTCTCCCGTTCCTAAGGTGCGGTTGCCTGAATTTTGA
       GTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCT
       TATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTG
       TTGTCTCATGACGATGCTTCGACTGC"""


    print('############################ old seq #########################')

    new_seqs = new_seqs[~new_seqs['accession'].isin(test.table['accession'])]  # ~ is the pd not in/!

    f = phyl_up.FilterUniqueAcc(test.config, test.table)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs
    new_seqs.loc[1, 'sseq'] = old_seq
    assert new_seqs.loc[1, 'sseq'] == old_seq
    new_seqs_before = new_seqs

    len_table_before = len(test.table)
    new_seqs = test.add_new_seqs(new_seqs)
    assert new_seqs.index.tolist() != new_seqs_before.index.tolist()
    # assert old_seq in [new_seqs['sseq']]
    assert new_seqs['sseq'].str.contains(old_seq).any()
    print(test.table.index)
    len_table_after = len(test.table)
    assert len_table_before < len_table_after, (len_table_before, len_table_after)
    before = len(new_seqs)

    assert test.table['sseq'].hasnans == False, test.table['sseq']

    print(new_seqs.index)
    f = phyl_up.FilterSeqIdent(test.config, test.table, test.status)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    print(f.del_table)
    del_tab = len(f.del_table)
    after = len(new_seqs)
    print(before, after, del_tab)
    assert before >= after
    assert del_tab >= 0
    assert after + del_tab == before, (after, del_tab, before)

    all_avail_data = test.table[test.table['status'].between(0, test.status, inclusive=True)]

    seq_table = all_avail_data['sseq']
    found = seq_table[seq_table.str.contains(old_seq)]
    print(found)
    print(found.count())
    assert found.count() == 1, found.count()


def test_no_similar():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"

    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

    new_seqs = test.extend()

    no_similar = """TCGAAACCTGCATAGCAGAACGACTTTTTTCTGTGAACATGTAAAAACAATTGGGTGTTCTAAG
       TATCGGGCTCTTGTTCGATTTCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCC
       TTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCC
       AAGGAAATATAAACTTAAGAAGGGCTTGTAAAAAATCCATGCATTGCCGTTCGCGGTGA
       TTGCATTGAAACTTGCTTCTTTATAATTCATAAACGACTCTCGGCAACGGATATCTCG
       GCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCC
       GTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTTGAGGGCACGT
       CTGCCTGGGCGTCACATATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTT
       TGAATGGGGACGGAGATTGGTCTCCCGTCCCCCCCTCCTAAGGTGCGGTTGCCTGAATTTTGA
       GTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCT
       TATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTG
       TTGTCTCATGACGATGCTTCGACTGCTCGAAACCTGCATAGCAGAACGACCTGTGAACATGT"""
    print('############################ no similar #########################')

    new_seqs = new_seqs[~new_seqs['accession'].isin(test.table['accession'])]  # ~ is the pd not in/!

    f = phyl_up.FilterUniqueAcc(test.config, test.table)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs
    new_seqs.loc[1, 'sseq'] = no_similar
    assert new_seqs.loc[1, 'sseq'] == no_similar
    new_seqs_before = new_seqs

    len_table_before = len(test.table)
    new_seqs = test.add_new_seqs(new_seqs)
    assert new_seqs.index.tolist() != new_seqs_before.index.tolist()
    # assert no_similar in [new_seqs['sseq']]
    assert new_seqs['sseq'].str.contains(no_similar).any()
    print(test.table.index)
    len_table_after = len(test.table)
    assert len_table_before < len_table_after, (len_table_before, len_table_after)
    before = len(new_seqs)

    assert test.table['sseq'].hasnans == False, test.table['sseq']

    print(new_seqs.index)
    f = phyl_up.FilterSeqIdent(test.config, test.table, test.status)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    print(f.del_table)
    del_tab = len(f.del_table)
    after = len(new_seqs)
    print(before, after, del_tab)
    assert before >= after
    assert del_tab >= 0
    assert after + del_tab == before, (after, del_tab, before)

    all_avail_data = test.table[test.table['status'].between(0, test.status, inclusive=True)]

    seq_table = all_avail_data['sseq']
    found = seq_table[seq_table.str.contains(no_similar)]
    print(found)
    print(found.count())
    assert found.count() == 1, found.count()

def test_not_add_identical():
    workdir = "tests/output/test_runs"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
    mattype = "fasta"
    configfi = "data/localblast.config"


    conf = config.ConfigObj(configfi, workdir, interactive=False)
    conf.threshold = 2
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf, mrca=18794)

    new_seqs = test.extend()

    new_ident_key = 'KM592514.1'
    new_ident = """TCGAAACCTGCATAGCAGAACGACCTGTGAACATGTAAAACAATTGGGTGTTCTAAGTATCGGGCTCTTGTCCGATTCCTAGGATGCCATGTTGACGTGCGTCTTTGGCAAGCCCCTTGGGTGTCTAAGGACGTCACGTCGACGCAACAACAAACCCCCGGCACGGCATGTGCCAAGGAAATATAAACTTAAGAAGGGCTTGTTCCATGCATTGCCGTTCGTGGTGACTGCATTGAAACTTGCTTCTCTATAATTAATAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAACCTTTTGGTTGAGGGCACGTCTGCCTGGGCGTCACACATCGCGTCGCCCCCATCACACCTCTTGACGGGGATGTTTGAATGGGGACGAAGATTGGTCTCCTGTTCCTAAGGTGCGGTTGGCTGAATTTTGAGTCCTCTTCGACGGACGCACGATTAGTGGTGGTTGACAAGACCTTCTTATCGAGTTGTGTGTTCCAAGAAGTAAGGAATATCTCTTTAACGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGC"""

    assert test.table['sseq'].str.contains(new_ident).any()

    print('############################ new ident #########################')

    new_seqs = new_seqs[~new_seqs['accession'].isin(test.table['accession'])]  # ~ is the pd not in/!

    f = phyl_up.FilterUniqueAcc(test.config, test.table)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs
    new_seqs.loc[1, 'sseq'] = new_ident
    new_seqs.loc[1, 'accession'] = new_ident_key

    assert new_seqs.loc[1, 'sseq'] == new_ident

    assert new_seqs['sseq'].str.contains(new_ident).any()

    assert new_seqs[new_seqs['sseq'].str.contains(new_ident)].count() > 0
    assert new_seqs[new_seqs['accession'].str.contains(new_ident_key)].count() > 0

    f = phyl_up.FilterSeqIdent(test.config, test.table, test.status)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    assert new_seqs['sseq'].str.contains(new_ident) == 0
    assert new_seqs['accession'].str.contains(new_ident_key) == 0

    all_avail_data = test.table[test.table['status'].between(0, test.status, inclusive=True)]

    seq_table = all_avail_data['sseq']
    found = seq_table[seq_table.str.contains(new_ident)]
    print(found)
    print(found.count())
    assert found.count() == 1, found.count()


