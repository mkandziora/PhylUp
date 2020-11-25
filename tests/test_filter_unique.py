import os
import pandas as pd
from distutils.dir_util import copy_tree

from PhylUp import phyl_up, config

def test_filter_unique():
    print('test_filter_unique')
    workdir = "tests/output/test_runs_unique"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
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
    conf.threshold = 2
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)

    new_seqs = test.extend()
    before = len(new_seqs)

    f = phyl_up.FilterUniqueAcc(test.config, test.table)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    del_tab = len(f.del_table)
    after = len(new_seqs)
    print(before, after, del_tab)
    assert before > after
    assert del_tab > 0
    assert after + del_tab == before

def test_filters_to_longest():
    print('test_filter_unique')
    workdir = "tests/output/test_runs_unique"
    trfn = "data/tiny_test_example/test.tre"
    schema_trf = "newick"
    id_to_spn = "data/tiny_test_example/test_nicespl.csv"
    seqaln = "data/tiny_test_example/test.fas"
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
    conf.threshold = 2
    conf.blast_folder = os.path.abspath("./data/blast_for_tests")
    test = phyl_up.PhylogeneticUpdater(id_to_spn, seqaln, mattype, trfn, schema_trf, conf)

    new_seqs = test.extend()

    # make own longest selection
    accs = new_seqs.drop_duplicates(subset=['accession'], keep='first')

    len_supposed = len(accs)

    select = pd.DataFrame(columns=['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date", 'accession',
                                   'pident', 'evalue', 'bitscore', 'sseq', 'title'])
    for acc in accs.accession:
        acc_df = new_seqs[new_seqs.accession == acc]
        lenseq_acc_df = acc_df['sseq'].apply(len)
        idxmax_val = lenseq_acc_df.idxmax()
        select = select.append(acc_df.loc[idxmax_val])

    assert len(select) == len_supposed, (len(select), len_supposed)

    # run filter from phylup
    before = len(new_seqs)

    f = phyl_up.FilterUniqueAcc(test.config, test.table)
    f.filter(new_seqs)
    new_seqs = f.upd_new_seqs

    assert new_seqs.accession.isin(select.accession).all()

