"""
PhylUp: automatically update alignments.
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

All rights reserved. No warranty, explicit or implicit, provided. No distribution or modification of code allowed.
All classes and methods will be distributed under an open license in the near future.

Package to automatically update alignments and phylogenies using local sequences or a  local Genbank database.

Parts are inspired by the program physcraper developed by me and Emily Jane McTavish.
"""


import os
import subprocess
import sys
import pandas as pd

from copy import deepcopy
from dendropy import DnaCharacterMatrix, Tree

import numpy as np
import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser

from . import cd
from . import debug

def truncate_papara_aln(aln):
    """
    Used to split papara alingment into new and old sequences. Needed by EPA-ng.

    EPA-ng has functionality but that one requires extended and old aln to have same amount of bases...

    :return:
    """
    # # split papara.extended into new and olds seqs (1 arg = old, 2 = all)
    # # split does not work, as original aln_papara often shorter than extended:
    # subprocess.call(["epa-ng", "--split", "aln_papara.phy", "papara_alignment.extended"])
    debug('truncate papara aln')
    len_aln = len(aln.taxon_namespace)
    # read the data file in as a list
    ref_msa_fn = open(os.path.abspath('papara_alignment.phylip_trim'), "r")
    ref_msa = ref_msa_fn.readlines()
    ref_msa_fn.close()
    query_msa = deepcopy(ref_msa)

    del query_msa[1:(len_aln+1)]
    fout = open("new_seqs_papara.phylip", "w")
    fout.writelines(query_msa)
    fout.close()

    del ref_msa[len_aln+1:]
    fout = open("old_seqs_papara.phylip", "w")
    fout.writelines(ref_msa)
    fout.close()

    tmpaln = DnaCharacterMatrix.get(path='old_seqs_papara.phylip', schema='phylip')
    tmpaln.write(path="old_seqs.fasta", schema='fasta')
    tmpaln = DnaCharacterMatrix.get(path='new_seqs_papara.phylip', schema='phylip')
    tmpaln.write(path="new_seqs.fasta", schema='fasta')


def write_papara_alnfile(aln, workdir):
    """This writes out aln files for papara (except query sequences).
    Papara needs phylip format for the alignment.
    """
    aln.write(path=os.path.join(workdir, "aln_papara.phy"), schema="phylip")


def write_papara_trefile(tre, workdir):
    """This writes out tree files for papara.
    """
    print('write papara files')
    tre.resolve_polytomies()
    tre.deroot()
    tmptre = tre.as_string(schema="newick", unquoted_underscores=True, suppress_rooting=True)
    tmptre = tmptre.replace(":0.0;", ";")
    tmptre = tmptre.replace("'", "_")
    fi = open(os.path.join(workdir, "papara_tre.tre"), "w")
    fi.write(tmptre)
    fi.close()


def write_aln(aln, workdir, alnpath="updt_aln.fasta", alnschema="fasta"):
    """
    Outputs dendropy alignment  as file.

    :param aln: dendropy aln object
    :param workdir: working directory
    :param alnpath: alignemnt file name
    :param alnschema: format for alignment in file
    :return:
    """
    aln.write(path=os.path.join(workdir, alnpath),
              schema=alnschema)


def write_tre(tre, workdir, treepath="updt_tre.tre", treeschema="newick"):
    """
    Outputs dendropy tree  as file.

    :param tre: dendropy tree object
    :param workdir: working directory
    :param treepath: tree file name
    :param treeschema: format for tree in file
    :return:
    """
    tre.write(path=os.path.join(workdir, treepath), schema=treeschema, unquoted_underscores=True, suppress_rooting=True)


def replace_uid_with_name(file_path, table, type):
    """
    translate unique ids into species names in file.

    :param file_path:
    :param table:
    :param type: aln or tree, defines which type is being relabeled.
    :return:
    """
    # print(os.getcwd())
    name_list = []
    with open(file_path, "r") as label_new:
        labelled = label_new.read()
        with open("{}_relabel".format(file_path), "wt") as fout:
            present = table[table['status'] >= 0]
            for idx in present.index:
                if 'concat_id' in present.columns:
                    id = present.loc[idx, 'concat_id']
                    split_name = 'concat_{}'.format(int(id))
                else:
                    split_name = present.loc[idx, 'accession'].split('.')[0]
                if split_name not in name_list:
                    if split_name in labelled:
                        debug(split_name)
                        spn = present.loc[idx, 'ncbi_txn'].replace(" ", "_").replace("-", "_").replace("'", "")
                        if type == 'tree':
                            labelled = replace_in_tree(idx, labelled, present, split_name, spn)
                        else:
                            labelled = replace_in_aln(idx, labelled, present, split_name, spn)
                name_list.append(split_name)

            fout.write(labelled)


def replace_in_tree(idx, labelled, present, split_name, spn):
    try:
        labelled = labelled.replace("{}:".format(split_name), "{}_{}:".format(spn, split_name))
    except AttributeError:
        labelled = labelled.replace("{}:".format(split_name),
                                    "{}_{}:".format(present.loc[idx, 'ncbi_txid'],
                                                    split_name))
    return labelled


def replace_in_aln(idx, labelled, present, split_name, spn):

    try:
        labelled = labelled.replace(">{}".format(split_name),
                                    ">{}_{}".format(spn, split_name))
    except AttributeError:
        labelled = labelled.replace(">{}".format(split_name),
                                    ">{}_{}".format(present.loc[idx, 'ncbi_txid'],
                                                    split_name))
    return labelled


def check_align(aln):
    """
    Checks if the alignment is aligned, and if yes it returns the length of the alignment.

    :param aln: alignment in dendropy format.
    :return:
    """
    i = 0
    seqlen = len(aln[i])
    while seqlen == 0:
        i = i + 1
        seqlen = len(aln[i])
    for tax in aln:
        if len(aln[tax]) != seqlen:
            sys.stderr.write("Alignment is not aligned.")
            return
    return seqlen


def resolve_polytomies(tre, workdir):
    """
    Randomly resolves polytomies in provided dendropy tree.

    :param tre:
    :return:
    """
    print('resolve_polytomies')
    tre.resolve_polytomies()
    tre.deroot()
    tre.as_string(schema='newick')
    tre_fn = os.path.join(workdir, "papara_tre.tre")
    with open(tre_fn, "w") as tre_file:
        tre_file.write("{}".format(tre.as_string(schema='newick', suppress_rooting=True)))


def run_modeltest(aln_fn, workdir, model, partition=None):
    """
    Run modeltest to get the best substitution model for the alignments.

    Note: if number of threads is supplied it often breaks, .e.g. '-p {}'.format(2)

    :param aln_fn:
    :param workdir:
    :param model:
    :param partition:
    :return:
    """
    aln_fn = os.path.abspath(os.path.join(workdir, aln_fn))

    if not os.path.exists('{}.out'.format(aln_fn)):
        if os.path.exists('{}.ckp'.format(aln_fn)):
            os.remove('{}.ckp'.format(aln_fn))
        if partition is None:
            subprocess.run(['modeltest-ng', '-i', aln_fn,  '--force'], shell=False)
        else:
            subprocess.run(['modeltest-ng', '-i', aln_fn, '--partition', partition, '--force'], shell=False)

    best_model = find_best_model(aln_fn, model, partition)
    return best_model


def find_best_model(model_outfn, model='AICc', partition=None):
    """
    Read output of modeltest-ng and get best model.

    :param model_outfn: Output file from modeltest-ng
    :param model: information criterium for model selection
    :param partition: partition file
    :return:
    """
    print('get best substitution model')
    assert model in ['AICc', 'AIC', 'BIC']
    # fn = '{}.part.{}'.format(model_outfn, model.lower())
    model_outfn = '{}.out'.format(model_outfn)
    with open(model_outfn) as f:
        datafile = f.readlines()
    found_model = False
    all_models = list()
    if partition is None:
        for line in datafile:
            if model in line:
                found_model = True
            if found_model is True:
                if 'raxml-ng --msa ' in line:
                    all_models = line.split(' ')[-1][:-1]
                    return all_models
    else:
        all_models = 'DNA'
    return all_models


def estimate_number_threads_raxml(workdir, aln_fn, model):
    """
    Get the optimal number of threads for the raxml-ng tree search.

    :param workdir: working directory
    :param aln_fn: filename of alignment
    :param model: substitution model
    :return: optimal number of threads
    """
    print('estimate_number_threads_raxml')
    with cd(workdir):
        subprocess.run(['raxml-ng-mpi', '--parse', '--msa', aln_fn,
                        '--prefix', 'numthreads', '--model', model, '--redo'])
        with open('numthreads.raxml.log') as f:
            datafile = f.readlines()
        for line in datafile:
            if 'Recommended number of threads / MPI processes:' in line:
                num_threads = line.split(': ')[1][0]
                return num_threads


def build_table_from_file(id_to_spn, config, downtorank=None):
    """
    Build self.table from file. This is the object that basically holds all needed information.

    :param id_to_spn: file that translates tip names to species names. Must be accepted by ncbi
    :param config: configuration object
    :param downtorank: name of rank if sequences shall be filtered to a higher rank
    :return: table with all sequence information available
    """
    debug("Build table with information about sequences and taxa.")
    ncbi_parser = ncbi_data_parser.Parser(names_file=config.ncbi_parser_names_fn,
                                          nodes_file=config.ncbi_parser_nodes_fn,
                                          interactive=False)
    debug(id_to_spn)
    table = get_txid_for_name_from_file(id_to_spn, ncbi_parser)  # builds name id link
    table['status'] = 0
    try:
        table['date'] = pd.Timestamp.strptime('01/01/00', '%d/%m/%y')
    except NotImplementedError:
        table['date'] = pd.to_datetime('01/01/00', format='%d/%m/%y')
    table['sseq'] = None
    table['status'] = table['status'].astype(int)
    table['ncbi_txid'] = table['ncbi_txid'].astype(int)
    #if downtorank is not None:
        #table['original_ncbi_txid'] = table['ncbi_txid'].astype(int)
        #table['ncbi_txid'] = np.vectorize(ncbi_parser.get_downtorank_id)(table['original_ncbi_txid'])
        #table['ncbi_txid'] = np.vectorize(ncbi_parser.get_downtorank_id)(table['original_ncbi_txid'], downtorank)
    #else:
    #    table['ncbi_txid'] = table['ncbi_txid'].astype(int)
    return table


def read_in_tree(tre_fn, tr_schema=None, taxon_namespace=None):
    """
    Reads in phylogeny from file.

    :param tre_fn: tree filename
    :param tr_schema: format of tree
    :param taxon_namespace: needs taxon_namespace from aln so that they match
    :return: tre als dendropy object
    """
    if tr_schema is None:
        tr_schema = 'newick'
    if taxon_namespace:
        tre = Tree.get(path=tre_fn, schema=tr_schema,
                       taxon_namespace=taxon_namespace, preserve_underscores=True)
        assert tre.taxon_namespace == taxon_namespace
    else:
        tre = Tree.get(path=tre_fn, schema=tr_schema, preserve_underscores=True)
    # print(tre.as_string('newick'))
    return tre


def read_in_aln(aln_fn, aln_schema):
    """
    Reads text alignment in as dendropy object.

    :return: alignment as dendroppy object
    """
    aln = DnaCharacterMatrix.get(path=aln_fn, schema=aln_schema)
    assert aln.taxon_namespace
    return aln


def get_txid_for_name_from_file(tipname_id_fn, ncbi_parser):
    """
    Get taxon id for tipname, speciesname files using ncbi translation.

    Only used for option add unpublished.

    :param tipname_id_fn: file with tipname and corresponding species name
    :param ncbi_parser: ncbiTAXONparser instance
    :return:
    """
    print(tipname_id_fn)
    columns = ['accession', 'ncbi_txn']
    name_id = pd.read_csv(tipname_id_fn, names=columns,  sep=",", header=None)
    name_id['ncbi_txid'] = name_id['ncbi_txn'].apply(ncbi_parser.get_id_from_name)
    name_id['ncbi_txid'].isnull().values.any() == False
    return name_id


def add_seq_to_table(aln, table):
    """
    Puts input sequences from alignment into the pandas table.

    :return:
    """
    queried_taxa = []
    for taxon, seq in aln.items():
        seq = seq.symbols_as_string().replace('?', '')
        seq = seq.replace('-', '')
        for index in table.index:
            tip_name = table.loc[index, 'accession']
            if taxon.label in tip_name:
                table.at[index, 'sseq'] = seq
                queried_taxa.append(tip_name)
    assert tip_name in queried_taxa, (tip_name, queried_taxa)
    return table
