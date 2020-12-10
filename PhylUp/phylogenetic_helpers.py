"""
PhylUp: phylogenetic alignment building with custom taxon sampling
Copyright (C) 2020  Martha Kandziora
martha.kandziora@mailbox.org

Package to automatically update alignments and phylogenies using local sequences or a local Genbank database
while controlling for the number of sequences per OTU.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.


"""


import os
import subprocess
import sys
import pandas as pd
import numpy as np
from copy import deepcopy
from dendropy import DnaCharacterMatrix, Tree

import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser

from . import cd
from . import debug
from . import suppress_stdout


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


def make_mafft_aln(aln, workdir):
    """ Write fasta for mafft alignment"""
    write_papara_alnfile(aln, workdir)

    aln = DnaCharacterMatrix.get(path=os.path.join(workdir, "aln_papara.phy"), schema='phylip')
    aln.write(path=os.path.join(workdir, "aln_papara.fasta"), schema='fasta')

    with open(os.path.join(workdir, "aln_papara.fasta")) as aln_fasta:
        existing = aln_fasta.read()
        with open(os.path.join(workdir, "new_seqs.fasta")) as new:
            new_seqs = new.read()
            with open(os.path.join(workdir, "mafft.fasta"), "wt") as fout:
                fout.write(existing)
                fout.write(new_seqs)


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


def replace_uid_with_name(file_path, table, matrix_type):
    """
    translate unique ids into species names in file.

    :param file_path:
    :param table:
    :param matrix_type: aln or tree, defines which type is being relabeled.
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
                        spn = present.loc[idx, 'ncbi_txn'].replace(" ", "_").replace("-", "_").replace("'", "")
                        if matrix_type == 'tree':
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
        labelled = labelled.replace(">{}\n".format(split_name),
                                    ">{}_{}\n".format(spn, split_name))
    except AttributeError:
        labelled = labelled.replace(">{}\n".format(split_name),
                                    ">{}_{}\n".format(present.loc[idx, 'ncbi_txid'],
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
    # while seqlen == 0:  #useless, seqlen never == 0
    #     i = i + 1
    #     seqlen = len(aln[i])
    for tax in aln:
        if len(aln[tax]) != seqlen:
            sys.stderr.write("Alignment is not aligned: {} vs {} .\n".format(len(aln[tax]), seqlen))
            return
    return seqlen

def run_papara():
    """
    Runs papara and adds new sequences to the alignment.

    :return:
    """
    with suppress_stdout():
        subprocess.check_call(["papara_static_x86_64", "-t", "papara_tre.tre", "-s", "aln_papara.phy",
                               #  "-j", "{}".format(self.config.num_threads),  # FIXME: only works when papara is compiled.
                               "-q", "new_seqs.fasta", "-n", 'phylip'], shell=False)

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
    sys.stdout.write('Estimate number of threads raxml')
    with cd(workdir):
        subprocess.run(['raxml-ng', '--parse', '--msa', aln_fn,
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
    ncbi_parser = ncbi_data_parser.Parser(names_file=config.ncbi_parser_names_fn,
                                          nodes_file=config.ncbi_parser_nodes_fn,
                                          interactive=False)
    sys.stdout.write("Build table with information about sequences and taxa.\n")
    if os.path.exists(os.path.join(config.workdir, 'spn_input_ncbiid.txt')):
        columns = ['accession', 'ncbi_txn', 'ncbi_txid']
        table = pd.read_csv(os.path.join(config.workdir, 'spn_input_ncbiid.txt'), names=columns,  sep=' ', header=0 )

    else:
        debug(id_to_spn)
        table = get_txid_for_name_from_file(id_to_spn, ncbi_parser)  # builds name id link
        table.to_csv(os.path.join(config.workdir, 'spn_input_ncbiid.txt'), sep=' ', header=True, index=False)
    table = table.drop_duplicates(subset='accession', keep='first')  # drop duplicated entries from file
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
    columns = ['accession', 'ncbi_txn']
    name_id = pd.read_csv(tipname_id_fn, names=columns,  sep=",", header=None)
    name_id['ncbi_txid'] = name_id['ncbi_txn'].apply(ncbi_parser.get_id_from_name)
    assert name_id['ncbi_txid'].isnull().values.any() == False, ('not all ncbi taxon ids assigned.',
                                                                 name_id['ncbi_txid'].isnull().values.any())
    return name_id


def add_seq_to_table(aln, table):
    """
    Puts input sequences from alignment into the pandas table.

    :return:build_table_from_file
    """
    queried_taxa = []
    for taxon, seq in aln.items():
        seq = seq.symbols_as_string().replace('?', '')
        seq = seq.replace('-', '')
        contains_string = table['accession'].str.contains(taxon.label)
        if contains_string.any():
            table.at[table['accession'] == taxon.label, 'sseq'] = seq
            queried_taxa.append(taxon.label)
        assert taxon.label in queried_taxa, (taxon.label, queried_taxa, 'often missing in spn translation table')

    table['sseq'].replace('', np.nan, inplace=True)
    table.dropna(subset=['sseq'], inplace=True)
    return table


def make_preferred_taxon_list(other_runs, fn, overlap_complete=True):
    """

    :param other_runs: dictionary with locus name, folderpath to runs that are part of preferred assesment.
    :param fn: path to the overlap file
    :param overlap_complete: either all must overlap, or intersections of diff combined
    :return:
    """
    preferred_taxa_all = dict()
    for locus in other_runs.keys():
        preferred_taxa = pd.read_csv(other_runs[locus], names=['taxon_name', 'ncbi_id'])
        unique_taxa = set(preferred_taxa.ncbi_id.to_list())
        preferred_taxa_all[locus] = unique_taxa
    key_list = list(preferred_taxa_all.keys())
    for i in range(0, len(key_list)-1):
        set_1 = preferred_taxa_all[key_list[i]]
        set2 = preferred_taxa_all[key_list[i+1]]
        overlap2 = set_1.intersection(set2)
        if overlap_complete == True:
            print('complete')
            if i == 0:
                overlap_all = overlap2
            else:
                set_1 = overlap_all
                set2 = overlap2
                overlap_all = set_1.intersection(set2)
        else:
            overlap_all = overlap_all + overlap2

        with open(fn, mode='wt') as f:
            for item in overlap_all:
                f.write("%s\n" % item)