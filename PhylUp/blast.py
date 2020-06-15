"""
PhylUp: automatically update alignments.
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

All rights reserved. No warranty, explicit or implicit, provided. No distribution or modification of code allowed.
All classes and methods will be distributed under an open license in the near future.

Package to automatically update alignments and phylogenies using local sequences or a local Genbank database
while controlling for the number of sequences per OTU.

Parts are inspired by the program physcraper developed by me and Emily Jane McTavish.
"""

# sequence retrieval methods

import os
import sys
import datetime
import subprocess
from tempfile import NamedTemporaryFile
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser  # is the ncbi data parser class and associated functions
from . import cd
from . import suppress_stdout
from . import debug
# import line_profiler

"""
blastdbcmd outfmt options:
https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastdbcmd_application_opti/
"""

def get_acc_from_blast(query_string):
    """
    Get the accession number from a blast query using local blast.

    :param query_string: string from blast search where acc and gi are encoded in.
    :return: genbank accession
    """
    # Note: To get acc is more difficult now, as new seqs not always have gi number, then query changes
    query_split = query_string.split("|")
    if len(query_split) > 3:
        gb_acc = query_split[3]
    elif len(query_split) == 3:
        gb_acc = query_split[1]
    else:
        gb_acc = query_split[0]
    return gb_acc


def get_taxid_from_acc(gb_acc, blastdb, workdir):
    """
    Use the blastdb to get the taxon id from a queried gb acc.

    Sometimes there are more than a single id, as database has merged redundant seqs.

    :param gb_acc: Genbank accession number
    :param workdir: working directory
    :param blastdb: path to blast db
    :return: list of taxon ids associated with the GB id - there are multiple because of the merging of redundant data
    """
    # print('get taxid from acc')
    if not os.path.exists(os.path.join(workdir, 'tmp')):
        os.mkdir(os.path.join(workdir, 'tmp'))
    fn = os.path.join(workdir, 'tmp/tmp_search.csv')
    fn_open = open(fn, "w+")
    fn_open.write("{}\n".format(gb_acc))
    fn_open.close()
    # cmd1 = "blastdbcmd -db {}/nt_v5 -entry_batch {} -outfmt %T -out {}/tmp/tax_id_{}.csv".format(blastdb, fn, workdir,
    #                                                                                             gb_acc)
    cmd1 = "blastdbcmd -db {} -entry_batch {} -outfmt %T -out {}/tmp/tax_id_{}.csv".format(blastdb, fn, workdir,
                                                                                                 gb_acc)
    fn_out = "{}/tmp/tax_id_{}.csv".format(workdir, gb_acc)
    if not os.path.exists(fn_out):
        with suppress_stdout():
            os.system(cmd1)
    f = open(fn_out)
    tax_id_l = []
    for line in iter(f):
        line = line.rstrip().lstrip()
        tax_id_l.append(int(line))
    f.close()
    return tax_id_l


# todo: not used

def get_taxid_from_acc_stdout(gb_acc, blastdb, workdir):
    """
    Use the blastdb to get the taxon id from a queried gb acc.

    Sometimes there are more than a single id, as database has merged redundant seqs.

    :param gb_acc: Genbank accession number
    :param workdir: working directory
    :param blastdb: path to blast db
    :return: list of taxon ids associated with the GB id - there are multiple because of the merging of redundant data
    """
    # print('get_taxid_from_acc_stdout')
    if not os.path.exists(os.path.join(workdir, "tmp")):
        os.mkdir(os.path.join(workdir, "tmp"))
    fn = os.path.join(workdir, "tmp/tmp_search.csv")
    fn_open = open(fn, "w+")
    fn_open.write("{}\n".format(gb_acc))
    fn_open.close()
    #cmd1 = "blastdbcmd -db {}/nt_v5 -entry_batch {} -outfmt %T -out -".format(blastdb, fn)
    cmd1 = "blastdbcmd -db {} -entry_batch {} -outfmt %T -out -".format(blastdb, fn)
    result = subprocess.check_output(cmd1, shell=True).decode(sys.stdout.encoding)  # todo maybe replace with .run
    taxid_l = result.split('\n')
    taxid_l = list(filter(None, taxid_l))
    return taxid_l


def write_local_blast_files(workdir, seq_id, seq, db=False, fn=None):
    """
    Writes the files needed to run a local blast search.
    Output will be read by run_blast_query()

    :param workdir: working directory
    :param seq_id: unique sequence identifier: Genbank accession, localID
    :param seq: sequence to write
    :param db: if True, will be written to database file instead
    :param fn: file name
    :return: files with sequences written to it in fasta format
    """
    # print("writing blast files")
    if not os.path.exists(os.path.join(workdir, "tmp")):
        os.mkdir(os.path.join(workdir, "tmp"))
    if db:
        if fn is None:
            fnw = os.path.join(workdir, "tmp/filter_seq_db")
        else:
            fnw = os.path.join(workdir, "tmp/{}_db".format(fn))
        fi_o = open(fnw, "a+")
    else:
        fnw = os.path.join(workdir, "tmp/query_seq.fas")
        fi_o = open(fnw, "w")
    fi_o.write(">{}\n".format(seq_id))
    fi_o.write("{}\n".format(str(seq).replace("-", "")))
    fi_o.close()


def make_local_blastdb(workdir, db, taxid, path_to_db=None):
    """
    Adds sequences into a new blast database, which then can be used to blast aln seq against it
    and adds sequences that were found to be similar to input.
    If this option is used, it queries against local database first and only in "2" round
    it goes back to blasting against GenBank.

    :param workdir:  where to write the files
    :param db: string that defines, if local, genbank or filter
    :param taxid: taxon id/name of file
    :param path_to_db: path where db is stored
    :return: writes local blast databases for the local sequences
    """
    # print("make_local_blastdb")
    if db == "unpublished" or db == 'own':
        taxid_map = os.path.abspath(taxid)
        print('Make local blast database from: {}'.format(path_to_db))
        localfiles = os.listdir(path_to_db)

        if db == 'unpublished':
            path_db = os.path.join(workdir, 'tmp/unpublished_seq_db')
        else:
            path_db = os.path.join(path_to_db, 'db/own_seq_db')
        if os.path.exists(path_db):
            os.remove(path_db)
        db_abspath = os.path.abspath(path_db)
        fn = '{}_seq'.format(db)

        for index, item in enumerate(localfiles):
            item = str(item)
            if item.startswith(".~"):  # remove those files from list
                localfiles[index] = None

        localfiles = list(filter(None, localfiles))
        for filename in localfiles:
            filepath = os.path.join(path_to_db, filename)
            open_file = open(filepath)
            content = open_file.readlines()
            content = [x.strip() for x in content]
            count = 0
            gb_id_l = content[::2]
            gb_id_l = list(filter(None, gb_id_l))
            seq_l = content[1::2]
            seq_l = list(filter(None, seq_l))
            # note: a file with multiple seqs can be read in as well
            assert len(gb_id_l) == len(seq_l), (len(gb_id_l), len(seq_l), gb_id_l, seq_l)

            for i in range(0, len(gb_id_l)):
                key = gb_id_l[i].replace(">", "")
                count = count + 1
                seq = seq_l[i]
                write_local_blast_files(workdir, key, seq, db=True, fn=fn)

        with cd(path_to_db):
            cmd1 = "makeblastdb -in {} -dbtype nucl -taxid_map {} -parse_seqids".format(db_abspath, taxid_map)
            with suppress_stdout():
                os.system(cmd1)
    elif db == "filterrun":
        with cd(workdir):
            cmd1 = "makeblastdb -in ./tmp/filter_seq_db -dbtype nucl -parse_seqids -taxid {}".format(int(taxid))
            with suppress_stdout():
                os.system(cmd1)


def get_full_seq_stdout(gb_acc, blast_seq, workdir, blastdb, db):
    """
    Get full sequence from gb_acc that was retrieved via blast.

    Currently only used for local searches.
    Genbank database sequences are retrieving them in batch mode, which is hopefully faster.

    :param gb_acc: unique sequence identifier (often genbank accession number)
    :param blast_seq: sequence retrived by blast
    :param workdir: working directory
    :param blastdb: location of file/db with full seq
    :param db: type of db - local, filter, Genbank
    :return: full sequence, the whole submitted sequence, not only the part that matched the blast query sequence
    """
    # print("get full seq stdout")
    if db is not "Genbank":  # no need to make a db first (it already exists), we just open it and get full seq
        seq_set = False
        seq, seq_set = get_seq_from_file(gb_acc, blastdb, seq_set)
        if seq_set is False:
            seq, seq_set = get_seq_from_file(gb_acc, os.path.join(workdir, 'tmp/query_seq.fas'), seq_set)
        assert seq_set is True
    else:
        # for the Genbank db query it runs using stdout
        fn = os.path.join(workdir, "tmp/tmp_search.csv")
        fn_open = open(fn, "w+")
        fn_open.write("{}\n".format(gb_acc.split(".")[0]))
        fn_open.close()
        db_path = "{}".format(blastdb)
        #cmd1 = "blastdbcmd -db {}/nt_v5  -entry_batch {} -outfmt %f -out -".format(db_path, fn)
        cmd1 = "blastdbcmd -db {}  -entry_batch {} -outfmt %f -out -".format(db_path, fn)
        result = subprocess.check_output(cmd1, shell=True).decode(sys.stdout.encoding)
        seqn = result.split('\n')[1:]
        seperator = ''
        seq = seperator.join(seqn)
        acc_str = result.split('\n')[0]
        assert gb_acc in acc_str, (gb_acc, acc_str)
        # print(gb_acc)
    full_seq = check_directionality(blast_seq, seq)
    return full_seq

# todo get full seq tmp and stdout both not used, I still use the one using the file...


def get_full_seq_tmp(gb_acc, blast_seq, workdir, blastdb, db):
    """
    Get full sequence from gb_acc that was retrieved via blast.

    Currently only used for local searches,

    Genbank database sequences are retrieving them in batch mode, which is hopefully faster.

    :param gb_acc: unique sequence identifier (often genbank accession number)
    :param blast_seq: sequence retrived by blast
    :param workdir: working directory
    :param blastdb: location of file/db with full seq
    :param db: type of db - local, filter, Genbank
    :return: full sequence, the whole submitted sequence, not only the part that matched the blast query sequence
    """
    # print("get_full_seq_tmp")
    if db is not "Genbank":  # no need to make a db first (it already exists), we just open it and get full seq
        seq_set = False  # refactor: put seq_set var into function
        seq, seq_set = get_seq_from_file(gb_acc, blastdb, seq_set)
        if seq_set is False:
            seq, seq_set = get_seq_from_file(gb_acc, os.path.join(workdir, 'tmp/query_seq.fas'), seq_set)
        assert seq_set is True
    else:
        # for the Genbank db query it runs using stdout
        fn = os.path.join(workdir, "tmp/tmp_search.csv")
        fn_open = open(fn, "w+")
        fn_open.write("{}\n".format(gb_acc.split(".")[0]))
        fn_open.close()
        db_path = "{}".format(blastdb)
        with NamedTemporaryFile() as f:
            subprocess.check_call(
                ["blastdbcmd", "-db", "{}".format(db_path), "-entry_batch", "{}".format(fn), "-outfmt", "%f",
                 "-out", "-"], stdout=f)
            f.seek(0)
            result = f.read().decode(sys.stdout.encoding)
        seqn = result.split('\n')[1:]
        seperator = ''
        seq = seperator.join(seqn)
        acc_str = result.split('\n')[0]
        assert gb_acc in acc_str, (gb_acc, acc_str)
    full_seq = check_directionality(blast_seq, seq)
    return full_seq


def get_full_seq(gb_acc, blast_seq, workdir, blastdb, db):
    """
    Get full sequence from gb_acc that was retrieved via blast.

    Currently only used for local searches,
    Genbank database sequences are retrieving them in batch mode, which is hopefully faster.

    :param gb_acc: unique sequence identifier (often genbank accession number)
    :param blast_seq: sequence retrived by blast
    :param workdir: working directory
    :param blastdb: location of file/db with full seq
    :param db: type of db - local, filter, Genbank
    :return: full sequence, the whole submitted sequence, not only the part that matched the blast query sequence
    """
    # print("get full seq")
    if db is not "Genbank":  # no need to make a db first (it already exists), we just open it and get full seq
        seq_set = False
        seq, seq_set = get_seq_from_file(gb_acc, blastdb, seq_set)
        if seq_set is False:
            seq, seq_set = get_seq_from_file(gb_acc, '{}/tmp/query_seq.fas'.format(workdir), seq_set)
        assert seq_set is True
    else:
        if not os.path.exists(os.path.join(workdir, "tmp")):
            os.mkdir(os.path.join(workdir, "tmp"))
        full_seq_fn = os.path.join(workdir, "tmp/full_seq_{}.fasta".format(gb_acc))
        if os.path.exists(full_seq_fn):
            if not os.stat(full_seq_fn).st_size > 0:
                os.remove(full_seq_fn)
        if not os.path.exists(full_seq_fn):
            # print('get full seq blast query')
            fn = os.path.join(workdir, "tmp/tmp_search.csv")
            fn_open = open(fn, "w+")
            fn_open.write("{}\n".format(gb_acc.split(".")[0]))
            fn_open.close()
            db_path = "{}".format(blastdb)
            #cmd1 = "blastdbcmd -db {}/nt_v5  -entry_batch {} -outfmt %f -out {}/tmp/full_seq_{}.fasta".format(db_path, fn, workdir, gb_acc)
            cmd1 = "blastdbcmd -db {} -entry_batch {} -outfmt %f -out {}/tmp/full_seq_{}.fasta".format(db_path, fn, workdir, gb_acc)
            # print(cmd1)
            with suppress_stdout():
                os.system(cmd1)
            #DEVNULL = open(os.devnull, 'wb')
            #subprocess.call(cmd1, shell=False, stderr=DEVNULL)
            assert os.stat(full_seq_fn).st_size > 0, ('file {} has no content'.format(full_seq_fn))
        # read in file to get full seqseq
        fn_seq = os.path.join(workdir, "tmp/full_seq_{}.fasta".format(gb_acc))
        f = open(fn_seq)
        seq = ""
        for line in iter(f):
            line = line.rstrip().lstrip()
            if line[0] != ">":
                seq += line
            elif line[0] == ">":
                assert gb_acc in line, (gb_acc, line)
        f.close()
    full_seq = check_directionality(blast_seq, seq)
    return full_seq


def get_seq_from_file(gb_acc, fn, seq_set):
    """
    Read full seq out of file.

    :param gb_acc: gb acc associated with id
    :param fn: file name where seq is in
    :param seq_set: T/F - seq found?
    :return: seq and seq_set
    """
    found = False
    seq = None
    with open(fn) as f:
        for i, line in enumerate(f):
            first_item = line.split('\t')[0].rstrip().lstrip()[1:]
            if found:
                seq = line.rstrip().lstrip()
                seq = seq.upper()
                seq_set = True
                break
            elif gb_acc == first_item:
                found = True
    return seq, seq_set


def check_directionality(blast_seq, seq):
    """
    Check direction of sequence as seq is sometimes not stored in Genbank with the correct directionality.

    :param blast_seq: seq that was found by blast
    :param seq: Genbank full seq as saved
    :return: full seq in correct directionality
    """
    assert blast_seq is not "", blast_seq
    assert seq is not "", seq
    orig = Seq(seq, generic_dna)  # generic_dna and Alphabet will be deprecated soon. I've not yet found a solution.
    dna_comp = orig.complement()
    dna_rcomp = orig.reverse_complement()
    dna_r = orig[::-1]
    full_seq = str()
    seq_replace = blast_seq.replace("-", "")
    if seq_replace in orig:
        full_seq = seq
    elif seq_replace in dna_r:
        full_seq = dna_r
    elif seq_replace in dna_comp:
        full_seq = dna_comp
    elif seq_replace in dna_rcomp:
        full_seq = dna_rcomp
    assert seq_replace in full_seq, (seq_replace, full_seq, seq)
    full_seq = str(full_seq)
    assert type(full_seq) == str, (type(full_seq))
    return full_seq


def get_new_seqs(query_seq, taxon, db_name, config, mrca=None):
    """
    Produces the pandas df with the new sequences. Main function of this class.

    :param query_seq: blast query sequence for search
    :param taxon: taxon name (=gb_acc, ncbi id); used for fn
    :param db_path: path to blast db
    :param db_name: name of blast db
    :param config: config obj
    :param mrca: ncbi id to limit blast results
    :return: returns the new sequences found via blast
    """
    # print('get new seqs')
    run_blast_query(query_seq, taxon, db_name, config, mrca)
    new_blastseqs = read_blast_query_pandas(taxon, config, db_name)  # pandas implementation of read_blast_query()
    assert new_blastseqs["ncbi_txn"].isnull() is not None, (new_blastseqs["ncbi_txn"])
    return new_blastseqs


def run_blast_query(query_seq, taxon, db_name, config, mrca=None):
    """
    Contains the cmds used to run a local blast query.

    :param query_seq: blast query sequence for search
    :param taxon: taxon name (=gb_acc, ncbi id); used for fn
    :param db_name: name of blast db
    :param config: config obj
    :param mrca: optional, new blast+ can limit results to list of taxids - not yet implemented in code
    :return: runs local blast query and writes it to file
    """
    # print("run_blast_query")
    #errf = open("./{}/blast_stderr".format(config.workdir), 'a')

    assert taxon not in [None, 'nan', 'NA', 'na']
    taxon = str(taxon)
    if len(taxon.split('.')) > 1:
        taxon = taxon.split('.')[0]
    if db_name == "unpublished":  # Run a local blast search if the data is unpublished or for filtering.
        query_output_fn = os.path.join(config.workdir, "tmp/unpublished_query_result.txt")
        input_fn = os.path.join(config.blast_folder, "./query_seq.fas")
        db = 'unpublished_seq_db'
    elif db_name == "filterrun":  # Run a local blast search if the data is unpublished or for filtering.
        query_output_fn = os.path.join(config.workdir, "tmp/{}.txt".format(taxon))
        input_fn = os.path.join(config.workdir, "tmp/{}_tobeblasted.fas".format(taxon))
        db = 'filter_seq_db'
    elif db_name == "Genbank":
        query_output_fn = os.path.join(config.blast_folder, "./{}.txt".format(taxon))
        input_fn = os.path.join(config.blast_folder, "./{}_tobeblasted.fas".format(taxon))
    else:
        print('DOES THIS EVER HAPPEN?')
        print(db_name)
        sys.exit(-33)
    query_output_fn = os.path.abspath(query_output_fn)
    input_fn = os.path.abspath(input_fn)

    toblast = open(input_fn, "w")
    toblast.write(">{}\n".format(taxon))
    toblast.write("{}\n".format(query_seq))
    toblast.close()

    # print('run blastn')
    outfmt = " -outfmt '6 sseqid staxids sscinames pident evalue bitscore sseq salltitles sallseqid'"

    DEVNULL = open(os.devnull, 'wb')
    if db_name == "Genbank":
        with cd(os.path.abspath(config.blastdb_path)):
            # TODO MK: blast+ v. 2.8 code - then we can limit search to taxids: -taxids self.mrca_ncbi
            #  !no mrca information here avail! - needs also some form of taxonomy db - unsure which
            #blastcmd = "blastn -query " + input_fn + " -db {}/nt_v5 -out ".format(os.path.abspath(config.blastdb)) + \
            #           query_output_fn + " {} -num_threads {}".format(outfmt, config.num_threads) + \
            #           " -max_target_seqs {} -max_hsps {}".format(config.hitlist_size, config.hitlist_size)
            blastcmd = "blastn -query " + input_fn + " -db {} -out ".format(os.path.abspath(config.blastdb)) + \
                       query_output_fn + " {} -num_threads {}".format(outfmt, config.num_threads) + \
                       " -max_target_seqs {} -max_hsps {}".format(config.hitlist_size, config.hitlist_size)
            systemcmd = "{} | tee ./{}/blast_stderr, shell = False".format(blastcmd, config.workdir)
                        # "-evalue {} - taxids {}".format(config.evalue, mrca)
            # needs to run from within the folder:
            with suppress_stdout():
                if not os.path.isfile(query_output_fn):
                    #subprocess.call(blastcmd, shell=False, stderr=DEVNULL)
                    os.system(blastcmd)
                    # print(blastcmd)
                elif not os.stat(query_output_fn).st_size > 0:
                    #subprocess.call(blastcmd, shell=False, stderr=DEVNULL)
                    os.system(blastcmd)

    else:
        # print("run against local data")
        with cd(os.path.join(config.workdir, "tmp")):
            blastcmd = "blastn -query {} -db {} -out ".format(input_fn, db) + query_output_fn + \
                       " {} -num_threads {}".format(outfmt, config.num_threads) + \
                       " -max_target_seqs {} -max_hsps {}".format(config.hitlist_size_unpublished, config.hitlist_size_unpublished)
            with suppress_stdout():
                os.system(blastcmd)
            #subprocess.call(blastcmd, shell=False, stderr=DEVNULL)
            # todo: produces blastn taxdb warning, taxids here not needed and not part of taxdb anyways as local search.
            #  I keep it to make it coherent for reading in results.


def read_blast_query_pandas(blast_fn, config, db_name):
    """
    Implementation to read in results of local blast searches.

    :param blast_fn: filename to read in; often taxon or gb_acc
    :param config: config object
    :param db_name: name that specifies if filterrun, unpublished or normal search
    :return: updated self.new_seqs and self.data.gb_dict dictionaries
    """
    debug('read_blast_query_pandas')
    blast_fn = str(blast_fn)
    if len(blast_fn.split('.')) > 1:
        blast_fn = blast_fn.split('.')[0]
    if db_name == "filterrun":  # Run a local blast search if the data is unpublished or for filtering.
        query_output_fn = os.path.join(config.workdir, "tmp/{}.txt".format(blast_fn))
    elif db_name == 'unpublished':
        query_output_fn = os.path.join(config.workdir, "tmp/unpublished_query_result.txt")
    else:
        query_output_fn = os.path.join(config.blast_folder, "./{}.txt".format(blast_fn))
    query_output_fn = os.path.abspath(query_output_fn)
    colnames = ['accession;gi', 'ncbi_txid', 'ncbi_txn', 'pident', 'evalue', 'bitscore', 'sseq', 'title', 'accession']
    data = pd.read_csv(query_output_fn, names=colnames, sep="\t", header=None)

    if data['ncbi_txid'].isnull().values.any() == True:
        colnames = ['accession;gi', 'ncbi_txid']
        tax_id_map = pd.read_csv(os.path.join(config.workdir, 'tmp/unpublished_txid.txt'), names=colnames, sep=" ",
                                 header=None)
        before = len(data)
        for idx in data.index:
            if data.loc[idx, 'ncbi_txid'] == None:
                acc = data.loc[idx, 'accession;gi']
                map_id = tax_id_map[tax_id_map['accession;gi'] == acc, 'ncbi_txid']
                data.loc[idx, 'ncbi_txid'] = map_id
        after=len(data)
        assert after == before, (before, after)

    if not db_name == 'unpublished':
        assert len(data) > 0, data
    redundant = data[data['accession'].str.contains(';')]
    non_redundant_redundant = get_non_redundant_data(config, redundant)  # new data frame with split value columns
    non_redundant = data[data['accession'].str.contains(';') == False]
    non_redundant['accession'] = non_redundant['accession'].apply(get_acc_from_blast)  # todo: throws pandas warning about sort, but if added, its breaking
    assert len(redundant) + len(non_redundant) == len(data)
    new_seqs = pd.concat([non_redundant_redundant, non_redundant], sort=True, ignore_index=True)


    assert new_seqs['ncbi_txid'].isnull().values.any() == False, ('Not all new seqs have a taxon id. Often related to error in input of unpublished.')
    new_seqs['sseq'] = new_seqs['sseq'].str.replace("-", "")
    new_seqs['date'] = datetime.datetime.strptime('01/01/00', '%d/%m/%y')
    # todo this could be made faster by running it on the redundant/non_redundant data first
    new_seqs = wrapper_get_fullseq(config, new_seqs, db_name)
    return new_seqs


def get_non_redundant_data(config, redundant):
    """

    Get per taxonid the information into a pandas dataframe.

    Structure of redundant data (tab-delimited):
    gi|429489218|gb|JX895383.1|	1268578;1269244	Senecio provincialis;Senecio provincialis x Senecio lagascanus
    95.879	0.0	1177
    AACAAGGTTTCCGTCCAAGAAGTAAGGAATATCTCTTTAATGACCCTAAAGTGTTGTCTCATGACGATGCTTCGACTGCGACCCCAGGTCAGGCGGGATTAGCATATCAAT
    Senecio provincialis clone JC4662-14 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S
    ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene,
    partial sequence<>Senecio provincialis clone JC4659-14 18S ribosomal RNA gene, partial sequence;
    internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence;
    and 28S ribosomal RNA gene, partial sequence<>Senecio provincialis x Senecio lagascanus clone JC3622-8 18S
    ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal
    transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence
    gi|429489218|gb|JX895383.1|;gi|429489247|gb|JX895412.1|;gi|429489255|gb|JX895420.1|


    :param config:
    :param redundant:
    :return:
    """
    debug('get_non_redundant_data')
    colnames = ['accession;gi', 'ncbi_txid', 'ncbi_txn', 'pident', 'evalue', 'bitscore', 'sseq', 'title', 'accession']
    non_redundant_redundant = pd.DataFrame(columns=colnames)
    queried_acc = set()  # used to test if gb_acc was added before  aka query_dict in physcraper

    acc_column = redundant["accession"].str.split(";", n=-1, expand=True)  # todo: try to expand in row instead of column
    # todo: think about making a single df first where unique acc corresponds to taxid and same index val
    #  as in redundant, then use that whole df to make new non_redundant_redundant df. - discussion moritz
    for idx in acc_column.index:
        gb_acc = get_acc_from_blast(redundant.loc[idx, 'accession;gi'])
        all_taxids = get_taxid_from_acc(gb_acc, config.blastdb, config.workdir)

        # ##################################
        # #todo: moritz suggestions
        # all_accs = redundant.loc[idx, 'accession'].str.split(';')  # get list of all accs
        # for i in len(all_taxids):
        #     if txid not in used_txids:
        #         # get all relevant data - create new line for those
        #
        #         # add txid to used_taxids
        #
        # # todo: to make it faster get indices of txids
        # ##################################

        all_taxids_set = set(all_taxids)
        found_taxids = set()
        t = acc_column.loc[idx]
        shortened = t.dropna()
        ncbi_parser = ncbi_data_parser.Parser(
            names_file=config.ncbi_parser_names_fn,
            nodes_file=config.ncbi_parser_nodes_fn)
        while len(found_taxids) < len(set(all_taxids_set)):
            for i in range(len(shortened)):
                # if we found all taxon_ids present in the initial list, stop looking for more
                if len(found_taxids) == len(all_taxids_set):
                    break  # this leaves the for loop, no need to iterate through if we have all ids
                tax_id = all_taxids[i]
                if tax_id in found_taxids:
                    continue  # this goes back to beginning of for loop
                gb_acc = get_acc_from_blast(shortened.loc[i])
                spn = ncbi_parser.get_name_from_id(tax_id)  # todo could be done later before creating table - for efficiency do later
                split_df = redundant.loc[[idx]]
                split_df['accession'] = gb_acc
                split_df['ncbi_txid'] = tax_id
                split_df['ncbi_txn'] = spn
                # non_redundant_redundant = non_redundant_redundant.append(split_df, ignore_index=True, sort=True)
                non_redundant_redundant = pd.concat([non_redundant_redundant, split_df], ignore_index=True, sort=True)
                found_taxids.add(tax_id)
                queried_acc.add(gb_acc)
    if len(redundant) > 0:
        assert len(non_redundant_redundant) >= len(redundant), \
            (len(non_redundant_redundant), len(redundant), non_redundant_redundant['accession'],
             redundant['accession;gi'])
    return non_redundant_redundant


def wrapper_get_fullseq(config, new_blast_seq_dict, db):
    """
    Wrapper function to get the full sequences from acc

    :param config: config data
    :param new_blast_seq_dict: new_seqs dataframe
    :param db: name of the database
    :return: updated blast seq dict, with long seq
    """
    if db == 'Genbank':
        blastdb = config.blastdb
    elif db == 'unpublished':
        blastdb = os.path.join(config.workdir, 'tmp/unpublished_seq_db')
    else:
        blastdb = os.path.join(config.workdir, 'tmp/filter_seq_db')
    for idx in new_blast_seq_dict.index:
        gb_acc = new_blast_seq_dict.loc[idx, "accession"]
        if len(gb_acc.split(".")) >= 2 or db == 'unpublished':  # Do not add sequences that are not in Genbank accession format, e.g. PDB
            # replace sequence
            full_seq = get_full_seq(gb_acc, new_blast_seq_dict.loc[idx, "sseq"], config.workdir, blastdb, db)
            # new_blast_seq_dict.loc[idx, "sseq"] = full_seq
            new_blast_seq_dict.at[idx, "sseq"] = full_seq  # note: set_value is much faster than .loc method, but deprecated, try at
    return new_blast_seq_dict
