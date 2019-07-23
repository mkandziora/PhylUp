# class sequence and id retrieval
import os
import datetime
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from . import ncbi_data_parser  # is the ncbi data parser class and associated functions
from . import cd


def get_acc_from_blast(query_string):
    """
    Get the accession number from a blast query using local blast.

    """
    # Note: get acc is more difficult now, as new seqs not always have gi number, then query changes
    if len(query_string.split("|")) >= 3:
        gb_acc = query_string.split("|")[3]
    else:
        gb_acc = query_string.split("|")[0]
    # assert len(gb_acc.split(".")) >= 2, (len(gb_acc.split(".")), gb_acc)
    return gb_acc


# todo unused
def get_gi_from_blast(query_string):
    """
    Get the gi number from a blast query using local blast. If not available return None.

    Note: get acc is more difficult now, as new seqs not always have gi number, then query changes

    :param query_string:
    :return:
    """
    if len(query_string.split("|")) >= 3:
        gb_id = query_string.split("|")[1]
    assert len(gb_id.split(".")) < 2, (len(gb_id.split(".")), gb_id)
    return gb_id


def get_taxid_from_acc(gb_acc, blastdb, workdir):
    """
    Use the blastdb to get the taxon id from a queried gb acc.

    Sometimes there are more than a single id, as database has merged redundant seqs.

    :param gb_acc: Genbank accession number
    :param workdir: working directory
    :param blastdb: path to blast db
    :return: list of taxon ids associated with the GB id - there are multiple because of the merging of redundant data
    """
    if not os.path.exists("{}/tmp".format(workdir)):
        os.mkdir("{}/tmp".format(workdir))
    fn = "{}/tmp/tmp_search.csv".format(workdir)
    fn_open = open(fn, "w+")
    fn_open.write("{}\n".format(gb_acc))
    fn_open.close()

    cmd1 = "blastdbcmd -db {}/nt_v5 -entry_batch {} -outfmt %T -out {}/tmp/tax_id_{}.csv".format(blastdb, fn, workdir,
                                                                                                 gb_acc)
    os.system(cmd1)
    f = open("{}/tmp/tax_id_{}.csv".format(workdir, gb_acc))
    tax_id_l = []
    for line in iter(f):
        line = line.rstrip().lstrip()
        tax_id_l.append(int(line))
    f.close()
    return tax_id_l


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
    if not os.path.exists("{}/".format(workdir)):
        os.makedirs("{}/tmp/".format(workdir))
    if db:
        if fn is None:
            fnw = "{}/tmp/filter_seq_db".format(workdir)
        else:
            fnw = "{}/tmp/{}_db".format(workdir, fn)
        fi_o = open(fnw, "a")
    else:
        fnw = "{}/tmp/query_seq.fas".format(workdir, fn)
        fi_o = open(fnw, "w")
    fi_o.write(">{}\n".format(seq_id))
    fi_o.write("{}\n".format(str(seq).replace("-", "")))
    fi_o.close()


def make_local_blastdb(workdir, db, path_to_db=None):
    """
    Adds sequences into a  new blast database, which then can be used to blast aln seq against it
    and adds sequences that were found to be similar to input.
    If this option is used, it queries against local database first and only in "2" round
    it goes back to blasting against GenBank.

    :param workdir:  where to write the files
    :param db: string that defines, if local, genbank or filter
    :param path_to_db: path where db is stored
    :return: writes local blast databases for the local sequences
    """
    print("make_local_blastdb")

    if db == "local":
        # self.path_to_local_seq = path_to_local_seq
        localfiles = os.listdir(output_db_path)
        for index, item in enumerate(localfiles):
            item = str(item)
            if item.startswith(".~"):  # remove those files from list
                localfiles[index] = None
        localfiles = filter(None, localfiles)
        for filename in localfiles:
            filepath = "{}/{}".format(output_db_path, filename)
            open_file = open(filepath)
            content = open_file.readlines()
            content = [x.strip() for x in content]
            content = filter(None, content)  # fastest
            count = 0
            gb_id_l = content[::2]
            seq_l = content[1::2]
            # in current setup 1 seq per file, but this is written in a way,
            # that a file with multiple seqs can be read in as well
            os.remove('{}/tmp/filter_seq_db'.format(workdir))
            len_gb = len(gb_id_l)
            for i in range(0, len_gb):
                key = gb_id_l[i].replace(">", "")
                count = count + 1
                seq = seq_l[i]
                write_local_blast_files(workdir, key, seq, db=True, fn="local_seq")
        with cd(output_db_path):
            cmd1 = "makeblastdb -in ./tmp/local_seq_db -dbtype nucl"
            os.system(cmd1)
    elif db == "filterrun":
        with cd(workdir):
            cmd1 = "makeblastdb -in ./tmp/filter_seq_db -dbtype nucl"
            os.system(cmd1)


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
        if not os.path.exists("{}/tmp".format(workdir)):
            os.mkdir("{}/tmp".format(workdir))
        if not os.path.exists("{}/tmp/full_seq_{}.fasta".format(workdir, gb_acc)):
            # print('get full seq blast query')
            fn = "{}/tmp/tmp_search.csv".format(workdir)
            fn_open = open(fn, "w+")
            fn_open.write("{}\n".format(gb_acc.split(".")[0]))
            fn_open.close()
            db_path = "{}".format(blastdb)
            cmd1 = "blastdbcmd -db {}/nt_v5  -entry_batch {} " \
                   "-outfmt %f -out {}/tmp/full_seq_{}.fasta".format(db_path, fn, workdir, gb_acc)
            os.system(cmd1)
        # else:
        #     print('get full seq blast query done earlier')

        # read in file to get full seq
        fn_seq = "{}/tmp/full_seq_{}.fasta".format(workdir, gb_acc)
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
            if found:
                seq = line.rstrip().lstrip()
                seq = seq.upper()
                seq_set = True
                break
            elif gb_acc in line:
                found = True
    return seq, seq_set


def check_directionality(blast_seq, seq):
    """
    Check direction of sequence as seq is sometimes not stored in Genbank with the correct directionality.

    :param blast_seq: seq that was found by blast
    :param seq: Genbank full seq as saved
    :return: full seq in correct directionality
    """
    orig = Seq(seq, generic_dna)
    dna_comp = orig.complement()
    dna_rcomp = orig.reverse_complement()
    dna_r = orig[::-1]
    full_seq = str()
    if blast_seq.replace("-", "") in orig:
        full_seq = seq
    elif blast_seq.replace("-", "") in dna_r:
        full_seq = dna_r
    elif blast_seq.replace("-", "") in dna_comp:
        full_seq = dna_comp
    elif blast_seq.replace("-", "") in dna_rcomp:
        full_seq = dna_rcomp
    assert blast_seq.replace("-", "") in full_seq, (blast_seq.replace("-", ""), full_seq, seq)
    full_seq = str(full_seq)
    assert type(full_seq) == str, (type(full_seq))
    return full_seq


def get_new_seqs(query_seq, taxon, db_path, db_name, config):
    """
    Produces the pandas df with the new sequences. Main function of this class.

    :param query_seq: blast query sequence for search
    :param taxon: taxon name (=gb_acc, ncbi id); used for fn
    :param db_path: path to blast db
    :param db_name: name of blast db
    :param config: config obj
    :return:
    """
    run_blast_query(query_seq, taxon, db_path, db_name, config)
    new_blastseqs = read_blast_query_pandas(taxon, config, db_name)  # pandas implementation of read_blast_query()
    assert new_blastseqs["ncbi_txn"].isnull() is not None, (new_blastseqs["ncbi_txn"])
    return new_blastseqs


def run_blast_query(query_seq, taxon, db_path, db_name, config):
    """
    Contains the cmds used to run a local blast query.

    :param query_seq:
    :param taxon:
    :param db_path:
    :param db_name:
    :param config:
    :return: runs local blast query and writes it to file
    """
    # print("run_blast_query")
    assert taxon not in [None, 'nan', 'NA', 'na']

    taxon = str(taxon)
    if len(taxon.split('.')) > 1:
        taxon = taxon.split('.')[0]
    db_path = os.path.abspath(db_path)
    if db_name == "unpublished":  # Run a local blast search if the data is unpublished or for filtering.
        query_output_fn = "{}/tmp/unpublished_query_result.txt".format(config.workdir)
        input_fn = "{}/blast/query_seq.fas".format(config.workdir)
        db = 'unpublished_seq_db'
    elif db_name == "filterrun":  # Run a local blast search if the data is unpublished or for filtering.
        query_output_fn = "{}/tmp/{}.txt".format(config.workdir, taxon)
        input_fn = "{}/tmp/{}_tobeblasted.fas".format(config.workdir, taxon)
        db = 'filter_seq_db'
    else:
        query_output_fn = "{}/blast/{}.txt".format(config.workdir, taxon)
        input_fn = "{}/blast/{}_tobeblasted.fas".format(config.workdir, taxon)
    query_output_fn = os.path.abspath(query_output_fn)
    input_fn = os.path.abspath(input_fn)

    toblast = open(input_fn, "w")
    toblast.write(">{}\n".format(taxon))
    toblast.write("{}\n".format(query_seq))
    toblast.close()

    outfmt = " -outfmt '6 sseqid staxids sscinames pident evalue bitscore sseq salltitles sallseqid'"

    if db_name == "Genbank":
        with cd(db_path):
            # this format (6) allows to get the taxonomic information at the same time
            # outfmt = " -outfmt 5"  # format for xml file type
            # TODO MK: blast+ v. 2.8 code - then we can limit search to taxids: -taxids self.mrca_ncbi
            blastcmd = "blastn -query " + input_fn + \
                       " -db {}/nt_v5 -out ".format(db_path) + query_output_fn + \
                       " {} -num_threads {}".format(outfmt, config.num_threads) + \
                       " -max_target_seqs {} -max_hsps {}".format(config.hitlist_size,
                                                                  config.hitlist_size)
            # needs to run from within the folder:
            if not os.path.isfile(query_output_fn):
                os.system(blastcmd)
            elif not os.stat(query_output_fn).st_size > 0:
                os.system(blastcmd)
    else:
        print("run against local data")
        with cd("{}/tmp/".format(config.workdir)):
            blastcmd = "blastn -query {} -db {} -out ".format(input_fn, db) + query_output_fn + \
                       " {} -num_threads {}".format(outfmt, config.num_threads) + \
                       " -max_target_seqs {} -max_hsps {}".format(config.hitlist_size, config.hitlist_size)
            os.system(blastcmd)
            # todo: produces blastn taxdb warning, taxids here not needed and not part of taxdb anyways as local search.
            # I keep it to make it coherent for reading in results.


def read_blast_query(taxon, config, db_name):
    """
    Implementation to read in results of local blast searches.

    :param taxon:


def wrapper_get_fullseq(config, new_blast_seq_dict, db):
    """

    :param config:
    :param new_blast_seq_dict:
    :param db:
    :return:
    """
    # print(db)
    if db == 'Genbank':
        blastdb = config.blastdb
    else:
        blastdb = '{}/tmp/filter_seq_db'.format(config.workdir)
    for idx in new_blast_seq_dict.index:
        gb_acc = new_blast_seq_dict.loc[idx, "accession"]
        if len(gb_acc.split(".")) >= 2:  # Do not add sequences that are not in Genbank accession format, e.g. PDB
            # replace sequence
            full_seq = get_full_seq(gb_acc, new_blast_seq_dict.loc[idx, "sseq"], config.workdir, blastdb, db)
            new_blast_seq_dict.loc[idx, "sseq"] = full_seq
    return new_blast_seq_dict


# ###########################################################################################


def get_full_seq_batch(list_gb_acc, blast_db, workdir):
    """
    Get full sequences using a batch query. Is likely faster than doing it one by one.

    Implemented currently only to use with the Genbank database, not for local databases.

    :param blast_db: path of db
    :param list_gb_acc: list with Genbank accession numbers, where we want to get the full sequence from the db
    :param workdir: working directory
    :return: dictionary with full sequences
    """
    print("get_full_seq_batch")

    if not os.path.exists("{}/tmp".format(workdir)):
        os.mkdir("{}/tmp".format(workdir))
    fn = "{}/tmp/tmp_search.csv".format(workdir)
    fn_open = open(fn, "w+")
    for gb_acc in list_gb_acc:
        gb_nv = gb_acc.split(".")[0]
        fn_open.write("{}\n".format(gb_acc))

        fn_open.write("{}\n".format(gb_nv))
    fn_open.close()
    if len(list_gb_acc.keys()) >= 2:
        cmd1 = "blastdbcmd -db {} -entry_batch {} -outfmt %f " \
               "-out {}/tmp/all_full_seqs.fasta".format(blast_db, fn, workdir)
    else:
        cmd1 = "blastdbcmd -db {} -entry {} -outfmt %f " \
               "-out {}/tmp/all_full_seqs.fasta".format(blast_db, list_gb_acc.keys()[0], workdir)
    os.system(cmd1)
    # read in file to get full seq
    fn = "{}/tmp/all_full_seqs.fasta".format(workdir)
    # assert that every gb_acc is found in file
    for item in list_gb_acc:
        found = False
        df = open(fn)
        for line in df:
            if item in line:
                found = True
        assert found is True, item
    # find correct strand for each gb_acc
    f = open(fn)
    seq = ""
    full_seqs = {}  # dictionary to be filled with full seqs that are in the correct direction
    gb_acc_l_intern = set()
    count = 0
    for line in iter(f):
        line = line.rstrip().lstrip()
        if line[0] != ">":
            seq += line
        elif line[0] == ">":
            count += 1
            print(count)
            if seq != "":
                # print("get following seqs")
                full_seqs = find_correct_strand_batch(full_seqs, gb_acc_l_intern, list_gb_acc, seq)
                print(len(full_seqs.keys()))
                # seq = ""  # seems to not be needed
            # now get new gb_acc for this line
            splitline = line.split(">")
            gb_acc_l_intern = set()  # reassign for new list
            del splitline[0]  # removes first element of list, which is an empty string bc of split ">"
            for item in splitline:
                gb_acc = get_acc_from_blast(item)
                if gb_acc in list_gb_acc:
                    gb_acc_l_intern.add(gb_acc)
            seq = ""  # make new empty string for next seq
    f.close()
    # get information for last for loop cycle
    full_seqs = find_correct_strand_batch(full_seqs, gb_acc_l_intern, list_gb_acc, seq)
    assert set(full_seqs.keys()) == set(list_gb_acc.keys()), (
        "missing in query_dict:", [x for x in set(full_seqs.keys()) if x not in set(list_gb_acc.keys())],
        "missing in full seqs:", [x for x in set(list_gb_acc.keys()) if x not in set(full_seqs.keys())])
    return full_seqs


def find_correct_strand_batch(full_seqs, gb_acc_l_intern, list_gb_acc, seq):
    """
    Find the right direction for the sequence to be added.

    :param full_seqs: is a dictionary to where the full seqs are added
    :param gb_acc_l_intern: is a list of all the gb_acc that are merged into one sequence
    :param list_gb_acc: dictionary with all gb_acc and short_blast_seqs for which we want to get the longer sequences
    :param seq: is the longer sequence from the blast db
    :return: dictionary with the full sequences in the correct direction  is returned
    """
    # print("find_correct_strand")
    assert seq != ""
    for gb_acc in gb_acc_l_intern:
        assert gb_acc is not None
        blast_seq = list_gb_acc[gb_acc]
        orig = Seq(seq, generic_dna)
        dna_comp = orig.complement()
        dna_rcomp = orig.reverse_complement()
        dna_r = orig[::-1]
        full_seq = str()
        if blast_seq.replace("-", "") in orig:
            full_seq = seq
        elif blast_seq.replace("-", "") in dna_r:
            full_seq = dna_r
        elif blast_seq.replace("-", "") in dna_comp:
            full_seq = dna_comp
        elif blast_seq.replace("-", "") in dna_rcomp:
            full_seq = dna_rcomp
        assert blast_seq.replace("-", "") in full_seq, (blast_seq.replace("-", ""), full_seq, seq, gb_acc)
        full_seq = str(full_seq)
        assert type(full_seq) == str, (type(full_seq))
        full_seqs[gb_acc] = full_seq
    return full_seqs


def read_blast_query(blast_fn, config, db_name):
    """
    Implementation to read in results of local blast searches.

    :param blast_fn: file to read
    :param config: config object
    :param db_name:
    :return: updated self.new_seqs and self.data.gb_dict dictionaries
    """
    print("read_blast_query")
    print(blast_fn)
    ncbi_parser = ncbi_data_parser.Parser(
        names_file=config.ncbi_parser_names_fn,
        nodes_file=config.ncbi_parser_nodes_fn)
    blast_fn = str(blast_fn)
    if len(blast_fn.split('.')) > 1:
        blast_fn = blast_fn.split('.')[0]

    if db_name == "filterrun":  # Run a local blast search if the data is unpublished or for filtering.
        query_output_fn = "{}/tmp/{}.txt".format(config.workdir, blast_fn)
    else:
        query_output_fn = "{}/blast/{}.txt".format(config.workdir, blast_fn)
    query_output_fn = os.path.abspath(query_output_fn)
    columns = ['ncbi_txn', 'ncbi_txid', 'status', "date", 'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title']
    new_blast_seq_dict = pd.DataFrame(columns=columns)
    queried_acc = set()  # used to test if gb_acc was added before
    with open(query_output_fn, mode="r") as infile:
        for lin in infile:
            bitscore, evalue, gb_acc, pident, sallseqid, salltitles, \
                sscinames, sseq, staxids, stitle = get_blastval_from_line(db_name, lin, ncbi_parser, blast_fn)
            # NOTE: sometimes there are seq which are identical & are combined in the local blast db...
            # Get all of them! (they can be of a different taxon id = get redundant seq info)
            found_taxids = set()
            found_spn = set()

            # get additional info only for seq that pass the eval
            # TODO: e filter should maybe not be here...but makes it much faster
            if evalue < float(config.e_value_thresh):
                # print(sallseqid)
                if len(sallseqid.split(";")) > 1:
                    # FOR MERGED SEQS
                    sallseqid_l, salltitles_l, sscinames_l, staxids_l = split_multiple_tolist(sallseqid, salltitles,
                                                                                              sscinames, staxids)
                    # make sure you have the correct infos
                    tax_id_l = get_taxid_from_acc(gb_acc, config.blastdb, config.workdir)
                    for item in tax_id_l:
                        assert str(item) in staxids_l
                    print(tax_id_l)
                    # this while loop is here to speed up the process of finding the correct information
                    count = 0
                    stop_while = False
                    while len(found_taxids) < len(set(staxids_l)):  # as long as i have not found all taxids for the seq
                        count += 1
                        if stop_while:
                            break
                        if count == 5:
                            break  # too many tries to find correct number of redundant taxa
                        elif count == 1:
                            id_before = 0
                            for i in range(0, len(sallseqid_l)):
                                # if we found all taxon_ids present in the initial list, stop looking for more
                                if len(found_taxids) == len(staxids_l):
                                    break
                                gb_acc = get_acc_from_blast(sallseqid_l[i])

                                # if gb acc was already read in before stop the for loop
                                if gb_acc in queried_acc:
                                    stop_while = True
                                    break

                                stitle = salltitles_l[i]
                                # if both var are the same, we do not need to search GB for taxon info
                                staxids = tax_id_l[i]
                                qtaxid = int(get_taxid_from_acc(gb_acc, config.blastdb, config.workdir)[0])
                                id_now = qtaxid
                                assert str(qtaxid) in staxids_l, (str(qtaxid), staxids_l)
                                assert qtaxid in tax_id_l, (qtaxid, tax_id_l)
                                assert str(staxids) in staxids_l, (staxids, staxids_l)

                                # sometimes if multiple seqs are merged,
                                # we lack information about which taxon is which gb_acc...
                                # test it here:
                                # if we have same number of gb_acc and taxon id go ahead as usual
                                if len(sallseqid_l) == len(staxids_l):
                                    sscinames = sscinames_l[i]
                                # if only one taxon id present, all are from same taxon
                                elif len(staxids_l) == 1:
                                    sscinames = sscinames_l[0]
                                    qtaxid = staxids_l[0]
                                # if not the first item and id different from before: get name
                                elif i != 0 and id_before != id_now:
                                    sscinames = ncbi_parser.get_name_from_id(qtaxid)
                                elif i == 0:  # for first item in redundant data, always get info
                                    sscinames = ncbi_parser.get_name_from_id(qtaxid)
                                else:  # if id_before and id_now were the same, we do not need to add same seq again
                                    continue
                                # next vars are used to stop loop if all taxids were found
                                found_taxids.add(staxids)
                                found_spn.add(sscinames)
                                id_before = qtaxid
                                # add info to new_blast_seq_dict
                                if gb_acc not in queried_acc:
                                    queried_acc.add(gb_acc)
                                    new_blast_seq_dict = new_blast_seq_dict.append(
                                        {'ncbi_txn': sscinames, 'ncbi_txid': abs(int(staxids)), 'status': -100,
                                         "date": datetime.datetime.strptime('01/01/00', '%d/%m/%y'),
                                         "accession": gb_acc, 'pident': float(pident), 'evalue': float(evalue),
                                         'bitscore': float(bitscore), 'sseq': sseq, 'title': stitle}, ignore_index=True)
                else:
                    staxids = int(staxids)
                    if gb_acc not in queried_acc:
                        queried_acc.add(gb_acc)
                        new_blast_seq_dict = new_blast_seq_dict.append(
                            {'ncbi_txn': sscinames, 'ncbi_txid': abs(int(staxids)), 'status': -100,
                             "date": datetime.datetime.strptime('01/01/00', '%d/%m/%y'), "accession": gb_acc,
                             'pident': float(pident), 'evalue': float(evalue),
                             'bitscore': float(bitscore), 'sseq': sseq, 'title': stitle}, ignore_index=True)
    # add data which was not added before and that passes the evalue thresh
    new_blast_seq_dict = wrapper_get_fullseq(config, new_blast_seq_dict, db_name)
    print('new_blast_dict')
    print(len(new_blast_seq_dict))
    print(new_blast_seq_dict[['ncbi_txn', 'ncbi_txid', 'status', 'accession']])
    return new_blast_seq_dict


def get_blastval_from_line(db_name, lin, ncbi_parser, taxon):
    """
    Method used in read_blast_query() that splits the lin of the result file as needed
    :param db_name: name that specifies if filterrun, unpublished or normal search
    :param lin: line in file
    :param ncbi_parser: ncbi parser class to translate names, id's,...
    :param taxon:
    :return:
    """
    sseqid, staxids, sscinames, pident, evalue, bitscore, sseq, salltitles, sallseqid = lin.strip().split('\t')
    # get acc is more difficult now, as new seqs not always have gi number, then query changes
    if db_name == 'filterrun':
        staxids = str(taxon)
        sscinames = ncbi_parser.get_name_from_id(staxids)
        stitle = None
        salltitles = None
        sallseqid = sseqid
    else:
        sscinames = sscinames.replace(" ", "_").replace("/", "_")
        stitle = salltitles

    gb_acc = get_acc_from_blast(sseqid)
    sseq = sseq.replace("-", "")
    pident = float(pident)
    evalue = float(evalue)
    bitscore = float(bitscore)
    return bitscore, evalue, gb_acc, pident, sallseqid, salltitles, sscinames, sseq, staxids, stitle


def split_multiple_tolist(sallseqid, salltitles, sscinames, staxids):
    """
    For redundant data this will split the string with multiple information into a list.

    :param sallseqid:
    :param salltitles:
    :param sscinames:
    :param staxids:
    :return:
    """
    staxids_l = staxids.split(";")
    sscinames_l = sscinames.split(";")
    sallseqid_l = sallseqid.split(";")
    salltitles_l = salltitles.split("<>")
    return sallseqid_l, salltitles_l, sscinames_l, staxids_l
