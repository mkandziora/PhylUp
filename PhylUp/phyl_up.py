"""
PhylUp: automatically update alignments.
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

All rights reserved. No warranty, explicit or implicit, provided. No distribution or modification of code allowed.
All classes and methods will be distributed under an open license in the near future.

Package to automatically update alignments and phylogenies using local sequences or a local Genbank database.

Parts are inspired by the program physcraper developed by me and Emily Jane McTavish.
"""

import sys
import numpy
import datetime
import pandas as pd
import numpy as np
import os
from dendropy import Tree, DnaCharacterMatrix
# from copy import deepcopy
import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser

from . import blast
from . import phylogen_updater
from . import write_msg_logfile
from . import phylogenetic_helpers

#import line_profiler

# TODO: write tests
# todo: write files to tmp... - i did - but then we cannot reuse them in different_level
#todo move some output to tmp folder and delete that


class PhylogeneticUpdater:
    def __init__(self, id_to_spn, aln, aln_schema, tre, tre_schema, config, mrca=None, blacklist=None):
        self.config = config
        self.workdir = self.config.workdir
        self.status = 0  # -1 = deleted, positive values = round, 0 = present at beginning
        self.aln_fn = aln
        self.aln_schema = aln_schema
        self.aln = DnaCharacterMatrix.get(path=self.aln_fn, schema=self.aln_schema)
        self.tre_fn = tre
        self.tr_schema = tre_schema
        self.tre = None
        self.blacklist = blacklist
        self.table = phylogenetic_helpers.build_table_from_file(id_to_spn, self.config, self.config.downtorank)
        phylogenetic_helpers.add_seq_to_table(self.aln, self.table)
        self.mrca = None
        self.set_mrca(mrca)

    def set_back_data_in_table(self):
        """
        Overwrites information in table for all blast date and status of newly added sequences.

        Needed for different_level == True

        :return:
        """
        if self.config.different_level == True and self.status == 0:
            try:
                self.table['date'] = pd.to_datetime('01/01/00', format='%d/%m/%y')
            except NotImplementedError:
                self.table['date'] = pd.Timestamp.strptime('01/01/00', '%d/%m/%y')
            self.table.at[self.table['status'] > 0, 'status'] = 0.5
            self.config.update_tree = False

    def set_mrca(self, mrca):
        """
        Set the mrca in the class for various input: None, one ncbi id, list of ncbi_ids.

        Note: list finds mrca of all ids, including ids that were not mentioned.
        Thus, not helpful for non-monophyletic groups!

        :param mrca: input mrca
        :return: sets self.mrca
        """
        if type(mrca) is int:
            self.mrca = mrca
        else:
            if type(mrca) is list:
                # print('get mrca of all input!')
                self.mrca = mrca
            else:
                print('get mrca of no input!')
                aln_taxids = set(self.table['ncbi_txid'].tolist())
                ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_parser_names_fn,
                                                      nodes_file=self.config.ncbi_parser_nodes_fn)
                self.mrca = ncbi_parser.get_mrca(taxid_set=aln_taxids)

    def extend_with_unpublished(self):
        """
        This will use a local database to search for new sequences to add.
        Mend for sequences that are not yet published.

        :return: new seqs from unpublished
        """
        sys.stdout.write('Find new sequences using a local database.\n')
        self.status += 1
        present_subset = self.table[self.table['status'] > -1]
        new_seqs_unpubl = pd.DataFrame(columns=['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date", 'accession',
                                                'pident', 'evalue', 'bitscore', 'sseq', 'title'])
        # convert name to file
        ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_parser_names_fn,
                                              nodes_file=self.config.ncbi_parser_nodes_fn)
        name_txid_unpublished = phylogenetic_helpers.get_txid_for_name_from_file(self.config.unpubl_names, ncbi_parser)
        if not os.path.exists(os.path.join(self.workdir, 'tmp')):
            os.mkdir(os.path.join(self.workdir, 'tmp'))
        name_txid_unpublished[['accession', 'ncbi_txid']].to_csv(os.path.join(self.workdir, 'tmp/unpublished_txid.txt'),
                                                                 sep=' ', header=False, index=False)
        blast.make_local_blastdb(self.config.workdir, db='unpublished',
                                 taxid=os.path.join(self.workdir, 'tmp/unpublished_txid.txt'),
                                 path_to_db=self.config.unpubl_data)
        # get new seqs from seqs in blast table seq
        for index in present_subset.index:
            tip_name = self.table.loc[index, 'accession']
            query_seq = self.table.loc[index, 'sseq']
            # Note: should be empty in later rounds...that is why it does not matter that new_seq_tax is reassigned
            new_seq_tax = blast.get_new_seqs(query_seq, tip_name, self.config.blastdb, "unpublished", self.config)
            print(new_seq_tax)
            # new_seqs_unpubl = new_seqs_unpubl.append(new_seq_tax, ignore_index=True)
            new_seqs_unpubl = pd.concat([new_seqs_unpubl, new_seq_tax], ignore_index=True, sort=True)
        new_seqs_unpubl = new_seqs_unpubl.drop(['accession;gi', 'ncbi_txid', 'ncbi_txn'], axis=1)
        new_seqs_unpubl = pd.merge(new_seqs_unpubl, name_txid_unpublished, on=['accession'], sort=True)
        new_seqs_unpubl['title'] = 'unpublished'
        new_seqs_unpubl['status'] = 0.5
        return new_seqs_unpubl

    def extend(self):
        """
        Gets new sequences from those existing once which have not yet been blasted before.

        :return:
        """
        # create list of indice of subset list
        sys.stdout.write("Find new sequences using the BLAST database.\n")
        # sets back date for rerun with different tribe or so
        self.set_back_data_in_table()
        self.status += 1
        #today = pd.Timestamp.today()
        #min_date_blast = today - pd.Timedelta(days=90)
        present_subset = self.table[self.table['status'] > -1]
        #present_subset = present_subset[(present_subset.date <= min_date_blast)]
        new_seqs_local = pd.DataFrame(columns=['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date", 'accession',
                                               'pident', 'evalue', 'bitscore', 'sseq', 'title'])
        # get new seqs from seqs in blast table seq
        msg = present_subset[['accession', 'ncbi_txn',  'status']].to_string()
        write_msg_logfile(msg, self.config.workdir)
        for index in present_subset.index:
            tip_name = self.table.loc[index, 'accession']
            # msg = tip_name
            # write_msg_logfile(msg, self.config.workdir)
            print('Blast: {}.'.format(tip_name))
            query_seq = self.table.loc[index, 'sseq']
            self.table.at[index, 'date'] = pd.Timestamp.today()  # this is the new version of pd.set_value(), sometimes it's iat
            new_seq_tax = blast.get_new_seqs(query_seq, tip_name, self.config.blastdb, "Genbank", self.config)
            # new_seqs_local = new_seqs_local.append(new_seq_tax, ignore_index=True)
            new_seqs_local = pd.concat([new_seqs_local, new_seq_tax], ignore_index=True, sort=True)
        assert len(new_seqs_local) > 0, new_seqs_local
        if self.blacklist is not None:
            assert type(self.blacklist) == list, (type(self.blacklist), self.blacklist)
            new_seqs_local = self.remove_blacklist_items(new_seqs_local)
        return new_seqs_local

    def remove_blacklist_items(self, new_seqs):
        """
        Removes accession numbers specified in the analysis file that shall not be added.

        :param new_seqs:  pandas dataframe with the new sequences retrieved earlier
        :return:
        """
        drop_boolean = np.where((new_seqs.accession.isin(self.blacklist)), True, False)
        new_seqs = new_seqs[drop_boolean != True]
        return new_seqs

    def add_new_seqs(self, new_seqs):
        """
        Adds new sequences to table, after removing those which were already present and returns indices.

        :param new_seqs:  pandas dataframe with the new sequences retrieved earlier
        :return: new_seqs without seqs that were added earlier
        """
        ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_parser_names_fn,
                                              nodes_file=self.config.ncbi_parser_nodes_fn)
        for index in new_seqs.index:
            taxid = new_seqs.loc[index, 'ncbi_txid']
            name = ncbi_parser.get_name_from_id(taxid)
            new_seqs.at[index, 'ncbi_txn'] = name

        new_seqs.ncbi_txn = new_seqs.ncbi_txn.str.replace(' ', '_')  # replace whitespace in name
        new_seqs.ncbi_txn = new_seqs.ncbi_txn.str.replace('.', '')

        subcols = new_seqs[['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date", 'sseq', 'accession']]
        subcols.at[:, 'status'] = self.status
        #self.table = self.table.append(subcols, ignore_index=True)
        self.table = pd.concat([self.table, subcols], ignore_index=True, sort=True)
        #assert self.table['status'].hasnans == False, self.table['status'].hasnans
        assert len(subcols) == len(self.table[self.table['status'] == self.status]), \
            (len(subcols), len(self.table[self.table['status'] == self.status]), subcols['accession'])
        new_seqs_upd_index = self.table[self.table['status'] == self.status]
        return new_seqs_upd_index

    def call_filter(self, new_seqs, aln):
        """
        Wrapper function where the different filters are being called

        :param new_seqs: pandas dataframe with the new sequences retrieved earlier
        :param aln: alignment as dendropy object
        :return:
        """
        msg = "Round of filters: {}\n".format(self.status)
        write_msg_logfile(msg, self.config.workdir)
        # everyone's filter
        new_seqs = new_seqs[~new_seqs['accession'].isin(self.table['accession'])]  # ~ is the pd not in/!
        print('Length of new seqs before filtering: {}'.format(len(new_seqs)))
        if len(new_seqs) > 0:
            print(new_seqs.columns)
            new_seqs = self.basic_filters(aln, self.mrca, new_seqs)
            # next filter need infos in table
            new_seqs = self.add_new_seqs(new_seqs)

        for idx in new_seqs.index:
            assert idx in self.table.index, (idx, self.table.index)
        # filter for seq identity process
        new_seqs = self.compare_filter(new_seqs)
        return new_seqs

    def compare_filter(self, new_seqs):
        """
        Function calls filters, that alter entries in the table.

        :param new_seqs: pandas dataframe with the new sequences retrieved/not filtered from earlier
        :return:
        """
        if len(new_seqs) > 0:
            for idx in new_seqs.index:
                assert idx in self.table.index, (idx, self.table.index)
            f = FilterSeqIdent(self.config, self.table, self.status)
            msg = "Time before Filter {}: {}.\n".format(f, datetime.datetime.now())
            write_msg_logfile(msg, self.config.workdir)
            before_filter = new_seqs
            f.filter(new_seqs)
            self.aln = f.aln
            new_seqs = f.upd_new_seqs  # assign new_seqs to last upd from filter before
            self.table = f.table
            del_seq = f.del_table

            # check that all del_tab entries are deleted in table
            assert any(
                x in del_seq['accession'].to_list() for x in self.table.loc[self.table['status'] >= 0, 'accession'].to_list()) == False, \
                ([i for i in del_seq['accession'].to_list() if i in self.table.loc[self.table['status'] >= 0, 'accession'].to_list()])
            check_df_index_unique(new_seqs)
            check_df_index_unique(self.table)
            check_filter_numbers(del_seq, new_seqs, before_filter)
            assert len(new_seqs) == len(self.table[self.table['status'] == self.status]), \
                (len(new_seqs), len(self.table[self.table['status'] == self.status]),
                 new_seqs['accession'], self.table[self.table['status'] == self.status]['accession'])

            # filter for number otu
            f = FilterNumberOtu(self.config, self.table, self.status)
            msg = "Time before Filter {}: {}.\n".format(f, datetime.datetime.now())
            write_msg_logfile(msg, self.config.workdir)
            f.filter(new_seqs, self.config.downtorank)
            self.table = f.table
            new_seqs = f.upd_new_seqs
        return new_seqs

    def basic_filters(self, aln, mrca, new_seqs):
        """
        Calls the basic filters, those which are independent of that had been added earlier.
        Sequences filtered here are not being added to table.

        :param aln: aln as dendropy object
        :param mrca: mrca of ingroup
        :param new_seqs: pandas dataframe with the new sequences retrieved earlier
        :return: subset of new seqs dataframe
        """
        orig_len = new_seqs
        columns = ['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date",
                   'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title']
        all_del = pd.DataFrame(columns=columns)
        remove_basics = [FilterUniqueAcc(self.config, self.table),
                         FilterBLASTThreshold(self.config),
                         FilterLength(self.config, aln),
                         FilterMRCA(self.config, mrca)
                         ]
        for f in remove_basics:
            msg = "Time before Filter {}: {}.\n".format(f, datetime.datetime.now())
            write_msg_logfile(msg, self.config.workdir)
            internal_new_seq = new_seqs
            f.filter(new_seqs)
            new_seqs = f.upd_new_seqs
            del_seq = f.del_table
            # all_del = all_del.append(del_seq, ignore_index=True)
            all_del = pd.concat([all_del, del_seq], ignore_index=True, sort=True)
            check_filter_numbers(del_seq, new_seqs, internal_new_seq)
        check_filter_numbers(all_del, new_seqs, orig_len)
        return new_seqs

    def run(self, status_end=None):
        """
        Main function of the class - Finds new sequences and adds them to the aln and tre.

        :return:
        """

        msg = "Begin update: {}.\n".format(datetime.datetime.now())
        write_msg_logfile(msg, self.config.workdir)

        columns = ['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date",
                   'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title']
        all_new_seqs = pd.DataFrame(columns=columns)

        self.call_input_cleaner()

        retrieved_seqs = 1
        if os.path.exists(os.path.join(self.workdir, 'all_new_seqs.updated')):
            all_new_seqs = pd.read_csv(os.path.join(self.workdir, 'all_new_seqs.updated'))
            self.table = pd.read_csv(os.path.join(self.workdir, 'table.updated'))
            self.aln = DnaCharacterMatrix.get(path=os.path.abspath(os.path.join(self.workdir, 'updt_aln.fasta')),
                                              schema='fasta')

            self.tre = Tree.get(path=os.path.abspath(os.path.join(self.workdir, 'updt_tre.tre')), schema='newick',
                                taxon_namespace=self.aln.taxon_namespace, preserve_underscores=True)
            msg = "New round of updating begins with mrca: {}.\n".format(self.mrca)
            write_msg_logfile(msg, self.config.workdir)
        if self.config.different_level == False and os.path.exists(os.path.join(self.workdir, 'all_new_seqs.updated')):
            next
        else:
            while retrieved_seqs > 0 and (status_end is None or self.status <= status_end):
                print(status_end is None or self.status <= status_end)
                msg = "Time before blast: {}.\n".format(datetime.datetime.now())
                write_msg_logfile(msg, self.config.workdir)
                if self.config.unpublished is True:
                    new_seqs = self.extend_with_unpublished()
                    self.config.unpublished = False
                    # print(new_seqs[['accession', 'ncbi_txid', 'ncbi_txn']])
                else:
                    new_seqs = self.extend()  # todo rename to find new seqs
                    msg = "Time after BLAST: {}.\n".format(datetime.datetime.now())
                    write_msg_logfile(msg, self.config.workdir)
                new_seqs = self.call_filter(new_seqs, self.aln)
                new_seqs = new_seqs[~new_seqs['accession'].isin(all_new_seqs['accession'])]  # ~ is the pd not in/!
                all_new_seqs = pd.concat([all_new_seqs, new_seqs], ignore_index=True, sort=True)
                print('Length of new seqs after filtering: {}'.format(len(new_seqs)))
                retrieved_seqs = len(new_seqs)
                msg = "Newly found seqs: {}.\n".format(len(new_seqs))
                write_msg_logfile(msg, self.config.workdir)
            msg = "Time before update tre/aln: {}.\n".format(datetime.datetime.now())
            write_msg_logfile(msg, self.config.workdir)
            msg = "Newly found seqs will be added to aln and tre: {}.\n".format(len(all_new_seqs))
            write_msg_logfile(msg, self.config.workdir)

            all_new_seqs.to_csv(os.path.join(self.workdir, 'all_new_seqs.updated'), index=False)
            if self.config.different_level:
                if os.path.exists(os.path.join(self.workdir, 'table.updated')):
                    os.rename(os.path.join(self.workdir, 'table.updated'),
                              os.path.join(self.workdir, "table_updated_tmp"))
            self.table.to_csv(os.path.join(self.workdir, 'table.updated'), index=False)
        self.update_aln()
        if self.tre is None:
            self.tre_fn = os.path.abspath(os.path.join(self.config.workdir, "updt_aln.fasta.tree"))
            self.tre = Tree.get(path=os.path.join(self.config.workdir, 'updt_tre.tre'),
                                schema="newick",
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
        self.update_tre()
        msg = "Time finished: {}.\n".format(datetime.datetime.now())
        write_msg_logfile(msg, self.config.workdir)

    def call_input_cleaner(self):
        """
        Calls InputCleaner class and updates input data.
        :return:
        """
        cleaner = phylogen_updater.InputCleaner(self.tre_fn, self.aln_fn, self.table, self.config, self.mrca)
        self.aln = cleaner.aln
        if self.tre_fn is not None:
            self.tre = cleaner.tre
        self.table = cleaner.table
        self.mrca = cleaner.mrca

    def update_aln(self):
        """
        Call the alignment updater class.

        :return:
        """
        # print(os.path.abspath(os.path.join(self.workdir, 'updt_aln.fasta')))
        aln = phylogenetic_helpers.read_in_aln(os.path.abspath(os.path.join(self.workdir, 'updt_aln.fasta')), self.aln_schema)
        if self.tre_fn is None:
            tre = None
        else:
            tre = Tree.get(path=os.path.abspath(os.path.join(self.workdir, 'updt_tre.tre')), schema='newick',
                           taxon_namespace=aln.taxon_namespace, preserve_underscores=True)
        aln_upd = phylogen_updater.AlnUpdater(self.config, aln, self.status, self.table, tre)
        self.aln = aln_upd.aln
        self.tre = aln_upd.tre
        msg = 'Updating of aln done.\n'
        write_msg_logfile(msg, self.config.workdir)

    def update_tre(self):
        """
        Call the tree updater class.

        :return:
        """
        phylogen_updater.TreeUpdater(self.config, self.aln, self.table, self.tre)
        msg = 'Updating of aln and tre done.\n'
        write_msg_logfile(msg, self.config.workdir)


# ################################################################################################################
class Filter(object):
    """
    Super class for the different Filter classes.
    """
    def __init__(self, config):
        self.config = config
        columns = ['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date",
                   'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title']
        self.upd_new_seqs = pd.DataFrame(columns=columns)
        self.del_table = pd.DataFrame(columns=columns)
        self.ncbi_parser = None

    def assert_after_filter(self, new_seqs):
        """
        Different assert statements to ensure everything works correctly.

        :param new_seqs: pandas new seqs
        :return:
        """
        assert all(x in new_seqs['accession'].tolist() for x in self.upd_new_seqs['accession'].tolist())

        assert all(x in new_seqs['accession'].tolist() for x in self.del_table['accession'].tolist())

        assert any(x in self.upd_new_seqs['accession'].tolist() for x in self.del_table['accession'].tolist()) == False, \
            ([i for i in self.upd_new_seqs['accession'].tolist() if i in self.del_table['accession'].tolist()])

        assert self.upd_new_seqs.index.is_unique is True, self.upd_new_seqs.index
        assert self.del_table.index.is_unique is True, self.upd_new_seqs.index

        assert len(self.del_table) + len(self.upd_new_seqs) == len(new_seqs), \
            (len(self.del_table), len(self.upd_new_seqs), len(new_seqs), [x in new_seqs['accession'].tolist() for x
                                                                          in self.del_table['accession'].tolist()],
             [x in new_seqs['accession'].tolist() for x in self.upd_new_seqs['accession'].tolist()])
        assert self.upd_new_seqs.index.is_unique is True, self.upd_new_seqs.index


def assert_new_seqs_table(new_seqs, table, status):
    """
    Different assert statements to ensure everything works correctly.

    :param new_seqs: pandas new seqs
    :param table: pandas table with all information
    :param status: round of updating
    :return:
    """
    assert len(new_seqs) == len(table[table['status'] == status]), \
        (len(new_seqs), len(table[table['status'] == status]), new_seqs['accession'])
    assert new_seqs.index.is_unique is True, new_seqs.index
    assert new_seqs['ncbi_txid'].hasnans == False, new_seqs[['ncbi_txid', 'accession']].hasnans
    for idx in new_seqs.index:
        assert idx in table.index, (idx, table.index)


def calculate_mean_sd(bitscores):
    """
    Calculates standard deviation and mean of scores which are used as a measure of sequence differentiation
    for a given taxon, is being used to select a random representative of a taxon later.

    :param bitscores:
    """
    len_bits = len(bitscores)
    sum_bits = bitscores.sum()
    bit_sd = float(numpy.std(bitscores))
    mean_hsp_bits = float(sum_bits / len_bits)
    mean_sd = {"mean": mean_hsp_bits, "sd": bit_sd}
    return mean_sd


class FilterNumberOtu(Filter):
    """
    Filter new sequences to the number defined as threshold. Either via local blast or select longest
    """
    def __init__(self, config, table, status):
        super().__init__(config)
        self.table = table
        assert self.table.index.is_unique is True, self.table.index
        self.status = status
        assert len(self.upd_new_seqs) == 0, self.upd_new_seqs


    def initialize(self, config):
        self.ncbi_parser = ncbi_data_parser.Parser(names_file=config.ncbi_parser_names_fn,
                                                   nodes_file=config.ncbi_parser_nodes_fn)

    def filter(self, new_seqs, downtorank=None):
        assert_new_seqs_table(new_seqs, self.table, self.status)
        print("filter FilterNumberOtu")
        # get all seqs and ids from aln and before for next filter
        present = self.table[self.table['status'] > -1]
        added_before = present[present['status'].astype(int) < self.status]

        if downtorank is None:
            tax_ids_newseqs = new_seqs['ncbi_txid']
            txids_added = added_before['ncbi_txid']
        else:
            new_seqs = self.set_downtorank(new_seqs, downtorank)  # todo: this is likely doubled now, since i chnages now ncbi id to rank and original_ncbi_id for the real one...
            new_seqs_downto = self.set_downtorank(added_before, downtorank)
            tax_ids_newseqs = new_seqs['downtorank']
            txids_added = new_seqs_downto['downtorank']

        for txid in set(tax_ids_newseqs):
            # generate taxon_subsets
            # print('filter per otu')
            ns_txid_df = new_seqs[tax_ids_newseqs == txid]
            os_txid_df = added_before[txids_added == txid]
            os_txid_df = os_txid_df[os_txid_df['status'] >= 0]
            if downtorank is not None and downtorank not in ['species', 'subspecies', 'variety']:
                # figure out which lower rank are already available
                avail_txid = os_txid_df['ncbi_txid'].tolist()
                # filter ns_txid_df to remove lower ranks already available
                ns_txid_df = ns_txid_df[~ns_txid_df['ncbi_txid'].isin(avail_txid)]
            # if we still want to add loop through
            if len(os_txid_df) < self.config.threshold:
                if len(ns_txid_df) + len(os_txid_df) > self.config.threshold:  # filter
                    if len(os_txid_df) == 0:  # new taxa, select random seq for blast
                        # print('taxa new')
                        if self.config.filtertype == 'blast':
                            filtered_acc = self.wrapper_filter_blast_otu(ns_txid_df, os_txid_df, 'accession', txid)  # returns only add column
                            filtered = new_seqs[new_seqs['accession'].isin(filtered_acc)]
                            assert len(filtered_acc) == len(filtered), (len(filtered_acc), len(filtered),
                                                                        [i for i in filtered_acc if i in filtered['accession'].to_list()])
                        elif self.config.filtertype == 'length':
                            print('filter otu by length')
                            filtered = self.select_seq_by_length(ns_txid_df, os_txid_df)
                        else:
                            # print('should not happen')
                            sys.exit(2)
                    else:  # select seq from aln
                        # print('taxa present already')
                        if self.config.filtertype == 'blast':
                            filtered_acc = self.wrapper_filter_blast_otu(ns_txid_df, os_txid_df, 'tip_name', txid)
                            filtered = new_seqs[new_seqs['accession'].isin(filtered_acc)]
                            assert len(filtered_acc) == len(filtered), (len(filtered_acc), len(filtered))
                        elif self.config.filtertype == 'length':
                            print('filter otu by length')
                            filtered = self.select_seq_by_length(ns_txid_df, os_txid_df)
                        else:
                            # print('should not happen')
                            sys.exit(2)
                elif len(ns_txid_df) + len(os_txid_df) <= self.config.threshold:  # filter
                    # print('add all')
                    filtered = ns_txid_df
                self.upd_new_seqs = pd.concat([self.upd_new_seqs, filtered], ignore_index=True, sort=True)  # if this is further below, old filtered entry will be added
            elif len(os_txid_df) > self.config.threshold:
                if not self.config.downtorank:
                    # print('sample size to big')
                    sys.exit(-3)
            else:
                # print('sample size correct - nothing to add')
                filtered = pd.DataFrame()
            check_df_index_unique(self.upd_new_seqs)
        not_selected = list(set(new_seqs['accession'].values) - set(self.upd_new_seqs['accession'].values))
        del_tab = new_seqs[new_seqs['accession'].isin(not_selected)]
        self.del_table = del_tab
        if len(not_selected) > 0:
            self.table.at[self.table['accession'].isin(not_selected), 'status'] = -1
            self.table.at[self.table['accession'].isin(not_selected), 'status_note'] = 'too many seqs of same tax_id'
        check_filter_numbers(not_selected, self.upd_new_seqs, new_seqs)

        msg = "Filter FilterNumberOtu reduced the new seqs from {} to {}.\n".format(len(new_seqs),
                                                                                    len(self.upd_new_seqs))
        write_msg_logfile(msg, self.config.workdir)

    def set_downtorank(self, new_seqs, downtorank):
        """
        Method is used to get the corresponding rank if higher than lowest otu by ncbi.

        :param new_seqs:
        :param downtorank:
        :return:
        """
        # print("set_downtorank")
        new_seqs.at[:, 'downtorank'] = 0
        if self.ncbi_parser is None:
            self.initialize(self.config)
        for txid in set(new_seqs['ncbi_txid']):
            downtorank_id = self.ncbi_parser.get_downtorank_id(txid, downtorank)
            new_seqs.at[new_seqs.ncbi_txid == txid, 'downtorank'] = downtorank_id
        return new_seqs


    def select_seq_by_length(self, filter_dict, exist_dict):
        """
        This is another mode to filter the sequences, if there are more than the threshold amount available.
        This one selects new sequences by length instead of by score values. It is selected by "selectby='length'".

        :param filter_dict:
        :param exist_dict:
        :return: filtered sequences
        """
        # print("select_seq_by_length")
        len_seqs_new = filter_dict['sseq'].apply(len)
        amnt_missing_seq = self.config.threshold - len(exist_dict)
        if len(len_seqs_new) > amnt_missing_seq:
            select = pd.DataFrame(columns=['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date", 'accession',
                                           'pident', 'evalue', 'bitscore', 'sseq', 'title'])
            for i in range(1, amnt_missing_seq+1):
                idxmax_val = len_seqs_new.idxmax()
                select = select.append(filter_dict.loc[idxmax_val])
                filter_dict = filter_dict.drop([idxmax_val])
                len_seqs_new = len_seqs_new.drop([idxmax_val])
            assert len(select) == amnt_missing_seq, (len(select), amnt_missing_seq)
        else:
            print('DO WE EVER GET HERE? length filter')
            select = len_seqs_new
            assert len(select) <= amnt_missing_seq
        return select


    def wrapper_filter_blast_otu(self, filter_dict, exist_dict, columnname, txid):
        """
        get subsample of new seqs for otu_sample by selecting random new seqs of taxa that are completely new or for
        random existing seq.

        :param filter_dict: table with all new seqs
        :param exist_dict: table with all available/added seqs
        :param columnname: depending on columnname, filter_dict is being altered to remove the one seqs that is used
                            for the blasting
        :param txid: ncbi taxon id
        :return:
        """
        # print('wrapper_filter_blast_otu')
        rndm = filter_dict.sample(1)
        idx = rndm.index.values.tolist()[0]
        query_seq = rndm.loc[idx, 'sseq']
        query_seq_id = rndm.loc[idx, 'accession']
        blast.write_local_blast_files(self.config.workdir, query_seq_id, query_seq, db=False, fn=txid)
        # write local db and create it
        if columnname == 'tip_name':
            for_db = filter_dict
        else:
            for_db = filter_dict.drop([idx])
        # filter_seq_db has to be deleted for every new filterun
        if os.path.exists(os.path.join(self.config.workdir, 'tmp/filter_seq_db')):
            os.remove(os.path.join(self.config.workdir, 'tmp/filter_seq_db'))
        for idx in for_db.index:
            seq_id = for_db.loc[idx, 'accession']
            seq = for_db.loc[idx, 'sseq']
            # print('write local blast files')
            blast.write_local_blast_files(self.config.workdir, seq_id, seq, db=True)
        blast.make_local_blastdb(self.config.workdir, db='filterrun', taxid=txid)
        filter_blast_seqs = blast.get_new_seqs(query_seq, txid, os.path.join(self.config.workdir, 'tmp'),
                                               'filterrun', self.config)
        subfilter_blast_seqs = remove_mean_sd_nonfit(filter_blast_seqs)
        filtered_seqs = self.select_number_missing(subfilter_blast_seqs, exist_dict)
        filtered_acc = set(list(filtered_seqs['accession']))
        return filtered_acc


    def select_number_missing(self, subfilter_blast_seqs, table_otu_dict):
        """
        Selects random sequences according to how many are missing.

        :param subfilter_blast_seqs:
        :param table_otu_dict:
        :return:
        """
        amnt_missing_seq = self.config.threshold - len(table_otu_dict)
        if amnt_missing_seq < len(subfilter_blast_seqs):
            # print('subfilter-sample')
            filtered_seqs = subfilter_blast_seqs.sample(amnt_missing_seq)
        elif amnt_missing_seq >= len(subfilter_blast_seqs):
            # print('subfilter-all')
            filtered_seqs = subfilter_blast_seqs
        else:
            print('should not happen')
            sys.exit(-2)
        return filtered_seqs


def remove_mean_sd_nonfit(filter_blast_seqs):
    """
    Remove blast results that are outside of the bitscore range.

    :param filter_blast_seqs:
    :return:
    """
    # remove non fitting entries - too far from blast query seq....
    mean_sd = calculate_mean_sd(filter_blast_seqs['bitscore'])
    subfilter_blast_seqs = pd.DataFrame()
    for idx in filter_blast_seqs.index:
        seq_bitscore = filter_blast_seqs.loc[idx, 'bitscore']
        if (seq_bitscore >= mean_sd['mean'] - mean_sd['sd']) & \
                (seq_bitscore <= mean_sd['mean'] + mean_sd['sd']):
            subfilter_blast_seqs = subfilter_blast_seqs.append(filter_blast_seqs.loc[idx])
        # else:
        #     print('filter thresh too large')
    return subfilter_blast_seqs

#todo table probably not needed anymore
class FilterSeqIdent(Filter):
    """
    Filters new sequences that are identical (seq and taxon id) to something already in the data.
    """
    def __init__(self, config, table, status):
        super().__init__(config)
        self.status = status
        self.table = table
        assert self.table.index.is_unique is True, self.table.index
        assert self.table['ncbi_txid'].hasnans == False, self.table['ncbi_txid'].hasnans
        assert self.table['sseq'].hasnans == False, self.table['sseq']


    def filter(self, new_seqs):
        assert_new_seqs_table(new_seqs, self.table, self.status)
        print("filter FilterSeqIdent rm duplicate")
        #self.upd_new_seqs = deepcopy(new_seqs)
        self.upd_new_seqs = new_seqs

        # drop all sequences that are identical (seq and ncbi_txid) from new_seqs
        dupl_df = self.upd_new_seqs[self.upd_new_seqs.duplicated(['ncbi_txid', 'sseq'], keep='first')]
        for ind in dupl_df.index:
            to_del = new_seqs.loc[ind]
            self.del_table = self.del_table.append(to_del)
            self.upd_new_seqs = self.upd_new_seqs.drop([ind])

        for idx in set(self.upd_new_seqs.index):
            # generate seq and compare it
            txid_compare = new_seqs.loc[idx, 'ncbi_txid']
            seq_compare = new_seqs.loc[idx, 'sseq']

            # go through new_seqs to check if there is more than one identical (1 is the one we compare here)
            same_seqs = self.upd_new_seqs[(self.upd_new_seqs.sseq.str.contains(seq_compare)) & (self.upd_new_seqs['ncbi_txid'] == txid_compare)]
            if len(same_seqs.index) > 1:
                print("identical seq new")
                to_del = new_seqs.loc[idx]
                print(new_seqs.loc[idx, "accession"])

                self.del_table = self.del_table.append(to_del)
                self.upd_new_seqs = self.upd_new_seqs.drop([idx])
            else:
                same_table = self.table[self.table['status'].between(-1, self.status, inclusive=False)]
                same_old = same_table[(same_table.sseq.str.contains(seq_compare)) & (same_table['ncbi_txid'] == int(txid_compare))]

                # go through avail seq in table to check if there is one identical
                if len(same_old.index) > 0:
                    print("identical old")
                    print(new_seqs.loc[idx, "accession"])

                    to_del = new_seqs.loc[idx]
                    self.del_table = self.del_table.append(to_del)
                    self.upd_new_seqs = self.upd_new_seqs.drop([idx])
                    assert idx not in self.upd_new_seqs.index
        for index in self.del_table.index:
            self.table.at[index, 'status'] = -1
            self.table.loc[index, 'status_note'] = 'same seq ident'
            # assert self.table.loc[index, 'status'] == -1
        self.assert_after_filter(new_seqs)
        assert self.table['sseq'].hasnans == False, self.table['sseq']
        msg = "Filter FilterSeqIdent reduced the new seqs from {} to {}.\n".format(len(new_seqs),
                                                                                   len(self.upd_new_seqs))
        write_msg_logfile(msg, self.config.workdir)
        self.aln = self.remove_existing_for_longer()

    def remove_existing_for_longer(self):
        print("remove_existing_for_longer")
        # check if one of the new seqs is identical but longer than old
        existing_old = self.table[self.table['status'].between(-1, self.status, inclusive=False)]
        aln = DnaCharacterMatrix.get(path=os.path.abspath(os.path.join(self.config.workdir, "updt_aln.fasta")), schema='fasta')
        counter = 0
        for idx in existing_old.index:
            # generate seq and compare it
            txid_compare = existing_old.loc[idx, 'ncbi_txid']
            seq_compare = existing_old.loc[idx, 'sseq']
            # use those to compare to new seqs:
            same_new = self.upd_new_seqs[(self.upd_new_seqs.sseq.str.contains(seq_compare)) & (
                            self.upd_new_seqs['ncbi_txid'] == int(txid_compare))]
            # go through avail seq in table to check if there is one identical
            if len(same_new.index) > 0:
                print("old seq found in new")
                self.table.at[idx, 'status'] = -1
                self.table.at[idx, 'status_note'] = 'new seq found - identical but longer'
                tip_name = self.table.loc[idx, "accession"]
                for tax, seq in aln.items():
                    if tip_name == tax.label:
                        print("remove from aln")
                        aln.remove_sequences([tax])
                aln.write(path="updt_aln.fasta", schema="fasta")
                counter += 1

                # assert that seq is removed:
                count = 0
                for tax, seq in aln.items():
                    seq = seq.symbols_as_string().replace('-', '')
                    if seq == seq_compare:
                        print('found seq')
                        count = 1
                assert count == 0
        msg = "Filter FilterSeqIdent removed {} sequence(s) from existing alignment due to longer" \
              " but identical sequence in new seqs.\n".format(counter)
        write_msg_logfile(msg, self.config.workdir)
        return aln


def get_taxid_using_table(idx_ident_new_seqs, new_seqs):
    """
    Get the taxon ids using the index and the corresponding table.
    Used for entrys which have the same sequence as the one we compare in main function.

    :param idx_ident_new_seqs: index of sequence we want to get the taxon id from.
    :param new_seqs: table with new sequences
    :return:
    """
    try:
        txid_new_seqident = set(new_seqs.loc[idx_ident_new_seqs, 'ncbi_txid'].to_list())
    except AttributeError:  # series vs something else issue
        txid_new_seqident = set(new_seqs.loc[idx_ident_new_seqs, 'ncbi_txid'].tolist())
    return txid_new_seqident


class FilterMRCA(Filter):
    """
    Filters sequences that are not part of the defined mrca.
    """
    def __init__(self, config, mrca):
        super().__init__(config)
        self.mrca = mrca

    def initialize(self, config):
        self.ncbi_parser = ncbi_data_parser.Parser(names_file=config.ncbi_parser_names_fn,
                                                   nodes_file=config.ncbi_parser_nodes_fn)


    def filter(self, new_seqs):
        print("FilterMRCA")
        if self.ncbi_parser is None:
            self.initialize(self.config)
        to_del = pd.DataFrame()
        for tax_id in set(new_seqs['ncbi_txid'].astype(int)):
            tax_id = int(tax_id)
            mrca_tx = self.ncbi_parser.match_id_to_mrca(tax_id, self.mrca)
            select_tf = new_seqs['ncbi_txid'].astype(int) == tax_id
            if type(self.mrca) == int:
                # TODO: make self.mrca to be a set - single id is type int
                print("DO I EVER GET HERE - MRCA IS INT")
                if mrca_tx == self.mrca:
                    to_add = new_seqs[select_tf]
                    assert len(to_add) != 0, len(to_add)
                    self.upd_new_seqs = pd.concat([self.upd_new_seqs, to_add], axis=0, ignore_index=True, sort=True)
                else:
                    # print('out of mrca')
                    to_del = new_seqs[select_tf]
                    to_del.at[:, 'status'] = 'deleted - mrca'
                    self.del_table = pd.concat([self.del_table, to_del], axis=0, ignore_index=True, sort=True)
            elif type(self.mrca) == set:
                # print(mrca_tx in self.mrca, mrca_tx, self.mrca)
                if mrca_tx in self.mrca:
                    to_add = new_seqs[select_tf]
                    assert len(to_add) != 0, len(to_add)
                    self.upd_new_seqs = pd.concat([self.upd_new_seqs, to_add], axis=0, ignore_index=True, sort=True)
                else:
                    # print('out of mrca')
                    to_del = new_seqs[select_tf]
                    to_del.at[:, 'status'] = 'deleted - mrca'
                    self.del_table = pd.concat([self.del_table, to_del], axis=0, ignore_index=True, sort=True)
            else:
                sys.stderr.write('MRCA FILTER DOES NOT BEHAVE AS EXPECTED!')
                sys.exit(-1000)

        check_df_index_unique(self.upd_new_seqs)
        check_filter_numbers(self.del_table, self.upd_new_seqs, new_seqs)

        msg = "Filter FilterMRCA reduced the new seqs from {} to {}.\n".format(len(new_seqs), len(self.upd_new_seqs))
        write_msg_logfile(msg, self.config.workdir)
        if len(to_del) > 0:
            to_del.to_csv(os.path.join(self.config.workdir, 'wrong_mrca.csv'), mode='a')


class FilterBLASTThreshold(Filter):
    """
    Removes sequences that do not pass the similarity filter (blast threshold).
    """
    def __init__(self, config):
        super().__init__(config)


    def filter(self, new_seqs):
        print("FilterBLASTThreshold")
        tf_eval = new_seqs['evalue'].astype(float) < float(self.config.e_value_thresh)
        upd_new_seqs = new_seqs[tf_eval == True]
        deltab = new_seqs[tf_eval != True]
        deltab['status_note'] = 'deleted - evalue'

        self.del_table = deltab
        self.upd_new_seqs = upd_new_seqs

        check_df_index_unique(self.upd_new_seqs)
        check_filter_numbers(deltab, upd_new_seqs, new_seqs)

        msg = "Filter FilterBLASTThreshold reduced the new seqs from {} to {}.\n".format(len(new_seqs),
                                                                                         len(upd_new_seqs))
        write_msg_logfile(msg, self.config.workdir)
        if len(deltab) > 0:
            deltab.to_csv(os.path.join(self.config.workdir, 'low_blast_threshold.csv'), mode='a')


class FilterUniqueAcc(Filter):
    """
    Removes sequences that were already there (have same accession number.)
    """
    def __init__(self, config, table):
        super().__init__(config)
        self.table = table


    def filter(self, new_seqs):
        print('FilterUniqueAcc')
        # delete things in table
        new_seqs_unique = new_seqs.drop_duplicates(subset=['accession'], keep='first')
        drop = new_seqs.drop(new_seqs_unique.index.tolist())

        self.upd_new_seqs = new_seqs_unique
        self.del_table = drop

        check_df_index_unique(self.upd_new_seqs)
        check_filter_numbers(drop, new_seqs_unique, new_seqs)

        msg = "FilterUniqueAcc has lowered the new seqs df from {} to {}.\n".format(len(new_seqs), len(new_seqs_unique))
        write_msg_logfile(msg, self.config.workdir)


def check_df_index_unique(df):
    assert df.index.is_unique == True, df.index


def check_filter_numbers(remove, add, before):
    assert len(remove) + len(add) == len(before), (len(remove), len(add), len(before))


class FilterLength(Filter):
    """
    Removes sequences that are too short or too long.
    """
    def __init__(self, config, aln):
        super().__init__(config)
        self.aln = aln

    # todo-maybe: rewrite to use table instead of aln - table not needed anywhere else...

    def avg_length_seqs_in_aln(self):
        """
        Get the average length of the sequences in the alignment.

        :return: average length
        """
        orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in self.aln]
        avg_seqlen = sum(orig_seqlen) / len(orig_seqlen)
        return avg_seqlen


    def filter(self, new_seqs):
        """
        Main function for class. Filters per length

        :param new_seqs: sequences to be filtered
        :return: filtered sequences
        """
        print('filter FilterLength')
        avg_len = self.avg_length_seqs_in_aln()
        del_seq = pd.DataFrame()
        filter_new_seqs = pd.DataFrame()
        # assert len(new_seqs.index) == len(set(new_seqs.index)), (len(new_seqs.index), len(set(new_seqs.index)))

        for idx in new_seqs.index:
            seq = str(new_seqs.loc[idx, 'sseq'])
            assert '-' not in seq
            # print(self.config.maxlen * avg_len, len(seq), self.config.minlen * avg_len)
            if self.config.maxlen * avg_len > len(seq) > self.config.minlen * avg_len:
                filter_new_seqs = filter_new_seqs.append(new_seqs.loc[idx], ignore_index=True)
            else:
                # print('too long/too short')
                del_seq = del_seq.append(new_seqs.loc[idx], ignore_index=True)
                # del_seq = pd.concat([del_seq, new_seqs.loc[idx]], ignore_index=True, sort=True)  # note:concat does not work, as it is df and series
        # important to not assign new value in the else above - it duplicates entries, thats why we have the next one
        if not del_seq.empty:
            len_seqs = del_seq['sseq'].apply(len)
            del_seq['status_note'] = 'seq len: ' + len_seqs.astype(str) + \
                                     ' ({}-{})'.format(self.config.maxlen * avg_len, self.config.minlen * avg_len)
            del_seq.to_csv(os.path.join(self.config.workdir, 'wrong_seq_length.csv'), mode='a')

        check_df_index_unique(self.upd_new_seqs)
        check_filter_numbers(del_seq, filter_new_seqs, new_seqs)

        self.upd_new_seqs = filter_new_seqs
        self.del_table = del_seq

        msg = "Filter FilterLength reduced the new seqs from {} to {}.\n".format(len(new_seqs), len(filter_new_seqs))
        write_msg_logfile(msg, self.config.workdir)
