import sys
import numpy
import datetime
import pandas as pd
import numpy as np
import os
from dendropy import DnaCharacterMatrix, Tree

from . import blast
from . import ncbi_data_parser
from . import aln_updater
from . import write_msg_logfile


# DONE TODO: write concat class
# DONE TODO: write mrca list for non-monophyletic
# TODO: add add_lower_taxa function (palms dataset) - work in progress
# TODO: write tests
# DONE todo: otu filter by length - work in progress
# todo: write files to tmp...

class Update_data:
    def __init__(self, id_to_spn, aln, aln_schema, tre, tre_schema, config, mrca=None, blacklist=None):
        self.date = pd.Timestamp.today()
        self.config = config
        self.workdir = self.config.workdir
        self.status = 0  # -100 = blast, -1 = deleted, positive values = round, 0 = present at beginning
        self.aln_fn = aln
        self.aln_schema = aln_schema
        self.aln = None
        self.tre_fn = tre
        self.tr_schema = tre_schema
        self.tre = None
        self.blacklist = blacklist
        self.table = self.build_table_from_file(id_to_spn)
        self.add_seq_to_table()
        self.mrca = None
        self.set_mrca(config, mrca)

    def set_back_data_in_table(self):
        if self.config.different_level == True and self.status == 0:
            self.table['date'] = pd.Timestamp.strptime('01/01/00', '%d/%m/%y')
            self.table.loc[self.table['status'] > 0.5, 'status'] = self.status
            self.config.update_tree = False

    def set_mrca(self, config, mrca):
        """
        Set the mrca in the class for various input: None, one ncbi id, list of ncbi_ids.

        Note: list finds mrca of all ids, including ids that were not mentioned.
        Thus, not helpful for non-monophyletic groups!

        :param config: configuration class
        :param mrca: input mrca
        :return: sets self.mrca
        """
        if type(mrca) is int:
            self.mrca = mrca
        else:
            if type(mrca) is list:
                print('get mrca of all input!')
                mrca_list = mrca
            else:
                print('get mrca of no input!')
                mrca_list = self.table['ncbi_txid'].to_list()

            # # not needed because of mrca_list
            # ncbi_parser = ncbi_data_parser.Parser(names_file=config.ncbi_parser_names_fn,
            #                                       nodes_file=config.ncbi_parser_nodes_fn)
            # mrca_aln = ncbi_parser.get_mrca(mrca_list)
            # self.mrca = mrca_aln
            self.mrca = mrca_list

    def build_table_from_file(self, id_to_spn):
        """
        Build self.table from file. This is the object that basically holds all needed information.

        :param id_to_spn: file that translates tipnames to speciesnames. Must be accepted by ncbi
        :return: table with all sequence information available
        """
        print("Build table with information about sequences and taxa.")
        table = self.get_txid_for_name_from_file(id_to_spn)  # builds name id link
        table['status'] = 0
        table['date'] = pd.Timestamp.strptime('01/01/00', '%d/%m/%y')
        table['sseq'] = None
        table['status'] = table['status'].astype(int)
        table['ncbi_txid'] = table['ncbi_txid'].astype(int)
        if self.config.downtorank is not None:
            table['original_ncbi_txid'] = table['ncbi_txid'].astype(int)
            ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_parser_names_fn,
                                                  nodes_file=self.config.ncbi_parser_nodes_fn)
            table['ncbi_txid'] = np.vectorize(ncbi_parser.get_downtorank_id)(table['original_ncbi_txid'],
                                                                             self.config.downtorank)
        else:
            table['ncbi_txid'] = table['ncbi_txid'].astype(int)
        return table

    def add_seq_to_table(self):
        """
        Puts input sequences into the pandas table.

        :return:
        """
        aln = self.read_in_aln()
        queried_taxa = []
        for index in self.table.index:
            tip_name = self.table.loc[index, 'accession']
            for taxon, seq in aln.items():
                if taxon.label == tip_name:
                    seq = seq.symbols_as_string().replace("-", "").replace("?", "")
                    self.table.loc[index, 'sseq'] = seq
                    queried_taxa.append(tip_name)
        assert tip_name in queried_taxa, (tip_name, queried_taxa)

    def read_in_tree(self, taxon_namespace):
        """
        Reads in phylogeny from file.

        :param taxon_namespace: needs taxon_namespace from aln so that they match
        :return: tre als dendropy object
        """
        if taxon_namespace:
            tre = Tree.get(path=self.tre_fn, schema=self.tr_schema,
                           taxon_namespace=taxon_namespace, preserve_underscores=True)
        assert tre.taxon_namespace == taxon_namespace
        print(tre.as_string('newick'))
        return tre

    def read_in_aln(self):
        """
        Reads text alignment in as dendropy object.

        :return: alignment as dendroppy object
        """
        aln = DnaCharacterMatrix.get(path=self.aln_fn, schema=self.aln_schema)
        assert aln.taxon_namespace
        return aln

    def extend_with_unpublished(self):
        """
        This will use a local database to search for new sequences to add.
        Mend for sequences that are not yet published.

        :return:
        """
        sys.stdout.write('Find new sequences using a local database.\n')
        self.status += 1
        present_subset = self.table[self.table['status'] > -1]
        new_seqs_unpubl = pd.DataFrame(columns=['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date", 'accession',
                                                'pident', 'evalue', 'bitscore', 'sseq', 'title'])
        blast.make_local_blastdb(self.config.workdir, db='unpublished', path_to_db=self.config.unpubl_data)
        # get new seqs from seqs in blast table seq
        queried_taxa = []
        for index in present_subset.index:  # should be empty in later rounds...that is why it does not matter that new_seqs_local is reassigned
            tip_name = self.table.loc[index, 'accession']
            query_seq = self.table.loc[index, 'sseq']
            new_seq_tax = blast.get_new_seqs(query_seq, tip_name, self.config.blastdb, "unpublished", self.config)
            queried_taxa.append(tip_name)
            assert tip_name in queried_taxa, (tip_name, queried_taxa)
            new_seqs_unpubl = new_seqs_unpubl.append(new_seq_tax, ignore_index=True)
        new_seqs_unpubl = new_seqs_unpubl.drop(['accession;gi', 'ncbi_txid', 'ncbi_txn'], axis=1)
        name_txid_unpublished = self.get_txid_for_name_from_file(self.config.unpubl_names)
        new_seqs_unpubl = pd.merge(new_seqs_unpubl, name_txid_unpublished, on=['accession'], sort=True)
        new_seqs_unpubl['title'] = 'unpublished'
        new_seqs_unpubl['status'] = self.status
        # print(new_seqs_unpubl[['accession', 'ncbi_txid', 'ncbi_txn']])
        return new_seqs_unpubl

    def get_txid_for_name_from_file(self, tipname_id_fn):
        """
        Get taxon id for tipname, speciesname files using ncbi translation

        :param tipname_id_fn: file with tipname and corresponding species name
        :return:
        """
        columns = ['accession', 'ncbi_txn']
        name_id = pd.read_csv(tipname_id_fn, names=columns,  sep=",", header=None)
        ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_parser_names_fn,
                                              nodes_file=self.config.ncbi_parser_nodes_fn)
        name_id['ncbi_txid'] = name_id['ncbi_txn'].apply(ncbi_parser.get_id_from_name)
        return name_id

    def extend(self, new_seqs=None):
        """
        Gets new sequences from thos existing once which have not yet been blasted before.
        :param new_seqs: pandas dataframe with the new sequences retrieved earlier
        :return:
        """
        # create list of indice of subset list
        sys.stdout.write("Find new sequences using the BLAST database.\n")
        # sets back date for rerun with different tribe or so
        self.set_back_data_in_table()
        self.status += 1
        present_subset = self.table['status'] > -1
        today = pd.Timestamp.today()
        min_date_blast = today - pd.Timedelta(days=90)
        present_subset = self.table[present_subset == True]
        present_subset = present_subset[(present_subset.date <= min_date_blast)]
        new_seqs_local = pd.DataFrame(columns=['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date", 'accession',
                                               'pident', 'evalue', 'bitscore', 'sseq', 'title'])
        # get new seqs from seqs in blast table seq
        queried_taxa = []
        msg = present_subset[['accession', 'ncbi_txn',  'status']].to_string()
        write_msg_logfile(msg, self.config.workdir)
        for index in present_subset.index:  # should be empty in later rounds...that is why it does not matter that new_seqs_local is reassigned
            tip_name = self.table.loc[index, 'accession']
            msg = tip_name
            write_msg_logfile(msg, self.config.workdir)
            print('Blast: {}.'.format(tip_name))
            query_seq = self.table.loc[index, 'sseq']
            bef = self.table.iloc[index]['date']
            self.table.at[index, 'date'] = today  # this is the new version of pd.set_value(), sometimes it's iat
            aft = self.table.iloc[index]['date']
            assert bef != aft, (bef, aft)
            new_seq_tax = blast.get_new_seqs(query_seq, tip_name, self.config.blastdb, "Genbank", self.config)
            queried_taxa.append(tip_name)
            assert tip_name in queried_taxa, (tip_name, queried_taxa)
            new_seqs_local = new_seqs_local.append(new_seq_tax, ignore_index=True)
        assert len(new_seqs_local) > 0, new_seqs_local
        if self.blacklist is not None:
            new_seqs_local = self.remove_blacklist_items(new_seqs_local)
        return new_seqs_local

    def remove_blacklist_items(self, new_seqs):
        """
        Removes accession numbers specified in the config file that shall not be added
        :param new_seqs:  pandas dataframe with the new sequences retrieved earlier
        :return:
        """
        drop_boolean = np.where((new_seqs.accession.isin(self.blacklist)), True, False)
        new_seqs = new_seqs[drop_boolean != True]
        return new_seqs

    def add_new_seqs(self, new_seqs):
        """
        Adds new sequences to table, after removing those which were already present

        :param new_seqs:  pandas dataframe with the new sequences retrieved earlier
        :return: new_seqs without seqs that were added earlier
        """
        subcols = new_seqs[['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date", 'sseq', 'accession']]
        #subcols = subcols[~subcols['accession'].isin(self.table['accession'])]  # ~ is the pd not in/!  #todo move to unique filter - should not be needed anymore here
        subcols.loc[:, 'status'] = self.status

        assert subcols['status'].hasnans == False, subcols['status']
        self.table = self.table.append(subcols, ignore_index=True)
        assert self.table['status'].hasnans == False, self.table['status'].hasnans
        assert len(subcols) == len(self.table[self.table['status'] == self.status]), \
            (len(subcols), len(self.table[self.table['status'] == self.status]), subcols['accession'])
        msg = "add_new_seqs has lowered the new seqs df from {} to {}.\n".format(len(new_seqs), len(subcols))
        write_msg_logfile(msg, self.config.workdir)
        return subcols

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
        new_seqs = self.basic_filters(aln, self.mrca, new_seqs)
        # next filter need infos in table
        new_seqs = self.add_new_seqs(new_seqs)
        # filter for seq identity process
        new_seqs = self.compare_filter(new_seqs)
        return new_seqs

    def compare_filter(self, new_seqs):
        """
        This filter alteres entries in the table.

        :param new_seqs: pandas dataframe with the new sequences retrieved/not filtered from earlier
        :return:
        """
        if len(new_seqs) > 0:
            # #############
            present = self.table[self.table['status'] > -1]
            old_seqs = present[present['status'].astype(int) < self.status]
            seqs_to_comp = old_seqs
            # #######################
            f = FilterSeqIdent(self.config, seqs_to_comp, self.table, self.status)
            msg = "Time before Filter {}: {}.\n".format(f, datetime.datetime.now())
            write_msg_logfile(msg, self.config.workdir)
            before_filter = len(new_seqs)
            f.filter(new_seqs)
            new_seqs = f.upd_new_seqs  # assign new_seqs to last upd from filter before
            self.table = f.table
            del_seq = f.del_table
            assert len(del_seq) + len(new_seqs) == before_filter, (
                len(del_seq), len(new_seqs), before_filter)
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
        orig_len = len(new_seqs)
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
            internal_new_seq = len(new_seqs)
            f.filter(new_seqs)
            new_seqs = f.upd_new_seqs
            del_seq = f.del_table
            all_del = all_del.append(del_seq, ignore_index=True)
            assert len(del_seq) + len(new_seqs) == internal_new_seq, (
                len(del_seq), len(new_seqs), internal_new_seq)
        assert len(all_del) + len(new_seqs) == orig_len, (len(all_del), len(new_seqs), orig_len, all_del, new_seqs)
        all_del.to_csv('{}/deleted_early.csv'.format(self.config.workdir), mode='a')

        return new_seqs

    def run(self):
        """
        Main function of the class - Finds new sequences and adds them to the aln and tre.

        :return:
        """
        msg = "Begin update: {}.\n".format(datetime.datetime.now())
        write_msg_logfile(msg, self.config.workdir)

        columns = ['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date",
                   'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title']
        all_new_seqs = pd.DataFrame(columns=columns)

        aln = self.read_in_aln()
        tre = self.read_in_tree(aln.taxon_namespace)
        cleaner = aln_updater.InputCleaner(tre, aln, self.tre_fn, self.aln_fn, self.table, self.config, self.mrca)
        self.aln = cleaner.aln
        self.tre = cleaner.tre
        self.table = cleaner.table
        self.mrca = cleaner.mrca
        retrieved_seqs = 1
        if os.path.exists('{}/all_new_seqs.updated'.format(self.workdir)):
            # all_new_seqs = pd.read_csv('{}/all_new_seqs.updated'.format(self.workdir))
            self.table = pd.read_csv('{}/table.updated'.format(self.workdir))
        if self.config.different_level == False and os.path.exists('{}/all_new_seqs.updated'.format(self.workdir)):
            next
        else:
            while retrieved_seqs > 0:
                msg = "Time before blast: {}.\n".format(datetime.datetime.now())
                write_msg_logfile(msg, self.config.workdir)
                if self.config.unpublished is True:
                    new_seqs = self.extend_with_unpublished()
                    # print(new_seqs[['accession', 'ncbi_txid', 'ncbi_txn']])
                else:
                    new_seqs = self.extend()  # todo rename to find new seqs
                    msg = "Time after BLAST: {}.\n".format(datetime.datetime.now())
                    write_msg_logfile(msg, self.config.workdir)
                print('Length of new seqs before filtering: {}'.format(len(new_seqs)))
                new_seqs = self.call_filter(new_seqs, aln)
                new_seqs = new_seqs[~new_seqs['accession'].isin(all_new_seqs['accession'])]  # ~ is the pd not in/!
                # all_new_seqs = all_new_seqs.append(new_seqs, ignore_index=True)
                # all_new_seqs = new_seqs  # i replaces the one above for the different filter setting - TODO can probable all be removed - assert, etc
                print('Length of new seqs after filtering: {}'.format(len(new_seqs)))
                retrieved_seqs = len(new_seqs)
                msg = "Newly found seqs: {}.\n".format(len(new_seqs))
                write_msg_logfile(msg, self.config.workdir)
                # following assert is probably only important for debugging and all_new_seqs_table not used otherwise
                # all_new_seqs_table = self.table[self.table['status'] > 0]
                # assert len(all_new_seqs_table) == len(all_new_seqs), (len(all_new_seqs_table), len(all_new_seqs))
                if self.config.unpublished is True:
                    self.config.unpublished = False
                    # retrieved_seqs = 0
            msg = "Time before update tre/aln: {}.\n".format(datetime.datetime.now())
            write_msg_logfile(msg, self.config.workdir)
            msg = "Newly found seqs will be added to aln and tre: {}.\n".format(len(all_new_seqs))
            write_msg_logfile(msg, self.config.workdir)

            all_new_seqs.to_csv('{}/all_new_seqs.updated'.format(self.workdir))
            self.table.to_csv('{}/table.updated'.format(self.workdir))
        print('go to update aln and tre')
        self.update_aln_tre()
        self.table.to_csv('{}/seq_table.csv'.format(self.config.workdir))
        self.table.to_pickle('{}/seq_table.pickle'.format(self.config.workdir))
        msg = "Time finished: {}.\n".format(datetime.datetime.now())
        write_msg_logfile(msg, self.config.workdir)

    def update_aln_tre(self):
        print("update aln")
        aln = self.read_in_aln()
        tre = self.read_in_tree(aln.taxon_namespace)
        aln_upd = aln_updater.PhyAlnUpdater(tre, aln, self.status, self.table, self.config)
        aln_upd.update_data()
        aln_upd.write_labelled('updt_tre.tre')
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


def calculate_mean_sd(bitscores):
    """Calculates standard deviation and mean of scores which are used as a measure of sequence differentiation
    for a given taxon, being used to select a random representative of a taxon later.

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
        self.status = status

    def initialize(self, config):
        self.ncbi_parser = ncbi_data_parser.Parser(names_file=config.ncbi_parser_names_fn,
                                                   nodes_file=config.ncbi_parser_nodes_fn)

    def filter(self, new_seqs, downtorank=None):
        print("filter FilterNumberOtu")
        assert len(new_seqs) == len(self.table[self.table['status'] == self.status]), \
            (len(new_seqs), len(self.table[self.table['status'] == self.status]), new_seqs['accession'])

        # get all seqs and ids from aln and before for next filter
        # #############
        present = self.table[self.table['status'] > -1]
        added_before = present[present['status'].astype(int) < self.status]
        # #######################
        # present = self.table['status'] != -1
        # old_status = self.table['status'].astype(int) < self.status
        # filter_table = self.table[present & old_status]
        if downtorank is None:
            tax_ids = new_seqs['ncbi_txid']  # doubles with next line
            ns_otu_df = new_seqs['ncbi_txid']
            table_otu_df = added_before['ncbi_txid']
        else:
            new_seqs = self.set_downtorank(new_seqs, downtorank)  # todo: this is likely doubled now, since i chnages now ncbi id to rank and original_ncbi_id for the real one...
            table = self.set_downtorank(added_before, downtorank)
            tax_ids = new_seqs['downtorank']
            ns_otu_df = new_seqs['downtorank']
            table_otu_df = table['downtorank']
        self.table[self.table['status'] == self.status]['accession'].to_csv('debug_table_samestatus', header=True)
        assert len(new_seqs) == len(self.table[self.table['status'] == self.status]), \
            (len(new_seqs), len(self.table[self.table['status'] == self.status]),
             new_seqs['accession'], self.table[self.table['status'] == self.status]['accession'])

        for txid in set(tax_ids):
            # generate otu_subset
            # print('filter per otu')
            ns_otu_dict = new_seqs[ns_otu_df == txid]
            table_otu_dict = added_before[table_otu_df == txid]
            table_otu_dict = table_otu_dict[~table_otu_dict['status'].isin([-1, 0])]
            if len(table_otu_dict) < self.config.threshold:  # if we still want to add - might be obsolete
                if len(ns_otu_dict) + len(table_otu_dict) > self.config.threshold:  # filter
                    if len(table_otu_dict) == 0:  # new taxa, select random seq for blast
                        # print('taxa new')
                        if self.config.filtertype == 'blast':
                            filtered_acc = self.wrapper_filter_blast_otu(ns_otu_dict, table_otu_dict, 'accession', txid)
                            filtered = new_seqs[new_seqs['accession'].isin(filtered_acc)]
                        elif self.config.filtertype == 'length':
                            print('filter otu by length')
                            filtered = self.select_seq_by_length(ns_otu_dict, table_otu_dict)
                        else:
                            print('should not happen')
                            sys.exit(2)
                    else:  # select seq from aln
                        # print('taxa present already')
                        if self.config.filtertype == 'blast':
                            filtered_acc = self.wrapper_filter_blast_otu(ns_otu_dict, table_otu_dict, 'tip_name', txid)
                            filtered = new_seqs[new_seqs['accession'].isin(filtered_acc)]
                        elif self.config.filtertype == 'length':
                            print('filter otu by length')
                            filtered = self.select_seq_by_length(ns_otu_dict, table_otu_dict)
                        else:
                            print('should not happen')
                            sys.exit(2)
                elif len(ns_otu_dict) + len(table_otu_dict) <= self.config.threshold:  # filter
                    # add all
                    # print('add all')
                    filtered = ns_otu_dict
                    # print(len(filtered))
                self.upd_new_seqs = self.upd_new_seqs.append(filtered)
            elif len(table_otu_dict) > self.config.threshold:
                print('sample size to big')
                sys.exit(-3)
            else:
                print('sample size correct')
        print('len upd new seqs: ')
        print(len(self.upd_new_seqs))
        not_selected = list(set(new_seqs['accession'].values) - set(self.upd_new_seqs['accession'].values))
        selected = list(set(self.upd_new_seqs.index))
        assert len(selected) + len(not_selected) == len(new_seqs), (len(selected), len(not_selected), len(new_seqs))
        del_tab = new_seqs[new_seqs['accession'].isin(not_selected)]
        del_tab['status'] = 'deleted - number otu'
        self.del_table = del_tab
        assert len(self.del_table) + len(self.upd_new_seqs) == len(new_seqs), \
            (len(self.del_table), len(self.upd_new_seqs), len(new_seqs))
        if len(not_selected) > 0:
            self.table.loc[self.table['accession'].isin(not_selected), 'status'] = -1
            self.table.loc[self.table['accession'].isin(not_selected), 'status_note'] = 'too many seqs of same tax_id'

        # find why new seqs and self.table show different acc:
        added = self.table[self.table['status'] == self.status]
        assert len(added) == len(selected), \
            (len(added), len(selected), self.table[['accession', 'status']], self.upd_new_seqs[['accession', 'status']])
        assert set(added['accession'].values) == set(self.upd_new_seqs['accession'].values)
        msg = "Filter FilterNumberOtu has lowered the new seqs df from {} to {}.\n".format(len(new_seqs), len(self.upd_new_seqs))
        write_msg_logfile(msg, self.config.workdir)

    def set_downtorank(self, new_seqs, downtorank):
        """
        Method is used to get the corresponding rank if higher than lowest otu by ncbi.
        
        :param new_seqs:
        :param downtorank:
        :return:
        """
        # print("set_downtorank")
        new_seqs.loc[:, 'downtorank'] = 0
        if self.ncbi_parser is None:
            self.initialize(self.config)
        for txid in set(new_seqs['ncbi_txid']):
            downtorank_id = self.ncbi_parser.get_downtorank_id(txid, downtorank)
            new_seqs.loc[new_seqs.ncbi_txid == txid, 'downtorank'] = downtorank_id
        return new_seqs

    def select_seq_by_length(self, filter_dict, exist_dict):
        """

        This is another mode to filter the sequences, if there are more than the threshold amount available.
        This one selects new sequences by length instead of by score values. It is selected by "selectby='length'".

         :param filter_dict:
        :param exist_dict:
        :return: filtered sequences
        """
        print("select_seq_by_length")

        #len_seqs_new = filter_dict[len(filter_dict['sseq'])]
        len_seqs_new = filter_dict['sseq'].apply(len)
        # print(filter_dict['sseq'])
        print(len_seqs_new)
        # print(type(len_seqs_new))
        amnt_missing_seq = self.config.threshold - len(exist_dict)
        print('amount missing data:', amnt_missing_seq)
        if len(len_seqs_new) > amnt_missing_seq:
            select = pd.DataFrame(columns=['ncbi_txn', 'ncbi_txid', 'status', 'status_note', "date", 'accession',
                                                'pident', 'evalue', 'bitscore', 'sseq', 'title'])
            for i in range(1,amnt_missing_seq+1):
                print(i)
                idxmax_val = len_seqs_new.idxmax()
                max_val = len_seqs_new.max()
                # print(max_val, idxmax_val)
                select = select.append(filter_dict.loc[idxmax_val])
                filter_dict = filter_dict.drop([idxmax_val])
                len_seqs_new = len_seqs_new.drop([idxmax_val])
            assert len(select) == amnt_missing_seq, (len(select), amnt_missing_seq)
        else:
            print('DO WE EVER GET HERE? length filter')
            select = len_seqs_new
            assert len(select) <= amnt_missing_seq
        print(select)
        return select

    def wrapper_filter_blast_otu(self, filter_dict, exist_dict, columnname, txid):
        """
        get subsample of new seqs for otu_sample by selecting random new seqs of taxa that are completely new or for
        random existing seq.

        :param filter_dict:
        :param exist_dict:
        :param columnname:
        :param txid:
        :return:
        """
        print('wrapper_filter_blast_otu')
        rndm = filter_dict.sample(1)
        idx = rndm.index.values.tolist()[0]
        query_seq = rndm.loc[idx, 'sseq']
        query_seq_id = rndm.loc[idx, 'accession']
        blast.write_local_blast_files(self.config.workdir, query_seq_id, query_seq, db=False,
                                      fn=txid)
        # write local db and create it
        if columnname == 'tip_name':
            for_db = filter_dict
        else:
            for_db = filter_dict.drop([idx])
        # filter_seq_db has to be deleted for every new filterun
        if os.path.exists('{}/tmp/filter_seq_db'.format(self.config.workdir)):
            os.remove('{}/tmp/filter_seq_db'.format(self.config.workdir))
        for idx in for_db.index:
            seq_id = for_db.loc[idx, 'accession']
            seq = for_db.loc[idx, 'sseq']
            # print('write local blast files')
            blast.write_local_blast_files(self.config.workdir, seq_id, seq, db=True)
        blast.make_local_blastdb(self.config.workdir, db='filterrun')
        filter_blast_seqs = blast.get_new_seqs(query_seq, txid, '{}/tmp'.format(self.config.workdir),
                                               'filterrun', self.config)
        subfilter_blast_seqs = remove_mean_sd_nonfit(filter_blast_seqs)
        filtered_seqs = self.select_number_missing(subfilter_blast_seqs, exist_dict)
        filtered_acc = list(filtered_seqs['accession'])
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


class FilterSeqIdent(Filter):
    """
    Filters new sequences that are identical (seq and taxon id) to something already in the data.
    """
    def __init__(self, config, old_seqs, table, status):
        super().__init__(config)
        self.old_seqs = old_seqs
        self.table = table
        self.status = status

    def filter(self, new_seqs):
        assert len(new_seqs) == len(self.table[self.table['status'] == self.status]), \
            (len(new_seqs), len(self.table[self.table['status'] == self.status]), new_seqs['accession'])

        print("filter FilterSeqIdent")
        to_add = pd.DataFrame()
        new_seqs_seqs = new_seqs['sseq']

        for idx in set(new_seqs.index):
            # generate seq and compare lsit
            new_seqs_seqs_drop = new_seqs_seqs.drop([idx])
            txid_new = new_seqs.loc[idx, 'ncbi_txid']

            old_seqs_seq = self.old_seqs['sseq']
            seq = new_seqs.loc[idx, 'sseq']
            txid_old = self.old_seqs.loc[self.old_seqs['sseq'] == seq].index

            if seq not in old_seqs_seq:  # if we still want to add - might be obsolete
                if seq not in new_seqs_seqs_drop:  # seq never ident
                    to_add = to_add.append(new_seqs.loc[idx])
                else:  # seq ident to new
                    txid_new_seqident = new_seqs_seqs_drop.loc[new_seqs_seqs_drop['sseq'] == seq].index
                    if txid_new_seqident == txid_old:
                        to_del = new_seqs.loc[idx]
                        self.del_table = self.del_table.append(to_del)
                        self.table.loc[idx, 'status'] = -1
                        self.table.loc[idx, 'status_note'] = 'same seq ident'
                    else:
                        to_add = to_add.append(new_seqs.loc[idx])
            else:
                print('sequence ident matches old, same txid?')
                if txid_new == txid_old:
                    to_del = new_seqs.loc[idx]
                    self.del_table = self.del_table.append(to_del)
                    self.table.loc[idx, 'status'] = -1
                    self.table.loc[idx, 'status_note'] = 'same seq ident'
                else:
                    to_add = to_add.append(new_seqs.loc[idx])
        assert len(self.del_table) + len(to_add) == len(new_seqs), (len(self.del_table), len(to_add), len(new_seqs))

        self.upd_new_seqs = to_add
        assert self.upd_new_seqs.index.is_unique is True, self.upd_new_seqs.index

        msg = "Filter FilterSeqIdent has lowered the new seqs df from {} to {}.\n".format(len(new_seqs), len(to_add))
        write_msg_logfile(msg, self.config.workdir)


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
        for tax_id in set(new_seqs['ncbi_txid'].astype(int)):
            tax_id = int(tax_id)
            # print(tax_id, self.mrca)
            mrca_tx = self.ncbi_parser.match_id_to_mrca(tax_id, self.mrca)
            select_TF = new_seqs['ncbi_txid'].astype(int) == tax_id
            if mrca_tx == self.mrca:
                to_add = new_seqs[select_TF]
                assert len(to_add) != 0, len(to_add)
                self.upd_new_seqs = pd.concat([self.upd_new_seqs, to_add], axis=0, sort=True)
            else:
                # print('out of mrca')
                to_del = new_seqs[select_TF]
                to_del.loc[:, 'status'] = 'deleted - mrca'
                self.del_table = pd.concat([self.del_table, to_del], axis=0, sort=True)
        assert len(self.del_table) + len(self.upd_new_seqs) == len(new_seqs), \
            (len(self.del_table), len(self.upd_new_seqs), len(new_seqs))
        assert self.upd_new_seqs.index.is_unique is True, self.upd_new_seqs.index

        msg = "Filter FilterMRCA has lowered the new seqs df from {} to {}.\n".format(len(new_seqs),
                                                                                      len(self.upd_new_seqs))
        write_msg_logfile(msg, self.config.workdir)
        if len(to_del) > 0:
            to_del.to_csv('{}/wrong_mrca.csv'.format(self.config.workdir), mode='a')


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
        deltab['status'] = 'deleted - evalue'
        assert len(deltab) + len(upd_new_seqs) == len(new_seqs), (len(deltab), len(upd_new_seqs), len(new_seqs))

        self.del_table = deltab
        self.upd_new_seqs = upd_new_seqs
        assert self.upd_new_seqs.index.is_unique == True, self.upd_new_seqs.index

        msg = "Filter FilterBLASTThreshold has lowered the new seqs df from {} to {}.\n".format(len(new_seqs),
                                                                                                len(upd_new_seqs))
        write_msg_logfile(msg, self.config.workdir)
        if len(deltab) > 0:
            deltab.to_csv('{}/low_blast_threshold.csv'.format(self.config.workdir), mode='a')


class FilterUniqueAcc(Filter):
    """
    Removes sequences that were already there (have same accession number.)
    """
    def __init__(self, config, table):
        super().__init__(config)
        self.table = table

    def filter(self, new_seqs):  # todo: rewrite: this assumes thats upd_new-seqs has information in before, that is not true anymore
        print('FilterUniqueAcc')
        # delete things in table
        new_seqs_drop = new_seqs.drop_duplicates(subset=['accession'])
        assert len(self.upd_new_seqs) == 0, self.upd_new_seqs  # if this does not break, remove next statement
        new_seqs_unique = self.upd_new_seqs.drop_duplicates(subset=['accession'], keep='first')
        # print(len(new_seqs_drop) > 0, len(new_seqs_unique) > 0)
        if len(new_seqs_drop) > 0 and len(new_seqs_unique) > 0:
            print('if you never see this, delete it...')
            # Select all duplicate rows based on one column
            del_seq = new_seqs[new_seqs.duplicated(['accession'])]
            new_seqs_unique_fn = new_seqs_unique
        else:
            new_seqs_unique_fn = new_seqs_drop
            dropped = new_seqs.duplicated(subset=['accession'])
            del_seq = pd.concat([self.del_table, new_seqs[dropped == True]], sort=True)
        self.upd_new_seqs = new_seqs_unique_fn
        self.del_table = del_seq
        assert len(del_seq) + len(new_seqs_unique_fn) == len(new_seqs), \
            (len(del_seq), len(new_seqs_unique_fn), len(new_seqs))
        assert self.upd_new_seqs.index.is_unique == True, self.upd_new_seqs.index
        msg = "Filter FilterUniqueAcc has lowered the new seqs df from {} to {}.\n".format(len(new_seqs),
                                                                                           len(new_seqs_unique_fn))
        write_msg_logfile(msg, self.config.workdir)


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
        assert len(new_seqs.index) == len(set(new_seqs.index)), (len(new_seqs.index), len(set(new_seqs.index)))

        for idx in new_seqs.index:
            seq = str(new_seqs.loc[idx, 'sseq'])
            assert '-' not in seq
            if self.config.maxlen * avg_len > len(seq) > self.config.minlen * avg_len:
                filter_new_seqs = filter_new_seqs.append(new_seqs.loc[idx], ignore_index=True)
            else:
                print('too long/too short')
                del_seq = del_seq.append(new_seqs.loc[idx], ignore_index=True)
        # important to not assign new value in the else above - it duplicates entries, thats why we have the next one
        if len(del_seq) > 0:
            len_seqs = del_seq['sseq'].apply(len)
            print(len_seqs)
            del_seq.loc[:, 'status'] = 'deleted - seq len wrong: {}'.format(len(del_seq.loc[:, 'sseq']))
            # TODO: seq len is not recalculated.... maybe this fixes it
            del_seq.to_csv('{}/seq_too_short-long.csv'.format(self.config.workdir), mode='a')
            print(some)

        assert self.upd_new_seqs.index.is_unique == True, self.upd_new_seqs.index  # used to find errors -
        assert len(del_seq) + len(filter_new_seqs) == len(new_seqs), (len(del_seq), len(filter_new_seqs), len(new_seqs))
        self.upd_new_seqs = filter_new_seqs
        self.del_table = del_seq
        msg = "Filter FilterLength has lowered the new seqs df from {} to {}.\n".format(len(new_seqs), len(filter_new_seqs))
        write_msg_logfile(msg, self.config.workdir)
