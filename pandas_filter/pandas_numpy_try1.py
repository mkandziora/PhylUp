import sys
import numpy
import datetime
import pandas as pd
import numpy as np
import os
from ete3 import NCBITaxa
from dendropy import DnaCharacterMatrix, Tree
# from copy import deepcopy


from . import blast
from . import ncbi_data_parser
from . import aln_updater


def clean_inputlabel(tipname):
    """
    Cleans original tipnames in input data - so that python does not do funny things with it.

    :param tipname: original input from file
    :return: cleaned label
    """
    clean_lab = tipname.replace("-", "")
    clean_lab = clean_lab.replace(" ", "")
    clean_lab = clean_lab.replace("_", "")
    clean_lab = clean_lab.replace("'", "")
    clean_lab = clean_lab.replace("/", "")
    return clean_lab

# TODO: add add_lower_taxa function (palms dataset)
# TODO: add gb_id_filename naming scheme for reuse of blast files
# TODO: add unublished seq add
# TODO: write tests


class Update_data:
    def __init__(self, id_to_spn, aln, aln_schema, tre, tre_schema, config, mrca=None, blacklist=None):
        self.otu_counter = 1
        self.date = datetime.date.today()
        self.config = config
        self.workdir = self.config.workdir
        self.status = 0
        self.aln_fn = aln
        self.aln_schema = aln_schema
        self.aln = None
        self.tr_fn = tre
        self.tr_schema = tre_schema
        self.tre = None
        self.blacklist = blacklist

        self.table = self.build_table_from_file(id_to_spn)
        self.add_seq_to_table()

        self.mrca = None
        self.set_mrca(config, mrca)

        self.filtertypes = None

    def set_mrca(self, config, mrca):
        '''
        Set the mrca in the class for various input: None, one ncbi id, list of ncbi_ids.

        Note list finds mrca of all ids, including ids that were not mentioned.
        Thus, not helpful for non-monophyletic groups!

        :param config: configuration class
        :param mrca: input mrca
        :return: sets self.mrca
        '''
        if type(mrca) is int:
            self.mrca = mrca
        else:
            if type(mrca) is list:
                print('get mrca of all input!')
                mrca_list = mrca
            else:
                print('get mrca of no input!')
                mrca_list = self.table['ncbi_txid'].to_list()
            ncbi_parser = ncbi_data_parser.Parser(names_file=config.ncbi_parser_names_fn,
                                                  nodes_file=config.ncbi_parser_nodes_fn)
            mrca_aln = ncbi_parser.get_mrca(mrca_list)
            self.mrca = mrca_aln

    def build_table_from_file(self, id_to_spn):
        print("build table")
        otu_dict = pd.DataFrame(
            columns=['ncbi_txn', 'ncbi_txid', 'org_sp_name', 'tip_name', 'status', "date", 'sseq', 'index'])
        with open(id_to_spn, mode="r") as infile:
            for lin in infile:
                # clean input data
                tipname, species = lin.strip().split(",")
                clean_lab = clean_inputlabel(tipname)
                assert tipname not in otu_dict, ("standardized label ('{}') of `{}` "
                                                 "already exists".format(clean_lab, tipname))
                species = species.replace("_", " ")
                # get ncbi infos
                ncbi = NCBITaxa()
                name2taxid = ncbi.get_name_translator([species])
                # print(len(name2taxid.items()))
                if len(name2taxid.items()) == 1:
                    ncbiid = name2taxid[species][0]
                elif len(name2taxid.items()) > 1:
                    print("write method")
                    sys.exit(-1)
                else:
                    print('Name not known by ncbi')
                    sys.exit(-1)
                ncbi_spn = species
                status = 0
                otu_dict = otu_dict.append({'ncbi_txn': ncbi_spn, 'ncbi_txid': ncbiid, 'org_sp_name': ncbi_spn,
                                            'tip_name': tipname, 'status': status,
                                            "date": datetime.datetime.strptime('01/01/00', '%d/%m/%y'),
                                            'sseq': None, "index": self.otu_counter}, ignore_index=True)
                self.otu_counter += 1
        return otu_dict

    def add_seq_to_table(self):
        '''
        Puts input sequences into the pandas table.

        :return:
        '''
        aln = self.read_in_aln()
        queried_taxa = []
        for index in self.table.index:
            tip_name = self.table.loc[index, 'tip_name']
            for taxon, seq in aln.items():
                if taxon.label == tip_name:
                    seq = seq.symbols_as_string().replace("-", "").replace("?", "")
                    self.table.loc[index, 'sseq'] = seq
                    queried_taxa.append(tip_name)
        assert tip_name in queried_taxa, (tip_name, queried_taxa)

    def read_in_tree(self, taxon_namespace):
        tre = Tree.get(path=self.tr_fn,
                       schema=self.tr_schema,
                       taxon_namespace=taxon_namespace,
                       preserve_underscores=True)
        return tre

    def read_in_aln(self):
        '''
        Reads text alignment in as dendropy object.

        :return:
        '''
        # with open(self.aln_fn, "r") as fin:
        #     filedata = fin.read()
        # # Write the file out again
        # with open("{}/orig_inputaln.fasta".format(self.workdir), "w") as aln_file:
        #     aln_file.write(filedata)
        # filedata = filedata.replace("?", "-")
        # self.aln_fn = "{}/updt_aln.fasta".format(self.workdir)
        # with open(self.aln_fn, "w") as aln_file:
        #     aln_file.write(filedata)
        # # use replaced aln as input
        aln = DnaCharacterMatrix.get(path=self.aln_fn, schema=self.aln_schema)
        assert aln.taxon_namespace
        return aln

    def extend(self, new_seqs=None, date=None):
        """
        """
        # todo: write using status, currently everything is blasted again

        # create list of indice of subset list
        sys.stdout.write("func extend")
        self.status += 1
        present_subset_df = self.table['status'] != "excluded"
        today = datetime.date.today()
        min_date_blast = today - datetime.timedelta(days=90)
        present_subset_df = self.table[present_subset_df == True]
        present_subset_df = present_subset_df[(present_subset_df.date <= min_date_blast)]

        present_subset = present_subset_df

        print(len(present_subset))

        new_seqs_local = pd.DataFrame(columns=['ncbi_txn', 'ncbi_txid', 'org_sp_name', 'tip_name', 'status', 'status_note', "date",
                                               'index', 'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title'])

        queried_taxa = []
        for index in present_subset.index:
            print('blast table seq')

            tip_name = self.table.loc[index, 'tip_name']
            self.aln = self.read_in_aln()
            query_seq = self.table.loc[index, 'sseq']
            self.table.loc[index, 'date'] = today  # TODO: is here the blast date the issue???
            new_seq_tax = blast.get_new_seqs(query_seq, tip_name, self.config.blastdb, "Genbank", self.config)
            queried_taxa.append(tip_name)
            assert tip_name in queried_taxa, (tip_name, queried_taxa)
            new_seqs_local = new_seqs_local.append(new_seq_tax, ignore_index=True)

        if new_seqs is not None:
            print(len(new_seqs))
            print(new_seqs)
            print('add new seqs to blast table')
            # add to present_subset the sequences which have been added from former new_seqs searches
            new_seqs_subset = new_seqs[(new_seqs.date <= min_date_blast)]
            print(new_seqs_subset.date)

            print(len(new_seqs_subset))
            for index in new_seqs_subset.index:
                tip_name = new_seqs_subset.loc[index, 'accession']
                self.aln = self.read_in_aln()

                query_seq = new_seqs_subset.loc[index, 'sseq']
                new_seqs_subset.loc[index, 'date'] = today  # TODO: is here the blast date issue???
                new_seq_tax = blast.get_new_seqs(query_seq, tip_name, self.config.blastdb, "Genbank", self.config)

                queried_taxa.append(tip_name)

                assert tip_name in queried_taxa, (tip_name, queried_taxa)
                new_seqs_local = new_seqs_local.append(new_seq_tax, ignore_index=True)

        assert len(new_seqs_local) > 0, new_seqs_local
        new_seqs_local = pd.concat([new_seqs_local, new_seqs], sort=True)
        if self.blacklist is not None:
            new_seqs_local = self.remove_blacklist_items(new_seqs_local)
        return new_seqs_local

    def remove_blacklist_items(self, new_seqs_local):
        remove_TF = np.where((new_seqs_local.accession.isin(self.blacklist)), True, False)
        new_seqs_local = new_seqs_local[remove_TF != True]
        return new_seqs_local

    def add_new_seqs(self, new_seqs, status=0):
        subcols = new_seqs[['ncbi_txn', 'ncbi_txid', 'org_sp_name', 'tip_name', 'status', 'status_note', "date", 'sseq', 'index']]
        subcols['tip_name'] = new_seqs['accession']

        self.table = self.table.append(subcols, ignore_index=True)

    def call_filter(self, new_seqs, aln, mrca):
        print('call filter')
        print(len(new_seqs))
        orig_len = len(new_seqs)
        all_del = pd.DataFrame(
            columns=['ncbi_txn', 'ncbi_txid', 'org_sp_name', 'tip_name', 'status', 'status_note', "date", 'index',
                     'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title'])

        # everyones filter
        self.filtertypes = [FilterUniqueAcc(self.config),
                            FilterBLASTThreshold(self.config),
                            FilterLength(self.config, aln),
                            FilterMRCA(self.config, mrca),
                            FilterSeqIdent(self.config, aln)
                            ]

        for f in self.filtertypes:
            print('length new seqs before individual filters')
            print(len(new_seqs))
            print(f)
            internal_new_seq = len(new_seqs)
            f.filter(new_seqs)
            new_seqs = f.upd_new_seqs
            print(len(new_seqs))

            del_seq = f.del_table
            all_del = all_del.append(del_seq, ignore_index=True)
            assert len(del_seq) + len(new_seqs) == internal_new_seq, (
                len(del_seq), len(new_seqs), internal_new_seq)
        assert len(all_del) + len(new_seqs) == orig_len, (len(all_del), len(new_seqs), orig_len, all_del, new_seqs)

        new_seqs = f.upd_new_seqs  # assign new_seqs to last upd from filter before

        f = FilterNumberOtu(self.config, self.table, aln)
        f.filter(new_seqs, self.config.downtorank)
        new_seqs = f.upd_new_seqs
        print('new seqs after number otu')
        print(len(new_seqs))
        all_del = f.del_table
        return new_seqs, all_del

    def run(self):
        all_new_seqs = pd.DataFrame(
            columns=['ncbi_txn', 'ncbi_txid', 'org_sp_name', 'tip_name', 'status', 'status_note', "date", 'index',
                     'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title'])
        aln = self.read_in_aln()
        tre = self.read_in_tree(aln.taxon_namespace)

        cleaner = aln_updater.InputCleaner(tre, aln, self.aln_fn, self.table, self.config, self.mrca)

        self.aln = cleaner.aln
        self.tre = cleaner.tre
        self.table = cleaner.table
        self.mrca = cleaner.mrca

        retrieved_seqs = 1
        new_seqs = None
        while retrieved_seqs > 0:
            new_seqs = self.extend(new_seqs)  # todo rename to find new seqs
            print(' len new seqs before filter - round after')
            print(self.status)

            print(len(new_seqs))
            new_seqs, all_del = self.call_filter(new_seqs, aln, self.mrca)
            print('len new seqs after filter - round after')
            print(len(new_seqs))

            retrieved_seqs = len(new_seqs)

            # todo add loop to only add new gb_ids.
            all_new_seqs = all_new_seqs.append(new_seqs, ignore_index=True)
        print('last filter')
        all_new_seqs, all_del = self.call_filter(all_new_seqs, self.aln, self.mrca)
        print(len(all_new_seqs))
        self.add_new_seqs(all_new_seqs)
        self.update_aln_tre(all_new_seqs[['accession', 'sseq']])
        self.table.to_csv('{}/seq_table.csv'.format(self.config.workdir))
        self.table.to_pickle('{}/seq_table.pickle'.format(self.config.workdir))

    def update_aln_tre(self, new_seq_acc_seq):
        print("update aln")
        aln = self.read_in_aln()
        tre = self.read_in_tree(aln.taxon_namespace)

        aln_upd = aln_updater.PhyAlnUpdater(tre, aln, self.status, new_seq_acc_seq, self.table, self.config)
        aln_upd.update_data()
        aln_upd.write_labelled()


# ################################################################################################################
class Filter(object):
    def __init__(self, config):
        self.config = config
        self.upd_new_seqs = pd.DataFrame(
            columns=['ncbi_txn', 'ncbi_txid', 'org_sp_name', 'tip_name', 'status', 'status_note', "date", 'index',
                     'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title'])
        self.del_table = pd.DataFrame(
            columns=['ncbi_txn', 'ncbi_txid', 'org_sp_name', 'tip_name', 'status', 'status_note', "date", 'index',
                     'accession', 'pident', 'evalue', 'bitscore', 'sseq', 'title'])
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

    def __init__(self, config, table, aln):
        super().__init__(config)
        self.table = table
        self.aln = aln

    def initialize(self, config):
        self.ncbi_parser = ncbi_data_parser.Parser(names_file=config.ncbi_parser_names_fn,
                                                   nodes_file=config.ncbi_parser_nodes_fn)

    # def build_table_from_dict(self, new_blast_seq_dict):
    #     print("build table dict")
    #
    #     otu_dict = pd.DataFrame(columns=['ncbi_txn', 'ncbi_txid', 'org_sp_name', 'tip_name', 'status', "date", 'index'])
    #
    #     print(new_blast_seq_dict.columns)
    #     subcols = new_blast_seq_dict[['ncbi_txn', 'ncbi_txid', 'status', "date"]]
    #     subcols['org_sp_name'] = None
    #     subcols['tip_name'] = None
    #     subcols['index'] = None
    #     print(subcols)
    #     for idx in subcols.index:
    #         print(new_blast_seq_dict.loc[idx, 'accession'])
    #         subcols.loc[idx, 'tip_name'] = new_blast_seq_dict.loc[idx, 'accession']
    #
    #         print(subcols.loc[idx, 'tip_name'])
    #         print(some)
    #         otu_dict = otu_dict.append(subcols.iloc[idx], ignore_index=True)
    #     return otu_dict

    def filter(self, new_seqs, downtorank=None):
        print("filter FilterNumberOtu")
        if downtorank is None:
            tax_ids = new_seqs['ncbi_txid']
            ns_otu_df = new_seqs['ncbi_txid']
            table_otu_df = self.table['ncbi_txid']
        else:
            new_seqs = self.set_downtorank(new_seqs, downtorank)
            table = self.set_downtorank(self.table, downtorank)
            tax_ids = new_seqs['downtorank']
            ns_otu_df = new_seqs['downtorank']
            table_otu_df = table['downtorank']

        for txid in set(tax_ids):
            # generate otu_subset
            print('filter per otu')
            ns_otu_dict = new_seqs[ns_otu_df == txid]
            table_otu_dict = self.table[table_otu_df == txid]
            if len(table_otu_dict) < self.config.threshold:  # if we still want to add - might be obsolete
                if len(ns_otu_dict) + len(table_otu_dict) > self.config.threshold:  # filter
                    if len(table_otu_dict) == 0:  # new taxa, select random seq for blast
                        print('taxa new')
                        print(ns_otu_dict[['accession', 'ncbi_txn', 'ncbi_txid']])
                        if self.config.filtertype == 'blast':
                            filtered = self.wrapper_filter_blast_otu(ns_otu_dict, table_otu_dict, 'accession', txid)
                        elif self.config.filtertype == 'length':
                            print('filter otu by length')
                            sys.exit(2)
                        else:
                            print('should not happen')
                            sys.exit(2)
                    else:  # select seq from aln
                        print('taxa present already')
                        if self.config.filtertype == 'blast':
                            filtered = self.wrapper_filter_blast_otu(ns_otu_dict, table_otu_dict, 'tip_name', txid)
                        elif self.config.filtertype == 'length':
                            print('filter otu by length')
                            sys.exit(2)
                        else:
                            print('should not happen')
                            sys.exit(2)

                else:  # add all
                    print('add all')
                    filtered = ns_otu_dict
                self.upd_new_seqs = pd.concat([self.upd_new_seqs, filtered], axis=0, ignore_index=True, sort=True)
            else:
                print('sample large enough')
        #subset = self.upd_new_seqs
        not_selected = list(set(new_seqs.index) - set(filtered.index))
        print(not_selected)
        # del_tab = new_seqs.ix[not_selected]  # is deprecated
        del_tab = new_seqs[new_seqs.index.isin(not_selected)]
        del_tab['status'] = 'deleted - number otu'
        self.del_table = del_tab

        # todo add assert

        lfd = "{}/logfile".format(self.config.workdir)
        with open(lfd, "a") as log:
            log.write("Filter FilterBLASTThreshold has lowered the new seqs df from {} to {}.\n".format(len(new_seqs),
                                                                                                        len(to_add)))

    def set_downtorank(self, new_seqs, downtorank):
        print("set_downtorank")
        new_seqs['downtorank'] = 0
        if self.ncbi_parser is None:
            self.initialize(self.config)
        for txid in set(new_seqs['ncbi_txid']):
            downtorank_id = self.ncbi_parser.get_downtorank_id(txid, downtorank)
            new_seqs.loc[new_seqs.ncbi_txid == txid, 'downtorank'] = downtorank_id
        return new_seqs

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
        query_seq_id = rndm.loc[idx, columnname]
        print('write query file')
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
            seq_id = for_db.loc[idx, columnname]
            seq = for_db.loc[idx, 'sseq']
            blast.write_local_blast_files(self.config.workdir, seq_id, seq, db=True, fn=txid)
        # ###########
        blast.make_local_blastdb(self.config.workdir, db='filterrun')
        filter_blast_seqs = blast.get_new_seqs(query_seq, txid,
                                               '{}/tmp'.format(self.config.workdir),
                                               'filterrun', self.config)
        filtered_seqs = self.remove_mean_sd_nonfit(filter_blast_seqs, exist_dict)
        # self.upd_new_seqs = pd.concat([filtered_seqs, self.upd_new_seqs], axis=0, ignore_index=True,
        #                               sort=True)
        return filtered_seqs


    def remove_mean_sd_nonfit(self, filter_blast_seqs, table_otu_dict):
        # remove non fitting entries - too far from blast query seq....
        mean_sd = calculate_mean_sd(filter_blast_seqs['bitscore'])
        subfilter_blast_seqs = pd.DataFrame()
        for idx in filter_blast_seqs.index:
            seq_bitscore = filter_blast_seqs.loc[idx, 'bitscore']
            if (seq_bitscore >= mean_sd['mean'] - mean_sd['sd']) & \
                    (seq_bitscore <= mean_sd['mean'] + mean_sd['sd']):
                subfilter_blast_seqs = subfilter_blast_seqs.append(filter_blast_seqs.loc[idx])
            else:
                print('filter thresh too large')
        amnt_missing_seq = self.config.threshold - len(table_otu_dict)
        if amnt_missing_seq < len(subfilter_blast_seqs):
            print('subfilter-sample')
            filtered_seqs = subfilter_blast_seqs.sample(amnt_missing_seq)
        elif amnt_missing_seq >= len(subfilter_blast_seqs):
            print('subfilter-all')
            filtered_seqs = subfilter_blast_seqs
        else:
            print('should not happen')
            sys.exit(-2)
        return filtered_seqs


class FilterSeqIdent(Filter):
    def __init__(self, config, aln):
        super().__init__(config)
        self.aln = aln

    def filter(self, new_seqs):
        # print("FilterSeqIdent")
        to_add = self.upd_new_seqs
        for idx in set(new_seqs.index):
            # generate seq and compare lsit
            seqs_list = [self.aln[tax].symbols_as_string().replace("-", "").replace("N", "") for tax in self.aln]
            seq = new_seqs.loc[idx, 'sseq']
            if seq not in seqs_list:  # if we still want to add - might be obsolete
                to_add = to_add.append(new_seqs.loc[idx], ignore_index=True)
            else:
                print('same sequence ident')
                to_del = new_seqs.loc[idx]
                self.del_table = self.del_table.append(to_del, ignore_index=True)
        self.del_table['status'] = 'deleted - same seq'
        self.upd_new_seqs = self.upd_new_seqs.append(to_add, ignore_index=True)
        lfd = "{}/logfile".format(self.config.workdir)
        with open(lfd, "a") as log:
            log.write("Filter FilterBLASTThreshold has lowered the new seqs df from {} to {}.\n".format(len(new_seqs),
                                                                                                        len(to_add)))

class FilterMRCA(Filter):
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

        for txid in set(new_seqs['ncbi_txid']):
            tax_id = int(txid)
            mrca_tx = self.ncbi_parser.match_id_to_mrca(tax_id, self.mrca)
            if mrca_tx == self.mrca:
                select_TF = new_seqs['ncbi_txid'] == tax_id
                to_add = new_seqs[select_TF]
                self.upd_new_seqs = pd.concat([self.upd_new_seqs, to_add], axis=0, ignore_index=True, sort=True)
            else:
                print('out of mrca')
                select_TF = new_seqs['ncbi_txid'] == tax_id
                to_del = new_seqs[select_TF]
                to_del['status'] = 'deleted - mrca'
                self.del_table = pd.concat([self.del_table, to_del], axis=0, ignore_index=True, sort=True)
        # todo add assert that rows of upd_new_seq not in del_table

        lfd = "{}/logfile".format(self.config.workdir)
        with open(lfd, "a") as log:
            log.write("Filter FilterBLASTThreshold has lowered the new seqs df from {} to {}.\n".format(len(new_seqs),
                                                                                                        len(to_add)))


class FilterBLASTThreshold(Filter):
    def __init__(self, config):
        super().__init__(config)

    def filter(self, new_seqs):
        print("FilterBLASTThreshold")
        TF_eval = new_seqs['evalue'] < float(self.config.e_value_thresh)
        upd_new_seqs = new_seqs[TF_eval == True]
        df = pd.concat([upd_new_seqs, self.upd_new_seqs], sort=True)
        new_seqs_unique_nd = df.drop_duplicates(subset=['accession'], keep='first')
        deltab = new_seqs[TF_eval != True]
        deltab['status'] = 'deleted - evalue'
        self.del_table = deltab
        self.upd_new_seqs = new_seqs_unique_nd
        lfd = "{}/logfile".format(self.config.workdir)
        with open(lfd, "a") as log:
            log.write("Filter FilterBLASTThreshold has lowered the new seqs df from {} to {}.\n".format(len(new_seqs),
                                                                                                   len(new_seqs_unique_nd)))


class FilterUniqueAcc(Filter):
    def __init__(self, config):
        super().__init__(config)

    def filter(self, new_seqs):
        # delete things in table
        new_seqs_drop = new_seqs.drop_duplicates(subset=['accession'])
        new_seqs_unique = self.upd_new_seqs.drop_duplicates(subset=['accession'], keep='first')
        if len(new_seqs_drop) > 0 and len(new_seqs_unique) > 0:
            df = pd.concat([new_seqs_drop, new_seqs_unique], sort=True)
            new_seqs_unique_nd = df.drop_duplicates(subset=['accession'], keep='first')  # inplace=True replaces df
            # Select all duplicate rows based on one column
            del_seq = df[df.duplicated(['accession'])]
            new_seqs_unique_fn = new_seqs_unique_nd
        else:
            new_seqs_unique_fn = new_seqs_drop
            dropped = new_seqs.duplicated(subset=['accession'])
            del_seq = pd.concat([self.del_table, new_seqs[dropped == True]], sort=True)
        self.upd_new_seqs = new_seqs_unique_fn
        self.del_table = del_seq
        assert len(del_seq) + len(new_seqs_unique_fn) == len(new_seqs), \
            (len(del_seq), len(new_seqs_unique_fn), len(new_seqs))
        lfd = "{}/logfile".format(self.config.workdir)
        with open(lfd, "a") as log:
            log.write("Filter FilterUniqueAcc has lowered the new seqs df from {} to {}.\n".format(len(new_seqs),
                                                                                                len(new_seqs_unique_fn)))


class FilterLength(Filter):
    def __init__(self, config, aln):
        super().__init__(config)
        self.aln = aln

    def avg_length_seqs_in_aln(self):

        orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in self.aln]
        avg_seqlen = sum(orig_seqlen) / len(orig_seqlen)
        return avg_seqlen

    def filter(self, new_seqs):
        avg_len = self.avg_length_seqs_in_aln()
        del_seq = pd.DataFrame()
        filter_new_seqs = pd.DataFrame()
        for idx in new_seqs.index:
            seq = str(new_seqs.loc[idx, 'sseq'])
            print(new_seqs.loc[idx, 'accession'])
            print(new_seqs.loc[idx])
            assert '-' not in seq
            if self.config.maxlen * avg_len > len(seq) > self.config.minlen * avg_len:
                filter_new_seqs = filter_new_seqs.append(new_seqs.loc[idx], ignore_index=True)
            else:
                print('too long/too short')
                del_seq = del_seq.append(new_seqs.loc[idx], ignore_index=True)
                del_seq.loc[idx, 'status'] = 'deleted - seq len wrong: {}'.format(len(seq))
                print(del_seq)
        assert len(del_seq) + len(filter_new_seqs) == len(new_seqs), (len(del_seq), len(filter_new_seqs), len(new_seqs))
        self.upd_new_seqs = pd.concat([self.upd_new_seqs, filter_new_seqs], sort=True)
        self.del_table = pd.concat([self.del_table, del_seq], sort=True)
        print(self.del_table)
        lfd = "{}/logfile".format(self.config.workdir)
        with open(lfd, "a") as log:
            log.write("Filter FilterLength has lowered the new seqs df from {} to {}.\n".format(len(new_seqs), len(filter_new_seqs)))
