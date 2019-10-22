"""uses ncbi databases to easily retrieve taxonomic information.

data is provided by ncbi - downloaded via db_updater

parts are altered/copied from https://github.com/zyxue/ncbitax2lin/blob/master/ncbitax2lin.py, I indicated those on top
of the methods. Those functions have a MIT license.


"""


import os
import sys
import pandas as pd
import logging
import re
import multiprocessing
import datetime

from . import debug, get_user_input

nodes = None
names = None


# extodo: make_lineage_table needs to be rerun, when nodes and names are updated: is now part of db_updater

def strip(str_):
    """ Strips of blank characters from string in pd dataframe.
    """
    return str_.strip()


def load_nodes(nodes_file):
    """ Loads nodes.dmp and converts it into a pandas.DataFrame.
    Contains the information about the taxonomic hierarchy of names.

    :param nodes_file: relative path to the files 'nodes'
    :return: pandas dataframe with the data
    """
    # print(nodes_file)
    assert os.path.exists(nodes_file), ("file `%s` does not exist. Make sure you downloaded the "
                                        "databases from ncbi." % nodes_file)
    col_names = ["tax_id",
                 "parent_tax_id",
                 "rank",
                 "embl_code",
                 "division_id",
                 "inherited_div_flag",
                 "genetic_code_id",
                 "inherited_GC__flag",
                 "mitochondrial_genetic_code_id",
                 "inherited_MGC_flag",
                 "GenBank_hidden_flag",
                 "hidden_subtree_root_flag",
                 "comments"]
    df = pd.read_csv(nodes_file, sep="|", header=None, index_col=False, names=col_names)
    # To get rid of flanking tab characters
    df["rank"] = df["rank"].apply(strip)
    df["embl_code"] = df["embl_code"].apply(strip)
    df["comments"] = df["comments"].apply(strip)
    return df


def load_names(names_file):
    """
    Loads names.dmp and converts it into a pandas.DataFrame.
    Includes only names which are accepted as scientific name by ncbi.

    :param names_file: relative path to the files 'names'
    :return: pandas dataframe with the data
    """
    assert os.path.exists(names_file), ("file `%s` does not exist. Make sure you downloaded the "
                                        "databases from ncbi." % names_file)
    df = make_clean_df(names_file)
    sci_df = df[df["name_class"] == "scientific name"]
    sci_df.reset_index(drop=True, inplace=True)
    return sci_df


def load_synonyms(names_file):
    """
    Loads names.dmp and converts it into a pandas.DataFrame.
        Includes only names which are viewed as synonym by ncbi.

    :param names_file: relative path to the files 'names'
    :return: pandas dataframe with the data
    """
    assert os.path.exists(names_file), ("file `%s` does not exist. Make sure you downloaded "
                                        "the databases from ncbi." % names_file)
    # print("load synonyms")
    df = make_clean_df(names_file)
    sci_df = df[df["name_class"] == "synonym"]
    sci_df.reset_index(drop=True, inplace=True)
    return sci_df


def make_clean_df(names_file):
    df = pd.read_csv(names_file, sep="|", header=None, index_col=False,
                     names=["tax_id", "name_txt", "unique_name", "name_class"])
    df["name_txt"] = df["name_txt"].apply(strip)
    df["unique_name"] = df["unique_name"].apply(strip)
    df["name_class"] = df["name_class"].apply(strip)
    return df


def extract_ncbi_data():
    os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./data/")
    os.system("gunzip -f -cd ./data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
    os.system("mv nodes.dmp ./data/")
    os.system("mv names.dmp ./data/")
    os.system("rm ./data/taxdump.tar.gz")


class Parser:
    """Reads in databases from ncbi to connect species names with the taxonomic identifier
    and the corresponding hierarchical information. It provides a fast way to get those information.

    Nodes includes the hierarchical information, names the scientific names and ID's.
    The files need to be updated regularly, best way to always do it when a new blast database was loaded.
    """
    def __init__(self, names_file, nodes_file, interactive=False):
        self.names_file = names_file
        self.nodes_file = nodes_file
        # self.initialize()
        if interactive is True:
            self._download_ncbi_parser()

    def _download_ncbi_parser(self):
        """Check if files are present and if they are up to date.
        If not files will be downloaded.
        """
        print(self.nodes_file)
        if not os.path.isfile(self.names_file):
            print("Do you want to download taxonomy databases from ncbi? Note: This is a US government website! "
                  "You agree to their terms. \n")
            x = get_user_input()
            if x == "yes":
                extract_ncbi_data()
            else:
                sys.stderr.write(
                    "You have no taxonomic database, which is needed to run PhyFilter. Restart and type 'yes'. \n")
                sys.exit(-10)
        else:
            download_date = os.path.getmtime(self.nodes_file)
            download_date = datetime.datetime.fromtimestamp(download_date)
            today = datetime.datetime.now()
            time_passed = (today - download_date).days
            if time_passed >= 30:
                print("Do you want to update taxonomy databases from ncbi? Note: This is a US government website! "
                      "You agree to their terms. \n")
                x = get_user_input()
                if x == "yes":
                    extract_ncbi_data()
                elif x == "no":
                    print("You did not agree to update data from ncbi. Old database files will be used.")
                else:
                    print("You did not type yes or no!")
        if x == 'yes':
            make_lineage_table()

    def initialize(self):
        """
         The data itself are not stored in __init__. Instead every time the function is loaded
        if it has not yet been run during a run.

        It was important when I used pickle. Might therefore not be needed anymore.
        """
        print("Initialize NODES and NAMES!!")
        global nodes
        nodes = load_nodes(self.nodes_file)
        global names
        names = load_names(self.names_file)
        global synonyms
        synonyms = load_synonyms(self.names_file)

    def get_mrca(self, taxid_set):
        """
        Find mrca for a set of taxon ids.

        Finds the most recent common ancestor of  all provided ncbi taxon ids.

        :param taxid_set: set of ncbi taxon ids
        :return: mrca
        """
        print('get mrca from taxid_set')
        if nodes is None:
            self.initialize()
        id_dict = dict()
        for item in taxid_set:
            id_list = list()
            parent_id = int(nodes[nodes["tax_id"] == item]["parent_tax_id"].values[0])
            while parent_id != 131567:
                id_list.append(parent_id)
                parent_id = int(nodes[nodes["tax_id"] == parent_id]["parent_tax_id"].values[0])
            id_dict[item] = id_list
        count = 0
        for taxid in next(iter(id_dict.values())):  # take random entry as comparison
            for id2 in id_dict:
                idl = id_dict[id2]
                if taxid in idl:
                    count += 1
                    if count == len(id_dict.keys()):
                        return taxid
                else:
                    continue

    def taxid_is_valid(self, tax_id):
        """
        Checks if input taxon id is known by ncbi.

        :param tax_id:  ncbi taxon id
        :return: True or False
        """
        if nodes is None:
            self.initialize()
        if type(tax_id) != int:
            # sys.stdout.write("WARNING: mrca_id {} is no integer. Will convert value to int\n".format(mrca_id))
            tax_id = int(tax_id)
        try:
            # can be any value, just checking if val is known
            int(nodes[nodes["tax_id"] == tax_id]["parent_tax_id"].values[0])
            return True
        except:
            return False

    def get_rank(self, tax_id):
        """
        Get rank for given ncbi tax id.

        :param tax_id: ncbi taxon id
        :return: rank
        """
        if nodes is None:
            self.initialize()
        rank = nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
        return rank

    def get_downtorank_id(self, tax_id, downtorank=None):
        """
        Recursive function to find the parent id of a taxon as defined by downtorank.

        Returns the taxonomic id of the specified rank to the taxon id provided,
        e.g. user gives a taxon id of a species and wants to find the corresponding family id.

        :param tax_id: ncbi taxon id
        :param downtorank: rank provided as string, e.g. 'species'. Rank must be known by ncbi,
        :return: rank id as defined by downtorank
        """
        # debug("get downtorank")
        # print(tax_id)
        if downtorank is None:
            sys.stderr.write("You try to get downtorank without supplying a rank.")
            return tax_id
        if downtorank == 'no rank':
            sys.stderr.write("Cannot provide an id of a given taxon id corresponding to 'no rank'.")
            sys.exit(-5)
        id_known = self.taxid_is_valid(tax_id)
        if id_known is False:
            sys.stderr.write("Taxon id is not known by nodes. Probably to new.\n.")
            sys.exit(-6)

        if nodes is None:
            self.initialize()
        if type(tax_id) != int:
            tax_id = int(tax_id)
        # following statement is to get id of taxa if taxa is higher ranked than specified
        rank = nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
        # next one is used to return id if rank of taxid was always higher than species
        # todo: should be working with any kind rank
        if rank != "species":
            if downtorank == "species":
                if (nodes[nodes["tax_id"] == tax_id]["rank"].values[0] != "varietas"
                        and nodes[nodes["tax_id"] == tax_id]["rank"].values[0] != "subspecies"):
                    return tax_id
        if nodes[nodes["tax_id"] == tax_id]["rank"].values[0] == downtorank:
            return tax_id
        elif nodes[nodes["tax_id"] == tax_id]["rank"].values[0] == "superkingdom":
            tax_id = 0
            return tax_id
        else:
            parent_id = int(nodes[nodes["tax_id"] == tax_id]["parent_tax_id"].values[0])
            return self.get_downtorank_id(parent_id, downtorank)

    def check_mrca_input(self, mrca_id):
        if type(mrca_id) == set:
            mrca_id_set = set()
            for item in mrca_id:
                id_known = self.taxid_is_valid(item)
                if id_known is True:
                    mrca_id_set.add(item)
                else:
                    sys.stderr.write("mrca id {}  is not known by nodes. Probably to new."
                                     " ID will be removed from mrca set.\n.".format(item))
        else:
            id_known = self.taxid_is_valid(mrca_id)
            if id_known is False:
                sys.stderr.write("mrca id {} is not known by nodes. Probably to new.\n.".format(mrca_id))
                sys.exit(-6)
            # do mrca format check
        if type(mrca_id) == set:
            if len(mrca_id) == 1:
                mrca_id = next(iter(mrca_id))
                mrca_id = int(mrca_id)
        elif type(mrca_id) != int:
            mrca_id = int(mrca_id)
        return mrca_id

    def match_id_to_mrca(self, tax_id, mrca_id):
        """
        Recursive function to find out if tax_id is part of mrca_id.

        :param tax_id: ncbi taxon id of otu that should be checked
        :param mrca_id: provided mrca ncbi taxon id
        :return: taxon id (= mrca id) if matches, if not 0
        """
        # todo rewrite to return true/false
        debug("match_id_to_mrca")
        if nodes is None:
            self.initialize()
        # do id testing
        id_known = self.taxid_is_valid(tax_id)
        if id_known is False:
            sys.stderr.write("Taxon id {} is not known by nodes. Probably to new.\n.".format(tax_id))
            return -500
        if type(tax_id) != int:
            tax_id = int(tax_id)
        self.check_mrca_input(mrca_id)

        # stop with weird tax ids...
        if tax_id in [81077, 28384, 131567, 1, 0]:  # other sequences/artificial sequences
            debug("artifical")
            tax_id = 0
            return tax_id

        rank_tax = nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
        if tax_id == mrca_id:
            debug("found right rank")
            return tax_id
        elif type(mrca_id) == set and tax_id in mrca_id:
            debug("found right rank")
            return tax_id
        elif rank_tax == "superkingdom":
            debug("superkingdom")
            tax_id = 0
            return tax_id
        else:
            debug("parent")
            parent_id = int(nodes[nodes["tax_id"] == tax_id]["parent_tax_id"].values[0])
            debug(parent_id)
            return self.match_id_to_mrca(parent_id, mrca_id)

    def get_name_from_id(self, tax_id):
        """
        Find the scientific name for a given ID.

        :param tax_id: ncbi taxon id to be checked
        :return: taxon name
        """
        # debug("get name from id")
        if names is None:
            self.initialize()
        if type(tax_id) != int:
            # sys.stdout.write("WARNING: mrca_id {} is no integer. Will convert value to int\n".format(mrca_id))
            tax_id = int(tax_id)
        try:
            if tax_id == 0:
                tax_name = "unidentified"
            else:
                tax_name = names[names["tax_id"] == tax_id]["name_txt"]
                tax_name = tax_name.values[0].replace(" ", "_")
                tax_name = tax_name.strip()
        except IndexError:
            sys.stdout.write("tax_id {} unknown by ncbi_parser files (names.dmp)\n".format(tax_id))
            tax_name = "unknown_{}".format(tax_id)
            if os.path.exists("ncbi_id_unknown.err"):
                fn = open("ncbi_id_unknown.err", "a")
                fn.write("{}".format(tax_id))
                fn.close()
            else:
                fn = open("ncbi_id_unknown.err", "w")
                fn.write("{}".format(tax_id))
                fn.close()
        return tax_name

    def get_id_from_name(self, tax_name):
        """ Find the ID for a given taxonomic name.
        """
        if names is None:
            self.initialize()
        org_tax = tax_name
        tax_name = tax_name.replace("_", " ")
        if len(tax_name.split(" ")) >= 2:
            if tax_name.split(" ")[1] == "sp.":
                tax_name = "{}".format(tax_name.split(" ")[0])
        try:
            tax_id = names[names["name_txt"] == tax_name]["tax_id"].values[0]
        except IndexError:
            debug(tax_name)
            debug(names[names["name_txt"] == tax_name])
            if len(tax_name.split(" ")) == 3:
                tax_name = "{} {}-{}".format(
                    tax_name.split(" ")[0],
                    tax_name.split(" ")[1],
                    tax_name.split(" ")[2])
                tax_id = names[names["name_txt"] == tax_name]["tax_id"].values[0]
                sys.stdout.write("tax_name {} unknown, modified to {} worked.\n".format(org_tax, tax_name))
            else:
                sys.stdout.write("Are you sure, its an accepted name and not a synonym: {}? "
                                 "I look in the synonym table now.\n".format(tax_name))
                tax_id = self.get_id_from_synonym(tax_name)
        tax_id = int(tax_id)
        return tax_id

    def get_id_from_synonym(self, tax_name):
        """ Find the ID for a given taxonomic name, which is not an accepted name.
        """
        if names is None:
            self.initialize()
        tax_name = tax_name.replace("_", " ")
        try:
            tax_id = synonyms[synonyms["name_txt"] == tax_name]["tax_id"].values[0]
        except IndexError:
            if len(tax_name.split(" ")) == 3:
                tax_name = "{} {}-{}".format(
                    tax_name.split(" ")[0],
                    tax_name.split(" ")[1],
                    tax_name.split(" ")[2])
                tax_id = names[names["name_txt"] == tax_name]["tax_id"].values[0]
            else:
                sys.stderr.write("ncbi taxon name unknown by parser files: {}, taxid set to 0.\n".format(tax_name))
                tax_id = 0
                if os.path.exists("ncbi_name_unknown.err"):
                    fn = open("ncbi_name_unknown.err", "a")
                    fn.write("{}".format(tax_name))
                    fn.close()
                else:
                    fn = open("ncbi_name_unknown.err", "w")
                    fn.write("{}".format(tax_name))
                    fn.close()
        return tax_id

    def get_lower_from_id(self, tax_id):
        """
        Find all lower taxon ids for given id.

        Takes ages, as it loops through almost all ids provided by ncbi.
        :param tax_id:
        :return:
        """
        if names is None:
            self.initialize()
        rank_list = ['tax_id', 'superkingdom', 'phylum', 'class', 'order', 'family', 'tribe', 'genus', 'species']
        rank_mrca = nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
        spn = self.get_name_from_id(tax_id)
        if rank_mrca in rank_list:  # corresponds to what is written into the lineage files.
            df = lineages_to_df(rank_mrca, spn)
        else:
            while rank_mrca not in rank_list:
                parent_id = int(nodes[nodes["tax_id"] == tax_id]["parent_tax_id"].values[0])
                rank_mrca = nodes[nodes["tax_id"] == parent_id]["rank"].values[0]
                if rank_mrca in rank_list:
                    df = lineages_to_df(rank_mrca, spn)
        return df['tax_id'].to_list()


def lineages_to_df(rank, name):
    """
    Returns pandas df subset of the taxon ids that correspond to the provided rank-name.
    E.g. 'family' and 'Asteraceae' will return all entries that exists for Asteraceae.

    :param rank: as known by ncbi taxonomy - not 'no rank'.
    :param name:
    :return:
    """
    if rank == 'no rank':
        sys.stderr.write("Cannot provide an id of a given taxon id corresponding to 'no rank'.")
        sys.exit(-5)
    if not os.path.exists(os.path.join('./data/lineages_cols.csv')):
        make_lineage_table()
    lin = pd.read_csv(os.path.join('./data/lineages_cols.csv'))
    subset = lin[lin[rank] == name]
    return subset


# copy from ncbitax2lin
def make_lineage_table():
    """
    Produces the table needed for get_lower_from_id().
    :return:
    """
    # data downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy/
    nodes_df = load_nodes('./data/nodes.dmp')
    names_df = load_names('./data/names.dmp')
    df = nodes_df.merge(names_df, on='tax_id')
    df = df[['tax_id', 'parent_tax_id', 'rank', 'name_txt']]
    df.reset_index(drop=True, inplace=True)
    # # log summary info about the dataframe
    # print('=' * 50)
    # df.info(verbose=True, memory_usage="deep")
    # print('=' * 50)

    global TAXONOMY_DICT
    TAXONOMY_DICT = dict(zip(df.tax_id.values, df.to_dict('records')))
    ncpus = multiprocessing.cpu_count()
    print(ncpus)
    print('map - takes time and a lot of memory... (multiprocessing)')
    pool = multiprocessing.Pool(ncpus)
    # take about 18G memory
    lineages_dd = pool.map(find_lineage, df.tax_id.values)
    pool.close()
    dd_for_df = dict(zip(range(len(lineages_dd)), lineages_dd))
    print('make pandas df')
    lineages_df = pd.DataFrame.from_dict(dd_for_df, orient='index')
    lineages_df.sort_values('tax_id', inplace=True)
    cols = ['tax_id',
            'superkingdom',
            'phylum',
            'class',
            'order',
            'family',
            'tribe',
            'genus',
            'species']
    lineages_df.to_csv(os.path.join('./data/lineages_cols.csv'), index=False, columns=cols)

    # zipping does not work, byte issue
    # util.backup_file(lineages_csv_output)
    # print('write zip')
    # lineages_csv_output = os.path.join('{0}.csv.gz'.format(output))
    # with open(lineages_csv_output, 'wb') as opf:
    #     # make sure the name and timestamp are not gzipped, (like gzip -n)
    #     opf_gz = gzip.GzipFile(
    #         filename='',        # empty string because fileobj is given
    #         mode='wb',          # wb doesn't seem to work sometimes
    #         compresslevel=9,
    #         fileobj=opf,
    #         mtime=0.   # an optional numeric timestamp, set to be deterministic
    #     )
    #     cols = ['tax_id',
    #             'superkingdom',
    #             'phylum',
    #             'class',
    #             'order',
    #             'family',
    #             'tribe'
    #             'genus',
    #             'species']
    #     other_cols = sorted([__ for __ in lineages_df.columns
    #                          if __ not in cols])
    #     output_cols = cols + other_cols
    #     output_cols_byte = []
    #     for item in output_cols:
    #         item = item.encode()
    #         output_cols_byte.append(item)
    #     lineages_df.to_csv(opf_gz, index=False, columns=output_cols_byte)
    #     opf_gz.close()


# copy from ncbitax2lin
def to_dict(lineage):
    """
    convert the lineage into a list of tuples in the form of

    [
        (tax_id1, rank1, name_txt1),
        (tax_id2, rank2, name_txt2),
        ...
    ]

    to a dict
    """
    dd = {}
    num_re = re.compile('[0-9]+')
    len_lineage = len(lineage)
    for k, __ in enumerate(lineage):
        tax_id, rank, name_txt = __
        # use the last rank as the tax_id, whatever it is, genus or species.
        if k == len_lineage - 1:
            dd['tax_id'] = tax_id

        # e.g. there could be multiple 'no rank'
        numbered_rank = rank
        while numbered_rank in dd:
            # print __, numbered_rank
            search = num_re.search(numbered_rank)
            if search is None:
                count = 1
            else:
                count = int(search.group()) + 1
            numbered_rank = '{0}{1}'.format(rank, count)
        dd[numbered_rank] = name_txt
    return dd


# copy from ncbitax2lin
def find_lineage(tax_id):
    """

    :param tax_id:
    :return:
    """
    if tax_id % 50000 == 0:
        logging.debug('working on tax_id: {0}'.format(tax_id))
    lineage = []
    while True:
        rec = TAXONOMY_DICT[tax_id]
        lineage.append((rec['tax_id'], rec['rank'], rec['name_txt']))
        tax_id = rec['parent_tax_id']
        if tax_id == 1:
            break
    # reverse results in lineage of Kingdom => species, this is helpful for
    # to_dict when there are multiple "no rank"s
    lineage.reverse()
    return to_dict(lineage)
