"""uses ncbi databases to easily retrieve taxonomic information.

parts are altered from https://github.com/zyxue/ncbitax2lin/blob/master/ncbitax2lin.py
"""

import os
import sys
import pandas as pd
from . import debug

debug("Current ncbi_parser version number: 07102019.0")

nodes = None
names = None


def strip(str_):
    """ Strips of blank characters from string in pd dataframe.
    """
    return str_.strip()


def load_nodes(nodes_file):
    """ Loads nodes.dmp and converts it into a pandas.DataFrame.
    Contains the information about the taxonomic hierarchy of names.
    """
    # print(nodes_file)
    assert os.path.exists(nodes_file), ("file `%s` does not exist. Make sure you downloaded the "
        "databases from ncbi." % nodes_file )
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
    """ Loads names.dmp and converts it into a pandas.DataFrame.
    Includes only names which are accepted as scientific name by ncbi.
    """
    assert os.path.exists(names_file), ("file `%s` does not exist. Make sure you downloaded the "
        "databases from ncbi." % names_file)
    df = pd.read_csv(names_file, sep="|", header=None, index_col=False,
                     names=["tax_id", "name_txt", "unique_name", "name_class"])
    df["name_txt"] = df["name_txt"].apply(strip)
    df["unique_name"] = df["unique_name"].apply(strip)
    df["name_class"] = df["name_class"].apply(strip)
    sci_df = df[df["name_class"] == "scientific name"]
    sci_df.reset_index(drop=True, inplace=True)
    return sci_df


def load_synonyms(names_file):
    """Loads names.dmp and converts it into a pandas.DataFrame.
        Includes only names which are viewed as synonym by ncbi.
    """
    assert os.path.exists(names_file), ("file `%s` does not exist. Make sure you downloaded "
                                        "the databases from ncbi." % names_file)
    # print("load synonyms")
    df = pd.read_csv(names_file, sep="|", header=None, index_col=False,
                     names=["tax_id", "name_txt", "unique_name", "name_class"])
    df["name_txt"] = df["name_txt"].apply(strip)
    df["unique_name"] = df["unique_name"].apply(strip)
    df["name_class"] = df["name_class"].apply(strip)
    sci_df = df[df["name_class"] == "synonym"]
    sci_df.reset_index(drop=True, inplace=True)
    return sci_df


class Parser:
    """Reads in databases from ncbi to connect species names with the taxonomic identifier
    and the corresponding hierarchical information. It provides a fast way to get those information.

    Nodes includes the hierarchical information, names the scientific names and ID's.
    The files need to be updated regularly, best way to always do it when a new blast database was loaded.
    """

    def __init__(self, names_file, nodes_file):
        self.names_file = names_file
        self.nodes_file = nodes_file
        # self.initialize()

    def initialize(self):
        """ The data itself are not stored in __init__. Instead every time the function is loaded
        if it has not yet been run during a run - was important for pickle.
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
        Find mrca for a set of taxa.

        :param taxid_set:
        :return:
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
        """check if input number is known by ncbi.
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
        """ Get rank for given ncbi tax id.
        """
        if nodes is None:
            self.initialize()
        rank = nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
        return rank

    def get_downtorank_id(self, tax_id, downtorank="species"):
        """ Recursive function to find the parent id of a taxon as defined by downtorank.
        """
        # debug("get downtorank")
        if nodes is None:
            self.initialize()
        if type(tax_id) != int:
            # sys.stdout.write("WARNING: mrca_id {} is no integer. Will convert value to int\n".format(mrca_id))
            tax_id = int(tax_id)
        # following statement is to get id of taxa if taxa is higher ranked than specified
        try:
            rank = nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
        except IndexError:
            rank = None
            return rank
        if rank != "species":
            if downtorank == "species":
                if (
                    nodes[nodes["tax_id"] == tax_id]["rank"].values[0] != "varietas"
                    and nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
                    != "subspecies"
                ):
                    return tax_id
        if nodes[nodes["tax_id"] == tax_id]["rank"].values[0] == downtorank:
            return tax_id
        elif nodes[nodes["tax_id"] == tax_id]["rank"].values[0] == "superkingdom":
            tax_id = 0
            return tax_id
        else:
            parent_id = int(nodes[nodes["tax_id"] == tax_id]["parent_tax_id"].values[0])
            return self.get_downtorank_id(parent_id, downtorank)

    def match_id_to_mrca(self, tax_id, mrca_id):
        """ Recursive function to find out if tax_id is part of mrca_id.
        """
        debug("match_id_to_mrca")
        if tax_id in [81077, 28384]:  # other sequences/artificial sequences
            debug("artifical")
            tax_id = 0
            return tax_id

        if tax_id == 131567 or tax_id == 1 or tax_id == 0:  # cellular organism
            tax_id = 0
            debug("cellular or 0")
            return tax_id
       
        if nodes is None:
            self.initialize()
        if type(mrca_id) == set:
            if len(mrca_id) == 1:
                mrca_id = next(iter(mrca_id))
            else:
                mrca_id = self.get_mrca(mrca_id)
            mrca_id = int(mrca_id)
        elif type(mrca_id) != int:
            # sys.stdout.write("WARNING: mrca_id {} is no integer. Will convert value to int\n".format(mrca_id))
            mrca_id = int(mrca_id)
        if type(tax_id) != int:
            # sys.stdout.write("WARNING: mrca_id {} is no integer. Will convert value to int\n".format(mrca_id))
            tax_id = int(tax_id)
        debug([tax_id, mrca_id])
        rank_mrca_id = nodes[nodes["tax_id"] == mrca_id]["rank"].values[0]
        rank_tax_id = nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
        debug([rank_mrca_id, rank_tax_id])
        if tax_id == mrca_id:
            debug("found right rank")
            return tax_id
        elif rank_tax_id == "superkingdom":
            debug("superkingdom")
            tax_id = 0
            return tax_id
        else:
            debug("parent")
            parent_id = int(nodes[nodes["tax_id"] == tax_id]["parent_tax_id"].values[0])
            debug(parent_id)
            return self.match_id_to_mrca(parent_id, mrca_id)

    def get_name_from_id(self, tax_id):
        """ Find the scientific name for a given ID.
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
