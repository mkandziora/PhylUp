"""
PhylUp: automatically update alignments.
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

All rights reserved. No warranty, explicit or implicit, provided. No distribution or modification of code allowed.
All classes and methods will be distributed under an open license in the near future.

Package to automatically update alignments and phylogenies using local sequences or a  local Genbank database.

Parts are inspired by the program physcraper developed by me and Emily Jane McTavish.
"""

# user settings class

import os
import sys
import configparser
import datetime

from . import db_updater
from . import debug
import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser


# TODO: make global blast folders that can be shared across runs


def is_number(s):
    """test if string can be coerced to float"""
    try:
        float(s)
        return True
    except ValueError:
        return False


class ConfigObj(object):
    def __init__(self, configfi, workdir, interactive=True):
        """
        Build a configuration class.

        During the initializing process the following self objects are generated:
            * self.workdir**: working directory
            * **self.e_value_thresh**: the defined threshold for the e-value during Blast searches, check out:
                https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
            * **self.hitlist_size**: the maximum number of sequences retrieved by a single blast search
            * **self.blast_loc**: defines which blasting method to use. Currently only local is implemented
                * if self.blast_loc == "local":
                    * self.blastdb: this defines the path to the local blast database
                    * self.ncbi_parser_nodes_fn: path to 'nodes.dmp' file, that contains the hierarchical information
                    * self.ncbi_parser_names_fn: path to 'names.dmp' file, that contains the different ID's
            * **self.num_threads**: number of cores to be used during a run
            * **self.minlen**: value from 0 to 1. Defines how much shorter new seq can be compared to input
            * **self.trim_perc**: value that determines how many seq need to be present before the beginning and end
                                    of alignment will be trimmed
            * **self.maxlen**: max length for values to add to aln
            * **self.different_level**: T/F, used to first filter for higher rank and then to a lower rank.
        * **self.filtertype**: either 'blast' or 'length'; filter number per otu by blast or choose longest
        * **self.backbone**: T/F; calculate complete tree or add to backbone
        * **self.update_tree**: T/F; update tree (T) or just update alignment
        * **self.threshold**: number of seq per otu
        * **self.downtorank**: filter number otu to specific rank
        * **self.unpublished**: T/F; add unpublished local sequences
        * **self.modeltest_criteria** BIC, AIC or AICc

        * **interactive**: T/F; checks if databases need to be updated
        * **self.logfile**: file location where some information during the run is logged
        * **self.fix_blast**: T/F; use same blast folder across runs

        :param configfi: a configuration file in a specific format. The file needs to have a heading of the format:
                        [blast] and then somewhere below that heading as string, e.g. e_value_thresh = value
        :param workdir: the working directory
        :param interactive: defaults to True, is used to interactively update the local blast databases
        """
        debug("Building config object\n")
        assert os.path.isfile(configfi), "file `%s` does not exists" % configfi

        self.workdir = workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        if not os.path.exists(os.path.join(self.workdir, "blast")):
            os.makedirs(os.path.join(self.workdir, "blast"))
        assert self.workdir

        config = configparser.ConfigParser()
        config.read_file(open(configfi))

        # general settings
        self.num_threads = int(config["general"]["num_threads"])
        # print(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))
        if os.environ.get('SLURM_JOB_CPUS_PER_NODE'):
            self.num_threads = int(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))
        sys.stdout.write('Workflow runs with {} threads.\n'.format(self.num_threads))

        self.mrca_input = config["general"]["mrca"]

        # unpublished settings
        self.unpublished = config["unpublished"]['unpublished']
        if self.unpublished == "True" or self.unpublished == "true":
            self.unpublished = True
            self.unpubl_data = config["unpublished"]['unpubl_data']
            self.unpubl_names = config["unpublished"]['unpubl_names']
            self.perpetual = config["unpublished"]['perpetual']
            self.blast_all = config["unpublished"]['blast_all']
            if self.blast_all == "True" or self.blast_all == "true":
                self.blast_all = True
            else:
                self.blast_all = False
        else:
            self.unpublished = False

        # read in blast settings
        self.blastdb = config["blast"]["localblastdb"]
        assert os.path.exists(self.blastdb), self.blastdb

        self.e_value_thresh = config["blast"]["e_value_thresh"]
        assert is_number(self.e_value_thresh), ("e_value_thresh is not defined as a number: {}.\n".format(self.e_value_thresh))

        self.hitlist_size = int(config["blast"]["hitlist_size"])
        assert is_number(self.hitlist_size), ("Hitlist size is not defined as a number: {}.\n".format(self.hitlist_size))

        self.fix_blast = config["blast"]["fix_blast_result_folder"]
        if self.fix_blast == "True" or self.fix_blast == "true":
            self.fix_blast = True
            self.blast_folder = os.path.abspath("./data/blast")
            if not os.path.exists(self.blast_folder):
                os.mkdir(self.blast_folder)
            sys.stdout.write(
                "You are using the same blast folder across runs ({}) - be careful. Make sure it is the same locus "
                "and that you did not change your blast settings.\n").format(self.blast_folder)
        else:
            self.fix_blast = False
            self.blast_folder = os.path.abspath(os.path.join(self.workdir, "blast"))

        # INTERNAL PhylUp SETTINGS
        self.ncbi_parser_nodes_fn = config["ncbi_parser"]["nodes_fn"]
        self.ncbi_parser_names_fn = config["ncbi_parser"]["names_fn"]
        # # TODO: test if this would make the whole thing faster, db building takes long, but likely blasting is faster
        # # make smaller blast database
        # name = self.make_db_from_taxid(self.mrca)
        # self.blastdb = '{}/{}_db'.format(self.blastdb, name)

        # read in alignment settings
        self.minlen = float(config["phylup"]["min_len"])
        assert 0 < self.minlen <= 1, ("Min len is not between 0 and 1: {}.\n".format(self.minlen))
        self.maxlen = float(config["phylup"]["max_len"])
        assert 1 < self.maxlen, ("Max len is not larger than 1: {}.\n".format(self.maxlen))
        self.trim_perc = float(config["phylup"]["trim_perc"])
        assert 0 < self.trim_perc < 1, ("Percentage for trimming is not between 0 and 1: {}.\n".format(self.trim_perc))

        # read in filter settings
        self.filtertype = config["filter"]["filtertype"]
        assert self.filtertype in ['blast', 'length'], ("self.filtertype `%s` "
                                                        "is not 'blast' or 'length'" % self.filtertype)
        self.threshold = int(config["filter"]["threshold"])
        self.downtorank = config["filter"]["downtorank"]

        self.different_level = config["filter"]["different_level"]
        if self.different_level == "True" or self.different_level == "true":
            self.different_level = True
        else:
            self.different_level = False
        assert self.different_level in [True, False], ("self.different_level `%s` "
                                                       "is not True or False" % self.different_level)

        # read in tree calculation settings
        self.update_tree = config["tree"]["update_tree"]
        if self.update_tree == "True" or self.update_tree == "true":
            self.update_tree = True
        else:
            self.update_tree = False
        self.backbone = config["tree"]["backbone"]
        if self.backbone == "True" or self.backbone == "true":
            self.backbone = True
        else:
            self.backbone = False

        self.modeltest_criteria = config['tree']['modeltest_criteria']
        assert self.modeltest_criteria in ['AIC', 'AICc', 'BIC'], ("self.modeltest_criteria `%s` "
                                                        "is not AIC, AICc or BIC" % self.modeltest_criteria)


        ####
        # check database status
        # print('Status for interactive is {}'.format(interactive))
        if interactive is None:
            interactive = sys.stdin.isatty()
        if interactive is False:
            sys.stdout.write("REMEMBER TO UPDATE THE NCBI DATABASES REGULARLY! ")
            download_date = os.path.getmtime("{}/nt_v5.60.nhr".format(self.blastdb))
            download_date = datetime.datetime.fromtimestamp(download_date)
            today = datetime.datetime.now()
            time_passed = (today - download_date).days
            sys.stdout.write("Looks like you last updated it {} days ago.\n".format(time_passed))
            sys.stdout.write("Run the file 'update_databases.py' from the data folder to automatically update. "
                             "Note, that you are accessing a US government website to do so.\n")
        if interactive is True:
            db_updater._download_localblastdb(self)
            ncbi_parser = ncbi_data_parser.Parser(names_file=self.ncbi_parser_names_fn,
                                                  nodes_file=self.ncbi_parser_nodes_fn)
            ncbi_parser._download_ncbi_parser()
            #db_updater._download_ncbi_parser(self)
        debug("check db file status?: {}".format(interactive))

        # ###########
        # internal settings
        self.logfile = os.path.join(self.workdir, "logfile")

    def make_db_from_taxid(self, taxid):
        """
        Generates a smaller Genbank database, that only contains sequences that are part of the mrca.

        :param taxid:
        :return:
        """
        ncbi_parser = ncbi_data_parser.Parser(names_file=self.ncbi_parser_names_fn,
                                              nodes_file=self.ncbi_parser_nodes_fn)
        if type(taxid) == list:
            mrca_id = ncbi_parser.get_mrca(taxid)
            taxid = mrca_id

        taxon_list = ncbi_parser.get_lower_from_id(taxid)
        taxonlist_fn = os.path.join(self.blastdb, "{}_idlist.txt".format(taxid))
        with open(taxonlist_fn) as fn:
            fn.write("\n".join(str(item) for item in taxon_list))

        cmd1 = "makeblastdb -dbtype nucl -parse_seqids -taxid_map {} " \
               "-out {}_db -title {}".format(taxonlist_fn, taxid, taxid)
        os.system(cmd1)
        return taxid
