# user settings class

import os
import sys
import configparser
from . import db_updater

_DEBUG_MK = 1


def debug(msg):
    """short debugging command
    """
    if _DEBUG_MK == 1:
        print(msg)


def is_number(s):
    """test if string can be coerced to float"""
    try:
        float(s)
        return True
    except ValueError:
        return False


class ConfigObj(object):
    """
    To build the class the following is needed:

      * **configfi**: a configuration file in a specific format, e.g. to read in self.e_value_thresh.

        The file needs to have a heading of the format: [blast] and then somewhere below that heading a string e_value_thresh = value

      * **interactive**: defaults to True, is used to interactively update the local blast databases

    During the initializing process the following self objects are generated:

      * **self.e_value_thresh**: the defined threshold for the e-value during Blast searches, check out: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
      * **self.hitlist_size**: the maximum number of sequences retrieved by a single blast search
      * **self.seq_len_perc**: value from 0 to 1. Defines how much shorter new seq can be compared to input
      * **self.trim_perc**: value that determines how many seq need to be present before the beginning and end of alignment will be trimmed
      * **self.maxlen**: max length for values to add to aln
      * **self.get_ncbi_taxonomy**: Path to sh file doing something...
      * **self.phylesystem_loc**: defines which phylesystem for OpenTree datastore is used - default: api, but can run on local version too.
      * **self.id_pickle**: path to pickle file
      * **self.email**: email address used for blast queries
      * **self.blast_loc**: defines which blasting method to use:

          * either web-query (=remote)
          * from a local blast database (=local)
      * **self.num_threads**: number of cores to be used during a run
      * **self.url_base**:

          * if blastloc == remote: it defines the url for the blast queries.
          * if blastloc == local: url_base = None
      * **self.unmapped**: used for OToL original tips that can not be assigned to a taxon

          * keep: keep the unmapped taxa and asign them to life
          * remove: remove the unmapped taxa from aln and tre
      * **self.delay**: defines when to reblast sequences in days
      * **self.add_lower_taxa**: T/F, enables to re-access formerly filtered seq by allowing them be passed into remove_identical. Used if we first filter for higher rank and then want to filter for a lower rank.
      * **optional self.objects**:

          * if blastloc == local:

              * self.blastdb: this defines the path to the local blast database
              * self.ncbi_parser_nodes_fn: path to 'nodes.dmp' file, that contains the hierarchical information
              * self.ncbi_parser_names_fn: path to 'names.dmp' file, that contains the different ID's
    """

    def __init__(self, configfi, workdir, interactive=None):
        debug(configfi)
        debug(os.path.isfile(configfi))
        assert os.path.isfile(configfi), "file `%s` does not exists" % configfi

        sys.stdout.write("Building config object\n")
        self.workdir = workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        if not os.path.exists("{}/blast".format(self.workdir)):
            os.makedirs("{}/blast".format(self.workdir))
        assert self.workdir
        config = configparser.ConfigParser()
        config.read_file(open(configfi))

        # read in blast settings
        self.email = config["blast"]["Entrez.email"]
        assert "@" in self.email, "your email `%s` does not have an @ sign" % self.email

        self.e_value_thresh = config["blast"]["e_value_thresh"]
        assert is_number(self.e_value_thresh), (
                "value `%s` does not exists" % self.e_value_thresh
        )
        self.hitlist_size = int(config["blast"]["hitlist_size"])
        assert is_number(self.hitlist_size), (
                "value `%s`is not a number" % self.e_value_thresh
        )
        self.blast_loc = config["blast"]["location"]
        assert self.blast_loc in ["local"], (
                "your blast location `%s` is not local" % self.blast_loc
        )
        if self.blast_loc == "local":
            self.blastdb = config["blast"]["localblastdb"]
            if not os.path.exists(self.blastdb):
                os.makedirs(self.blastdb)
            self.url_base = None
            self.ncbi_parser_nodes_fn = config["ncbi_parser"]["nodes_fn"]
            self.ncbi_parser_names_fn = config["ncbi_parser"]["names_fn"]
        if _DEBUG_MK:
            sys.stdout.write("{}\n".format(self.email))
            sys.stdout.write("{}\n".format(self.blast_loc))
            if self.blast_loc == "local":
                sys.stdout.write("local blast db {}\n".format(self.blastdb))
        self.num_threads = config["blast"].get("num_threads")
        debug("slurm threads")
        debug(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))
        if os.environ.get('SLURM_JOB_CPUS_PER_NODE'):
            self.num_threads = int(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))

        debug(self.num_threads)
        self.gb_id_filename = config["blast"].get("gb_id_filename", False)
        if self.gb_id_filename is not False:
            if self.gb_id_filename == "True" or self.gb_id_filename == "true":
                self.gb_id_filename = True
            else:
                self.gb_id_filename = False
        debug("shared blast folder? {}".format(self.gb_id_filename))
        self.delay = int(config["blast"]["delay"])
        assert is_number(self.delay), (
                "value `%s`is not a number" % self.delay
        )
        # #############
        # read in phy_filter settings
        self.minlen = float(config["phy_filter"]["min_len"])
        assert 0 < self.minlen <= 1, (
                "value `%s` is not between 0 and 1" % self.minlen
        )
        self.trim_perc = float(config["phy_filter"]["trim_perc"])
        assert 0 < self.trim_perc < 1, (
                "value `%s` is not between 0 and 1" % self.trim_perc
        )
        self.maxlen = float(config["phy_filter"]["max_len"])
        assert 1 < self.maxlen, (
                "value `%s` is not larger than 1" % self.maxlen
        )
        self.add_lower_taxa = config["phy_filter"]["add_lower_taxa"]
        if self.add_lower_taxa == "True" or self.add_lower_taxa == "true":
            self.add_lower_taxa = True
        else:
            self.add_lower_taxa = False
        assert self.add_lower_taxa in [True, False], (
                "self.add_lower_taxa `%s` is not True or False" % self.add_lower_taxa
        )

        self.filtertype = config["phy_filter"]["filtertype"]
        assert self.filtertype in ['blast', 'length'], (
                "self.filtertype `%s` is not 'blast' or 'length'" % self.filtertype
        )

        self.backbone = config["phy_filter"]["backbone"]
        if self.backbone == "True" or self.backbone == "true":
            self.backbone = True
        else:
            self.backbone = False

        if self.gb_id_filename == "True" or self.gb_id_filename == "true":
            self.gb_id_filename = True
        else:
            self.gb_id_filename = False

        self.threshold = int(config["phy_filter"]["threshold"])
        self.downtorank = config["phy_filter"]["downtorank"]

        # read in settings for internal phy_filter processes
        # default is api, but can run on local version of OpenTree datastore
        self.phylesystem_loc = config["phylesystem"]["location"]
        assert self.phylesystem_loc in ["local", "api"], \
            (
                "phylesystem location must be either local or api")
        # rewrites relative path to absolute path so that it behaves when changing dirs
        self.id_pickle = os.path.abspath(config["taxonomy"]["id_pickle"])  # TODO: unused self object
        ####
        # check database status
        print(interactive)
        # print(some)
        if interactive is None:
            interactive = sys.stdin.isatty()
            if interactive is False:
                sys.stdout.write("REMEMBER TO UPDATE THE NCBI DATABASES REGULARLY!!\n")
        if interactive is True:
            db_updater._download_ncbi_parser(self)
            db_updater._download_localblastdb(self)
        debug("check db file status?: {}".format(interactive))

        # ###########
        # internal settings
        self.logfile = "{}/logfile".format(self.workdir)

        # more settings
        self.settings = {}
        self.settings['blast_subdir'] = "{}/current_blast_run".format(self.workdir)


    def add_to_config(self, key, val):
        """Add more settings to config from analysis start file.
        """
        self.settings[key] = val

