# user settings class

import os
import sys
import configparser
from . import db_updater
from . import debug


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
            * **self.email**: email address used for blast queries
            * **self.e_value_thresh**: the defined threshold for the e-value during Blast searches, check out:
                https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
            * **self.hitlist_size**: the maximum number of sequences retrieved by a single blast search
            * **self.blast_loc**: defines which blasting method to use. Currently only local is implemented
                * if self.blast_loc == "local":
                    * self.blastdb: this defines the path to the local blast database
                    * self.ncbi_parser_nodes_fn: path to 'nodes.dmp' file, that contains the hierarchical information
                    * self.ncbi_parser_names_fn: path to 'names.dmp' file, that contains the different ID's


            * **self.num_threads**: number of cores to be used during a run
            * **self.gb_id_filename**: Set to True, if you want to reuse gb query files
            * **self.delay**: defines when to reblast sequences in days
            * **self.minlen**: value from 0 to 1. Defines how much shorter new seq can be compared to input
            * **self.trim_perc**: value that determines how many seq need to be present before the beginning and end
                                    of alignment will be trimmed
            * **self.maxlen**: max length for values to add to aln
            * **self.add_lower_taxa**: T/F, enables to re-access formerly filtered seq by allowing them be passed
                                       into remove_identical. Used if we first filter for higher rank and then want
                                       to filter for a lower rank.
        * **self.filtertype**: either 'blast' or 'length'; filter number per otu by blast or choose longest
        * **self.backbone**: T/F; calculate complete tree or add to backbone
        * **self.threshold**: number of seq per otu
        * **self.downtorank**: filter number otu to specific rank
        * **self.unpublished**: T/F; add unpublished local sequences
        * **self.logfile**: file where some things are logged


        :param configfi: a configuration file in a specific format. The file needs to have a heading of the format:
                        [blast] and then somewhere below that heading as string, e.g. e_value_thresh = value
        :param workdir: the working directory
        :param interactive: defaults to True, is used to interactively update the local blast databases
        """
        sys.stdout.write("Building config object\n")
        debug(configfi)
        assert os.path.isfile(configfi), "file `%s` does not exists" % configfi

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
        assert is_number(self.e_value_thresh), ("value `%s` does not exists" % self.e_value_thresh)

        self.hitlist_size = int(config["blast"]["hitlist_size"])
        assert is_number(self.hitlist_size), ("value `%s`is not a number" % self.e_value_thresh)

        self.blast_loc = config["blast"]["location"]
        assert self.blast_loc in ["local"], ("your blast location `%s` is not local" % self.blast_loc)
        if self.blast_loc == "local":
            self.blastdb = config["blast"]["localblastdb"]
            if not os.path.exists(self.blastdb):
                os.makedirs(self.blastdb)
            self.ncbi_parser_nodes_fn = config["ncbi_parser"]["nodes_fn"]
            self.ncbi_parser_names_fn = config["ncbi_parser"]["names_fn"]

        self.num_threads = int(config["blast"].get("num_threads"))

        debug("slurm threads")
        debug(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))
        if os.environ.get('SLURM_JOB_CPUS_PER_NODE'):
            self.num_threads = int(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))
        debug(self.num_threads)

        # # todo currently not used
        # self.gb_id_filename = config["blast"].get("gb_id_filename", False)
        # if self.gb_id_filename is not False:
        #     if self.gb_id_filename == "True" or self.gb_id_filename == "true":
        #         self.gb_id_filename = True
        #     else:
        #         self.gb_id_filename = False
        # debug("shared blast folder? {}".format(self.gb_id_filename))

        self.delay = int(config["blast"]["delay"])
        assert is_number(self.delay), ("value `%s`is not a number" % self.delay)

        # #############
        # read in phy_filter settings
        self.minlen = float(config["phy_filter"]["min_len"])
        assert 0 < self.minlen <= 1, ("value `%s` is not between 0 and 1" % self.minlen)
        self.trim_perc = float(config["phy_filter"]["trim_perc"])
        assert 0 < self.trim_perc < 1, ("value `%s` is not between 0 and 1" % self.trim_perc)
        self.maxlen = float(config["phy_filter"]["max_len"])
        assert 1 < self.maxlen, ("value `%s` is not larger than 1" % self.maxlen)

        self.different_level = config["phy_filter"]["different_level"]
        if self.different_level == "True" or self.different_level == "true":
            self.different_level = True
        else:
            self.different_level = False
        assert self.different_level in [True, False], ("self.different_level `%s` is not True or False" % self.different_level)

        # todo: currently unused - implement length
        self.filtertype = config["phy_filter"]["filtertype"]
        assert self.filtertype in ['blast', 'length'], ("self.filtertype `%s` is not 'blast' or 'length'" % self.filtertype)

        self.backbone = config["phy_filter"]["backbone"]
        if self.backbone == "True" or self.backbone == "true":
            self.backbone = True
        else:
            self.backbone = False
        self.update_tree = config["phy_filter"]["update_tree"]
        if self.update_tree == "True" or self.update_tree == "true":
            self.update_tree = True
        else:
            self.update_tree = False

        self.threshold = int(config["phy_filter"]["threshold"])
        self.downtorank = config["phy_filter"]["downtorank"]

        self.unpublished = config["blast"]['unpublished']
        if self.unpublished == "True" or self.unpublished == "true":
            self.unpublished = True
        else:
            self.unpublished = False
        if self.unpublished is True:
            self.unpubl_data = config["blast"]['unpubl_data']
            self.unpubl_names = config["blast"]['unpubl_names']
        else:
            self.unpubl_data = None
            self.unpubl_names = None

        ####
        # check database status

        #todo: check how to make it interactive again.
        print('Status for interactive is {}'.format(interactive))
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

    #     # more settings
    #     self.settings = {}
    #     self.settings['blast_subdir'] = "{}/current_blast_run".format(self.workdir)
    #
    #
    # def add_to_config(self, key, val):
    #     """Add more settings to config from analysis start file.
    #     """
    #     self.settings[key] = val
    #
