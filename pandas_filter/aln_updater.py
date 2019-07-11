# class for aln updating
#from __future__ import absolute_import

import os
import glob
import sys
import subprocess
import datetime
import re

import numpy as np

from . import ncbi_data_parser  # is the ncbi data parser class and associated functions
from . import cd

from dendropy import Tree, DnaCharacterMatrix  # , DataSet, datamodel

# next class is basically the ATT class of PhyScraper by EJ McTavish


class PhyAlnUpdater(object):
    """
    Class to update aln and tre with newly retrived sequences.

    Last step in the pipeline.
    """
    def __init__(self, tre, aln, status, new_seq_acc_seq, table, config):
        """
        :param tre: dendropy tree object
        :param aln: dendropy alignment object
        :param status: status, how many cycles since start...
        :param new_seq_acc_seq:
        :param table: pd df table with all the information about the seqs, old and new
        :param config: config object
        """
        self.aln = aln
        self.status = status
        self.table = table
        self.tre = tre
        self.config = config

        # todo is this needed???
        if self.status == 1:
            add_to_aln = self.table[self.table['status'] != 0]
        else:
            add_to_aln = self.table[self.table['status'] >= 1]
        self.new_seq_table = add_to_aln
        self.new_seqs = new_seq_acc_seq
        self.newseqs_file = "{}.fasta".format(str(datetime.date.today()))

    def update_data(self):
        """
        Key method of the class, cleans and updates tre and aln.

        :return: output files
        """
        if len(self.new_seq_table) > 0:
            self.write_papara_queryseqs()
            self.align_query_seqs()
            self.place_query_seqs()
            self.prune_short()
            self.write_papara_alnfile()
            self.trim()
            print('est tree')
            self.est_full_tree()  # for development speed up
            # self.calculate_final_tree()

            self.tre = Tree.get(path="{}/RAxML_bestTree.{}".format(self.config.workdir, str(datetime.date.today())),
                                schema="newick",
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
            os.rename("{}/RAxML_bestTree.{}".format(self.config.workdir, str(datetime.date.today())),
                      "{}/upd_tre".format(self.config.workdir))

    def align_query_seqs(self, papara_runname="extended"):
        """Runs papara on tre, aln and the new query sequences.

        :param papara_runname: possible file extension name for papara
        :return: writes out files after papara run/aligning seqs
        """
        cwd = os.getcwd()
        for filename in glob.glob('{}/papara*'.format(self.config.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.config.workdir, filename.split("/")[-1]))
        sys.stdout.write("aligning query sequences \n")
        self.write_papara_alnfile()
        self.write_papara_trefile()
        assert self.aln.taxon_namespace == self.tre.taxon_namespace

        os.chdir(self.config.workdir)  # todo Clean up dir moving
        try:
            subprocess.check_call(["papara_static_x86_64",
                                   "-t", "random_resolve.tre",
                                   "-s", "aln_papara.phy",
                                   #  "-j", "{}".format(self.config.num_threads),  # FIXME: Does not work on some machines
                                   "-q", self.newseqs_file,
                                   "-n", papara_runname])
            sys.stdout.write("Papara done \n")
        except subprocess.CalledProcessError as grepexc:
            print("error code", grepexc.returncode, grepexc.output)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                sys.stderr.write("Failed running papara. Is it installed?\n")
                sys.exit(-5)
            else:
                # Something else went wrong while trying to run ...
                raise
        path = "{}/papara_alignment.{}".format(os.getcwd(), papara_runname)
        assert os.path.exists(path), "{path} does not exists".format(path=path)
        os.chdir(cwd)

        fn = "{}/papara_alignment.{}".format(self.config.workdir, papara_runname)
        fn = os.path.abspath(fn)
        self.aln = DnaCharacterMatrix.get(path=fn, schema="phylip")
        self.aln.taxon_namespace.is_mutable = True  # Was too strict...
        lfd = "{}/logfile".format(self.config.workdir)
        with open(lfd, "a") as log:
            log.write("Following papara alignment, aln has {} seqs \n".format(len(self.aln)))

    def write_papara_queryseqs(self):
        """Writes out the query sequence file which is needed by papara.

        Is only used within func align_query_seqs.
        """
        print("write query seq")
        fi = open("{}/{}".format(self.config.workdir, self.newseqs_file), "w")
        for idx in self.new_seq_table.index:
            row = (self.new_seqs[self.new_seqs['accession'].str.contains(self.new_seq_table.loc[idx, 'tip_name'],
                                                                         regex=False)])
            if len(row) > 0:
                tip_name = self.new_seq_table.loc[idx, 'tip_name'].split('.')[0]
                fi.write(">{}\n".format(tip_name))
                idx_newseq = row.index
                # !! index is some form of list. only like this it will print it as str
                seq = self.new_seqs.loc[idx_newseq[0]]['sseq']
                fi.write("{}\n".format(seq))

    def write_papara_alnfile(self, ):
        """This writes out aln files for papara (except query sequences).
        Papara needs phylip format for the alignment.

        Is only used within func align_query_seqs.
        """
        alnfilename = "aln_papara.phy"
        self.aln.write(path="{}/{}".format(self.config.workdir, alnfilename), schema="phylip")

    def write_papara_trefile(self):
        """This writes out needed files for papara (except query sequences).
        Papara is finicky about trees and needs phylip format for the alignment.

        NOTE: names for tree and aln files should not be changed, as they are hardcoded in align_query_seqs().

        Is only used within func align_query_seqs.
        """
        print('write papara files')
        self.tre.resolve_polytomies()
        self.tre.deroot()
        tmptre = self.tre.as_string(schema="newick",
                                    unquoted_underscores=True,
                                    suppress_rooting=True)
        tmptre = tmptre.replace(":0.0;", ";")  # Papara is difficult about root
        tmptre = tmptre.replace("'", "_")

        filename = "random_resolve.tre"
        fi = open("{}/{}".format(self.config.workdir, filename), "w")
        fi.write(tmptre)
        fi.close()

    def prune_short(self):
        """Prunes sequences from alignment if they are shorter than specified in the config file.

        :return: prunes aln and tre
        """
        orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in self.aln]
        avg_seqlen = sum(orig_seqlen) / len(orig_seqlen)
        seq_len_cutoff = avg_seqlen * self.config.minlen
        prune = []
        aln_ids = set()
        for tax, seq in self.aln.items():
            aln_ids.add(tax.label)
            table = str.maketrans(dict.fromkeys('-?'))
            if len(seq.symbols_as_string().translate(table)) <= seq_len_cutoff:
                prune.append(tax)
        if prune:
            fi = open("{}/pruned_taxa".format(self.config.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in prune short "
                     "step due to sequence shorter than {}\n".format(seq_len_cutoff))
            for tax in prune:
                self.remove_taxa_aln_tre(tax.label)
                fi.write("{}, {}\n".format(tax.label, self.table[tax.label].get('tip_name')))
            fi.close()
        for tax in prune:
            self.table[tax.label]["status"] = "deleted in prune short"

    def remove_taxa_aln_tre(self, taxon_label):
        """Removes taxa from aln and tre and updates otu_dict,
        takes a single taxon_label as input.

        :param taxon_label: taxon_label from dendropy object - aln or phy
        :return: removes information/data from taxon_label
        """
        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        tax2 = self.tre.taxon_namespace.get_taxon(taxon_label)
        if tax:
            self.aln.remove_sequences([tax])
            self.aln.discard_sequences([tax])
            self.aln.taxon_namespace.remove_taxon_label(taxon_label)  # raises an error if label not found
            # the first prune does not remove it sometimes...
        if tax2:
            self.tre.prune_taxa([tax2])
            self.tre.prune_taxa_with_labels([taxon_label])
            self.tre.prune_taxa_with_labels([tax2])
            self.table[tax.label, 'status'] = "deleted, input conflict"
        else:
            self.table[tax.label, 'status'] = "deleted, updated table but was never in tre or aln!"

    def trim(self):
        """ It removes bases at the start and end of alignments, if they are represented by less than the value
        specified in the config file. E.g. 0.75 given in config.

        Means, that 75% of the sequences need to have a base present. Ensures, that not whole chromosomes get dragged in.
        It's cutting the ends of long sequences.
        """
        # print('in trim')
        taxon_missingness = self.config.trim_perc
        i = 0
        seqlen = len(self.aln[i])
        while seqlen == 0:
            i = i + 1
            seqlen = len(self.aln[i])
        for tax in self.aln:
            if len(self.aln[tax]) != seqlen:
                sys.stderr.write("can't trim un-aligned inputs, moving on")
                return
        start = 0
        stop = seqlen
        cutoff = len(self.aln) * taxon_missingness
        for i in range(seqlen):
            counts = {"?": 0, "-": 0}
            for tax in self.aln:
                call = self.aln[tax][i].label
                if call in ["?", "-"]:
                    counts[call] += 1
            if counts['?'] + counts['-'] <= cutoff:
                start = i
                break
        for i in range(seqlen, 0, -1):
            counts = {'?': 0, '-': 0}
            for tax in self.aln:
                call = self.aln[tax][i - 1].label
                if call in ['?', '-']:
                    counts[call] += 1
            if counts['?'] + counts['-'] <= cutoff:
                stop = i
                break
        # here alignment gets shortened to start:stop
        for taxon in self.aln:
            self.aln[taxon] = self.aln[taxon][start:stop]
        # make sure that tre is presented in aln
        self.check_tre_in_aln()

        # write to log
        lfd = "{}/logfile".format(self.config.workdir)
        with open(lfd, "a") as log:
            log.write("trimmed alignment ends to < {} missing taxa, "
                         "start {}, stop {}\n".format(taxon_missingness, start, stop))
        return

    def check_tre_in_aln(self):
        """ Makes sure that everything which is in tre is also found in aln.
        """
        print('check_tre_in_aln')
        aln_ids = set(taxon.label for taxon in self.aln)
        # todo assert does not work if i change aln name for file writing! (removing .1 from gb_acc
        # assert aln_ids.issubset(self.table['tip_name'].tolist()), ([x for x in aln_ids if x not in self.table['tip_name'].tolist()], self.table['tip_name'].tolist())
        treed_taxa = set(leaf.taxon.label for leaf in self.tre.leaf_nodes())
        treed_tax_not_in_aln = [tax for tax in treed_taxa if tax not in aln_ids]
        self.tre.prune_taxa(treed_tax_not_in_aln)

        assert treed_taxa.issubset(aln_ids), (treed_taxa, aln_ids)

    def place_query_seqs(self):
        """ Runs raxml on the tree, and the combined alignment including the new query seqs.
        Just for placement, to use as starting tree.
        """
        print("place query seq")
        if self.config.backbone is True:
            with cd(self.config.workdir):
                backbonetre = Tree.get(path="{}/backbone.tre".format(self.config.workdir),
                                       schema="newick",
                                       preserve_underscores=True)
                backbonetre.resolve_polytomies()
                backbonetre.write(path="random_resolve.tre", schema="newick", unquoted_underscores=True)

        if os.path.exists("RAxML_labelledTree.PLACE"):
            os.rename("RAxML_labelledTree.PLACE", "RAxML_labelledTreePLACE.tmp")
        sys.stdout.write("placing query sequences \n")
        with cd(self.config.workdir):
            try:
                print("try raxmlHPC-PTHREADS")
                subprocess.call(["raxmlHPC-PTHREADS",
                                 "-T", "{}".format(self.config.num_threads),
                                 "-m", "GTRCAT",
                                 "-f", "v",
                                 "-s", "papara_alignment.extended",
                                 "-t", "random_resolve.tre",
                                 "-n", "PLACE"])
                placetre = Tree.get(path="RAxML_labelledTree.PLACE",
                                    schema="newick",
                                    preserve_underscores=True)
            except:
                try:
                    sys.stderr.write("You do not have the raxmlHPC-PTHREADS installed, will fall down to slow version!")
                    subprocess.call(["raxmlHPC",
                                     "-m", "GTRCAT",
                                     "-f", "v",
                                     "-s", "papara_alignment.extended",
                                     "-t", "random_resolve.tre",
                                     "-n", "PLACE"])
                    placetre = Tree.get(path="RAxML_labelledTree.PLACE",
                                        schema="newick",
                                        preserve_underscores=True)
                except OSError as e:
                    if e.errno == os.errno.ENOENT:
                        sys.stderr.write("Failed running raxmlHPC. Is it installed?")
                        sys.exit(-6)
                    else:
                        # Something else went wrong while trying to run `wget`
                        raise
            placetre.resolve_polytomies()
            for taxon in placetre.taxon_namespace:
                if taxon.label.startswith("QUERY"):
                    taxon.label = taxon.label.replace("QUERY___", "")
            placetre.write(path="place_resolve.tre", schema="newick", unquoted_underscores=True)

    def est_full_tree(self, path="."):
        """Full raxml run from the placement tree as starting tree.
        The PTHREAD version is the faster one, hopefully people install it. if not it falls back to the normal raxml.
        the backbone options allows to fix the sceleton of the starting tree and just newly estimates the other parts.
        """
        print("est full tree")
        cwd = os.getcwd()
        os.chdir(self.config.workdir)

        date = str(datetime.date.today())  # Date of the run - may lag behind real date!
        for filename in glob.glob('{}/RAxML*'.format(self.config.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.config.workdir, filename.split("/")[-1]))
        try:
            num_threads = int(self.config.num_threads)
            if self.config.backbone is not True:
                print(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                       "-s", "{}/papara_alignment.extended".format(path),
                       "-t", "place_resolve.tre",
                       "-p", "1",
                       "-n", "{}".format(date)])
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-t", "place_resolve.tre",
                                 "-p", "1",
                                 "-n", "{}".format(date)])
            else:
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(num_threads), "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-r", "backbone.tre",
                                 "-p", "1",
                                 "-n", "{}".format(date)])
        except:
            sys.stderr.write("You do not have the raxmlHPC-PTHREADS installed, will fall down to slow version!")
            if self.config.backbone is not True:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-t", "place_resolve.tre",
                                 "-p", "1",
                                 "-n", "{}".format(date)])
            else:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                 "-s", "{}/papara_alignment.extended".format(path),
                                 "-r", "backbone.tre",
                                 "-p", "1",
                                 "-n", "{}".format(date)])
        os.chdir(cwd)

    def calculate_bootstrap(self):
        """Calculates bootstrap and consensus trees.

        -p: random seed
        -s: aln file
        -n: output fn
        -t: starting tree
        -b: bootstrap random seed
        -#: bootstrap stopping criteria
        -z: specifies file with multiple trees

        -f b: to make bipartition tree only
        """
        print("calculate bootstrap")
        cwd = os.getcwd()
        os.chdir(self.config.workdir)

        # check if job was started with mpi: this checks if actual several cores and nodes were allocated
        ntasks = os.environ.get('SLURM_NTASKS_PER_NODE')
        nnodes = os.environ.get("SLURM_JOB_NUM_NODES")
        date = str(datetime.date.today())  # Date of the run - may lag behind real date!
        aln_fn = os.path.abspath('papara_alignment.extended')

        mpi = False
        if nnodes is not None and ntasks is not None:
            mpi = True
        if mpi:
            print("run with mpi")
            env_var = int(nnodes) * int(ntasks)
            subprocess.call(["mpiexec", "-n", "{}".format(env_var), "raxmlHPC-MPI-AVX2",
                             # "raxmlHPC-PTHREADS", "-T", "{}".format(num_threads),
                             "-m", "GTRCAT",
                             "-s", "{}".format(aln_fn),
                             "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                             "-n", "{}".format(date)])
        else:
            try:
                subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                                 "-s", "{}".format(aln_fn),
                                 "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                                 "-n", "all{}".format(date)])
            except:
                subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                                 "-s", "{}".format(aln_fn),
                                 "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
                                 "-n", "all{}".format(date)])
        try:
            # strict consensus:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-J", "STRICT",
                             "-z", "RAxML_bootstrap.all{}".format(date),
                             "-n", "StrictCon{}".format(date)])
            # majority rule:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-J", "MR",
                             "-z", "RAxML_bootstrap.all{}".format(date),
                             "-n", "MR_{}".format(date)])
            # extended majority rule:
            subprocess.call(["raxmlHPC-PTHREADS", "-T", "{}".format(self.config.num_threads), "-m", "GTRCAT",
                             "-J", "MRE",
                             "-z", "RAxML_bootstrap.all{}".format(date),
                             "-n", "EMR{}".format(date)])
        except:
            sys.stderr.write("You do not have the raxmlHPC-PTHREADS installed, will fall down to slow version!")

            # strict consensus:
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-J", "STRICT",
                             "-z", "RAxML_bootstrap.all{}".format(date),
                             "-n", "StrictCon{}".format(date)])
            # majority rule:
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-J", "MR",
                             "-z", "RAxML_bootstrap.all{}".format(date),
                             "-n", "MR_{}".format(date)])
            # extended majority rule:
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-J", "MRE",
                             "-z", "RAxML_bootstrap.all{}".format(date),
                             "-n", "EMR{}".format(date)])
        os.chdir(cwd)

    def calculate_final_tree(self):
        """Calculates the final tree using a trimmed alignment.

        :return: final PS data
        """
        print("calculate final tree")
        self.write_files(treepath="physcraper_final_notrim.tre", alnpath="physcraper_final_notrim.fas")
        self.prune_short()
        self.trim()
        self.write_files(treepath="physcraper_final_trim.tre", alnpath="physcraper_final_trim.fas")

        if os.path.exists("[]/previous_run".format(self.config.workdir)):
            self.est_full_tree(path="previous_run")
        else:
            self.est_full_tree()
        self.calculate_bootstrap()

        # todo reassign self tre and aln to updated data

    def write_files(self, treepath="PhyFilter.tre", treeschema="newick", alnpath="PhyFilter.fas", alnschema="fasta"):
        """Outputs dendropy tre and aln as file
        """
        print("write_files")
        self.tre.write(path="{}/{}".format(self.config.workdir, treepath),
                       schema=treeschema, unquoted_underscores=True)
        self.aln.write(path="{}/{}".format(self.config.workdir, alnpath),
                       schema=alnschema)

    # todo check pandas version works
    def write_labelled(self, treepath=None, alnpath=None):
        """ Output tree and alignment with human readable labels.

        :param treepath: optional: full file name (including path) for phylogeny
        :param alnpath:  optional: full file name (including path) for alignment
        :return: writes out labelled phylogeny and alignment to file
        """
        print("write labelled files")
        if treepath == None:
            treepath = "{}/RAxML_bestTree.{}".format(self.config.workdir, str(datetime.date.today()))
        tmp_newick = self.tre.as_string(schema="newick")
        tmp_tre = Tree.get(data=tmp_newick,
                           schema="newick",
                           preserve_underscores=True)
        new_names = set()
        for taxon in tmp_tre.taxon_namespace:
            new_label = '{}_{}'.format(self.table['ncbi_txn'],
                                       self.table['tip_name'])  # tipname is either tipname or gi number!
            taxon.label = new_label
            new_names.add(new_label)
        tmp_tre.write(path=treepath,
                      schema="newick",
                      unquoted_underscores=True,
                      suppress_edge_lengths=False)

        if alnpath is not None:
            tmp_fasta = self.aln.as_string(schema="fasta")
            tmp_aln = DnaCharacterMatrix.get(data=tmp_fasta, schema="fasta",
                                             taxon_namespace=tmp_tre.taxon_namespace)
            tmp_aln.write(path=alnpath,
                          schema="fasta")


class InputCleaner(object):
    """
    This is the input class, that cleans the data before running.
    """
    def __init__(self, tre, aln, aln_fn, table, config_obj, mrca=None):
        print('Input Cleaner')
        self.tre = tre
        self.aln = aln
        self.config = config_obj
        self.table = table
        self.aln = self.write_clean_aln(aln_fn)

        self._reconcile()
        self._reconcile_names()

        if not os.path.exists(self.config.workdir):
            os.makedirs(self.config.workdir)

        # ncbi parser contains information about spn, tax_id, and ranks
        self.ncbi_parser = ncbi_data_parser.Parser(names_file=config_obj.ncbi_parser_names_fn,
                                                   nodes_file=config_obj.ncbi_parser_nodes_fn)
        # set mrca id
        assert type(mrca) in [int, list] or mrca is None, ("mrca is not an integer, list or None: {}".format(mrca))
        mrca = self.format_mrca_set(mrca)
        assert type(mrca) in {list, set, int}, ("your ingroup_mrca '%s' is not an integer/list/set." % mrca)
        self.mrca = mrca

    def format_mrca_set(self, mrca):
        """takes the input and makes a set of valid mrca ids.
        """
        print('format_mrca_set')
        print(mrca)
        if type(mrca) is int:
            valid = self.ncbi_parser.taxid_is_valid(mrca)
            if valid:
                print(type({mrca}))
                return mrca
            else:
                sys.stderr("Your mrca id is not known by ncbi: {}".format(mrca))
                return None
        elif type(mrca) in {list, set}:
            mrca_set = set()
            for item in mrca:
                valid = self.ncbi_parser.taxid_is_valid(item)
                if valid:
                    mrca_set.add(item)
                else:
                    sys.stderr("Your mrca id is not known by ncbi: {}".format(item))
            mrca = self.ncbi_parser.get_mrca(taxid_set=mrca_set)
            return mrca
        elif mrca is None:
            print('mrca is None')
            aln_taxids = set(self.table['ncbi_txid'].tolist())
            mrca = self.ncbi_parser.get_mrca(taxid_set=aln_taxids)
            return mrca
        else:
            sys.stderr.write("Method 'format_mrca_set' does not behave as expected")

    def _reconcile(self):
        """Taxa that are only found in the tree, or only in the alignment are deleted.

        This checks that the tree "original labels" from phylesystem
        align with those found in the alignment.
        """
        print("reconcile")
        treed_taxa = set(leaf.taxon.label for leaf in self.tre.leaf_nodes())
        aln_tax = set(taxon.label for taxon in self.aln)
        prune = treed_taxa ^ aln_tax
        if len(prune) > 1:
            errmf = 'RECONCILIATION Some of the taxa in the tree are not in the alignment or vice versa' \
                    ' and will be deleted. Missing "{}"\n'
            errm = errmf.format('", "'.join(prune))
            sys.stderr.write(errm)

        del_aln = []
        del_tre = []
        for taxon in prune:
            assert (taxon in aln_tax) or (taxon in treed_taxa)
            if taxon in aln_tax:
                del_aln.append(taxon)
            if taxon in treed_taxa:
                del_tre.append(taxon)
        # self.aln.remove_sequences(del_aln)  # raises error if taxon not in aln
        self.aln.discard_sequences(del_aln)
        self.tre.prune_taxa(del_tre)

        # remove associated taxa from namespace
        if del_aln != []:
            print('remove from aln del')
            for item in del_aln:
                self.aln.taxon_namespace.remove_taxon_label(item)  # can only be removed one by one
        if del_tre != []:
            print('remove from tre del')
            for item in del_tre:
                self.aln.taxon_namespace.remove_taxon_label(item)  # can only be removed one by one
        prune_list = list(prune)
        if len(prune) > 0:
            TF = np.where((self.table.tip_name.isin(prune_list)), True, False)
            self.table.loc[TF, 'status'] = "excluded"
            self.table.loc[TF, 'status_note'] = "deleted in reconciliation"
            print('upd_table')
            print(self.table['status'])
        assert self.aln.taxon_namespace == self.tre.taxon_namespace

    def _reconcile_names(self):
        """It rewrites some tip names, which kept being an issue when it starts with a number at the beginning.
        Then somehow a n was added to the tip names. Only used for own input data

        :return: replaced tip names
        """
        tipnames = self.table['tip_name'].tolist()
        for tax in self.aln.taxon_namespace:
            if tax.label in tipnames:
                pass
            else:
                found_label = 0
                match = re.match("'n[0-9]{1,3}", tax.label)
                newname = ""
                if match:
                    newname = tax.label[2:]
                    newname = newname[:-1]
                for idx in self.table.index:
                    original = self.table.loc[idx, "tip_name"]
                    if original == tax.label or original == newname:
                        tax.label = self.table.loc[idx, "accession"]
                        found_label = 1
                if found_label == 0:
                    sys.stderr.write("could not match tiplabel {} or {} to an OTU\n".format(tax.label, newname))

    def write_clean_aln(self, aln_fn):
        with open(aln_fn, "r") as fin:
            filedata = fin.read()
        # Write the file out again
        with open("{}/orig_inputaln.fasta".format(self.config.workdir), "w") as aln_file:
            aln_file.write(filedata)
        filedata = filedata.replace("?", "-")
        upd_aln_fn = "{}/updt_aln.fasta".format(self.config.workdir)
        with open(upd_aln_fn, "w") as aln_file:
            aln_file.write(filedata)
        # use replaced aln as input
        aln = DnaCharacterMatrix.get(path=upd_aln_fn, schema='fasta', taxon_namespace=self.tre.taxon_namespace)
        return aln

# #################################

# todo unused
def aln_to_seqs(aln, unalign):
    """Turn a dendropy alignment object into a dict, with otu_id and seqs.
    """
    print("aln_to_seqs query seq")
    tax_seq_dict = {}
    for tax, seq in aln.items():
        if unalign:
            table = str.maketrans(dict.fromkeys('-?'))
            seq = seq.symbols_as_string().translate(table)
            # TODO: maybe below is the correct syntax
            # self.data.aln[tax].symbols_as_string().replace("-", "").replace("N", "")
        else:
            seq = seq.symbols_as_string()
        tax_seq_dict[tax.label] = seq
    return tax_seq_dict