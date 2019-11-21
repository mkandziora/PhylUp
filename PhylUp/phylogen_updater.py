"""
PhylUp: automatically update alignments and phylogenies.
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

Package to automatically update alignments and phylogenies using local and Genbank datasets

Parts of the code are inspired by the program physcraper developed by me and Emily Jane McTavish.

All classes and methods are distributed under the following license.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import re
import os
import glob
import sys
import subprocess
import random
import shutil
import numpy as np
from dendropy import Tree, DnaCharacterMatrix

import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser  # it is the ncbi data parser class and associated functions

from . import cd
from . import write_msg_logfile
from . import suppress_stdout
from . import phylogenetic_helpers


class AlnUpdater(object):
    """
    Class to update aln with newly retrieved sequences.
    """
    def __init__(self, config, aln, status, table, tre=None):
        """
        :param tre: dendropy tree object
        :param aln: dendropy alignment object
        :param status: status, how many cycles since start...
        :param table: pd df table with all the information about the seqs, old and new
        :param config: config object
        """
        self.aln = aln
        self.status = status
        self.table = table
        self.tre = tre
        self.config = config

        if self.tre is None:  # generate random tree, e.g. from modeltest
            self.tre_fn = "updt_aln.fasta.tree"
            best_subst_model = phylogenetic_helpers.run_modeltest('updt_aln.fasta',
                                                                  self.config.workdir, self.config.modeltest_criteria)
            self.tre = Tree.get(path=os.path.join(self.config.workdir, self.tre_fn),
                                schema="newick",
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
            phylogenetic_helpers.write_papara_trefile(self.tre, self.config.workdir)
        else:
            phylogenetic_helpers.write_aln(self.aln, self.config.workdir)
            phylogenetic_helpers.write_tre(self.tre, self.config.workdir)
        self.new_seq_table = self.table[self.table['status'] > 0]  # gets all new seqs (status>0)
        self.newseqs_file = "new_seqs.fasta"
        self.update_data()

    def update_data(self):
        """
        Key method of the class, cleans and updates tre and aln.

        :return: output files
        """
        print('update data')
        print(len(self.new_seq_table))
        if len(self.new_seq_table) > 0:
            self.write_papara_queryseqs()
            self.align_query_seqs()
            phylogenetic_helpers.write_aln(self.aln, self.config.workdir)
            phylogenetic_helpers.write_tre(self.tre, self.config.workdir)
            self.prune_short()
            phylogenetic_helpers.write_aln(self.aln, self.config.workdir)
            phylogenetic_helpers.write_tre(self.tre, self.config.workdir)
            self.aln = self.trim(os.path.abspath(os.path.join(self.config.workdir, 'papara_alignment.extended')),
                                 format_aln='phylip')
            phylogenetic_helpers.write_aln(self.aln, self.config.workdir)
            phylogenetic_helpers.write_tre(self.tre, self.config.workdir)
        self.write_labelled('updt_aln.fasta')

    def align_query_seqs(self, papara_runname="extended"):
        """Runs papara to add new sequences to the alignment.

        :param papara_runname:
        :return: writes new extended alignment
        """
        sys.stdout.write("aligning query sequences \n")
        for filename in glob.glob('{}/papara*'.format(self.config.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.config.workdir, filename.split("/")[-1]))
        phylogenetic_helpers.write_papara_alnfile(self.aln, self.config.workdir)
        if self.tre is not None:
            phylogenetic_helpers.write_papara_trefile(self.tre, self.config.workdir)
            assert self.aln.taxon_namespace == self.tre.taxon_namespace
        with cd(self.config.workdir):
            self.call_papara(papara_runname)
            path = os.path.join(os.getcwd(), "papara_alignment.{}".format(papara_runname))
            assert os.path.exists(path), "{} does not exists".format(path)
        self.aln = self.trim()
        msg = "Following papara alignment, aln has {} seqs \n".format(len(self.aln))
        write_msg_logfile(msg, self.config.workdir)

    def call_papara(self, papara_runname="extended"):
        """
        Calls external papara and returns error is program is not installed properly.

        :param papara_runname:  possible file extension name for papara
        :return:
        """
        try:
            with suppress_stdout():
                subprocess.check_call(["papara_static_x86_64",
                                       "-t", "random_resolve.tre",
                                       "-s", "aln_papara.phy",
                                       #  "-j", "{}".format(self.config.num_threads),  # FIXME: only works when papara is compiled.
                                       "-q", self.newseqs_file,
                                       "-n", papara_runname])
        except subprocess.CalledProcessError as grepexc:
            print("error code", grepexc.returncode, grepexc.output)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                sys.stderr.write("Failed running papara. Is it installed?\n")
                sys.exit(-5)
            else:
                # Something else went wrong while trying to run ...
                raise

    # todo?: move to papara_helper_class (next 3 methods)
    def write_papara_queryseqs(self):
        """Writes out the query sequence file which is needed by papara.

        Is only used within func align_query_seqs.
        """
        print("write query seq")
        fi = open(os.path.join(self.config.workdir, self.newseqs_file), "w")
        for idx in self.new_seq_table.index:
            tip_name = self.new_seq_table.loc[idx, 'accession'].split('.')[0]
            fi.write(">{}\n".format(tip_name))
            # !! index is some form of list. only like this it will print it as str
            seq = self.new_seq_table.loc[idx]['sseq']
            fi.write("{}\n".format(seq))

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
            msg = "Taxa pruned from tree and alignment in prune short " \
                  "step due to sequence shorter than {}:\n".format(seq_len_cutoff)
            write_msg_logfile(msg, self.config.workdir)
            for tax in prune:
                self.remove_taxa_aln_tre(tax.label)
                msg = '{}, '.format(tax.label)
                write_msg_logfile(msg, self.config.workdir)
            msg = '\n'
            write_msg_logfile(msg, self.config.workdir)
        for tax in prune:
            self.table.loc[self.table['accession'] == tax.label, "status"] = -1
            self.table.loc[self.table['accession'] == tax.label, "status_note"] = "deleted in prune short"

    def remove_taxa_aln_tre(self, taxon_label):
        """Removes taxa from aln and tre, takes a single taxon_label as input.

        :param taxon_label: taxon_label from dendropy object - aln or phy
        :return: removes data for taxon_label
        """
        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        tax2 = self.tre.taxon_namespace.get_taxon(taxon_label)

        self.aln.remove_sequences([tax])
        self.aln.discard_sequences([tax])
        self.aln.taxon_namespace.remove_taxon_label(taxon_label)  # raises an error if label not found

        self.tre.prune_taxa([tax2])
        self.tre.prune_taxa_with_labels([taxon_label])
        self.tre.prune_taxa_with_labels([tax2])
        self.table[tax.label, 'status'] = -1
        self.table[tax.label, 'status_note'] = "deleted, input conflict"

    def trim(self, aln_fn=False, format_aln=None):
        """ It removes bases at the start and end of alignments, if they are represented by less than the value
        specified in the config file (trim_perc). E.g. 0.75 given in config.

        Means, that 75% of the sequences need to have a base present. This ensures, that not whole chromosomes
        get dragged in by cutting the ends of long sequences.
        """
        # print('in trim')
        print(aln_fn, format_aln)
        if aln_fn:
            aln = DnaCharacterMatrix.get(path=aln_fn, schema=format_aln)
        else:
            aln = self.aln
        seqlen = phylogenetic_helpers.check_align(aln)
        start = 0
        stop = seqlen
        cutoff = len(aln) * self.config.trim_perc
        for i in range(seqlen):
            counts = {"?": 0, "-": 0}
            for tax in aln:
                call = aln[tax][i].label
                if call in ["?", "-"]:
                    counts[call] += 1
            if counts['?'] + counts['-'] <= cutoff:
                start = i
                break
        for i in range(seqlen, 0, -1):
            counts = {'?': 0, '-': 0}
            for tax in aln:
                call = aln[tax][i - 1].label
                if call in ['?', '-']:
                    counts[call] += 1
            if counts['?'] + counts['-'] <= cutoff:
                stop = i
                break
        # shorten aln to start:stop
        for taxon in aln:
            aln[taxon] = aln[taxon][start:stop]
        msg = "Trimmed beginning and end of alignment to > {}. " \
              "Start: {}, stop: {}.\n".format(self.config.trim_perc, start, stop)
        write_msg_logfile(msg, self.config.workdir)

        if not aln_fn:
            aln_fn = os.path.abspath(os.path.join(self.config.workdir, 'updt_aln.fasta'))
        else:
            aln_fn = '{}_trim'.format(aln_fn)
            aln.write(path=os.path.join(self.config.workdir, 'updt_aln.fasta'), schema='fasta')
        aln.write(path="{}".format(aln_fn), schema='fasta')
        return aln

    def write_labelled(self, alnpath):
        """ Output tree and alignment with human readable labels.

        :param alnpath:  optional: full file name (including path) for alignment
        :return: writes out labelled phylogeny and alignment to file
        """
        print("write labelled files")
        phylogenetic_helpers.replace_uid_with_name(os.path.join(self.config.workdir, alnpath), self.table)


class PhyUpdater(object):
    """
    This updates the phylogeny with the new alignment before running.
    """
    def __init__(self, config, aln, table, tre=None):
        """
        :param tre: dendropy tree object
        :param aln: dendropy alignment object
        :param table: pd df table with all the information about the seqs, old and new
        :param config: config object
        """
        self.aln = aln
        self.table = table
        self.tre = tre
        if self.tre is None:  # generate random tree, e.g. from modeltest
            self.tre = "updt_aln.fasta.tree"
        self.config = config
        self.update_phyl()

    def update_phyl(self):
        self.place_query_seqs_epa()
        self.check_tre_in_aln()
        self.calculate_final_tree()  # comment out for development speed up
        if self.config.update_tree is True:
            self.tre = Tree.get(path=os.path.join(self.config.workdir, "fulltree.raxml.bestTree"),
                                schema="newick",
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
            fn = os.path.join(self.config.workdir, 'papara_alignment.extended')
            self.aln = DnaCharacterMatrix.get(path=fn, schema="phylip")
            self.aln.write(path=os.path.join(self.config.workdir, "updt_aln.fasta"), schema='fasta')
            shutil.copy(os.path.join(self.config.workdir, "fulltree.raxml.bestTree"),
                        os.path.join(self.config.workdir, "updt_tre.tre"))
            sys.stdout.write('Updating of aln and tre done.\n')

    def check_tre_in_aln(self):
        """ Makes sure that everything which is in tre is also found in aln.
        """
        aln_ids = set(taxon.label for taxon in self.aln)
        # todo assert does not work if i change aln name for file writing! (removing .1 from gb_acc)
        # assert aln_ids.issubset(self.table['tip_name'].tolist()), \
        #     ([x for x in aln_ids if x not in self.table['tip_name'].tolist()], self.table['tip_name'].tolist())
        treed_taxa = set(leaf.taxon.label for leaf in self.tre.leaf_nodes())
        treed_tax_not_in_aln = [tax for tax in treed_taxa if tax not in aln_ids]
        self.tre.prune_taxa(treed_tax_not_in_aln)
        assert treed_taxa.issubset(aln_ids), (treed_taxa, aln_ids)
        assert aln_ids.issubset(treed_taxa), (treed_taxa, aln_ids)

    def place_query_seqs_epa(self):
        """ Runs epa-ng on the tree, and the combined alignment including the new query seqs.
        Just for placement, to use as starting tree.
        """
        print("place query seq")
        if self.config.backbone is True:
            with cd(self.config.workdir):
                backbonetre = Tree.get(path=os.path.join(self.config.workdir, "backbone.tre"), schema="newick",
                                       preserve_underscores=True)
                phylogenetic_helpers.resolve_polytomies(backbonetre, self.config.workdir)
        phylogenetic_helpers.write_papara_trefile(self.tre, self.config.workdir)
        phylogenetic_helpers.resolve_polytomies(self.tre, self.config.workdir)

        with cd(self.config.workdir):
            # get model for epa-ng
            phylogenetic_helpers.truncate_papara_aln(self.aln)
            with suppress_stdout():
                subprocess.call(["epa-ng", "--ref-msa", "old_seqs.fasta", "--tree", "random_resolve.tre",
                                 "--query",  "new_seqs.fasta", "--outdir", "./", "--model", 'GTR+G',
                                 '--redo'])
                # make newick tre
                subprocess.call(["gappa", "examine", "graft",
                                 "--jplace-path", "epa_result.jplace",
                                 "--allow-file-overwriting"])
            placetre = Tree.get(path="epa_result.newick", schema="newick")
            placetre.write(path="place_resolve.tre", schema="newick", unquoted_underscores=True, suppress_rooting=True)
            # phylogenetic_helpers.evaluate_raxmlng_alignment()

    def est_full_tree_ng(self, aln_fn, best_subst_model, num_threads):
        """Full raxml run from the placement tree as starting tree.

        the backbone options allows to fix the skeleton of the starting tree and just newly estimates the other parts.
        """
        print("est full tree ng")
        with cd(self.config.workdir):
            seed = str(random.randint(1, 21))

            subprocess.call(["raxml-ng-mpi", '--check', '--msa', "{}".format(aln_fn),
                             '--model', best_subst_model, '--prefix', 'check'])
            if self.config.update_tree is True:
                if self.config.backbone is not True:
                    print('tree')
                    subprocess.call(["raxml-ng-mpi",
                                     "--threads", "{}".format(num_threads),
                                     "--msa", "{}".format(aln_fn),
                                     # '--model', "GTR+G",
                                     "--tree", "epa_result.newick",
                                     "--seed", "{}".format(seed),
                                     "--prefix", "fulltree"])
                else:
                    subprocess.call(["raxml-ng-mpi", "--threads", "{}".format(num_threads),
                                     "--msa {}".format(aln_fn),
                                     "--constraint-tree", "backbone.tre",
                                     "--seed", "{}".format(seed),
                                     "--prefix", "fulltree"])
            else:
                todo = 'To update the data run the following command in your working directory.'
                if self.config.backbone is not True:
                    cmd1 = "raxml-ng-mpi --threads {} --msa {} --tree epa_result.newick --seed {}" \
                           " --prefix fulltree".format(num_threads, aln_fn, seed)
                else:
                    cmd1 = "raxml-ng-mpi --threads {} --msa {} --constraint-tree backbone.tre --seed {}" \
                           " --prefix fulltree".format(num_threads, aln_fn, seed)
                print(todo)
                print(cmd1)

    def calculate_bootstrap_ng(self, best_subst_model, num_threads, aln_fn='updt_aln.fasta'):
        """Calculates bootstrap and consensus trees.
        """
        print("calculate bootstrap")
        with cd(self.config.workdir):
            # check if job was started with mpi: this checks if actual several cores and nodes were allocated
            ntasks = os.environ.get('SLURM_NTASKS_PER_NODE')
            ntasks = os.environ.get('SLURM_JOB_CPUS_PER_NODE')
            nnodes = os.environ.get("SLURM_JOB_NUM_NODES")
            seed = str(random.randint(1, 21))
            mpi = False
            if self.config.update_tree is True:
                if nnodes is not None and ntasks is not None:
                    mpi = True
                    env_var = int(nnodes) * int(ntasks)
                if mpi:  # todo add "mpiexec", "-n", "{}".format(env_var)
                    # try:
                    #     print("run with mpi")
                    #     subprocess.call(["mpiexec", "-n", "{}".format(env_var), "raxml-ng-mpi", '--all', "--msa",
                    #                      "{}".format(aln_fn), '--model', "GTR+G", '--bs-trees', 'autoMRE',
                    #                      '--seed', seed, "--threads", "{}".format(str(self.config.num_threads)),
                    # except:
                    subprocess.call(["raxml-ng-mpi", '--all', "--msa", "{}".format(aln_fn), '--model', best_subst_model,
                                     '--bs-trees',  # 'tbe',  #
                                     'autoMRE', '--seed', seed, "--threads", "{}".format(num_threads),
                                     "--prefix", "fulltree"])
                else:
                    print('else')
                    subprocess.call(["raxml-ng-mpi", '--all', "--msa", "{}".format(aln_fn), '--model', best_subst_model,
                                     '--bs-trees',  # 'fbp,tbe',
                                     'autoMRE', '--seed', seed, "--threads", "{}".format(num_threads),
                                     "--prefix", "fulltree"])
                subprocess.call(["raxml-ng-mpi", '--consense', 'MRE', '--tree', 'fulltree.raxml.bootstraps',
                                 "--prefix", 'consMRE'])
                subprocess.call(["raxml-ng-mpi", '--consense', 'STRICT', '--tree', 'fulltree.raxml.bootstraps',
                                 "--prefix", 'consSTRICT'])
                subprocess.call(["raxml-ng-mpi", '--consense', 'MR', '--tree', 'fulltree.raxml.bootstraps',
                                 "--prefix", 'consMR'])
            else:
                todo = 'To update the data run the following command in your working directory.'
                cmd1 = "raxml-ng-mpi  --all --threads {} --msa {} --model {} --bs-trees autoMRE --seed {}" \
                       " --prefix fulltree".format(num_threads, aln_fn, best_subst_model, seed)
                cmd2 = "raxml-ng-mpi --consense MRE --tree fulltree.raxml.bootstraps --prefix consMRE"
                cmd3 = "raxml-ng-mpi --consense STRICT --tree fulltree.raxml.bootstraps --prefix consSTRICT"
                cmd4 = "raxml-ng-mpi --consense MR --tree fulltree.raxml.bootstraps --prefix consMR"

                print(todo)
                print(cmd1)
                print(cmd2)
                print(cmd3)
                print(cmd4)

    def calculate_final_tree(self):
        """Calculates the final tree using a trimmed alignment.

        :return: final data
        """
        print("calculate final tree")
        if self.config.backbone is True:
            aln_fn = 'papara_alignment.extended'
        else:
            aln_fn = 'updt_aln.fasta'

        best_subst_model = phylogenetic_helpers.run_modeltest(aln_fn, self.config.workdir, self.config.modeltest_criteria)
        num_threads = phylogenetic_helpers.estimate_number_threads_raxml(self.config.workdir,
                                                                         "updt_aln.fasta", best_subst_model)

        if self.config.backbone is True:
            self.est_full_tree_ng(aln_fn, best_subst_model, num_threads)
        else:
            self.calculate_bootstrap_ng(best_subst_model, num_threads, aln_fn)

    def write_labelled(self, treepath):
        """ Output tree with human readable labels.

        :param treepath: full file name (including path) for phylogeny
        :return: writes out labelled phylogeny and alignment to file
        """
        print("write labelled files")
        if treepath is None:
            treepath = os.path.join(self.config.workdir, "fulltree.raxml.bestTree")
        else:
            treepath = os.path.join(self.config.workdir, treepath)
        phylogenetic_helpers.replace_uid_with_name(treepath, self.table)


class InputCleaner(object):
    """
    This is the input class, that cleans the data before updating the phylogeny.
    """
    def __init__(self, tre_fn, aln_fn, table, config_obj, mrca=None):
        print('Clean the input data: {}, {}.'.format(tre_fn, aln_fn))
        self.config = config_obj
        if not os.path.exists(self.config.workdir):
            os.makedirs(self.config.workdir)
        # ncbi parser contains information about spn, tax_id, and ranks
        self.ncbi_parser = ncbi_data_parser.Parser(names_file=config_obj.ncbi_parser_names_fn,
                                                   nodes_file=config_obj.ncbi_parser_nodes_fn)
        # set mrca id
        assert type(mrca) in [int, list] or mrca is None, ("mrca is not an integer, list or None: {}".format(mrca))
        self.mrca = self.format_mrca_set(mrca)

        self.table = table
        self.aln = self.write_clean_aln(aln_fn)
        if tre_fn is not None:
            self.tre = self.write_clean_tre(tre_fn)
            if self.config.different_level is False:
                self._reconcile()  # turned of for different level, as tre is not updated between runs, aln is.
        self._reconcile_names()

    def format_mrca_set(self, mrca):
        """
        Takes the input and makes a set of valid mrca ids.

        :param mrca:
        :return:
        """
        print('Format mrca: {}'.format(mrca))
        if type(mrca) is int:
            valid = self.ncbi_parser.taxid_is_valid(mrca)
            mrca = {mrca}
            if not valid:
                sys.stderr.write("Your mrca id is not known by ncbi: {}".format(mrca))
        elif type(mrca) in {list, set}:
            mrca_set = set()
            for item in mrca:
                valid = self.ncbi_parser.taxid_is_valid(item)
                if valid:
                    mrca_set.add(item)
                else:
                    sys.stderr("Your mrca id is not known by ncbi: {}".format(item))
            mrca = mrca_set
            # mrca = self.ncbi_parser.get_mrca(taxid_set=mrca_set)
        elif mrca is None:
            # print('mrca is None - find mrca of all input data.')
            aln_taxids = set(self.table['ncbi_txid'].tolist())
            mrca = {self.ncbi_parser.get_mrca(taxid_set=aln_taxids)}
        else:
            sys.stderr.write("Method 'format_mrca_set' does not behave as expected")
        assert type(mrca) in {list, set, int}, ("your ingroup_mrca '%s' is not an integer/list/set." % mrca)
        return mrca

    def _reconcile(self):
        """Taxa that are only found in the tree, or only in the alignment are deleted.
        """
        assert type(self.aln) == DnaCharacterMatrix

        print("Reconcile taxa between input data")
        treed_taxa = set(leaf.taxon for leaf in self.tre.leaf_nodes())
        aln_tax = set(taxon for taxon in self.aln)
        prune = treed_taxa ^ aln_tax

        del_aln = []
        del_tre = []
        for taxon in prune:
            assert (taxon in aln_tax) or (taxon in treed_taxa)
            if taxon in aln_tax:
                del_aln.append(taxon)
            if taxon in treed_taxa:
                del_tre.append(taxon)
        self.aln.discard_sequences(del_aln)
        self.tre.prune_taxa(del_tre)

        # remove associated taxa from namespace
        if del_aln != []:
            print('Remove data from aln.')
            for item in del_aln:
                self.aln.taxon_namespace.remove_taxon_label(item.label)  # can only be removed one by one
                msg = '{} is only in the alignment, not in the tree and will be deleted.\n'.format(item.label)
                write_msg_logfile(msg, self.config.workdir)
        if del_tre != []:
            print('Remove data from tre.')
            for item in del_tre:
                self.tre.taxon_namespace.remove_taxon_label(item.label)  # can only be removed one by one
                msg = '{} is only in the tree, not in the alignment and will be deleted.\n'.format(item.label)
                write_msg_logfile(msg, self.config.workdir)
        prune_list = list(prune)
        if len(prune) > 0:
            true_false = np.where((self.table.accession.isin(prune_list)), True, False)
            self.table.loc[true_false, 'status'] = -1
            self.table.loc[true_false, 'status_note'] = "deleted in reconciliation"
        assert self.aln.taxon_namespace == self.tre.taxon_namespace

    def _reconcile_names(self):
        """It rewrites some tip names, which kept being an issue when it starts with a number at the beginning.
        Then somehow a n was added to the tip names. Only used for own input data

        :return: replaced tip names
        """
        tipnames = self.table['accession'].tolist()
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
                    original = self.table.loc[idx, "accession"]
                    if original == tax.label or original == newname:
                        tax.label = self.table.loc[idx, "accession"]
                        found_label = 1
                if found_label == 0:
                    sys.stderr.write("could not match tip label {} or {} to "
                                     "any ncbi taxon name\n".format(tax.label, newname))

    def write_clean_aln(self, aln_fn):
        """
        Write out original and cleaned alignemnt (? converted to -, no whitespaces)

        :param aln_fn: filename of alignment
        :return:
        """
        with open(aln_fn, "r") as fin:
            filedata = fin.read()
        # Write the file out again
        with open(os.path.join(self.config.workdir, "orig_inputaln.fasta"), "w") as aln_file:
            aln_file.write(filedata)
        filedata = filedata.replace("?", "-")
        filedata = filedata.replace(" ", "_")

        upd_aln_fn = os.path.join(self.config.workdir, "updt_aln.fasta")
        with open(upd_aln_fn, "w") as aln_file:
            aln_file.write(filedata)
        # use replaced aln as input
        aln = DnaCharacterMatrix.get(path=upd_aln_fn, schema='fasta')  # , taxon_namespace=self.tre.taxon_namespace)
        return aln

    def write_clean_tre(self, tre_fn):
        """
        Write out original and cleaned tre (no whitespaces)

        :param tre_fn: tree file name
        :return:
        """
        if tre_fn is not None:
            with open(tre_fn, "r") as fin:
                filedata = fin.read()
            # Write the file out again
            with open(os.path.join(self.config.workdir, "orig_tre.tre"), "w") as tre_file:
                tre_file.write(filedata)
            filedata = filedata.replace(" ", "_")

            upd_tre_fn = os.path.join(self.config.workdir, "updt_tre.tre")
            with open(upd_tre_fn, "w") as tre_file:
                tre_file.write(filedata)
            # use replaced aln as input
            tre = Tree.get(path=upd_tre_fn, schema='newick', taxon_namespace=self.aln.taxon_namespace,
                           preserve_underscores=True)
            return tre
