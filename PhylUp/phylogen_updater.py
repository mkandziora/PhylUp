"""
PhylUp: automatically update alignments.
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

All rights reserved. No warranty, explicit or implicit, provided. No distribution or modification of code allowed.
All classes and methods will be distributed under an open license in the near future.

Package to automatically update alignments and phylogenies using local sequences or a  local Genbank database.

Parts are inspired by the program physcraper developed by me and Emily Jane McTavish.
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
                                preserve_underscores=True)
            phylogenetic_helpers.write_papara_trefile(self.tre, self.config.workdir)
        else:
            phylogenetic_helpers.write_aln(self.aln, self.config.workdir)
            phylogenetic_helpers.write_tre(self.tre, self.config.workdir)
        self.new_seq_table = self.table[self.table['status'] >= 1]  # gets all new seqs (status>0.5)
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
            self.delete_short_seqs()
            self.write_papara_queryseqs()
            self.add_query_seqs_to_aln()
            phylogenetic_helpers.write_aln(self.aln, self.config.workdir)
            phylogenetic_helpers.write_tre(self.tre, self.config.workdir)
            # self.aln = self.trim(os.path.abspath(os.path.join(self.config.workdir, 'papara_alignment.phylip')),
            #                      format_aln='phylip')
            # phylogenetic_helpers.write_aln(self.aln, self.config.workdir)
            # phylogenetic_helpers.write_tre(self.tre, self.config.workdir)
        self.write_labelled('updt_aln.fasta')

    def add_query_seqs_to_aln(self):
        """Runs papara to add new sequences to the alignment.

        :return: writes new extended alignment
        """
        sys.stdout.write("aligning query sequences \n")
        for filename in glob.glob('{}/papara*'.format(self.config.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.config.workdir, filename.split("/")[-1]))
        phylogenetic_helpers.write_papara_alnfile(self.aln, self.config.workdir)
        if self.tre is not None:
            phylogenetic_helpers.write_papara_trefile(self.tre, self.config.workdir)
        with cd(self.config.workdir):
            self.run_papara()
            path = os.path.join(os.getcwd(), "papara_alignment.phylip".format())
            assert os.path.exists(path), "{} does not exists".format(path)
        aln_old = self.aln
        self.aln = DnaCharacterMatrix.get(path=os.path.join(self.config.workdir, "papara_alignment.phylip"),
                                          schema='phylip')
        self.aln = self.trim(os.path.join(self.config.workdir, 'papara_alignment.phylip'), 'phylip')
        with cd(self.config.workdir):
            phylogenetic_helpers.truncate_papara_aln(aln_old)
        msg = "Following papara alignment, aln has {} seqs \n".format(len(self.aln))
        write_msg_logfile(msg, self.config.workdir)

    def run_papara(self):
        """
        Runs papara and adds new sequences to the alignment.

        :return:
        """
        with suppress_stdout():
            subprocess.check_call(["papara_static_x86_64", "-t", "papara_tre.tre", "-s", "aln_papara.phy",
                                   #  "-j", "{}".format(self.config.num_threads),  # FIXME: only works when papara is compiled.
                                   "-q", self.newseqs_file, "-n", 'phylip'])

    def write_papara_queryseqs(self):
        """Writes out the query sequence file which is needed by papara.

        Is only used within func add_query_seqs_to_aln.
        """
        print("write query seq")
        fi = open(os.path.join(self.config.workdir, self.newseqs_file), "w")
        for idx in self.new_seq_table.index:
            tip_name = self.new_seq_table.loc[idx, 'accession'].split('.')[0]
            fi.write(">{}\n".format(tip_name))
            # !! index is some form of list. only like this it will print it as str
            seq = self.new_seq_table.loc[idx]['sseq']
            fi.write("{}\n".format(seq))

    def delete_short_seqs(self):
        """Delete sequences from alignment if they are shorter than specified in the config file.

        :return: pruned aln and tre
        """
        orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("?", "")) for tax in self.aln]
        avg_seqlen = sum(orig_seqlen) / len(orig_seqlen)
        seq_len_cutoff = avg_seqlen * self.config.minlen
        delete_seqs = []
        aln_ids = set()
        for tax, seq in self.aln.items():
            aln_ids.add(tax.label)
            # table = str.maketrans(dict.fromkeys('-?'))
            # if len(seq.symbols_as_string().translate(table)) <= seq_len_cutoff:
            if len(seq.symbols_as_string().replace("-", "").replace("?", "")) <= seq_len_cutoff:
                delete_seqs.append(tax)
        if delete_seqs:
            msg = "Taxa deleted from tree and alignment in delete short sequences" \
                  "as sequences are shorter than %.0f:\n" % seq_len_cutoff
            write_msg_logfile(msg, self.config.workdir)
            for tax in delete_seqs:
                self.remove_taxa_aln_tre(tax.label)
                msg = '{}, '.format(tax.label)
                write_msg_logfile(msg, self.config.workdir)
                self.table.at[self.table['accession'] == tax.label, "status"] = -1
                self.table.at[self.table['accession'] == tax.label, "status_note"] = "deleted - too short"
            msg = '\n'
            write_msg_logfile(msg, self.config.workdir)
        self.aln.write(path=os.path.join(self.config.workdir, 'updt_aln.fasta'), schema='fasta')
        self.tre.write(path=os.path.join(self.config.workdir, 'updt_tre.tre'), schema='newick',
                       unquoted_underscores=True, suppress_rooting=True)
        self.table.to_csv(os.path.join(self.config.workdir, 'table.updated'), index=False)

    def remove_taxa_aln_tre(self, taxon_label):
        """Removes taxon from aln and tre, takes a single taxon_label as input.

        :param taxon_label: taxon_label from dendropy object - aln or phy
        :return:
        """
        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        tax2 = self.tre.taxon_namespace.get_taxon(taxon_label)

        self.aln.remove_sequences([tax])
        self.aln.discard_sequences([tax])
        self.aln.taxon_namespace.remove_taxon_label(taxon_label)  # raises an error if label not found

        self.tre.prune_taxa([tax2])
        self.tre.prune_taxa_with_labels([taxon_label])
        self.tre.prune_taxa_with_labels([tax2])

    def trim(self, aln_fn=False, format_aln=None):
        """
        Removes bases at the start and end of alignment, if they are represented by less than the value
        specified in the config file (trim_perc).

        Means, that trim_perc (e.g. 0.75 = 75%) of the sequences need to have a base present.
        This ensures, that not whole chromosomes get dragged in by cutting the ends of long sequences.
        """
        if aln_fn:
            aln = DnaCharacterMatrix.get(path=aln_fn, schema=format_aln)
        else:
            aln = self.aln
        if aln_fn:
            aln_fn_trim = '{}_untrimmed'.format(aln_fn)
            aln.write(path="{}".format(aln_fn_trim), schema=format_aln)
        aln.write(path=os.path.join(self.config.workdir, 'updt_aln_untrimmed.fasta'), schema='fasta')
        seqlen = phylogenetic_helpers.check_align(aln)
        start = 0
        stop = seqlen
        cutoff = len(aln) * self.config.trim_perc
        range_start = range(seqlen)
        start = get_base_position(aln, cutoff, range_start, start)
        range_end = range(seqlen - 1, 0, -1)
        stop = get_base_position(aln, cutoff, range_end, stop)
        # shorten aln to start:stop
        for taxon in aln:
            aln[taxon] = aln[taxon][start:stop]
        msg = "Original aln was trimmed from 1-{} to {}-{}.\n".format(seqlen, start, stop)
        write_msg_logfile(msg, self.config.workdir)

        if aln_fn:
            aln_fn_trim = '{}_trim'.format(aln_fn)
            aln.write(path="{}".format(aln_fn_trim), schema=format_aln)
        aln.write(path=os.path.join(self.config.workdir, 'updt_aln.fasta'), schema='fasta')
        return aln

    def write_labelled(self, alnpath):
        """ Output tree and alignment with human readable labels.

        :param alnpath:  optional: full file name (including path) for alignment
        :return: writes out labelled phylogeny and alignment to file
        """
        # print("write labelled aln")
        phylogenetic_helpers.replace_uid_with_name(os.path.join(self.config.workdir, alnpath), self.table, 'aln')


def get_base_position(aln, cutoff, range_val, baseposition):
    """
    Find baseposition in aln where sufficient data is available.

    :param aln: alignment as dendropy object
    :param cutoff: min value that needs to be present in aln
    :param range_val: range to test for - different for start and stop
    :param baseposition: position for cutting
    :return: baseposition
    """
    for i in range_val:
        counts = 0
        for tax in aln:
            base = aln[tax][i].label
            if base in ['?', '-']:
                counts += 1
        if counts <= cutoff:
            baseposition = i
            break
    return baseposition


class TreeUpdater(object):
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
            self.tre = Tree.get(path=os.path.join(self.config.workdir, "updt_aln.fasta.tree"), schema="newick",
                                preserve_underscores=True, taxon_namespace=self.aln.taxon_namespace)

        self.config = config
        self.update_phyl()

    def update_phyl(self):
        """
        Main function of the class. Does the updating of the phylogeny.
        :return:
        """
        print("update phylogeny")
        self.place_query_seqs_epa()
        self.check_tre_in_aln()
        self.calculate_final_tree()  # comment out for development speed up
        if self.config.update_tree is True:
            # self.tre = Tree.get(path=os.path.join(self.config.workdir, "fulltree.raxml.bestTree"),
            #                     schema="newick", preserve_underscores=True, taxon_namespace=self.aln.taxon_namespace)
            # fn = os.path.join(self.config.workdir, 'papara_alignment.phylip_trim', taxon_namespace=self.tre.taxon_namespace)
            # self.aln = DnaCharacterMatrix.get(path=fn, schema="phylip")
            # self.aln.write(path=os.path.join(self.config.workdir, "updt_aln.fasta"), schema='fasta')
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
        aln_ids = set(taxon.label for taxon in self.aln)
        treed_taxa = set(leaf.taxon.label for leaf in self.tre.leaf_nodes())
        assert treed_taxa.issubset(aln_ids), (len(treed_taxa), len(aln_ids))
        # assert aln_ids.issubset(treed_taxa),  (len(treed_taxa), len(aln_ids))

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
            with suppress_stdout():
                subprocess.call(["epa-ng", "--ref-msa", "old_seqs.fasta", "--tree", "papara_tre.tre",
                                 "--query",  "new_seqs.fasta", "--outdir", "./", "--model", 'GTR+G', '--redo'])
                # make newick tre
                subprocess.call(["gappa", "examine", "graft", "--jplace-path", "epa_result.jplace",
                                 "--allow-file-overwriting"])
            self.tre = Tree.get(path="epa_result.newick", schema="newick", preserve_underscores=True)
            #self.tre.write(path="place_resolve.tre", schema="newick", unquoted_underscores=True, suppress_rooting=True)
            self.tre.write(path="updt_tre.tre", schema="newick", unquoted_underscores=True, suppress_rooting=True)
            # phylogenetic_helpers.evaluate_raxmlng_alignment()

    def update_tree(self, aln_fn, best_subst_model, num_threads):
        """Estimate tree, using correct substitution model and a starting tree.

        The backbone options allows to fix the skeleton of the starting tree and just newly estimates the other parts.

        :param aln_fn: file name of input alignment
        :param best_subst_model: substitution model
        :param num_threads: number of threads
        :return:
        """
        print("est full tree ng")
        with cd(self.config.workdir):
            seed = str(random.randint(1, 21))

            subprocess.call(["raxml-ng-mpi", '--check', '--msa', "{}".format(aln_fn),
                             '--model', best_subst_model, '--prefix', 'check'])
            if self.config.update_tree is True:
                if self.config.backbone is not True:
                    # print('tree')
                    subprocess.call(["raxml-ng-mpi", "--threads", "{}".format(num_threads),
                                     "--msa", "{}".format(aln_fn), "--tree", "epa_result.newick",
                                     "--seed", "{}".format(seed), "--prefix", "fulltree"])
                else:
                    subprocess.call(["raxml-ng-mpi", "--threads", "{}".format(num_threads), "--msa {}".format(aln_fn),
                                     "--constraint-tree", "backbone.tre", "--seed", "{}".format(seed),
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
        aln_fn = 'updt_aln.fasta'

        best_subst_model = phylogenetic_helpers.run_modeltest(aln_fn, self.config.workdir, self.config.modeltest_criteria)
        num_threads = phylogenetic_helpers.estimate_number_threads_raxml(self.config.workdir, aln_fn, best_subst_model)

        if self.config.backbone is True:
            self.update_tree(aln_fn, best_subst_model, num_threads)
        else:
            self.calculate_bootstrap_ng(best_subst_model, num_threads, aln_fn)

    def write_labelled(self, treepath):
        """ Output tree with human readable labels.

        :param treepath: full file name (including path) for phylogeny
        :return: writes out labelled phylogeny and alignment to file
        """
        #print("write labelled files")
        if treepath is None:
            treepath = os.path.join(self.config.workdir, "fulltree.raxml.bestTree")
        else:
            treepath = os.path.join(self.config.workdir, treepath)
        phylogenetic_helpers.replace_uid_with_name(treepath, self.table, 'tree')


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
        assert type(self.aln) == DnaCharacterMatrix

        if tre_fn is not None:
            self.tre = self.write_clean_tre(tre_fn)
            if self.config.different_level is False:
                self.delete_missing()  # turned of for different level, as tre is not updated between runs, aln is.
        self.clean_inputname()

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

    def delete_missing(self):
        """ Remove taxa if only present in tree or aln.
        """
        print("Delete missing")
        treed_taxa = set(leaf.taxon for leaf in self.tre.leaf_nodes())
        aln_tax = set(taxon for taxon in self.aln)
        delete_tax = treed_taxa ^ aln_tax

        del_aln = []
        del_tre = []
        for taxon in delete_tax:
            assert (taxon in aln_tax) or (taxon in treed_taxa)
            if taxon in aln_tax:
                del_aln.append(taxon)
            if taxon in treed_taxa:
                del_tre.append(taxon)
        self.aln.discard_sequences(del_aln)
        self.tre.prune_taxa(del_tre)

        # remove associated taxa from namespace
        if del_aln != []:
            for item in del_aln:
                self.aln.taxon_namespace.remove_taxon_label(item.label)  # can only be removed one by one
                msg = '{} is only in the alignment, not in the tree and will be deleted.\n'.format(item.label)
                write_msg_logfile(msg, self.config.workdir)
        if del_tre != []:
            for item in del_tre:
                self.tre.taxon_namespace.remove_taxon_label(item.label)  # can only be removed one by one
                msg = '{} is only in the tree, not in the alignment and will be deleted.\n'.format(item.label)
                write_msg_logfile(msg, self.config.workdir)
        if len(delete_tax) > 0:
            true_false = np.where((self.table.accession.isin(list(delete_tax))), True, False)
            self.table.at[true_false, 'status'] = -1
            self.table.at[true_false, 'status_note'] = "deleted, was missing in aln or tre"
        assert self.aln.taxon_namespace == self.tre.taxon_namespace

    def clean_inputname(self):
        """It rewrites tip names if they start with a number at the beginning of the name.
        Python adds an 'n' to the name.

        :return: replaced tip names
        """
        tipnames = self.table['accession'].tolist()
        splitnames = []
        for item in tipnames:
            splititem = item.split('.')[0]
            splitnames.append(splititem)
        for tax in self.aln.taxon_namespace:
            if tax.label in splitnames:
                pass
            else:
                found_label = 0
                match = re.match("'n[0-9]{1,3}", tax.label)
                newname = ""
                if match:
                    newname = tax.label[2:]
                    newname = newname[:-1]
                for idx in self.table.index:
                    original = self.table.loc[idx, "accession"].split('.')[0]
                    # print(original, tax.label, newname)
                    if original == tax.label or original == newname:
                        #tax.label = self.table.loc[idx, "accession"].split('.')[0]
                        found_label = 1
                if found_label == 0:
                    sys.stderr.write("could not match tip label {} any ncbi taxon name\n".format(tax.label))

    def write_clean_aln(self, aln_fn):
        """
        Write out original and cleaned alignemnt (? converted to -, no whitespaces).

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
        aln = DnaCharacterMatrix.get(path=upd_aln_fn, schema='fasta')
        return aln

    def write_clean_tre(self, tre_fn):
        """
        Write out original and cleaned tre (no whitespaces).

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
            # use replaced tre as input
            tre = Tree.get(path=upd_tre_fn, schema='newick', taxon_namespace=self.aln.taxon_namespace,
                           preserve_underscores=True)
            return tre
