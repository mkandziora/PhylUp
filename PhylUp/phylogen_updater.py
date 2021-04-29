"""
PhylUp: phylogenetic alignment building with custom taxon sampling
Copyright (C) 2020  Martha Kandziora
martha.kandziora@mailbox.org

Package to automatically generate alignments (or update alignments and phylogenies)
using local sequences or a local Genbank database
while controlling for the number of sequences per OTU and taxonomic rank.

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
from dendropy import Tree, DnaCharacterMatrix

import numpy as np
import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser  # it is the ncbi data parser class and associated functions

from . import cd
from . import write_msg_logfile
from . import suppress_stdout
from . import phylogenetic_helpers
from . import debug


class AlnUpdater(object):
    """
    Class to update aln with newly retrieved sequences.
    """
    def __init__(self, config, aln, table, status_new=None, tre=None):
        """
        :param tre: dendropy tree object
        :param aln: dendropy alignment object
        :param table: pd df table with all the information about the seqs, old and new
        :param config: config object
        """
        self.aln = aln
        self.table = table
        self.tre = tre
        self.config = config
        self.newseqs_file = "new_seqs.fasta"
        self.tre_fn = "updt_aln.fasta.tree"

        if status_new is None:
            # gets here if not unpublished
            self.new_seq_table = self.table[self.table['status'] >= 1]  # gets all new seqs (status>0.5) -
        elif status_new == 0.5:
            self.new_seq_table = self.table[self.table['status'] > 0]  # gets all new seqs (status>0.5)
        else:
            self.new_seq_table = self.table[self.table['status'] >= status_new]  # gets all new seqs (status>0.5)
        self.update_data()

    def update_data(self):
        """
        Key method of the class, cleans and updates tre and aln.

        :return: output files
        """
        print('update aln data')
        print(len(self.new_seq_table.index))
        self.delete_short_seqs()

        if len(self.new_seq_table) > 0:

            self.write_papara_queryseqs()
            if len(self.aln) > 1:
                if self.tre is None:  # generate random tree, e.g. from modeltest
                    phylogenetic_helpers.run_modeltest('updt_aln.fasta', self.config.workdir,
                                                       self.config.modeltest_criteria)
                    self.tre = Tree.get(path=os.path.join(self.config.workdir, self.tre_fn),
                                        schema="newick", preserve_underscores=True)
                    phylogenetic_helpers.write_papara_trefile(self.tre, self.config.workdir)
                else:
                    phylogenetic_helpers.write_aln(self.aln, self.config.workdir)
                    phylogenetic_helpers.write_tre(self.tre, self.config.workdir)
                self.add_query_seqs_to_aln()
            else:
                self.add_queryseqs_to_singleseq()
                phylogenetic_helpers.run_modeltest('updt_aln.fasta', self.config.workdir,
                                                   self.config.modeltest_criteria)
                self.tre = Tree.get(path=os.path.join(self.config.workdir, 'updt_aln.fasta.tree'),
                                    schema="newick", preserve_underscores=True)
            phylogenetic_helpers.write_tre(self.tre, self.config.workdir)
        phylogenetic_helpers.write_aln(self.aln, self.config.workdir)

        # self.aln = self.trim(os.path.abspath(os.path.join(self.config.workdir, 'papara_alignment.phylip')),
        #                      format_aln='phylip')
        self.write_labelled('updt_aln.fasta')

    def add_queryseqs_to_singleseq(self):
        """
        Add query sequences to alignment using mafft.

        :return:
        """
        print('add query seq to single')
        phylogenetic_helpers.make_mafft_aln(self.aln, self.config.workdir)
        cmd_mafft = 'mafft --genafpair --leavegappyregion --maxiterate 16 --thread {} --reorder ' \
                    '{}/mafft.fasta > {}/mafft_align.fasta'.format(self.config.num_threads, self.config.workdir,
                                                                   self.config.workdir)
        subprocess.Popen(cmd_mafft, shell=True, stdout=subprocess.PIPE).stdout.read()
        self.aln = DnaCharacterMatrix.get(path=os.path.join(self.config.workdir, "mafft_align.fasta"), schema='fasta')
        self.aln = self.trim(os.path.join(self.config.workdir, 'mafft_align.fasta'), 'fasta')
        phylogenetic_helpers.write_aln(self.aln, self.config.workdir)

    def add_query_seqs_to_aln(self):
        """Runs papara to add new sequences to the alignment.

        :return: writes new extended alignment
        """
        sys.stdout.write("Aligning query sequences. \n")
        for filename in glob.glob('{}/papara*'.format(self.config.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.config.workdir, filename.split("/")[-1]))
        phylogenetic_helpers.write_papara_alnfile(self.aln, self.config.workdir)
        if self.tre is not None:
            phylogenetic_helpers.write_papara_trefile(self.tre, self.config.workdir)
        cwd = os.getcwd()

        try:
            os.chdir(self.config.workdir)
            # with cd(self.config.workdir):
            print('run papara')
            phylogenetic_helpers.run_papara()
            path = os.path.join(os.getcwd(), "papara_alignment.phylip".format())
            assert os.path.exists(path), "{} does not exists".format(path)
            os.chdir(cwd)
            aln_old = self.aln
            aln_path = os.path.join(self.config.workdir, "papara_alignment.phylip")
            self.aln = DnaCharacterMatrix.get(path=aln_path, schema='phylip')
            self.aln = self.trim(os.path.join(self.config.workdir, 'papara_alignment.phylip'), 'phylip')

            print('papara done and files loaded')

        except subprocess.CalledProcessError:
            os.chdir(cwd)
            print('Papara failed - using MAFFT now')
            aln_old = self.aln
            self.add_queryseqs_to_singleseq()
        # aln_old = self.aln
        # aln_path = os.path.join(self.config.workdir, "papara_alignment.phylip")
        with cd(self.config.workdir):
            phylogenetic_helpers.truncate_papara_aln(aln_old)
        phylogenetic_helpers.write_aln(self.aln, self.config.workdir)

        msg = "Following papara alignment, aln has {} seqs \n".format(len(self.aln))
        write_msg_logfile(msg, self.config.workdir)

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
            if len(seq.symbols_as_string().replace("-", "").replace("?", "")) <= seq_len_cutoff:
                delete_seqs.append(tax)
        if delete_seqs:
            msg = "Taxa deleted from tree and alignment in delete short sequences " \
                  "as sequences are shorter than {}.\n".format(seq_len_cutoff)
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
        # self.tre.write(path=os.path.join(self.config.workdir, 'updt_tre.tre'), schema='newick',
        #                unquoted_underscores=True, suppress_rooting=True)
        self.table.to_csv(os.path.join(self.config.workdir, 'table.updated'), index=False)

    def remove_taxa_aln_tre(self, taxon_label):
        """Removes taxon from aln and tre, takes a single taxon_label as input.

        :param taxon_label: taxon_label from dendropy object - aln or phy
        :return:
        """
        if self.tre is not None:
            # not sure why this function exist. None of them actually remove a tip.
            # tax2 = self.tre.taxon_namespace.get_taxon(taxon_label)
            # self.tre.prune_taxa([tax2])
            # self.tre.prune_taxa_with_labels([taxon_label])
            # self.tre.prune_taxa_with_labels([tax2])

            leaf_set = set(leaf.taxon for leaf in self.tre.leaf_nodes())
            for leaf in leaf_set:
                if leaf.label == taxon_label:
                    self.tre.prune_taxa([leaf])

        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        self.aln.remove_sequences([tax])
        self.aln.discard_sequences([tax])
        self.aln.taxon_namespace.remove_taxon_label(taxon_label)  # raises an error if label not found

    def trim(self, aln_fn=False, format_aln=None):
        """
        Removes bases at the start and end of alignment, if they are represented by less than the value
        specified in the config file (trim_perc).

        Means, that trim_perc (e.g. 0.75 = 75%) of the sequences need to have a base present.
        This ensures, that not whole chromosomes get dragged in by cutting the ends of long sequences.
        """
        print('trim alignment')
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
        self.config = config
        if self.tre is None:  # generate random tree, e.g. from modeltest
            self.tre = Tree.get(path=os.path.join(self.config.workdir, "updt_aln.fasta.tree"), schema="newick",
                                preserve_underscores=True, taxon_namespace=self.aln.taxon_namespace)

        self.update_phyl()

    def update_phyl(self):
        """
        Main function of the class. Does the updating of the phylogeny.
        :return:
        """
        print("update phylogeny")

        # this exists as from single seq there are no old seqs to place (no old_seqs.fasta) as no alignment is available
        if os.path.exists(os.path.join(self.config.workdir, "old_seqs.fasta")):
            self.place_query_seqs_epa()
            self.check_tre_in_aln()

        if self.config.update_tree is True:
            self.calculate_final_tree()  # comment out for development speed up
            # self.tre = Tree.get(path=os.path.join(self.config.workdir, "fulltree.raxml.bestTree"),
            #                     schema="newick", preserve_underscores=True, taxon_namespace=self.aln.taxon_namespace)
            # fn = os.path.join(self.config.workdir, 'papara_alignment.phylip_trim',
            #                   taxon_namespace=self.tre.taxon_namespace)
            # self.aln = DnaCharacterMatrix.get(path=fn, schema="phylip")
            # self.aln.write(path=os.path.join(self.config.workdir, "updt_aln.fasta"), schema='fasta')
            shutil.copy(os.path.join(self.config.workdir, "fulltree.raxml.bestTree"),
                        os.path.join(self.config.workdir, "updt_tre.tre"))
            self.write_labelled('updt_tre.tre')
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
        # assert aln_ids.issubset(treed_taxa), (len(treed_taxa), len(aln_ids))

    def place_query_seqs_epa(self):
        """ Runs epa-ng on the tree, and the combined alignment including the new query seqs.
        Just for placement, to use as starting tree.
        """
        print("place query seq")
        if len(self.tre.taxon_namespace) > 2:

            phylogenetic_helpers.write_papara_trefile(self.tre, self.config.workdir)

            # prepare tree
            self.tre.resolve_polytomies()
            self.tre.deroot()
            with open(os.path.join(self.config.workdir, "epa_tree.tre"), "w") as tre_file:
                tre_file.write("{}".format(self.tre.as_string(schema='newick', suppress_rooting=True)))

            with cd(self.config.workdir):
                # todo: get model for epa-ng
                with suppress_stdout():
                    subprocess.call(["epa-ng", "--ref-msa", "old_seqs.fasta", "--tree", "epa_tree.tre", "--query",
                                     "new_seqs.fasta", "--outdir", "./", "--model", 'GTR+G', '--redo'], shell=False)
                    # make newick tre
                    subprocess.call(["gappa", "examine", "graft", "--jplace-path", "epa_result.jplace",
                                     "--allow-file-overwriting"], shell=False)
                self.tre = Tree.get(path="epa_result.newick", schema="newick", preserve_underscores=True)
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
                             '--model', best_subst_model, '--prefix', 'check'], shell=False)
            if self.config.update_tree is True:
                if self.config.backbone is not True:
                    # print('tree')
                    subprocess.call(["raxml-ng-mpi", "--threads", "{}".format(num_threads), '--model', best_subst_model,
                                     "--msa", "{}".format(aln_fn), "--tree", "epa_result.newick",
                                     "--seed", "{}".format(seed), "--prefix", "fulltree"], shell=False)
                else:
                    print('backbone')
                    subprocess.call(["raxml-ng-mpi", "--threads", "{}".format(num_threads), '--model', best_subst_model,
                                     "--msa", aln_fn, "--tree-constraint", "backbone.tre",
                                     "--seed", "{}".format(seed), "--prefix", "fulltree"], shell=False)
            else:
                todo = 'To update the data run the following command in your working directory.'
                if self.config.backbone is not True:
                    cmd1 = "raxml-ng-mpi --threads {} --model {} --msa {} --tree epa_result.newick --seed {}" \
                           " --prefix fulltree".format(num_threads, best_subst_model, aln_fn, seed)
                else:
                    cmd1 = "raxml-ng-mpi --threads {} --model {} --msa {} --tree-constraint updt_tre.tre --seed {}" \
                           " --prefix fulltree".format(num_threads, best_subst_model, aln_fn, seed)
                sys.stdout.write(todo)
                sys.stdout.write(cmd1)

    def calculate_bootstrap_ng(self, best_subst_model, num_threads, aln_fn='updt_aln.fasta'):
        """Calculates bootstrap and consensus trees.

        :param best_subst_model:
        :param num_threads:
        :param aln_fn:
        :return:
        """
        sys.stdout.write("calculate bootstrap \n")
        with cd(self.config.workdir):
            # check if job was started with mpi: this checks if actual several cores and nodes were allocated
            ntasks = os.environ.get('SLURM_JOB_CPUS_PER_NODE')  # or ntasks = os.environ.get('SLURM_NTASKS_PER_NODE')
            nnodes = os.environ.get("SLURM_JOB_NUM_NODES")
            seed = str(random.randint(1, 21))
            mpi = False
            if self.config.update_tree is True:
                if nnodes is not None and ntasks is not None:
                    mpi = True
                if mpi:  # todo add "mpiexec", "-n", "{}".format(env_var)
                    # try:
                    #     print("run with mpi")
                    #     env_var = int(nnodes) * int(ntasks)
                    #     subprocess.call(["mpiexec", "-n", "{}".format(env_var), "raxml-ng-mpi", '--all', "--msa",
                    #                      "{}".format(aln_fn), '--model', "GTR+G", '--bs-trees', 'autoMRE',
                    #                      '--seed', seed, "--threads", "{}".format(str(self.config.num_threads)),
                    # except:
                    subprocess.call(["raxml-ng-mpi", '--all', "--msa", "{}".format(aln_fn), '--model', best_subst_model,
                                     '--bs-trees', 'autoMRE', '--seed', seed, "--threads", "{}".format(num_threads),
                                     "--prefix", "fulltree"], shell=False)  # 'tbe', #
                else:
                    debug('else')

                    subprocess.call(["raxml-ng-mpi", '--all', "--msa", "{}".format(aln_fn), '--model', best_subst_model,
                                     'autoMRE', '--seed', seed, "--threads", "{}".format(num_threads),
                                     "--prefix", "fulltree"], shell=False)  # 'tbe', #
                # subprocess.call(["raxml-ng-mpi", '--support', '--tree', 'fulltree.raxml.bestTree', '--bs-trees',
                #                 'fulltree.raxml.bootstraps', "--prefix", 'support'])
                print('make consensus')
                subprocess.call(["raxml-ng-mpi", '--consense', 'MRE', '--tree', 'fulltree.raxml.bootstraps',
                                 "--prefix", 'consMRE'], shell=False)
                subprocess.call(["raxml-ng-mpi", '--consense', 'STRICT', '--tree', 'fulltree.raxml.bootstraps',
                                 "--prefix", 'consSTRICT'], shell=False)
                subprocess.call(["raxml-ng-mpi", '--consense', 'MR', '--tree', 'fulltree.raxml.bootstraps',
                                 "--prefix", 'consMR'], shell=False)
            else:
                todo = 'To update the data run the following command in your working directory.'
                cmd1 = "raxml-ng-mpi  --all --threads {} --msa {} --model {} --bs-trees autoMRE --seed {}" \
                       " --prefix fulltree".format(num_threads, aln_fn, best_subst_model, seed)
                cmd2 = "raxml-ng-mpi --consense MRE --tree fulltree.raxml.bootstraps --prefix consMRE"
                cmd3 = "raxml-ng-mpi --consense STRICT --tree fulltree.raxml.bootstraps --prefix consSTRICT"
                cmd4 = "raxml-ng-mpi --consense MR --tree fulltree.raxml.bootstraps --prefix consMR"

                sys.stdout.write(todo)
                sys.stdout.write(cmd1)
                sys.stdout.write(cmd2)
                sys.stdout.write(cmd3)
                sys.stdout.write(cmd4)

                lfd = os.path.join(self.config.workdir, "logfile")
                with open(lfd, "a") as log:
                    log.write("{}\n".format(todo))
                    log.write("{}\n".format(cmd1))
                    log.write("{}\n".format(cmd2))
                    log.write("{}\n".format(cmd3))
                    log.write("{}\n".format(cmd4))

    def calculate_final_tree(self):
        """Calculates the final tree using a trimmed alignment.

        :return: final data
        """
        sys.stdout.write("calculate final tree")
        aln_fn = 'updt_aln.fasta'

        best_subst_model = phylogenetic_helpers.run_modeltest(aln_fn, self.config.workdir,
                                                              self.config.modeltest_criteria)
        num_threads = phylogenetic_helpers.estimate_number_threads_raxml(self.config.workdir, aln_fn, best_subst_model)

        if self.config.backbone is True:
            self.update_tree(aln_fn, best_subst_model, num_threads)
        else:
            self.calculate_bootstrap_ng(best_subst_model, num_threads, aln_fn)

    def write_labelled(self, treepath):
        """ Output tree with human readable labels.

        :param treepath: file name of tree
        :return: writes out labelled phylogeny to file
        """
        debug("write labelled files")
        if treepath is None:
            treepath = os.path.join(self.config.workdir, "fulltree.raxml.bestTree")
        else:
            treepath = os.path.join(self.config.workdir, treepath)
        phylogenetic_helpers.replace_uid_with_name(treepath, self.table, 'tree')


class InputCleaner(object):
    """
    This is the input class, that cleans the data before updating the phylogeny.
    """
    def __init__(self, tre_fn, tre_schema, aln_fn, aln_schema, table, config_obj):  # removed mrca
        sys.stdout.write('Clean the input data: {}, {}.'.format(tre_fn, aln_fn))
        self.config = config_obj
        if not os.path.exists(self.config.workdir):
            os.makedirs(self.config.workdir)
        # ncbi parser contains information about spn, tax_id, and ranks
        self.ncbi_parser = ncbi_data_parser.Parser(names_file=config_obj.ncbi_parser_names_fn,
                                                   nodes_file=config_obj.ncbi_parser_nodes_fn)
        # # set mrca id
        # assert type(mrca) in [int, list, set] or mrca is None, ("mrca is not an integer, list, set"
        #                                                         " or None: {}".format(mrca))
        # self.mrca = self.format_mrca_set(mrca)

        self.table = table
        self.aln = self.write_clean_aln(aln_fn, aln_schema)
        assert isinstance(self.aln, DnaCharacterMatrix), (type(self.aln))

        if tre_fn is not None:
            if tre_schema is None:
                tre_schema = 'newick'
            self.tre = self.write_clean_tre(tre_fn, tre_schema)
            if self.config.different_level is False:
                self.delete_missing()  # turned of for different level, as tre is not updated between runs, aln is.
            if self.config.backbone is True:
                phylogenetic_helpers.write_tre(self.tre, self.config.workdir, treepath="backbone.tre",
                                               treeschema="newick")
                with cd(self.config.workdir):
                    backbonetre = Tree.get(path="backbone.tre", schema="newick", preserve_underscores=True)
                backbonetre.resolve_polytomies()
                phylogenetic_helpers.write_tre(backbonetre, self.config.workdir, treepath="backbone.tre",
                                               treeschema="newick")
        self.clean_inputname()
        phylogenetic_helpers.write_aln(self.aln, self.config.workdir)
        if tre_fn is not None:
            phylogenetic_helpers.write_tre(self.tre, self.config.workdir)

    def format_mrca_set(self, mrca):
        """
        Takes the input and makes a set of valid mrca ids.

        :param mrca:
        :return:
        """
        assert type(mrca) in [int, list, set] or mrca is None, ("mrca is not an integer, list, set"
                                                                " or None: {}".format(mrca))
        mrca_name = self.ncbi_parser.get_name_from_id(list(mrca)[0])
        sys.stdout.write('Format mrca: {} - {}'.format(mrca, mrca_name))
        if isinstance(mrca, int):
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
            sys.exit(-3)
        # assert type(mrca) in {list, set, int}, ("your ingroup_mrca '%s' is not an integer/list/set." % mrca)
        return mrca

    def delete_missing(self):
        """ Remove taxa if only present in tree or aln.
        """
        debug("Delete missing taxa")
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
                self.table.loc[self.table['accession'] == item.label, 'status'] = -1
                self.table.loc[self.table['accession'] == item.label, 'status_note'] = 'not in the treer'
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

    # todo: does nothing? needed?
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
                # found_label = 0
                match = re.match("'n[0-9]{1,3}", tax.label)
                newname = ""
                if match:
                    newname = tax.label[2:]
                    newname = newname[:-1]
                for idx in self.table.index:
                    original = self.table.loc[idx, "accession"].split('.')[0]
                    if original == tax.label or original == newname:
                        tax.label = self.table.loc[idx, "accession"].split('.')[0]
                        found_label = 1
                if found_label == 0: # and self.table.loc[idx, "ncbi_txid"]:
                   sys.stderr.write("could not match tip label {} any ncbi taxon name\n".format(tax.label))

    def write_clean_aln(self, aln_fn, aln_schema):
        """
        Write out original and cleaned alignemnt (? converted to -, no whitespaces).

        :param aln_fn: filename of alignment
        :param aln_schema: format of alignment
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
        aln = DnaCharacterMatrix.get(path=upd_aln_fn, schema=aln_schema)

        # delete missing data only
        delete_seqs = []

        for tax in aln:
            seq = aln[tax].symbols_as_string()
            if seq == '-'*len(seq):
                delete_seqs.append(tax)
        if delete_seqs:
            print(delete_seqs)
            msg = "Taxa deleted from alignment: missing data only.\n"
            write_msg_logfile(msg, self.config.workdir)
            for tax in delete_seqs:
                aln.remove_sequences([tax])
                aln.discard_sequences([tax])
                aln.taxon_namespace.remove_taxon_label(tax.label)
                msg = '{}, '.format(tax.label)
                write_msg_logfile(msg, self.config.workdir)
                self.table.at[self.table['accession'] == tax.label, "status"] = -1
                self.table.at[self.table['accession'] == tax.label, "status_note"] = "deleted - missing data only"
            msg = '\n'
            write_msg_logfile(msg, self.config.workdir)
        aln.write(path=os.path.join(self.config.workdir, 'updt_aln.fasta'), schema='fasta')
        self.table.to_csv(os.path.join(self.config.workdir, 'table.updated'), index=False)
        return aln

    def write_clean_tre(self, tre_fn, tre_schema):
        """
        Write out original and cleaned tre (no whitespaces).

        :param tre_fn: tree file name
        :param tre_schema: schema of tree
        :return:
        """
        # if tre_fn is not None:
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
        tre = Tree.get(path=upd_tre_fn, schema=tre_schema, taxon_namespace=self.aln.taxon_namespace,
                       preserve_underscores=True)
        return tre
