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

from PhylUp import phyl_up, config
from PhylUp import phylogenetic_helpers
import os


def run_multiple(data, confs, end, overlap_folder=None, first_locus=False):
    """
    Main function to intitiate PhylUp.
    :param data: a dictionary containing a dictionary with the input per locus
    :param confs: a list containing the path to the different configuration file(s)
    :param end: a number defining how many blast rounds shall be done. Put 0 for a single run
    :param overlap_folder: provide a folder name, where data across loci is stored;
            only needed if multiple loci are being updated
    :param first_locus: Provide name of main locus (same as used in data),
           if preferred taxa shall be focus on samples from that locus
           (usually put the locus with most samples available here)
    """
    count = 0
    for conffi in confs:
        if len(confs) > 1:
            count += 1
        print('configuration file is:')
        print(conffi)

        for locus in data.keys():
            print(locus)
            files = data[locus]

            if len(confs) > 1:
                if count > 1:
                    if os.path.exists("{}/updt_tre.tre".format(files['workdir'])):
                        files['trfn'] = "{}/updt_tre.tre".format(files['workdir'])
                    files['seqaln'] = "{}/updt_aln.fasta".format(files['workdir'])
            conf = config.ConfigObj(conffi, files['workdir'], interactive=False)
            if conf.preferred_taxa is True:
                assert overlap_folder is not None
                if end is None or end >= 1:
                    end = 1

            test = phyl_up.PhylogeneticUpdater(files['idtospn'], files['seqaln'], files['mattype'], files['trfn'],
                                               files['schema_trf'], conf, ignore_acc_list=files['ignore_acc_list'])

            test.run(status_end=end)

        if conf.preferred_taxa is True:
            print('preferred taxa: select')
            found_taxa_list = {}
            for locus in data.keys():
                files = data[locus]
                found_taxa_list[locus] = "{}/found_taxa.csv".format(files['workdir'])
            mrca = test.set_mrca(test.config.mrca_input)
            if conf.preferred_taxa_fn is None:
                if not os.path.exists(overlap_folder):
                    os.mkdir(overlap_folder)
                conf.preferred_taxa_fn = os.path.join(overlap_folder, 'overlap{}.csv'.format(mrca))
                print(conf.preferred_taxa_fn)
                assert conf.preferred_taxa_fn is not None, conf.preferred_taxa_fn
                phylogenetic_helpers.make_preferred_taxon_list(found_taxa_list, conf.preferred_taxa_fn,
                                                               overlap_complete=True)
            else:
                assert os.path.exists(conf.preferred_taxa_fn), conf.preferred_taxa_fn

            for locus in data.keys():
                files = data[locus]

                if len(confs) > 1:
                    if count > 1:
                        if os.path.exists("{}/updt_tre.tre".format(files['workdir'])):

                            files['trfn'] = "{}/updt_tre.tre".format(files['workdir'])
                        files['seqaln'] = "{}/updt_aln.fasta".format(files['workdir'])

                conf = config.ConfigObj(conffi, files['workdir'], interactive=False)
                if conf.preferred_taxa_fn is None:
                    conf.preferred_taxa_fn = os.path.join(overlap_folder,  'overlap{}.csv'.format(mrca))
                    print(conf.preferred_taxa_fn)
                    assert conf.preferred_taxa_fn is not None, conf.preferred_taxa_fn

                test = phyl_up.PhylogeneticUpdater(files['idtospn'], files['seqaln'], files['mattype'], files['trfn'],
                                                   files['schema_trf'], conf, ignore_acc_list=files['ignore_acc_list'])
                test.run(status_end=end)

                # new implementation - restricts overlap to species from first locus.
                if first_locus is True:
                    tf = test.table['status'] >= 0
                    print(tf)
                    added_taxa = test.table[tf]
                    print(added_taxa['ncbi_txid'])
                    added_taxa['ncbi_txid'].to_csv(conf.preferred_taxa_fn, index=False)
                first_locus = False
