"""
PhylUp: phylogenetic alignment building with custom taxon sampling
Copyright (C) 2020  Martha Kandziora
martha.kandziora@mailbox.org

Package to automatically update alignments and phylogenies using local sequences or a local Genbank database
while controlling for the number of sequences per OTU.

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
    count = 0
    for item in confs:
        if len(confs) > 1:
            count += 1
        conffi = item
        print('configuration file is:')
        print(item)

        for locus in data.keys():
            print(locus)
            files = data[locus]

            if len(confs) > 1:
                if count > 1:
                    if os.path.exists("{}/updt_tre.tre".format(files['workdir'])):
                        files['trfn'] = "{}/updt_tre.tre".format(files['workdir'])
                    files['seqaln'] = "{}/updt_aln.fasta".format(files['workdir'])
            conf = config.ConfigObj(conffi, files['workdir'], interactive=False)
            if conf.preferred_taxa == True:
                assert overlap_folder != None
                if end == None or end >=1:
                    end = 1

            test = phyl_up.PhylogeneticUpdater(files['idtospn'], files['seqaln'], files['mattype'], files['trfn'],
                                               files['schema_trf'], conf, ignore_acc_list=files['ignore_acc_list'])

            test.run(status_end=end)

        if conf.preferred_taxa == True:
            print('preferred taxa: select')
            found_taxa_list = {}
            for locus in data.keys():
                files = data[locus]
                found_taxa_list[locus] = "{}/found_taxa.csv".format(files['workdir'])
            mrca = test.set_mrca(test.config.mrca_input)
            if conf.preferred_taxa_fn == None:
                conf.preferred_taxa_fn = os.path.join(overlap_folder, 'overlap.csv')
                if not os.path.exists(overlap_folder):
                    os.mkdir(overlap_folder)
                conf.preferred_taxa_fn = os.path.join(overlap_folder, 'overlap{}.csv'.format(mrca))
                assert conf.preferred_taxa_fn != None, conf.preferred_taxa_fn
                phylogenetic_helpers.make_preferred_taxon_list(found_taxa_list, conf.preferred_taxa_fn, overlap_complete=True)
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
                if conf.preferred_taxa_fn == None:
                    conf.preferred_taxa_fn = os.path.join(overlap_folder,  'overlap{}.csv'.format(mrca))
                    print(conf.preferred_taxa_fn)
                    assert conf.preferred_taxa_fn != None, conf.preferred_taxa_fn

                test = phyl_up.PhylogeneticUpdater(files['idtospn'], files['seqaln'], files['mattype'], files['trfn'],
                                                   files['schema_trf'], conf, ignore_acc_list=files['ignore_acc_list'])
                test.run(status_end=end)

                # new implementation - restricts overlap to species from first locus.
                if first_locus == True:
                    tf = test.table['status'] >= 0
                    print(tf)
                    added_taxa = test.table[tf]
                    print(added_taxa['ncbi_txid'])
                    added_taxa['ncbi_txid'].to_csv(conf.preferred_taxa_fn, index=False)
                first_locus = False
