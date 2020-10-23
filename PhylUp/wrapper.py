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

            if conf.preferred_taxa_fn == None:
                conf.preferred_taxa_fn = os.path.join(overlap_folder, 'overlap.csv')
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
