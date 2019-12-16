from Concat import concat, concat_config

workdir_its = 'Senecioneae_its'
workdir_ets = 'Senecioneae_ets'

workdir_comb = "remote/concat_nr"

genelist = {"ITS": workdir_its,
            "ETS": workdir_ets}

configfi = "./data/concat.config"  # configuration file


config = concat_config.ConcatConfigObj(configfi, workdir_comb)

conc = concat.Concat(config, workdir_comb)

conc = conc.run_phylup(genelistdict=genelist)

