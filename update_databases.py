from PhylUp import db_updater, config

configfi = "data/localblast.config"
workdir = 'test/output/updt_db'

conf = config.ConfigObj(configfi, workdir, interactive=False)

db_updater._download_localblastdb(conf)
db_updater._download_ncbi_parser(conf)