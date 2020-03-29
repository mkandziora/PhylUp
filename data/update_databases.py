from . import db_updater, config
from ncbiTAXONparser.ncbi_data_parser import _download_ncbi_parser

configfi = "data/localblast.config"
workdir = 'test/output/updt_db'

conf = config.ConfigObj(configfi, workdir, interactive=True)

db_updater._download_localblastdb(conf)
_download_ncbi_parser(conf)