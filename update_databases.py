from PhylUp import db_updater, config
import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser

configfi = "data/localblast.config"
workdir = 'test/output/updt_db'

conf = config.ConfigObj(configfi, workdir, interactive=True)

db_updater._download_localblastdb(conf)

ncbi_data_parser.Parser(names_file=config.ncbi_parser_names_fn,
                        nodes_file=config.ncbi_parser_nodes_fn)



ncbi_data_parser._download_ncbi_parser(conf)