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

# load and update Genbank
import os
import sys
import datetime
from . import get_user_input, cd
#from . import ncbi_data_parser
import ncbiTAXONparser.ncbi_data_parser as ncbi_data_parser

# Get the different needed databases.

def _download_localblastdb(config):
    """Check if files are present and if they are uptodate.
    If not files will be downloaded.
    """
    print(config.blastdb_path)
    if not os.path.isfile("{}/nt.22.nhr".format(config.blastdb_path)):
        print("Do you want to download the blast nt databases from ncbi? Note: "
              "This is a US government website! You agree to their terms. \n")
        x = get_user_input()
        if x == "yes":
            os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*' "
                      "-P {}/".format(config.blastdb_path))
            os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz' "
                      "-P {}/".format(config.blastdb_path))
            with cd(config.blastdb_path):
                print("update blast db")
                os.system("update_blastdb.pl nt")
                os.system("cat *.tar.gz | tar -xvzf - -i")
                os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                os.system("rm *.tar.gz*")
        else:
            sys.stderr.write("You have no nt database, which is needed to run PhylUp. Restart and type 'yes'. \n")
            sys.exit(-10)
    else:
        download_date = os.path.getmtime("{}/nt.22.nhr".format(config.blastdb_path))
        download_date = datetime.datetime.fromtimestamp(download_date)
        today = datetime.datetime.now()
        time_passed = (today - download_date).days
        if time_passed >= 90:
            print("""Your blast databases might not be up to date anymore. 
                  You downloaded them {} days ago. Do you want to update the blast databases from ncbi? 
                  Note: This is a US government website! You agree to their terms.\n""".format(time_passed))
            x = get_user_input()
            if x == "yes":
                with cd(config.blastdb_path):
                    os.system("update_blastdb.pl nt")
                    os.system("cat *.tar.gz | tar -xvzf - -i")
                    os.system("update_blastdb.pl taxdb")
                    os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                    os.system("rm *.tar.gz*")
            elif x == "no":
                print("You did not agree to update data from ncbi. Old database files will be used.")
            else:
                print("You did not type 'yes' or 'no'!")


def _download_ncbi_parser(config):
    """Check if files are present and if they are up to date.
    If not files will be downloaded.
    """

    if not os.path.isfile(config.ncbi_parser_nodes_fn):
        print("Do you want to download taxonomy databases from ncbi? Note: This is a US government website! "
              "You agree to their terms. \n")
        x = get_user_input()
        if x == "yes":
            os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./data/")
            os.system("gunzip -f -cd ./data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
            os.system("mv nodes.dmp ./data/")
            os.system("mv names.dmp ./data/")
            os.system("rm ./data/taxdump.tar.gz")
        else:
            sys.stderr.write("You have no taxonomic database, which is needed to run PhylUp. "
                             "Restart and type 'yes'. \n")
            sys.exit(-10)
    else:
        # x = 'no'
        download_date = os.path.getmtime(config.ncbi_parser_nodes_fn)
        download_date = datetime.datetime.fromtimestamp(download_date)
        today = datetime.datetime.now()
        time_passed = (today - download_date).days
        # print(time_passed)
        if time_passed >= 90:
            print("Do you want to update taxonomy databases from ncbi? Note: This is a US government website! "
                  "You agree to their terms. \n")
            x = get_user_input()
            if x == "yes":
                os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./data/")
                os.system("gunzip -f -cd ./data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                os.system("mv nodes.dmp ./data/")
                os.system("mv names.dmp ./data/")
                os.system("rm ./data/taxdump.tar.gz")
            elif x == "no":
                print("You did not agree to update data from ncbi. Old database files will be used.")
            else:
                print("You did not type yes or no!")
    # if x == 'yes':
    #     ncbi_data_parser.make_lineage_table()


# def _download_edirect():
#     os.system("sh -c '$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)'")
#     os.system("source ./install-edirect.sh")
