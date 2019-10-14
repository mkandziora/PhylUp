# load and update Genbank
import os
import sys
import datetime
from . import get_user_input, cd
from . import ncbi_data_parser

# Get the different needed databases.

def _download_localblastdb(config):
    """Check if files are present and if they are uptodate.
    If not files will be downloaded.
    """
    if not os.path.isfile("{}/nt.69.nhr".format(config.blastdb)):
        print("Do you want to download the blast nt databases from ncbi? Note: "
              "This is a US government website! You agree to their terms. \n")
        x = get_user_input()
        if x == "yes":
            os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*'"
                      "{}/".format(config.blastdb))
            os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'"
                      "{}/".format(config.blastdb))
            with cd(config.blastdb):
                print("update blast db")
                os.system("update_blastdb nt")
                os.system("cat *.tar.gz | tar -xvzf - -i")
                os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                os.system("rm *.tar.gz*")
        else:
            sys.stderr.write("You have no nt database, which is needed to run PhyFilter. Restart and type 'yes'. \n")
            sys.exit(-10)
    else:
        download_date = os.path.getmtime("{}/nt.60.nhr".format(config.blastdb))
        download_date = datetime.datetime.fromtimestamp(download_date)
        today = datetime.datetime.now()
        time_passed = (today - download_date).days
        if time_passed >= 90:
            print("""Your databases might not be uptodate anymore. 
                  You downloaded them {} days ago. Do you want to update the blast databases from ncbi? 
                  Note: This is a US government website! You agree to their terms.\n""".format(time_passed))
            x = get_user_input()
            if x == "yes":
                with cd(config.blastdb):
                    os.system("update_blastdb nt")
                    os.system("cat *.tar.gz | tar -xvzf - -i")
                    os.system("update_blastdb taxdb")
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
            sys.stderr.write("You have no taxonomic database, which is needed to run PhyFilter. Restart and type 'yes'. \n")
            sys.exit(-10)
    else:
        download_date = os.path.getmtime(config.ncbi_parser_nodes_fn)
        download_date = datetime.datetime.fromtimestamp(download_date)
        today = datetime.datetime.now()
        time_passed = (today - download_date).days
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
    if x == 'yes':
        ncbi_data_parser.make_lineage_table()

