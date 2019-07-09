# load and update Genbank
import os
import datetime
from . import get_user_input, cd




def _download_localblastdb(self):
    """Check if files are present and if they are uptodate.
    If not files will be downloaded.
    """
    print(self.blast_loc)
    if self.blast_loc == "local":
        # next line of codes exists to have interactive mode enabled while testing
        # this allows to not actually have a local ncbi database downloaded
        if not os.path.isfile("{}/empty_local_db_for_testing.nhr".format(self.blastdb)):
            if not os.path.isfile("{}/nt.69.nhr".format(self.blastdb)):
                print("Do you want to download the blast nt databases from ncbi? Note: "
                      "This is a US government website! You agree to their terms")
                x = get_user_input()
                if x == "yes":
                    os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*'"
                              "{}/".format(self.blastdb))
                    os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'"
                              "{}/".format(self.blastdb))
                    with cd(self.blastdb):
                        print("update blast db")
                        os.system("update_blastdb nt")
                        os.system("cat *.tar.gz | tar -xvzf - -i")
                        os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                        os.system("rm *.tar.gz*")
                elif x == "no":
                    print(
                        "You did not agree to download data from ncbi. Program will default to blast web-queries.")
                    self.blast_loc = "remote"
                else:
                    print("You did not type yes or no!")
            else:
                download_date = os.path.getmtime("{}/nt.60.nhr".format(self.blastdb))
                download_date = datetime.datetime.fromtimestamp(download_date)
                today = datetime.datetime.now()
                time_passed = (today - download_date).days
                if time_passed >= 90:
                    print("""Your databases might not be uptodate anymore. 
                          You downloaded them {} days ago. Do you want to update the blast databases from ncbi? 
                          Note: This is a US government website! You agree to their terms.""".format(time_passed))
                    x = get_user_input()
                    if x == "yes":
                        with cd(self.blastdb):
                            os.system("update_blastdb nt")
                            os.system("cat *.tar.gz | tar -xvzf - -i")
                            os.system("update_blastdb taxdb")
                            os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                            os.system("rm *.tar.gz*")
                    elif x == "no":
                        print("You did not agree to update data from ncbi. Old database files will be used.")
                    else:
                        print("You did not type 'yes' or 'no'!")

def _download_ncbi_parser(self):
    """Check if files are present and if they are up to date.
    If not files will be downloaded.
    """
    if self.blast_loc == "local":
        if not os.path.isfile(self.ncbi_parser_nodes_fn):
            print("Do you want to download taxonomy databases from ncbi? Note: This is a US government website! "
                  "You agree to their terms")
            x = get_user_input()
            if x == "yes":
                os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./tests/data/")
                os.system("gunzip -f -cd ./tests/data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                os.system("mv nodes.dmp ./tests/data/")
                os.system("mv names.dmp ./tests/data/")
                os.system("rm ./tests/data/taxdump.tar.gz")
            elif x == "no":
                print("You did not agree to download data from ncbi. Program will default to blast web-queries.")
                print("This is slow and crashes regularly!")
                self.blast_loc = "remote"
            else:
                print("You did not type yes or no!")
        else:
            download_date = os.path.getmtime(self.ncbi_parser_nodes_fn)
            download_date = datetime.datetime.fromtimestamp(download_date)
            today = datetime.datetime.now()
            time_passed = (today - download_date).days
            if time_passed >= 90:
                print("Do you want to update taxonomy databases from ncbi? Note: This is a US government website! "
                      "You agree to their terms")
                x = get_user_input()
                if x == "yes":
                    os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./tests/data/")
                    os.system("gunzip -f -cd ./tests/data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                    os.system("mv nodes.dmp ./tests/data/")
                    os.system("mv names.dmp ./tests/data/")
                elif x == "no":
                    print("You did not agree to update data from ncbi. Old database files will be used.")
                else:
                    print("You did not type yes or no!")