import os

os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./data/")
os.system("gunzip -f -cd ./data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
os.system("mv nodes.dmp ./data/")
os.system("mv names.dmp ./data/")
os.system("rm ./data/taxdump.tar.gz")