import os
from shutil import copyfile, rmtree


src = "./data/data_for_tests/updt_aln.fasta"
dst = "./tests/output/test_runs/updt_aln.fasta"
if os.path.exists("./tests/output/"):
    rmtree("./tests/output/")
    
os.mkdir("./tests/output/")
os.mkdir("./tests/output/test_runs/")
copyfile(src, dst)


os.system("wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./data/")
os.system("gunzip -f -cd ./data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
os.system("mv nodes.dmp ./data/")
os.system("mv names.dmp ./data/")
os.system("rm ./data/taxdump.tar.gz")
