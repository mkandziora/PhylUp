# Install Requirements:
# make folder where you want to install the needed phylogenetic methods

apt-get update
apt-get install -y libjson-perl
apt-get install -y mafft


###################################################
#install papara
mkdir PaPaRa
cd PaPaRa
wget https://sco.h-its.org/exelixis/resource/download/software/papara_nt-2.5-static_x86_64.tar.gz
gunzip  -cd papara_nt-2.5-static_x86_64.tar.gz | (tar xvf - )
echo export PATH="$PATH:$(pwd)" >> ~/.bashrc
cd ..


###################################################
# install EPA-NG
sudo apt-get install -y autotools-dev libtool flex bison cmake automake autoconf

mkdir EPA-ng
cd EPA-ng
wget https://github.com/Pbdas/epa-ng/archive/master.zip
unzip master.zip
cd epa-ng-master/
make
cd bin
echo export PATH="$PATH:$(pwd)" >> ~/.bashrc
cd ..
cd ..
cd ..


###################################################
# install gappa (for EPA)
mkdir gappa
cd gappa
wget https://github.com/lczech/gappa/archive/master.zip
unzip master.zip
cd gappa-master
make
cd ..
cd ..

###################################################
# raxml-ng
mkdir RAxML-ng
cd RAxML-ng
wget https://github.com/amkozlov/raxml-ng/releases/download/0.9.0/raxml-ng_v0.9.0_linux_x86_64_MPI.zip
unzip raxml-ng_v0.9.0_linux_x86_64_MPI.zip
sudo apt-get install -y libgmp3-dev

mkdir build && cd build
cmake -DUSE_MPI=ON ..
make
cd ..

#rm -r build
#mkdir build && cd build
#cmake ..
#make
#cd ..

cd bin
echo export PATH="$PATH:$(pwd)" >> ~/.bashrc

cd ..
cd ..


###################################################
# install modeltest-ng
mkdir modeltest-ng
cd modeltest-ng
git clone --recursive https://github.com/ddarriba/modeltest
cd modeltest
mkdir build
cd build
cmake ..
make
echo export PATH="$PATH:$(pwd)" >> ~/.bashrc

cd ..
cd ..
cd ..
cd ..

###################################################
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz
tar -xvf ncbi-blast-2.9.0+-x64-linux.tar.gz

wget 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
gunzip  -cd taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)
mv *.dmp ./data/


###################################################
source ~/.bashrc


##################################################

# install the Genbank databse
mkdir Genbank_database
cd Genbank_database
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.* # this downloads all nt-compressed files
cat *.tar.gz | tar -xvzf - -i # macOS tar does not support the -i flag, you need to use homebrew to brew install gnu-tar and replace the tar command by gtar
blastdbcmd -db nt -info # checks if it works
rm *.tar.gz*

# Install the taxonomy database:
wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz' # Download the taxdb archive
gunzip -cd taxdb.tar.gz | (tar xvf - ) # Install it in the BLASTDB directory
rm *.tar.gz*

cd ..
