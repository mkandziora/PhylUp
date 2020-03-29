# Install Requirements:
# make folder where you want to install the needed phylogenetic methods

apt-get update

#install papara
mkdir PaPaRa
cd PaPaRa
wget https://sco.h-its.org/exelixis/resource/download/software/papara_nt-2.5-static_x86_64.tar.gz
gunzip  -cd papara_nt-2.5-static_x86_64.tar.gz | (tar xvf - )
echo "export PATH=$PWD:\$PATH" >> ~/.bashrc
# source .bashrc
# papara
cd ..

# install EPA-NG
mkdir EPA-ng
cd EPA-ng
wget https://github.com/Pbdas/epa-ng/archive/master.zip
unzip master.zip
sudo apt-get install autotools-dev libtool flex bison cmake automake autoconf
cd epa-ng-master/
make
cd bin
echo "export PATH=$PWD:\$PATH" >> ~/.bashrc
cd ..
cd ..
cd ..

# install gappa (for EPA)
mkdir gappa
cd gappa
wget https://github.com/lczech/gappa/archive/master.zip
unzip master.zip
cd gappa-master
make
cd ..
cd ..

# raxml-ng
mkdir RAxML-ng
cd RAxML-ng
wget https://github.com/amkozlov/raxml-ng/releases/download/0.9.0/raxml-ng_v0.9.0_linux_x86_64_MPI.zip
unzip raxml-ng_v0.9.0_linux_x86_64_MPI.zip
sudo apt-get install flex bison libgmp3-dev

mkdir build && cd build
cmake -DUSE_MPI=ON ..
make
cd ..

rm -r build
mkdir build && cd build
cmake ..
make
cd ..

cd bin
echo "export PATH=$PWD:\$PATH" >> ~/.bashrctest

cd ..
cd ..

source ~/.bashrc
