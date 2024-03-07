#!/bin/bash

mkdir temp
make

# download and install slaMEM
curl -L -O https://github.com/fjdf/slaMEM/archive/master.tar.gz
tar -xf master.tar.gz
rm master.tar.gz
mv slaMEM-master slaMEM
cd slaMEM
make
cd ..

# download and install MAFFT
curl -L -O https://mafft.cbrc.jp/alignment/software/mafft-7.520-linux.tgz
tar xfvz mafft-7.520-linux.tgz
rm mafft-7.520-linux.tgz
mv mafft-linux mafft
