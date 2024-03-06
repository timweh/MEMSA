#!/bin/bash

mkdir temp

# download and install slaMEM
curl -L -O https://github.com/fjdf/slaMEM/archive/master.zip
unzip master.zip
rm master.zip
mv slaMEM-master slaMEM
cd slaMEM
make
cd ..

# download and install MAFFT
curl -L -O https://mafft.cbrc.jp/alignment/software/mafft-7.520-mac.zip?signed
unzip mafft-7.520-mac.zip
rm mafft-7.520-mac.zip
mv mafft-mac mafft
