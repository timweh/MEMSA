#!/bin/bash

mkdir temp
g++ -o memsa main.cpp

os=$(uname)
if [ $os == 'Linux' ]; then
    # download and install MAFFT
    curl -L -O https://mafft.cbrc.jp/alignment/software/mafft-7.520-linux.tgz
    tar xfvz mafft-7.520-linux.tgz
    rm mafft-7.520-linux.tgz
    mv mafft-linux64 mafft
    # download and install slaMEM
    curl -L -O https://github.com/fjdf/slaMEM/archive/master.tar.gz
    tar -xf master.tar.gz
    rm master.tar.gz
    mv slaMEM-master slaMEM
    cd slaMEM
    gcc -Wall -Wextra -Wunused -mpopcnt -Wuninitialized -O9 -fomit-frame-pointer -fcommon *.c -o slaMEM -lm
elif [ $os == 'Darwin' ]; then # MacOS
    # download and install MAFFT
    curl -L -O https://mafft.cbrc.jp/alignment/software/mafft-7.520-mac.zip?signed
    unzip mafft-7.520-mac.zip
    rm mafft-7.520-mac.zip
    mv mafft-mac mafft
    # download and install slaMEM
    curl -L -O https://github.com/fjdf/slaMEM/archive/master.zip
    unzip master.zip
    rm master.zip
    mv slaMEM-master slaMEM
    cd slaMEM
    make
else
    echo 'OS not supported'
    exit 1
fi
