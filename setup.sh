#!/usr/bin/env bash

set -e
set -o pipefail

cwd=`pwd`

mkdir -p downloads
mkdir -p bin

cd downloads

#wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz 
#tar xf bedtools-2.27.1.tar.gz 
#rm bedtools-2.27.1.tar.gz
#cd bedtools2
#make
#cp bin/bedtools bin/intersectBed ${cwd}/bin/
#cd ${cwd}

cd downloads
xpip=`which pip`
python3 ${xpip} install snakemake
cd ${cwd}
