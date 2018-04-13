#!/usr/bin/env bash

set -e
set -o pipefail

mkdir -p downloads
mkdir -p bin

cd downloads
cwd=`pwd`

#wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz 
#tar xf bedtools-2.27.1.tar.gz 
#rm bedtools-2.27.1.tar.gz
#cd bedtools2
#make
#cp bin/bedtools ../bin/intersectBed ${cwd}/bin/
#cd ${cwd}

xpip=`which pip`
python3 ${xpip} install snakemake
cd ${cwd}
