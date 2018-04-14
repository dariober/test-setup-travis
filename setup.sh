#!/usr/bin/env bash

VERSION=0.1.0

set -e
set -o pipefail
set -x

if [[ $1 == '-h' || $1 == '--help' ]]
then
cat <<EOF
DESCRIPTION
Install script

USAGE
bash setup.sh [--user]

Version $VERSION
EOF
exit 0
fi

if [[ $1 == '-v' || $1 == '--version' ]]
then
    echo $VERSION
    exit 0
fi

user=0
for x in "$@"
do
   if [[ "$x" == "--user" ]]
   then
       user=1
   fi
done

echo "VERSION"
echo "-------"
echo $PATH
which R
R --version
which Rscript
Rscript --version

which python3
python3 --version

which python3.6
python3.6 --version

which pip3
pip3 --version
echo "-------"

cwd=`pwd`
mkdir -p downloads
mkdir -p bin

# Snakemake
if [[ $user == 1 ]]
then
    pip3 install --user snakemake
else
    pip3 install snakemake
fi

# Samtools
cd ${cwd}/downloads
wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
tar xf samtools-1.8.tar.bz2
rm samtools-1.8.tar.bz2
cd samtools-1.8
./configure --prefix=`pwd`
make
make install
cp samtools ${cwd}/bin/
cd ${cwd}
bin/samtools --version

# bedtools
cd ${cwd}/downloads
wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz 
tar xf bedtools-2.27.1.tar.gz 
rm bedtools-2.27.1.tar.gz
cd bedtools2
make
cp bin/bedtools bin/intersectBed ${cwd}/bin/
cd ${cwd}
bin/bedtools --version

# R packages
cd ${cwd}
Rscript install/install_pkgs.R

