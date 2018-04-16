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

cwd=`pwd`
mkdir -p downloads
mkdir -p bin

# VEP
cd ${cwd}/downloads
export PERL_MM_USE_DEFAULT=1
wget https://github.com/Ensembl/ensembl-vep/archive/release/92.1.tar.gz
tar xf 92.1.tar.gz
rm 92.1.tar.gz
cd ensembl-vep-release-92.1
perl INSTALL.pl --NO_HTSLIB --AUTO a 
ln -sf `pwd`/vep ${cwd}/bin/
cd ${cwd}
bin/vep --help
bin/vep --database --vcf -i downloads/ensembl-vep-release-92.1/examples/homo_sapiens_GRCh38.vcf -o test.vep.vcf
rm test.vep.vcf*

#########
exit 0
#########

# R packages
cd ${cwd}
Rscript install/install_pkgs.R

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


