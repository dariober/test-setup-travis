#!/usr/bin/env bash

VERSION=0.1.0

set -e
set -o pipefail

# Parse arguments
# ===============

PG=`basename "$0"`
bin_dir=${HOME}/bin

# From https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -b|--bin_dir)
        bin_dir="$2"
        shift # past argument
        shift # past value
    ;;
    -v|--version)
        version=1
        shift # past argument
    ;;
    -h|--help)
        help=1
        shift # past argument
    ;;
    *)
        echo "Unknown option in $@"  
        exit 1
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

if [[ $help == 1 ]]
then
cat <<EOF
DESCRIPTION
Install script

-b|--bin_dir  Install missing programs here. This dir should writable and on
              your PATH. Default $bin_dir
-v|--version  Show version
-h|--help     Show help

Version $VERSION
EOF
exit 0
fi

if [[ $version == 1 ]]
then
    echo "$PG $VERSION"
    exit 0
fi
# End argument parsing
# ====================

set -x 

R --version

cwd=`pwd`
mkdir -p downloads

# GATK4
found=`command -v gatk` || true
if [[ -z $found ]]
then
    wget https://github.com/broadinstitute/gatk/releases/download/4.0.3.0/gatk-4.0.3.0.zip
    unzip gatk-4.0.3.0.zip
    rm gatk-4.0.3.0.zip
    ln -s `pwd`/gatk-4.0.3.0/gatk ${bin_dir}/
fi
command -v gatk
gatk Mutect2 --version

# Snakemake
found=`command -v snakemake` || true
if [[ -z $found ]]
then
    pip3 install --user 'snakemake==4.8.0'
fi
command -v snakemake
snakemake --version

# VEP
found=`command -v vep` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    export PERL_MM_USE_DEFAULT=1
    wget https://github.com/Ensembl/ensembl-vep/archive/release/92.1.tar.gz
    tar xf 92.1.tar.gz
    rm 92.1.tar.gz
    cd ensembl-vep-release-92.1
    perl INSTALL.pl --NO_HTSLIB --AUTO a 
    ln -s `pwd`/vep ${bin_dir}/
fi
command -v vep
vep --help

# Samtools
found=`command -v samtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
    tar xf samtools-1.8.tar.bz2
    rm samtools-1.8.tar.bz2
    cd samtools-1.8
    ./configure --prefix=`pwd`
    make
    make install
    cp samtools ${bin_dir}/
fi
command -v samtools
samtools --version

# bedtools
found=`command -v bedtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz 
    tar xf bedtools-2.27.1.tar.gz 
    rm bedtools-2.27.1.tar.gz
    cd bedtools2
    make
    cp bin/* ${bin_dir}/
fi
command -v bedtools
bedtools --version

exit 0

# R packages
cd ${cwd}
Rscript install/install_pkgs.R

