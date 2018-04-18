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
    shift # past argument
    ;;
esac
done

if [[ $help == 1 ]]
then
cat <<EOF
DESCRIPTION
Installer of several required programs. Only programs not found on PATH will be
installed. See code for programs and versions.

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

cwd=`pwd`
mkdir -p downloads
mkdir -p bin

# snp-pileup (facets)
# pindel (?)
# pysam <- Rewrite requesting package(s) in java?

# Picard 
cd ${cwd}/bin
rm -f picard.jar
wget https://github.com/broadinstitute/picard/releases/download/2.18.2/picard.jar
java -jar picard.jar MarkDuplicates --version || true 

# IGVTools
found=`command -v igvtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    wget http://data.broadinstitute.org/igv/projects/downloads/2.3/igvtools_2.3.98.zip
    unzip -q igvtools_2.3.98.zip
    rm igvtools_2.3.98.zip
    cp IGVTools/igvtools ${bin_dir}/
    cp IGVTools/igvtools.jar ${bin_dir}/
fi
command -v igvtools
igvtools help

# Manta
cd ${cwd}/bin
rm -f manta-1.3.2.centos6_x86_64.tar.bz2 manta-1.3.2.centos6_x86_64
wget https://github.com/Illumina/manta/releases/download/v1.3.2/manta-1.3.2.centos6_x86_64.tar.bz2
tar xf manta-1.3.2.centos6_x86_64.tar.bz2
rm manta-1.3.2.centos6_x86_64.tar.bz2
mv manta-1.3.2.centos6_x86_64 manta
manta/bin/configManta.py -h

# BBDuk
found=`command -v igvtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    wget --no-check-certificate https://sourceforge.net/projects/bbmap/files/BBMap_37.98.tar.gz
    tar xf BBMap_37.98.tar.gz
    rm BBMap_37.98.tar.gz
    ln -s bbmap/bbduk.sh ${bin_dir}/
fi
command -v bbduk.sh
bbduk.sh --help

# FastQC
found=`command -v fastqc` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    wget --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
    unzip -q -o fastqc_v0.11.7.zip
    rm fastqc_v0.11.7.zip
    ln -s `pwd`/FastQC/fastqc ${bin_dir}/
    chmod a+x `pwd`/FastQC/fastqc
fi
command -v fastqc
fastqc --version

# BWA
found=`command -v bwa` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz
    tar xf v0.7.17.tar.gz
    rm v0.7.17.tar.gz
    cd bwa-0.7.17
    make
    cp bwa ${bin_dir}/
fi
command -v bwa
bwa || true # bwa doesn't have a --version or --help option

# GATK4
found=`command -v gatk` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    wget https://github.com/broadinstitute/gatk/releases/download/4.0.3.0/gatk-4.0.3.0.zip
    unzip -q gatk-4.0.3.0.zip
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

# Python/Pandas
pip3 install --user pandas

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
    cp samtools ${bin_dir}/
fi
command -v samtools
samtools --version

# HTSLIB (tabix & bgzip)
found=`command -v tabix` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
    tar xf htslib-1.8.tar.bz2
    rm htslib-1.8.tar.bz2
    cd htslib-1.8
    ./configure --prefix=`pwd`
    make
    cp tabix ${bin_dir}/
    cp bgzip ${bin_dir}/
fi
command -v tabix 
tabix --version
command -v bgzip 
bgzip --version

# bcftools
found=`command -v bcftools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
    tar xf bcftools-1.8.tar.bz2
    rm bcftools-1.8.tar.bz2
    cd bcftools-1.8
    ./configure --prefix=`pwd`
    make
    cp bcftools ${bin_dir}/
fi
command -v bcftools 
bcftools --version

# BEDTOOLS
found=`command -v bedtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz 
    tar xf bedtools-2.27.1.tar.gz 
    rm bedtools-2.27.1.tar.gz
    cd bedtools2
    make -j 6
    cp bin/* ${bin_dir}/
fi
command -v bedtools
bedtools --version

# R packages
cd ${cwd}
Rscript install/install_pkgs.R
