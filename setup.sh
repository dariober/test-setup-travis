#!/usr/bin/env bash

VERSION=0.1.0

set -e
set -o pipefail

# Parse arguments
# ===============

PG=`basename "$0"`
bin_dir=${HOME}/bin
xlog=/dev/null

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
    -l|--log)
        xlog="$2"
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
-l|--log      Write to this log file where commands have been installed and 
              their version. Default /dev/null 
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
>$xlog

# pindel (?)
# pysam <- Rewrite requesting package(s) in java?

function install_htslib(){
    # Download and install htslib. Compiled stuff is in `pwd`/htslib 
    pushd .
    rm -f htslib-1.8.tar.bz2
    wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
    tar xf htslib-1.8.tar.bz2
    rm htslib-1.8.tar.bz2
    mv htslib-1.8 htslib
    cd htslib
    ./configure --prefix=`pwd`
    make -j 4
    make install
    popd 
}

# HTSLIB (tabix & bgzip)
found=`command -v tabix` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    install_htslib
    cp htslib/tabix ${bin_dir}/
    cp htslib/bgzip ${bin_dir}/
fi
command -v tabix | tee -a $xlog
tabix --version | tee -a $xlog
command -v bgzip | tee -a $xlog
bgzip --version | tee -a $xlog

# FACETS::snp-pileup
found=`command -v snp-pileup` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -rf facets
    git clone https://github.com/mskcc/facets.git
    cd facets/inst/extcode/
    install_htslib
    g++ -std=c++11 -I`pwd`/htslib/include snp-pileup.cpp \
        -L`pwd`/htslib/lib -lhts -Wl,-rpath=`pwd`/htslib/lib -o snp-pileup
    ln -sf `pwd`/snp-pileup ${bin_dir}/
    #cp snp-pileup ${bin_dir}/
    # rm -r htslib
fi
command -v snp-pileup | tee -a $xlog
snp-pileup --help | tee -a $xlog

# Picard 
found=`find ${cwd}/bin/ -name picard.jar`
if [[ -z $found ]]
then
    cd ${cwd}/bin
    rm -f picard.jar
    wget https://github.com/broadinstitute/picard/releases/download/2.18.2/picard.jar
fi
java -jar ${cwd}/bin/picard.jar MarkDuplicates --version || true | tee -a $xlog

# IGVTools
found=`command -v igvtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f igvtools_2.3.98.zip
    wget http://data.broadinstitute.org/igv/projects/downloads/2.3/igvtools_2.3.98.zip
    unzip -q igvtools_2.3.98.zip
    rm igvtools_2.3.98.zip
    cp IGVTools/igvtools ${bin_dir}/
    cp IGVTools/igvtools.jar ${bin_dir}/
fi
command -v igvtools | tee -a $xlog
igvtools help | tee -a $xlog

# Manta
found=`command -v ${cwd}/bin/manta/bin/configManta.py` || true
if [[ -z $found ]]
then
    cd ${cwd}/bin
    rm -rf manta-1.3.2.centos6_x86_64.tar.bz2 manta-1.3.2.centos6_x86_64
    wget https://github.com/Illumina/manta/releases/download/v1.3.2/manta-1.3.2.centos6_x86_64.tar.bz2
    tar xf manta-1.3.2.centos6_x86_64.tar.bz2
    rm manta-1.3.2.centos6_x86_64.tar.bz2
    mv manta-1.3.2.centos6_x86_64 manta
fi
${cwd}/bin/manta/bin/configManta.py -h | tee -a $xlog

# BBDuk
found=`command -v bbduk.sh` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f BBMap_37.98.tar.gz
    wget --no-check-certificate https://sourceforge.net/projects/bbmap/files/BBMap_37.98.tar.gz
    tar xf BBMap_37.98.tar.gz
    rm BBMap_37.98.tar.gz
    ln -s `pwd`/bbmap/bbduk.sh ${bin_dir}/
fi
command -v bbduk.sh | tee -a $xlog
bbduk.sh --help | tee -a $xlog

# FastQC
found=`command -v fastqc` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f fastqc_v0.11.7.zip
    wget --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
    unzip -q -o fastqc_v0.11.7.zip
    rm fastqc_v0.11.7.zip
    ln -s `pwd`/FastQC/fastqc ${bin_dir}/
    chmod a+x `pwd`/FastQC/fastqc
fi
command -v fastqc | tee -a $xlog
fastqc --version | tee -a $xlog

# BWA
found=`command -v bwa` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f v0.7.17.tar.gz
    wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz
    tar xf v0.7.17.tar.gz
    rm v0.7.17.tar.gz
    cd bwa-0.7.17
    make
    cp bwa ${bin_dir}/
fi
command -v bwa | tee -a $xlog
bwa || true | tee -a $xlog # bwa doesn't have a --version or --help option

# GATK4
found=`command -v gatk` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f gatk-4.0.3.0.zip
    wget https://github.com/broadinstitute/gatk/releases/download/4.0.3.0/gatk-4.0.3.0.zip
    unzip -q gatk-4.0.3.0.zip
    rm gatk-4.0.3.0.zip
    ln -s `pwd`/gatk-4.0.3.0/gatk ${bin_dir}/
fi
command -v gatk | tee -a $xlog
gatk Mutect2 --version | tee -a $xlog

# Snakemake
found=`command -v snakemake` || true
if [[ -z $found ]]
then
    pip3 install --user 'snakemake==4.8.0'
fi
command -v snakemake | tee -a $xlog
snakemake --version | tee -a $xlog

# Python/Pandas
pip3 install --user pandas

# VEP
found=`command -v vep` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    export PERL_MM_USE_DEFAULT=1
    rm -f 92.1.tar.gz
    wget https://github.com/Ensembl/ensembl-vep/archive/release/92.1.tar.gz
    tar xf 92.1.tar.gz
    rm 92.1.tar.gz
    cd ensembl-vep-release-92.1
    perl INSTALL.pl --NO_HTSLIB --AUTO a 
    ln -s `pwd`/vep ${bin_dir}/
fi
command -v vep | tee -a $xlog
vep --help | tee -a $xlog

# Samtools
found=`command -v samtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f samtools-1.8.tar.bz2
    wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
    tar xf samtools-1.8.tar.bz2
    rm samtools-1.8.tar.bz2
    cd samtools-1.8
    ./configure --prefix=`pwd`
    make
    cp samtools ${bin_dir}/
fi
command -v samtools | tee -a $xlog
samtools --version | tee -a $xlog


# bcftools
found=`command -v bcftools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f bcftools-1.8.tar.bz2
    wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
    tar xf bcftools-1.8.tar.bz2
    rm bcftools-1.8.tar.bz2
    cd bcftools-1.8
    ./configure --prefix=`pwd`
    make
    cp bcftools ${bin_dir}/
fi
command -v bcftools | tee -a $xlog
bcftools --version | tee -a $xlog

# BEDTOOLS
found=`command -v bedtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f bedtools-2.27.1.tar.gz
    wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz 
    tar xf bedtools-2.27.1.tar.gz 
    rm bedtools-2.27.1.tar.gz
    cd bedtools2
    make -j 6
    cp bin/* ${bin_dir}/
fi
command -v bedtools | tee -a $xlog
bedtools --version | tee -a $xlog

# R packages
cd ${cwd}
Rscript install/install_pkgs.R


