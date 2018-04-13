#!/usr/bin/env bash

VERSION=0.1.0

set -e
set -o pipefail

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

python3 --version
which python3
which python3.5

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

if [[ $user == 1 ]]
then
    pip3 install --user snakemake
else
    pip3 install snakemake
fi
cd ${cwd}
