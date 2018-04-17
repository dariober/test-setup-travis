#!/bin/bash

snakemake -n --printshellcmds \
    -s workflows/main.smk \
    --configfile test_data/config.yaml \
    --jobs 2
