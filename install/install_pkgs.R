#!/usr/bin/env Rscript

tryCatch({
        library(data.table)
    }, error = function(err) {
        install.packages('data.table', repos='http://cran.r-project.org')
    }
)

tryCatch({
        library(argparse)
    }, error = function(err) {
        install.packages('argparse', repos='http://cran.r-project.org')
    }
)

tryCatch({
        library(devtools)
    }, error = function(err) {
        install.packages('devtools', repos='http://cran.r-project.org')
    }
)

tryCatch({
        library(facets)
    }, error = function(err) {
        devtools::install_github("mskcc/facets")
    }
)
