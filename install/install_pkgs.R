#!/usr/bin/env Rscript

# userdir <- unlist(strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep))[1L]

install.packages('data.table', repos='http://cran.r-project.org')
install.packages('argparse', repos='http://cran.r-project.org')
install.packages('devtools', repos='http://cran.r-project.org')
devtools::install_github("mskcc/facets", build_vignettes = TRUE)
