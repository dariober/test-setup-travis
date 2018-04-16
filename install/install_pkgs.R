#!/usr/bin/env Rscript

install.packages('digest', repos='http://cran.r-project.org')
install.packages('data.table', repos='http://cran.r-project.org')
install.packages('argparse', repos='http://cran.r-project.org')
install.packages('devtools', repos='http://cran.r-project.org')
devtools::install_github("mskcc/facets")

