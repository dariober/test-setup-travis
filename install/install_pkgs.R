#!/usr/bin/env Rscript

# userdir <- unlist(strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep))[1L]

install.packages('data.table', repos='http://cran.r-project.org')
install.packages('argparse', repos='http://cran.r-project.org')
system('git clone https://github.com/mskcc/facets.git downloads/facets')
install.packages("downloads/facets", repos= NULL, type= "source", build_vignettes= TRUE)
