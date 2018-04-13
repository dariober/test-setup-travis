#!/usr/bin/env Rscript

install.packages('data.table', repos='http://cran.r-project.org')
install.packages('argparse', repos='http://cran.r-project.org')
install.packages('devtools', repos='http://cran.r-project.org')
devtools::install_github("mskcc/facets")

#system('git clone https://github.com/mskcc/pctGCdata.git downloads/pctGCdata')
#install.packages("downloads/pctGCdata", repos= NULL, type= "source")
#
#system('git clone https://github.com/mskcc/facets.git downloads/facets')
#install.packages("downloads/facets", repos= NULL, type= "source", build_vignettes= TRUE)
