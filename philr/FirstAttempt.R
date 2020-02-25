#Pipeline for processing BioLockJ output of Western gut data through philr
#Useful: https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.R
rm(list = ls()) #clear workspace

#Install philr
## try http:// if https:// URLs are not supported
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
#BiocManager::install("philr")
library("philr")

setwd(file.path('~','git','Western_gut'))
project = "RDP Western Gut"
#f_path = c('.', 'biolockj', 'rdp_2020Feb19', '04_BuildTaxaTables','output','rdp_2020Feb19_taxaCount_genus.tsv')
f_path = './biolockj/rdp_2020Feb19/04_BuildTaxaTables/output/rdp_2020Feb19_taxaCount_genus.tsv'

#load data
m = read.table(file.path(f_path), 
               sep="\t", 
               header=TRUE, 
               row.names = 1, 
               check.names = FALSE,
               stringsAsFactors=FALSE)

m = data.matrix(m + 1)

m_transform = philr(m, tree, sbp = NULL, part.weights = "uniform", ilr.weights = "uniform", return.all = FALSE)

