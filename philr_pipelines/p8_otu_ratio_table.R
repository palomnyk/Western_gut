# Author: Aaron Yerke
# Script for making ratio table of OTUS for ml purposes
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("phyloseq")
  BiocManager::install("philr")
  BiocManager::install("ape")
}
library("phyloseq")
library(philr); packageVersion("philr")
library("reshape2")
##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
setwd(file.path(home_dir))

##------------Import R objects and data preprocessing---------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects","otu_table.rds"))
otu_tab = readRDS(con)
close(con)

##--------------------Create Ratio table----------------------------##
otuRatios = data.frame(row.names = row.names(otu_tab))

for ( rn in 1:nrow(otu_tab)){
  ratios = vector(length = length(uniq_otus)^2 - length(uniq_otus))
  index = 1
  for (cn1 in 1:ncol(otu_tab)){
    for (cn2 in 1:ncol(otu_tab)){
      if (cn2 != cn1){
        ratios[index] = otu_tab[rn, cn1] / otu_tab[rn, cn2]
        names(ratios)[index] = paste0(colnames(otu_tab)[cn1], "/", colnames(otu_tab)[cn1])
        index = index + 1
      }# end if cn2 != cn1
    }#end cn2
  }#end cn1
  ratioDf = t(data.frame(ratios))
  colnames(ratioDf) = names(ratios)
  row.names(ratioDf) = row.names(otu_tab)[rn]
  
  otuRatios = rbind(otuRatios, ratioDf)
  # otuRatios = rbind( otuRatios, ratios)
}

write.table(otuRatios, 
            file = file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table.csv"),
            sep = ",")

saveRDS(otu_tab, file.path("philr_pipelines", "r_objects", "otu_ratio_table.rds"))
