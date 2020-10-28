# Author: Aaron Yerke
# Script for making ratio table of OTUS for ml purposes
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##
##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
setwd(file.path(home_dir))

##------------Import R objects and data preprocessing---------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects","otu_table.rds"))
otu_tab <- readRDS(con)
close(con)

otu_t <- otu_tab[1:5,1:4]

##--------------------Create Ratio table----------------------------##

simple_ratio_transform <- function(otu_df) {
  # Function for making ratio table of OTUS for ml purposes
  # The resulting table will look like: otu1/otu2, otu2/otu1
  otuRatios = data.frame(row.names = row.names(otu_df))
  for ( rn in 1:nrow(otu_df)){
    ratios = vector(length = ncol(otu_df)^2 - ncol(otu_df))
    index = 1
    for (cn1 in 1:ncol(otu_df)){
      for (cn2 in 1:ncol(otu_df)){
        if (cn2 != cn1){
          ratios[index] = otu_df[rn, cn1] / otu_df[rn, cn2]
          names(ratios)[index] = paste0(colnames(otu_df)[cn1], "/", colnames(otu_df)[cn2])
          index = index + 1
        }# end if cn2 != cn1
      }#end cn2
    }#end cn1
    ratioDf = t(data.frame(ratios))
    colnames(ratioDf) = names(ratios)
    row.names(ratioDf) = row.names(otu_df)[rn]
    
    otuRatios = rbind(otuRatios, ratioDf)
    # otuRatios = rbind( otuRatios, ratios)
  }
  return(otuRatios)
}

sum_ratio_transform <- function(otu_df) {
  # Function for making ratio table of OTUS for ml purposes
  # The resulting table will look like: otu1/otu2, otu2/otu1
  otuRatios = data.frame(row.names = row.names(otu_df))
  for ( rn in 1:nrow(otu_df)){
    ratios = vector(length = ncol(otu_df)^2 - ncol(otu_df))
    index = 1
    for (cn1 in 1:ncol(otu_df)){
      for (cn2 in 1:ncol(otu_df)){
        if (cn2 != cn1){
          ratios[index] = otu_df[rn, cn1] / otu_df[rn, cn2]
          names(ratios)[index] = paste0(colnames(otu_df)[cn1], "/", colnames(otu_df)[cn2])
          index = index + 1
        }# end if cn2 != cn1
      }#end cn2
    }#end cn1
    ratioDf = t(data.frame(ratios))
    colnames(ratioDf) = names(ratios)
    row.names(ratioDf) = row.names(otu_df)[rn]
    
    otuRatios = rbind(otuRatios, ratioDf)
    # otuRatios = rbind( otuRatios, ratios)
  }
  return(otuRatios)
}

my_ratios <- simple_ratio_transform(otu_t)

# write.table(otuRatios,
#             file = file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table.csv"),
#             sep = ",")
# 
# saveRDS(my_ratios, file.path("philr_pipelines", "r_objects", "otu_ratio_table.rds"))
