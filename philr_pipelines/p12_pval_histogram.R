# Author: Aaron Yerke
# Script for PERMANOVA R^2
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-Load Dependencies------------------------------------------------##
if (!requireNamespace("vegan", quietly = TRUE)){
  install.packages("vegan")
}
library(vegan)
##_Establish constants----------------------------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
setwd(file.path(home_dir))
##---------------------------Functions------------------------------##
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

##-Load data--------------------------------------------------------##
otu_ratio = read.table(file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table.csv"),
                       sep = ",",
                       header = TRUE)

philr_denovo = read.table(file.path(home_dir, "philr_pipelines", "tables", "denovo_ps_philr_transform.csv"),
                         sep = ",",
                         header = TRUE)

philr_ref = read.table(file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_philr_transform.csv"),
                         sep = ",",
                         header = TRUE)

metadata = read.table(file.path("philr_pipelines", "tables", "ps_sample_data_reduced.csv"), 
                      sep=",", 
                      header=TRUE)

my_datasets = list(otu_ratio, philr_denovo, philr_ref)
my_datasets = equal_num_columns_top_abund(my_datasets, 0.75)
my_ds_names = c("OTU ratio", "philr denovo tree", "philr ref tree")
##-Create histograms------------------------------------------------##

#set color palette
palette( c(rgb(0.99,0.1,0.1,0.3), rgb(0.1,0.9,0.1,0.3), rgb(0.1,0.1,0.9,0.3)) )

# pdf(file = file.path(output_dir, "pval_histogram.pdf"))
for(mta in 1:ncol(metadata)){
  my_meta = metadata[,mta]
  first_hist = FALSE #flag for adding to the original hist
  for( ds in 1:length(my_datasets)){
    my_table = my_datasets[ds][[1]]
    my_pvals = c()
    for( colmn in 1:ncol(my_table)){
      my_lm = lm(unlist(my_table[,colmn]) ~ my_meta)
      my_pval = anova(my_lm)$"Pr(>F)"[1]
      my_pvals = c(my_pvals, my_pval)
    }#for( colmn in 1:ncol(my_table)){
    hist(my_pvals, 
         add = first_hist, 
         col = ds,
         breaks = 100,
         main = names(metadata)[mta])
    print(ds)
    first_hist = TRUE #update flag
  }#for( ds in 1:length(my_datasets)){
  legend('topright', 
         legend = my_ds_names, 
         col = my_ds_names,
         pch = 15,
         cex=0.5)
}

# dev.off()


