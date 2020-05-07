# Author: Aaron Yerke
# This script is for exploring the B/V ratio and determining if it is the best ratio

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

##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
setwd(file.path(home_dir))

##------------Import R objects and data preprocessing----------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects","philr_transform.rds"))
ps.philr = readRDS(con)

con <- gzfile(file.path( "philr_pipelines", "r_objects","ps_philr_transform.rds"))
ps = readRDS(con)

tax = data.frame(tax_table(ps))

##----------------------calculate ratios----------------------------##
node_ratio_table = data.frame(row.names = )

sam_name = c()
ratio_nodes = c()

#iterate through each node and compare it to each other node
# for (i1 in 1:nrow(ps.philr)){
#   for (i2 in 1:nrow(ps.philr)){
#     print(paste(i1,i2))
#   }
#   
# }
plot_tree(ps, "sampledodge", nodeplotdefault, ladderize="left")

plot_tree(ps, "sampledodge", 
          ladderize="left",
          label.tips = "Genus",
          sizebase = 200,
          base.spacing = 0.00009,
          text.size = 2.75)

library('ape')

tl = data.frame(ps@phy_tree$tip.label)
e = data.frame(ps@phy_tree$edge
               )

plot_bar(ps, "Family", "Abundance", "Family", 
         facet_grid="Sample.Group~."
         )

n2<-extract.clade(ps@phy_tree,73)

plot_tree(n2, "sampledodge", 
          ladderize="left",
          # label.tips = "Genus", 
          sizebase = 200,
          base.spacing = 0.00009,
          text.size = 2.75)
