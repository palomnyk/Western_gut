# Author: Aaron Yerke
# Modified p6.4 to answer how OTUs compare to nodes with that OTU's name as highest voted branch
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
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
setwd(file.path(home_dir))
##-----------------------Functions----------------------------------##
onlyOne = function(namedVector){
  if(length(namedVector) == 1){
    return(names(namedVector)[1])
  }else {
    return("---")
  }
}

##------------Import R objects and data preprocessing---------------##
# con <- gzfile(file.path( "philr_pipelines", "r_objects","philr_transform.rds"))
# ps.philr = data.frame(readRDS(con))
# close(con)
# 
# con <- gzfile(file.path( "philr_pipelines", "r_objects","ps_philr_transform.rds"))
# ps = readRDS(con)
# close(con)

con <- gzfile(file.path( "philr_pipelines", "r_objects","ref_tree_philr_transform.rds"))
ps.philr = data.frame(readRDS(con))
close(con)

con <- gzfile(file.path( "philr_pipelines", "r_objects","ref_tree_ps_philr_transform.rds"))
ps = readRDS(con)
close(con)

# ps <- transform_sample_counts(ps, function(x) x+1)

philr_pb_ratio = ps.philr[,"n13"]

no = data.frame(ps@otu_table)
tax = data.frame(ps@tax_table)
otus = vector(length = ncol(no), mode = "character")
for (i in 1:ncol(no)){
  asv = names(no)[i]
  otus[i] = tolower(as.character(tax[asv, 5]))#this line sets the taxonomic level
}

uniq_otus = unique(otus)

# set up main data.frame (read.table allows for prenaming columns)
nodes_vs_otu = read.table(text = "",
                          colClasses = rep(x = "numeric", length(ps@phy_tree$node.label)),
                          col.names = ps@phy_tree$node.label)
# nodes_vs_otu = data.frame()

perfect_sep = rep(T, ncol(ps.philr))# vector(length = ncol(ps.philr))#tells whether node is perfect sep or not
node_vote_up = list()
node_vote_down = list()
for( m in 1:ncol(ps.philr)){
  votes <- name.balance(ps@phy_tree, ps@tax_table, names(ps.philr)[m], return.votes = c('up', 'down'))
  node_vote_up[m] = votes$up.votes$Genus
  node_vote_down[m] = votes$down.votes$Genus
  genus = paste0(onlyOne(votes$up.votes$Genus), "/", onlyOne(votes$down.votes$Genus))
  if( grepl( "---", genus, fixed = T) ){
    perfect_sep[[m]] = FALSE
  }
}

for(otu1 in 1:length(uniq_otus)){
  my_otu = uniq_otus[otu1]
  my_otus = otus == my_otu
  my_otus = replace(my_otus, is.na(my_otus), F)
  #need to add vectors so that they are the same length as nrow(ps.philr)
  my_otus = as.data.frame(no[,my_otus])
  my_otus = rowSums(my_otus, na.rm=T)
  
  my_corrs = vector(length = ncol(ps.philr))
  
  for (node1 in 1:ncol(ps.philr)){
    my_cor = cor(my_otus, unlist(ps.philr[,node1]), method = "kendall")
    my_corrs[node1] = my_cor
  }#end otu2 for loop
  nodes_vs_otu[otu1,] = my_corrs
}#end otu1 for loop

perf_df = nodes_vs_otu[,perfect_sep]
nperf_df = nodes_vs_otu[,!perfect_sep]

plot(perf_df, nperf_df,
     main = "OTU vs philr node kendall correlations")
# xlab = "perfect seperation nodes",
# ylab = "mixed nodes")

boxplot(unlist(as.vector(perf_df)) , unlist(as.vector(nperf_df)),
        main = "Kendall for perfect seperation and non perfect seperation vs OTU\n t-test p-value = 0.00227", 
        names = c("perfect sep", "blurry sep"), 
)
stripchart(unlist(as.vector(perf_df)), unlist(as.vector(nperf_df)), 
           method = "jitter",
           vertical=T, 
           add = T)
t.test(unlist(as.vector(perf_df)), unlist(as.vector(nperf_df)))




mean(as.vector(perf_df))

mean(nperf_df)
