# Author: Aaron Yerke
# This is a pipeline that was created to explore the philr transform.
# Many resources were used to create this pipeline:
#   Processing the sequences through DADA2, making trees with phangorn, and putting them into phyloseq objects:
#     This is by far the biggest source:
#     https://github.com/spholmes/F1000_workflow/blob/master/MicrobiomeWorkflow/MicrobiomeWorkflowII.Rmd
#     or (in article format)
#     https://f1000research.com/articles/5-1492/v2
#     
#   DADA2:
#     https://benjjneb.github.io/dada2/tutorial.html
#   
#   Using philr:
#     https://bioconductor.org/packages/release/bioc/html/philr.html
#     Code example:
#       https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.R

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("phyloseq")
  BiocManager::install("philr")
  BiocManager::install("ape")
}

library("phyloseq")
library("ape")
library(philr); packageVersion("philr")

##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
setwd(file.path(home_dir))

# ##------------------Import R objects and data-----------------------##
# con <- gzfile(file.path( "philr_pipelines", "r_objects", "phyloseq_obj.rds"))
# ps = readRDS(con)
# close(con)

##------------------Import R objects and data-----------------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects", "ref_tree_phyloseq_obj.rds"))
ps = readRDS(con)
close(con)

myMeta = read.table(file.path("fullMetadata.tsv"), 
                    sep="\t", 
                    header=TRUE, 
                    row.names = "run_accession", 
                    check.names = FALSE,
                    stringsAsFactors=FALSE)

##------------------------philr munging-----------------------------##
ps <-  filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
ps <-  filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
ps <- transform_sample_counts(ps, function(x) x+1)

##------------------------Test for reqs-----------------------------##
print("rooted tree?")
is.rooted(phy_tree(ps)) # Is the tree Rooted?
print('All multichotomies resolved?')
is.binary.tree(phy_tree(ps)) # All multichotomies resolved?

## ---- message=FALSE, warning=FALSE-----------------------------------------
phy_tree(ps) <- makeNodeLabel(phy_tree(ps), method="number", prefix='n')

name.balance(phy_tree(ps), tax_table(ps), 'n1')

##---------------------philr transform------------------------------##
ps.philr <- philr(ps@otu_table, ps@phy_tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')

# commented out to activate ref_tree_philr
saveRDS(ps, file.path("philr_pipelines", "r_objects", "denovo_ps_philr_transform.rds"))
saveRDS(ps.philr, file.path("philr_pipelines", "r_objects", "denovo_philr_transform.rds"))
write.table(ps.philr, 
            file.path("philr_pipelines", "tables", "denovo_ps_philr_transform.csv"),
            sep = ",")

# commented out to activate denovo tree
# saveRDS(ps, 
#         file.path("philr_pipelines", "r_objects", "ref_tree_ps_philr_transform.rds"))
# 
# write.table(ps.philr, 
#             file.path("philr_pipelines", "tables", "ref_tree_ps_philr_transform.csv"),
#             sep = ",")
# 
# write.table(metadata, 
#             file.path("philr_pipelines", "tables", "ps_sample_data.csv"),
#             sep = ",")
# 
# write.table(ps@otu_table, 
#             file.path("philr_pipelines", "tables", "asv_table.csv"),
#             sep = ",")