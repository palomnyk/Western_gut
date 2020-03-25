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

# ‘ape’, ‘dplyr’, ‘reshape2’, ‘plyr’
# .cran_packages <- c("ggplot2", "gridExtra")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# if (!requireNamespace("phangorn", quietly = TRUE))
#   install.packages("phangorn")  
# BiocManager::install("phyloseq")
# BiocManager::install("DECIPHER")

library("DECIPHER")
library("phangorn")
library("phyloseq")
#library("ape")

#home_dir = file.path('~','git','Western_gut')
home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
project = "RDP Western Gut"
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
# list.files(f_path)

setwd(file.path(home_dir))

con <- gzfile(file.path( "philr_pipelines","ForwardReads_DADA2.rds"))
seqtab = readRDS(con)

con <- gzfile(file.path( "philr_pipelines","ForwardReads_DADA2_alignment.rds"))
alignment <- readRDS(con)

con <- gzfile(file.path( "philr_pipelines","ForwardReads_DADA2_taxonomy.rds"))
taxTab <- readRDS(con)

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
print("phangorn completed")


myMeta = read.table(file.path('PRJEB28687.txt'), 
                    sep="\t", 
                    header=TRUE, 
                    row.names = "run_accession", 
                    check.names = FALSE,
                    stringsAsFactors=FALSE)

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(myMeta),
               tax_table(taxTab),
               phy_tree(fitGTR$tree))
# ps
print("Created ps")

# setwd(file.path(home_dir, "philr_pipelines"))
saveRDS(ps, file.path("philr_pipelines", "phyloseq_obj.rds"))