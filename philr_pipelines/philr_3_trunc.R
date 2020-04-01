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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("philr")
BiocManager::install("ape")

library("phyloseq")
library("ape")
library(philr); packageVersion("philr")

home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
project = "RDP Western Gut"
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
# list.files(f_path)

setwd(file.path(home_dir))

con <- gzfile(file.path( "philr_pipelines","phyloseq_obj.rds"))
ps = readRDS(con)

ps <-  filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
ps <-  filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
ps <- transform_sample_counts(ps, function(x) x+1)

#ps

print("rooted tree?")
is.rooted(phy_tree(ps)) # Is the tree Rooted?
print('All multichotomies resolved?')
is.binary.tree(phy_tree(ps)) # All multichotomies resolved?

## ---- message=FALSE, warning=FALSE-----------------------------------------
phy_tree(ps) <- makeNodeLabel(phy_tree(ps), method="number", prefix='n')

## --------------------------------------------------------------------------
name.balance(phy_tree(ps), tax_table(ps), 'n1')

## --------------------------------------------------------------------------
otu.table <- t(otu_table(ps))
tree <- phy_tree(ps)
print("accessing sample data")
metadata <- sample_data(ps)
print("accessing taxanomic data")
tax <- tax_table(ps)

# otu.table[1:2,1:2] # OTU Table
# tree # Phylogenetic Tree
# head(tax,2) # taxonomy table

## --------------------------------------------------------------------------
ps.philr <- philr(t(data.frame(otu.table)), tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')
ps.philr[1:5,1:5]

## --------------------------------------------------------------------------
ps.dist <- dist(ps.philr, method="euclidean")
ps.pcoa <- ordinate(ps, 'PCoA', distance=ps.dist)

catagoriesString = "Recruitment.Date	Sample.Date	Recruitment.Location	Researcher	Sub.Study	Sample.Month	Birth.Year	Age	Public.Housing	Medical.Assistance	Children.Free.Lunch	Highest.Education	Ethnicity	Religion	Birth.Location	Type.Birth.Location	Arrival.in.US	Years.in.US	Location.before.US	Type.location.before.US	Years.lived.in.Location.before.US	Tobacco.Use	Alcohol.Use	Height	Weight	Waist	BMI	BMI.Class	Medications	Breastfed	Years.Breastfed	Notes.Participants	Age.at.Arrival	Weight.M6	BMI.M6	Waist.M6	Measurement.Date.M6	Sample.Order	Sample.Group	Waist.Height.Ratio"
catagoriesString = unlist(strsplit(catagoriesString, split = "\t"))

setwd(file.path(output_dir, "philr_pcoa"))
for (catag in catagoriesString){
  print(catag)
  pdf(paste0("philr_PCOA_", catag, ".pdf"))
  plot_ordination(ps, ps.pcoa, color=trimws(catag))
  dev.off()
}
