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
project = "RDP Western Gut"
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
setwd(file.path(home_dir))

##---------------------Import R objects-----------------------------##
con <- gzfile(file.path( "philr_pipelines","phyloseq_obj.rds"))
ps = readRDS(con)

##------------------------philr munging-----------------------------##
ps <-  filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
ps <-  filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
ps <- transform_sample_counts(ps, function(x) x+1)

#ps
##------------------------Test for reqs-----------------------------##
print("rooted tree?")
is.rooted(phy_tree(ps)) # Is the tree Rooted?
print('All multichotomies resolved?')
is.binary.tree(phy_tree(ps)) # All multichotomies resolved?

## ---- message=FALSE, warning=FALSE-----------------------------------------
phy_tree(ps) <- makeNodeLabel(phy_tree(ps), method="number", prefix='n')

name.balance(phy_tree(ps), tax_table(ps), 'n1')

otu.table <- t(otu_table(ps))
tree <- phy_tree(ps)
print("accessing sample data")
metadata <- sample_data(ps)
print("accessing taxanomic data")
tax <- tax_table(ps)

# otu.table[1:2,1:2] # OTU Table
# tree # Phylogenetic Tree
# head(tax,2) # taxonomy table

##---------------------philr transform------------------------------##
ps.philr <- philr(t(data.frame(otu.table)), tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')
#ps.philr[1:5,1:5]

##---------------------PCOA transform and plots---------------------##
ps.dist <- dist(ps.philr, method="euclidean")
ps.pcoa <- ordinate(ps, 'PCoA', distance=ps.dist)

catagoriesString = "Recruitment.Date	Sample.Date	Recruitment.Location	Researcher	Sub.Study	Sample.Month	Birth.Year	Age	Public.Housing	Medical.Assistance	Children.Free.Lunch	Highest.Education	Ethnicity	Religion	Birth.Location	Type.Birth.Location	Arrival.in.US	Years.in.US	Location.before.US	Type.location.before.US	Years.lived.in.Location.before.US	Tobacco.Use	Alcohol.Use	Height	Weight	Waist	BMI	BMI.Class	Medications	Breastfed	Years.Breastfed	Notes.Participants	Age.at.Arrival	Weight.M6	BMI.M6	Waist.M6	Measurement.Date.M6	Sample.Order	Sample.Group	Waist.Height.Ratio"
catagoriesString = unlist(strsplit(catagoriesString, split = "\t"))

setwd(file.path(output_dir, "philr_pcoa"))

library("ggplot2")

theme_set(theme_classic())

#Shapes
mimicShapes = c(17, 17, 1, 1, 16, 16)

#Colors
mimicCols = c("gray", "maroon2", "palevioletred1", "deeppink2", "seagreen4", "seagreen4")
myTitle = paste("philr_PCOA_", "ethn", ".pdf", sep = "")

plt = plot_ordination(ps, ps.pcoa, color="Sample.Group", shape="Sample.Group") + 
  scale_colour_manual(values=mimicCols) + 
  scale_shape_manual(values = c(17, 1, 17, 16, 1, 16)) +
  stat_ellipse(show.legend=F, type="t", level=.6) +
  scale_x_reverse() + 
  geom_point(size=3)
ggsave(plt, filename = myTitle)

#this code still works, just  commenting so it doesn't run every time I test this
# for (catag in catagoriesString){
#   print(catag)
#   catag = trimws(catag)
#   myTitle = paste("philr_PCOA_", catag, ".pdf", sep = "")
#   plt = plot_ordination(ps, ps.pcoa, color=catag)
#   ggsave(plt, filename = myTitle)
# }

print("done with PCOA")

##----------------------------Boxplots------------------------------##
source(file.path(home_dir, "philr_pipelines", "b.p.ratio.r"))

source(file.path(home_dir, "philr_pipelines", "b.p.ratio.r"))

tax_names = row.names(otu.table)

for(i in 1:length(tax)){
  # print(tax[i])
  q = tolower(tax[i][[6]])
  # print(strsplit(as.character(tax[i]), " ")[2])
  if (!is.na(q)){
    if (q == "bacteroides"){
      # b = strsplit(as.character(tax[i]), " ")[2]
      print(tax[i])
    }
    if ( q == "prevotella" ) {
      # p = strsplit(tax[i], " ")[1]
      print(tax[i])
    }
  }

}
print(b)
print(p)
b = "AGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGCAGACGGGACTTTAAGTCAGCTGTGAAATTTTCCGGCTCAACCGGGAAACTGCAGTTGATACTGGCGTCCTTGAGTACGGTCGAGGCAGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCCTGCCAGACCGCAACTGACGTTC"
p = "AGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGTTTGTTAAGCGTGTTGTGAAATGTCGGGGCTCAACCTGGGCATTGCAGCGCGAACTGGCAGACTTGAGTGCGCGGGAAGTAGGCGGAATTCGTCGTGTAGCGGTGAAATGCTTAGATATGACGAAGAACTCCGATTGCGAAGGCAGCCTGCTGTAGCGTAACTGACGCTGA"

bacteroides_tax <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"
prevotella_tax <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella"

#Note: map0 is metadata, otu0 is otutable
plot.b.p.ratio.L(t(metadata, otu.table, b, p, "B/P_plot.pdf", num.clip.months=NULL, show.stats)
