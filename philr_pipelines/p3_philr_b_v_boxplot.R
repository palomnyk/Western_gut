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

##------------------Import R objects and data-----------------------##
con <- gzfile(file.path( "philr_pipelines","phyloseq_obj.rds"))
ps = readRDS(con)

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

library("ggplot2")

theme_set(theme_classic())

#Shapes
mimicShapes = c(17, 21, 17, 16, 21, 16)

#Colors
mimicCols = c("gray", "maroon2", "palevioletred1", "deeppink2", "seagreen4", "seagreen4")

##----------------------------Boxplots------------------------------##
source(file.path(home_dir, "philr_pipelines", "import_code", "b.p.ratio.r"))

tax_names = row.names(otu.table)

myOtu = as.data.frame(otu_table(ps))
myTax = na.omit(as.data.frame(tax_table(ps)))

myOtu_names = row.names(myOtu)

myOtu = myOtu[row.names(myOtu) %in% row.names(myMeta),]
myMeta = myMeta[row.names(myMeta) %in% row.names(myOtu),]

b_df = myTax[myTax[,"Genus"] == "Bacteroides",]
p_df = myTax[myTax[,"Genus"] == "Prevotella",]

b_names = row.names(b_df)
p_names = row.names(p_df)

b_vs_p = c()
for ( i in 1:length(row.names(myOtu))){
  b_val = sum(myOtu[i, b_names])
  p_val = sum(myOtu[i, p_names])
  # print(paste(b_val, p_val))
  b_vs_p = c(b_vs_p, b_val/p_val)
}

# scale_colour_manual(values=mimicCols) + 
#   scale_shape_manual(values = minicShapes)

# ggplot(x = as.factor(myMeta[,"Sample.Group"]), y =  b_vs_p, shape=mimicCols, colour=minicShapes)) + 
#   geom_boxplot() + geom_jitter(width=.3, shape=mimicCols, colour=minicShapes) +
#   ylab(ylab) + xlab("") + theme(axis.text.x = element_text(size=7)) +
#   ggtitle("B/P ratio")
# save_plot(plot=p, fn, useDingbats=FALSE, base_aspect_ratio = 1.3 )

plot(as.factor(myMeta[,"Sample.Group"]), log(b_vs_p), 
     xlab = "",
     ylab = "B/P ratio", 
     main = "Westernized gut log(bacteriodes/prevolta ratio)")
stripchart(as.factor(myMeta[,"Sample.Group"]) ~ log(b_vs_p), 
           vertical = TRUE, 
           method = "jitter",
           add = TRUE, 
           pch = mimicShapes, 
           col = mimicCols)

new_otu_names = vector(length = length(colnames(myOtu)))

for (i in 1:length(colnames(myOtu))){
  old_name = colnames(myOtu)[i]
  new_otu_names[i] = as.character(myTax[old_name, "Genus"])
}

colnames(myOtu) = new_otu_names

setwd(file.path(output_dir, "philr_pcoa"))

#Note: map0 is metadata, otu0 is otutable
# plot.b.p.ratio.all(myMeta, myOtu, "Bacteroides", "Prevotella", "BP_plot.pdf")

plot.b.p.ratio.L(myMeta, myOtu, "Bacteroides", "Prevotella", "BP_Lplot.pdf")
