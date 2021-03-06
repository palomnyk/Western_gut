# Author: Aaron Yerke
# This is a script to restart the pipeline that was created to explore the philr transform.
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

#try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
# BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("dada2")
BiocManager::install("phyloseq")
BiocManager::install("DECIPHER")
BiocManager::install("phangorn")
BiocManager::install("ggplot2")
BiocManager::install("gridExtra")
library("dada2")
library("DECIPHER")
library("phangorn")

home_dir = file.path('~','git','Western_gut')
output_dir = file.path('~','git','Western_gut', 'output')
project = "RDP Western Gut"
#f_path = c('.', 'biolockj', 'rdp_2020Feb19', '04_BuildTaxaTables','output','rdp_2020Feb19_taxaCount_genus.tsv')
f_path <- "~/git/Western_gut/sequences" # CHANGE ME to the directory containing the fastq files after unzipping.
#list.files(f_path)

setwd(file.path("~", "git","Western_gut", "philr_pipelines"))

con <- gzfile("ForwardReads_DADA2.rds")
seqtab = readRDS(con)

fastaRef <- file.path(home_dir, 'philr_pipelines', "taxonomy", "./rdp_train_set_16.fa.gz")
taxTab <- assignTaxonomy(seqtab, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))

seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
print("Alignment completed")

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
detach("package:dada2", unload=TRUE)
detach("package:DECIPHER", unload=TRUE)

print("phangorn completed")

print("Created ps")
library("phyloseq")

myMeta = read.table(file.path(home_dir,'metadata.txt'), 
                    sep="\t", 
                    header=TRUE, 
                    row.names = 1, 
                    check.names = FALSE,
                    stringsAsFactors=FALSE)

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sameple_data(myMeta),
               tax_table(taxTab),
               phy_tree(fitGTR$tree))
# ps

setwd(file.path("~", "git","Western_gut", "philr_pipelines"))
saveRDS(seqtab, "ForwardReads_DADA2.rds")

# Install philr
library(philr); packageVersion("philr")

print("philr installed")

ps <-  filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
ps <-  filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
ps <- transform_sample_counts(ps, function(x) x+1)

ps

G_tree <- phy_tree(ps)

## ---- message=FALSE, warning=FALSE-----------------------------------------
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

otu.table[1:2,1:2] # OTU Table
tree # Phylogenetic Tree
vershead(tax,2) # taxonomy table

## --------------------------------------------------------------------------
ps.philr <- philr(otu.table, tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')
ps.philr[1:5,1:5]

## --------------------------------------------------------------------------
ps.dist <- dist(ps.philr, method="euclidean")
ps.pcoa <- ordinate(ps, 'PCoA', distance=ps.dist)
plot_ordination(ps, ps.pcoa, color='SampleType') + geom_point(size=4)

## ---- message=FALSE, warning=FALSE-----------------------------------------
sample_data(ps)$human <- factor(get_variable(ps, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue"))

## ---- message=FALSE, warning=FALSE-----------------------------------------
library(glmnet); packageVersion('glmnet')
glmmod <- glmnet(ps.philr, sample_data(ps)$human, alpha=1, family="binomial")

## --------------------------------------------------------------------------
top.coords <- as.matrix(coefficients(glmmod, s=0.2526))
top.coords <- rownames(top.coords)[which(top.coords != 0)]
(top.coords <- top.coords[2:length(top.coords)]) # remove the intercept as a coordinate

## --------------------------------------------------------------------------
tc.names <- sapply(top.coords, function(x) name.balance(tree, tax, x))
tc.names

## --------------------------------------------------------------------------
votes <- name.balance(tree, tax, 'n730', return.votes = c('up', 'down'))

votes[[c('up.votes', 'Family')]] # Numerator at Family Level
votes[[c('down.votes', 'Family')]] # Denominator at Family Level

## ---- message=FALSE, warning=FALSE-----------------------------------------
library(ggtree); packageVersion("ggtree")
library(dplyr); packageVersion('dplyr')

## ---- message=FALSE, warning=FALSE-----------------------------------------
tc.nn <- name.to.nn(tree, top.coords)
tc.colors <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99')
p <- ggtree(tree, layout='fan') +
  geom_balance(node=tc.nn[1], fill=tc.colors[1], alpha=0.6) +
  geom_balance(node=tc.nn[2], fill=tc.colors[2], alpha=0.6) +
  geom_balance(node=tc.nn[3], fill=tc.colors[3], alpha=0.6) +
  geom_balance(node=tc.nn[4], fill=tc.colors[4], alpha=0.6) +
  geom_balance(node=tc.nn[5], fill=tc.colors[5], alpha=0.6)
p <- annotate_balance(tree, 'n16', p=p, labels = c('n16+', 'n16-'),
                      offset.text=0.15, bar=FALSE)
annotate_balance(tree, 'n730', p=p, labels = c('n730+', 'n730-'),
                 offset.text=0.15, bar=FALSE)

## --------------------------------------------------------------------------
ps.philr.long <- convert_to_long(ps.philr, get_variable(ps, 'human')) %>%
  filter(coord %in% top.coords)

library("ggplot2")
library("gridExtra")

ggplot(ps.philr.long, aes(x=labels, y=value)) +
  geom_boxplot(fill='lightgrey') +
  facet_grid(.~coord, scales='free_x') +
  xlab('Human') + ylab('Balance Value') +
  theme_bw()

## ---- message=FALSE, warning=FALSE-----------------------------------------
library(tidyr); packageVersion('tidyr')

ps.philr.long %>%
  rename(Human=labels) %>%
  filter(coord %in% c('n16', 'n730')) %>%
  spread(coord, value) %>%
  ggplot(aes(x=n16, y=n730, color=Human)) +
  geom_point(size=4) +
  xlab(tc.names['n16']) + ylab(tc.names['n730']) +
  theme_bw()