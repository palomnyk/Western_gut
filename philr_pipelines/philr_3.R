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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("philr")
BiocManager::install("ape")

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
source(file.path('~','git', "Western_gut", "philr_pipelines", "b.p.ratio.r"))

bacteroides <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"
prevotella <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella"

#Note: map0 is metadata, otu0 is otutable
plot.b.p.ratio.L(t(metadata), otu.table, "Bacteroides", "Prevotella", "B/P_plot.pdf", num.clip.months=NULL, show.stats)


## ---- message=FALSE, warning=FALSE-----------------------------------------
# sample_data(ps)$human <- factor(get_variable(ps, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue"))

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

