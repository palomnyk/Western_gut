# Author: Aaron Yerke
# Comparison of philr and OTU table using glmod
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-Load Dependencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("phyloseq")
  BiocManager::install("philr")
  BiocManager::install("ape")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
}
library("RColorBrewer")
library("phyloseq")
library(philr); packageVersion("philr")
library(ggtree); packageVersion("ggtree")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion('dplyr')

##_Establish constants----------------------------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
setwd(file.path(home_dir))

##-Load data--------------------------------------------------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects","philr_transform.rds"))
ps.philr = readRDS(con)
close(con)

con <- gzfile(file.path( "philr_pipelines", "r_objects","ps_philr_transform.rds"))
ps = readRDS(con)
close(con)

otu_ratio = read.table(file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table.csv"),
                       sep = ",",
                       header = TRUE)

philr_trans = read.table(file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_philr_transform.csv"),
                         sep = ",",
                         header = TRUE)
philr_trans = philr_trans - min(philr_trans)

otu_tab = read.table(file.path(home_dir, "philr_pipelines", "tables", "otu_table.csv"),
                     sep = ",",
                     header = TRUE)

asv_table = read.table(file.path(home_dir, "philr_pipelines", "tables", "asv_table.csv"),
                       sep = ",",
                       header = TRUE)

metadata = read.table(file.path("philr_pipelines", "tables", "ps_sample_data_reduced.csv"), 
                      sep=",", 
                      header=TRUE)

my_datasets = list(asv_table, otu_tab, otu_ratio, philr_trans, cbind(otu_ratio, philr_trans), cbind(otu_tab, otu_ratio))
my_datasets = equal_num_columns_top_abund(my_datasets, 0.7)
my_ds_names = c("ASV", "OTU", "OTU ratio", "PhilR", "OTU + philR", "OTU + OTU ratio")

##-Load Dependencies------------------------------------------------##
##-Load Dependencies------------------------------------------------##
library(glmnet); packageVersion('glmnet')
glmmod <- glmnet::glmnet(as.matrix(philr_trans), metadata$Sample.Group, alpha=0.3, family="multinomial")

coef <- coefficients(glmmod, s=0.16)
top_coords <- data.frame(row.names = colnames(philr_trans))

for (n in names(coef)){
# according to https://slowkow.com/notes/sparse-matrix/ i contains the row number
  my_coef = as.data.frame(as.matrix(coef[n][[1]]))
  my_coef = my_coef[2:nrow(my_coef),]
  top_coords[n] = my_coef
  # tc <- rownames(top.coords)[which(top.coords != 0)]
}
# top.coords <- as.matrix(coefficients(glmmod, s=0.16)[m])
# top.coords <- rownames(top.coords)[which(top.coords != 0)]
# (top.coords <- top.coords[2:length(top.coords)]) # remove the intercept as a coordinate

top_coords <- top_coords[ rowSums(top_coords != 0) > 0, ]

tc_names <- sapply(row.names(top_coords), function(x) philr::name.balance(ps@phy_tree, ps@tax_table, x))

# votes <- name.balance(ps@phy_tree, ps@tax_table, 'n2', return.votes = c('up', 'down'))
# votes[[c('up.votes', 'Family')]] # Numerator at Family Level
# votes[[c('down.votes', 'Family')]] # Denominator at Family Level

tc.nn <- philr::name.to.nn(ps@phy_tree, row.names(top_coords))
tc.colors <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', "#fb9a33", "#12a02c")

my_colors <- c(rgb(0.99,0.1,0.1,0.2), 
               rgb(0.1,0.9,0.1,0.2), 
               rgb(0.1,0.1,0.9,0.2),
               rgb(0.3,0.6,0.1,0.2), 
               rgb(0.1,0.5,0.6,0.2), 
               rgb(0.8,0.1,0.2,0.2))


my_col <- c("red", "green", "blue", "yellow")

relevent_nodes <- c()

for (n in ps@phy_tree$node.label){
  if (n %in% row.names(top_coords)){
    relevent_nodes <- c(relevent_nodes, n)
  }else{
    relevent_nodes <- c(relevent_nodes, "")
  }
}

relevent_nodes <- c(ps@phy_tree$tip.label, relevent_nodes)

# set color palette
palette( RColorBrewer::brewer.pal(as.integer(ncol(top_coords)),"Accent") )

pdf(file = file.path(output_dir, "philr_tree_with_label.pdf"))
p <- ggtree::ggtree(ps@phy_tree, layout='rectangular') # geom_tiplab("genus") + geom_nodelab(geom = "label" = row.names(top_coords))
for (i1 in 1:nrow(top_coords)){#
  print(row.names(top_coords)[i1])
  for (i2 in 1:ncol(top_coords)){
    if (top_coords[i1,i2] != 0){
      print(palette()[i2])
      p <- p + ggtree::geom_balance(node=tc.nn[i1],
                                    color = "grey",
                                    fill = palette()[i2],
                                    alpha=0.1)
      # p <- philr::annotate_balance(ps@phy_tree,
      #                       row.names(top_coords)[i1],
      #                       p=p,
      #                       labels = c(paste0(row.names(top_coords)[i1],'+'),
      #                                  paste0(row.names(top_coords)[i1],'-')),
      #                       # offset.text=0.15,
      #                       bar=FALSE)
    }
  }
}
p <- p + ggtree::geom_text(aes(label = relevent_nodes), hjust=-.3)
p
dev.off()


ps.philr.long <- convert_to_long(ps.philr, get_variable(ps, 'gen2')) %>%
  filter(coord %in% top.coords)


ggplot(ps.philr.long, aes(x=labels, y=value)) +
  geom_boxplot(fill='lightgrey') +
  facet_grid(.~coord, scales='free_x') +
  xlab('Differentiate 2nd gen from from others') + ylab('Balance Value') +
  theme_bw()
