# Author: Aaron Yerke
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
library(ggtree); packageVersion("ggtree")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion('dplyr')

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

# sg = as.factor(ps@sam_data$Sample.Group)
# sample_data(GP)$human <- factor(get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue"))
sample_data(ps)$gen2 <- factor(get_variable(ps, "Sample.Group") %in% c("Karen2nd", "Hmong2nd"))


library(glmnet); packageVersion('glmnet')
glmmod <- glmnet(ps.philr, sample_data(ps)$gen2, alpha=0.3, family="binomial")
top.coords <- as.matrix(coefficients(glmmod, s=0.16))
top.coords <- rownames(top.coords)[which(top.coords != 0)]
(top.coords <- top.coords[2:length(top.coords)]) # remove the intercept as a coordinate

tc.names <- sapply(top.coords, function(x) name.balance(ps@phy_tree, ps@tax_table, x))

votes <- name.balance(ps@phy_tree, ps@tax_table, 'n2', return.votes = c('up', 'down'))
votes[[c('up.votes', 'Family')]] # Numerator at Family Level
votes[[c('down.votes', 'Family')]] # Denominator at Family Level


tc.nn <- name.to.nn(ps@phy_tree, top.coords)
tc.colors <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99')
p <- ggtree(ps@phy_tree, layout='rectangular') +
  geom_text(aes(label=node), hjust=-.3) +# geom_tiplab("genus") +
  geom_balance(node=tc.nn[1], fill=tc.colors[1], alpha=0.6) +
  geom_balance(node=tc.nn[2], fill=tc.colors[2], alpha=0.6) +
  geom_balance(node=tc.nn[3], fill=tc.colors[3], alpha=0.6) +
  geom_balance(node=tc.nn[4], fill=tc.colors[4], alpha=0.6) +
  geom_balance(node=tc.nn[5], fill=tc.colors[5], alpha=0.6)
p = annotate_balance(ps@phy_tree, top.coords[1], p=p, labels = c('n3+', 'n3-'),
                      offset.text=0.15, bar=TRUE)
p = annotate_balance(ps@phy_tree, 'n9', p=p, labels = c('n9+', 'n9-'),
                 offset.text=0.1, bar=T)
p


ps.philr.long <- convert_to_long(ps.philr, get_variable(ps, 'gen2')) %>%
  filter(coord %in% top.coords)


ggplot(ps.philr.long, aes(x=labels, y=value)) +
  geom_boxplot(fill='lightgrey') +
  facet_grid(.~coord, scales='free_x') +
  xlab('Differentiate 2nd gen from from others') + ylab('Balance Value') +
  theme_bw()
