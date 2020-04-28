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
con <- gzfile(file.path( "philr_pipelines", "r_objects", "phyloseq_obj.rds"))
ps = readRDS(con)

# con <- gzfile(file.path( "philr_pipelines","phyloseq_obj.rds"))
# ps = readRDS(con)

##------------Import R objects and data preprocessing----------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects","philr_transform.rds"))
ps.philr = readRDS(con)

con <- gzfile(file.path( "philr_pipelines", "r_objects","ps_philr_transform.rds"))
ps = readRDS(con)


##---------------------PCOA transform and plots---------------------##
ps.dist <- dist(ps.philr, method="euclidean")


library(vegan)
an = anosim(ps.dist, grouping = as.factor(sample_data(ps)[,"Sample.Group"]), permutations = 999, distance = "bray", strata = NULL,
            parallel = getOption("mc.cores"))
print(an)
#  Unfortunately anosim does not have philr as a distance matrix
# https://www.rdocumentation.org/packages/vegan/versions/2.3-5/topics/vegdist
