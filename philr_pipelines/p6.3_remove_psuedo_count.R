# Author: Aaron Yerke
# Modified p4.1 to answer this question: Find nodes with perfect separation of any taxa at any level
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
library("reshape2")
##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
setwd(file.path(home_dir))

##------------Import R objects and data preprocessing---------------##
# con <- gzfile(file.path( "philr_pipelines", "r_objects","philr_transform.rds"))
# ps.philr = data.frame(readRDS(con))
# close(con)
# 
# con <- gzfile(file.path( "philr_pipelines", "r_objects","ps_philr_transform.rds"))
# ps = readRDS(con)
# close(con)

con <- gzfile(file.path( "philr_pipelines", "r_objects","ref_tree_philr_transform.rds"))
ps.philr = data.frame(readRDS(con))
close(con)

con <- gzfile(file.path( "philr_pipelines", "r_objects","ref_tree_ps_philr_transform.rds"))
ps = readRDS(con)
close(con)

metadata = data.frame(ps@sam_data)
metadata = metadata[, !apply(is.na(metadata), 2, all)]#remove columns where all are NA
metadata[c(1:3)] <- list(NULL)
drop = c("sample_accession","secondary_sample_accession","experiment_accession",             
         "tax_id","scientific_name","instrument_model","library_layout",
         "experiment_title","sample_title", "study_title", "run_alias",
         "fastq_ftp" ,"fastq_galaxy","submitted_ftp"                    
         ,"submitted_galaxy","sra_ftp", "sra_galaxy", "study_accession")
metadata = metadata[ , !(names(metadata) %in% drop)]
##----------------------t-test w metadata----------------------------##

philr_pb_ratio = ps.philr[,"n13"]

no = data.frame(ps@otu_table)
tax = data.frame(ps@tax_table)
otus = vector(length = ncol(no), mode = "character")
for (i in 1:ncol(no)){
  asv = names(no)[i]
  otus[i] = tolower(as.character(tax[asv, 5]))
}

##----------------------calculate ratios----------------------------##
# bact = otus == "bacteroides"
# bact = replace(bact, is.na(bact), F)
# prev = otus == "prevotella"
# prev = replace(prev, is.na(prev), F)

bact = otus == "bacteroidaceae"
bact = replace(bact, is.na(bact), F)
prev = otus == "prevotellaceae"
prev = replace(prev, is.na(prev), F)

class_pb_ratio= c(length = nrow(no))

for (i in 1:nrow(no)){
  class_pb_ratio[i] = sum(no[i, bact])/ sum(no[i, prev])
}

print(paste("sanity check length:", length(class_pb_ratio) == length(philr_pb_ratio)))
print(paste("sanity check order: ", all.equal(row.names(no), row.names(ps.philr) )))
metadata = data.frame(ps@sam_data)
metadata = metadata[, !apply(is.na(metadata), 2, all)]#remove columns where all are NA
metadata[c(1:3)] <- list(NULL)
drop = c("sample_accession","secondary_sample_accession","experiment_accession",             
         "tax_id","scientific_name","instrument_model","library_layout",
         "experiment_title","sample_title", "study_title", "run_alias",
         "fastq_ftp" ,"fastq_galaxy","submitted_ftp"                    
         ,"submitted_galaxy","sra_ftp", "sra_galaxy", "study_accession")
metadata = metadata[ , !(names(metadata) %in% drop)]
##----------------------t-test w metadata----------------------------##

philr_pb_ratio = ps.philr[,"n13"]

no = data.frame(ps@otu_table)
tax = data.frame(ps@tax_table)
otus = vector(length = ncol(no), mode = "character")
for (i in 1:ncol(no)){
  asv = names(no)[i]
  otus[i] = tolower(as.character(tax[asv, 5]))
}

##----------------------calculate ratios----------------------------##
# bact = otus == "bacteroides"
# bact = replace(bact, is.na(bact), F)
# prev = otus == "prevotella"
# prev = replace(prev, is.na(prev), F)

bact = otus == "bacteroidaceae"
bact = replace(bact, is.na(bact), F)
prev = otus == "prevotellaceae"
prev = replace(prev, is.na(prev), F)

#looking only at only bact and prev
targ_cols = Reduce( "|", list(prev, bact))
rs = rowSums(ps@otu_table[targ_cols])

philr_mat = matrix(unlist(ps.philr))
View(table(philr_mat))
View(table(ps@otu_table))
View(table(no[, bact]))
View(table(ps.philr[,13]))
View(table(no[, prev]))


#for making plot with lines hidden
drops_class = no[targ_cols] == 1
drops_philr = ps.philr[,"n13"] == -0.169192981810488

drops_all = Reduce( "|", list(drops_class, drops_philr))

plot( ps.philr[,"n13"] ~ log(class_pb_ratio),
      main = paste("Philr vs DADA2 Classifier bacteroidaceae/prevotellaceae"
      )
)
plot( ps.philr[,"n13"] ~ class_pb_ratio,
      main = paste("Philr vs DADA2 Classifier bacteroidaceae/prevotellacea"
      )
)
plot( ps.philr[,"n13"][!drops_all] ~ log(class_pb_ratio[!drops_all]),
      main = paste("Philr vs DADA2 Classifier bacteroidaceae/prevotellaceae\n", 
                   "with pseudo count subtracted out")
)
plot( ps.philr[,"n13"][!drops_all] ~ class_pb_ratio[!drops_all],
      main = paste("Philr vs DADA2 Classifier bacteroidaceae/prevotellaceae\n", 
                   "with pseudo count subtracted out")
)