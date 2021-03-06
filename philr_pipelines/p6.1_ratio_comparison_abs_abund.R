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

# ps <- transform_sample_counts(ps, function(x) x+1)

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

# hist(philr_pb_ratio, breaks = 50)
# barplot(philr_pb_ratio, 1:length((philr_pb_ratio)), col = as.factor(metadata[,"Sample.Group"]))
# 
# hist(class_pb_ratio, breaks = 50)
# 
# hist(log(philr_pb_ratio), breaks = 50)
# 
# hist(log(class_pb_ratio), breaks = 50)

#absolute abundance with counts -hidden points
rs = rowSums(ps@otu_table[])

max_rs = max(rs)
min_rs = min(rs)

r_2 = c()
percent_progress = c()

percentage_step = 0.005
step = (max_rs - min_rs) * percentage_step

loop = 1

# pdf(file = file.path(output_dir, "philr_vs_Class_all_abs_abund_bact_prev.pdf"))
for (i in seq(min_rs, max_rs, by = step)){
  print(i)
  rm = rs > i
  if(all(rm == F)){stop("all values false")}
  
  my_philr_rat = philr_pb_ratio[rm]
  my_class_rat = class_pb_ratio[rm]
  
  my_lm = lm(my_philr_rat ~ log(my_class_rat ))
  # if( is.na( summary(my_lm)$adj.r.squared )){stop("all values false")}
  myCorr = cor.test( my_philr_rat, log(my_class_rat),
                     method = "kendall")
  r_2 = c(r_2, myCorr$p.value)
  
  # r_2 = c(r_2, summary(my_lm)$adj.r.squared)
  percent_progress = c(percent_progress, loop * percentage_step)
  
  plot(my_philr_rat ~ log(my_class_rat),
       main = paste0("Philr vs DADA2 Classifier bacteroidaceae/prevotellaceae\n",
                     "Hiding b/p all abs abund < ", i, " of ", max_rs, 
                     "\nadj_r^2: ", round(summary(my_lm)$adj.r.squared, 5)),
       # col = rm,
  )
  loop = loop + 1
}
plot( percent_progress ~ r_2, 
      main = paste(bact, "/", prev, "\n",
                   "Max r^2:", round(max(r_2, na.rm = T), 5), 
                   "at", percent_progress[match(max(r_2), r_2)][1],
                   "of max abundance hidden"),
)
# dev.off()


#looking only at only bact and prev
targ_cols = Reduce( "|", list(prev, bact))
rs = rowSums(ps@otu_table[targ_cols])

max_rs = max(rs)
min_rs = min(rs)

r_2 = c()
percent_progress = c()

percentage_step = 0.01
step = (max_rs - min_rs) * percentage_step

loop = 1

pdf(file = file.path(output_dir, "philr_vs_Class_abs_abund_specific_bact_prev.pdf"))
for (i in seq(min_rs, max_rs, by = step)){
  print(i)
  rm = rs > i
  if(all(rm == F)){stop("all values false")}

  my_philr_rat = philr_pb_ratio[rm]
  my_class_rat = class_pb_ratio[rm]

  my_lm = lm(my_philr_rat ~ log(my_class_rat ))
  
  myCorr = cor.test( my_philr_rat, log(my_class_rat),
                     method = "kendall")
  r_2 = c(r_2, myCorr$p.value)
  # r_2 = c(r_2, summary(my_lm)$adj.r.squared)
  percent_progress = c(percent_progress,loop * percentage_step)

  plot(my_philr_rat ~ log(my_class_rat),
       main = paste0("Philr vs DADA2 Classifier bacteroidaceae/prevotellaceae\n",
                     "Hiding specific b/p abs abund < ", i, " of ", max_rs, 
                     "\nadj_r^2: ", round(summary(my_lm)$adj.r.squared, 5),
                     "\nKendall corr: ", myCorr$p.value,
                     ),
       # col = rm,
  )
  loop = loop + 1
}
plot( percent_progress ~ r_2, 
     main = paste("Max r^2:", max(r_2), 
                  "at", percent_progress[match(max(fr_2), r_2)][1] * 100, "percent of max abundance hidden"),
)
dev.off()

philr_mat = matrix(unlist(ps.philr))
# View(table(philr_mat))
# View(table(ps@otu_table))
# View(table(no[, bact]))
# View(table(ps.philr[,13]))
# View(table(no[, prev]))


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