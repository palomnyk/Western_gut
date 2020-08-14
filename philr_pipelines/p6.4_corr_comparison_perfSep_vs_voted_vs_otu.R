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
##-----------------------Functions----------------------------------##
onlyOne = function(namedVector){
  if(length(namedVector) == 1){
    return(names(namedVector)[1])
  }else {
    return("---")
  }
}

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

philr_pb_ratio = ps.philr[,"n13"]

no = data.frame(ps@otu_table)
tax = data.frame(ps@tax_table)
otus = vector(length = ncol(no), mode = "character")
for (i in 1:ncol(no)){
  asv = names(no)[i]
  otus[i] = tolower(as.character(tax[asv, 5]))#this line sets the taxonomic level
}

uniq_otus = unique(otus)
# 
# pdf(file = file.path(output_dir, "philr_v_class", "philr_vs_Class_abs_abund_bact_prev.pdf"))

# set up 
nodes_vs_otu = read.table(text = "",
                          colClasses = rep(x = "numeric", length(ps@phy_tree$node.label)),
                          col.names = ps@phy_tree$node.label)

nodes_status = c()

for(otu1 in 1:length(uniq_otus)){
  
  other_otus = uniq_otus[-c(1:otu1)]
  
  if (length(other_otus) > 0){
    my_otu = uniq_otus[otu1]
    print(paste("otu1: ", my_otu))
    
    my_otus1 = otus == my_otu
    my_otus1 = replace(my_otus1, is.na(my_otus1), F)
    
    for (otu2 in 1:length(other_otus)){
      print(length(other_otus))
      other_otu = other_otus[otu2]
      print(paste("otu2: ", other_otu))
      my_otus2 = otus == other_otu
      my_otus2 = replace(my_otus2, is.na(my_otus2), F)
    
      class_pb_ratio = c(length = nrow(no))
      for (i in 1:nrow(no)){
        class_pb_ratio[i] = sum(no[i, my_otus1])/ sum(no[i, my_otus2])
      }
      
      # print(paste("sanity check length: length(class_pb_ratio) == length(philr_pb_ratio)"))
      stopifnot(length(class_pb_ratio) == length(philr_pb_ratio), local = TRUE)
      # print(paste("sanity check order:  all.equal(row.names(no), row.names(ps.philr) )"))
      stopifnot(all.equal(row.names(no), row.names(ps.philr)), local = TRUE)
      
      
      #looking only at only target bugs
      targ_cols = Reduce( "|", list(my_otus1, my_otus2))
      rs = rowSums(ps@otu_table[targ_cols])
      
      max_rs = max(rs)
      min_rs = min(rs)
      
      r_2 = c()
      percent_progress = c()
      
      percentage_step = 0.005
      step = (max_rs - min_rs) * percentage_step
      
      loop = 1
      
      # pdf(file = file.path(output_dir, "philr_vs_Class_abs_abund_specific_bact_all_bact.pdf"))
      for (i in seq(min_rs, max_rs, by = step)){
        print( paste( "current i:", i))
        rm = rs > i
        if(all(rm == F)){
          print("all values false")
          break
        }
        my_philr_rat = philr_pb_ratio[rm]
        my_class_rat = class_pb_ratio[rm]
        
        my_lm = lm(my_philr_rat ~ log(my_class_rat ))
        if( is.na( summary(my_lm)$adj.r.squared )){
          print("all values false")
          break
          }
        r_2 = c(r_2, summary(my_lm)$adj.r.squared)
        percent_progress = c(percent_progress,loop * percentage_step)
        
        # plot(my_philr_rat ~ log(my_class_rat),
        #      main = paste0("Philr vs DADA2 Classifier bacteroidaceae/prevotellaceae\n",
        #                    "Hiding specific b/p abs abund < ", i, " of ", max_rs, 
        #                    "\nadj_r^2: ", round(summary(my_lm)$adj.r.squared, 5)),
        #      # col = rm,
        # )
        loop = loop + 1
      }#end loop through abund
      print("end abund loop")
      plot( percent_progress ~ r_2, 
            main = paste(my_otu, "/", other_otu, "\n",
                         "Max r^2:", max(r_2), 
                         "at", percent_progress[match(max(r_2), r_2)][1] * 100, "percent of max abundance hidden"),
      )
      # dev.off()
    }#end otu2 for loop
  }#end if if (length(other_otus) > 0)
}#end otu1 for loop


