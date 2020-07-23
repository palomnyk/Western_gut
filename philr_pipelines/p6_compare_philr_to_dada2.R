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
le = ncol(ps.philr)# * ncol(metadata)#number of rows in df
philrColumns <- vector(length = le)
philrIndex <- vector(length = le)
philrAnnot = vector(length = le)

#taxonomy
phylum = vector(length = le)
class = vector(length = le)
order = vector(length = le)
family = vector(length = le)
genus = vector(length = le)

index <- 1

for( m in 1:ncol(ps.philr))
{
  votes <- name.balance(ps@phy_tree, ps@tax_table, names(ps.philr)[m], return.votes = c('up', 'down'))

  philrColumns[index] <- names(ps.philr)[m]
  philrIndex[index] <- m
  philrAnnot[index] = name.balance(ps@phy_tree, ps@tax_table, names(ps.philr)[m])
  
  genus[index] = paste0(onlyOne(votes$up.votes$Genus), "/", onlyOne(votes$down.votes$Genus))
  family[index] = paste0(onlyOne(votes$up.votes$Family), "/", onlyOne(votes$down.votes$Family))
  order[index] = paste0(onlyOne(votes$up.votes$Order), "/", onlyOne(votes$down.votes$Order))
  class[index] = paste0(onlyOne(votes$up.votes$Class), "/", onlyOne(votes$down.votes$Class))
  phylum[index] = paste0(onlyOne(votes$up.votes$Phylum), "/", onlyOne(votes$down.votes$Phylum))
  
  index <- index + 1
}

print(philrColumns[index - 1])

pValAdj <- p.adjust( pValues, method = "BH" )
# dFrame <- data.frame(pValAdj,philrColumns,philrIndex,metaNames,metaIndex,philrIndex, philrAnnot, genus, family, order, class, phylum)
# dFrame <- data.frame(pValAdj,philrIndex,metaNames, genus, family, order, class, phylum)
# dFrame <- data.frame(pValAdj,philrIndex,metaNames,metaIndex,philrIndex, philrAnnot, genus)
dFrame <- data.frame(philrIndex, genus, family, order, class, phylum)

dFrame <- dFrame [order(dFrame$pValAdj),]

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

hist(philr_pb_ratio, breaks = 50)
barplot(philr_pb_ratio, 1:length((philr_pb_ratio)), col = as.factor(metadata[,"Sample.Group"]))

hist(class_pb_ratio, breaks = 50)

hist(log(philr_pb_ratio), breaks = 50)

hist(log(class_pb_ratio), breaks = 50)

plot(philr_pb_ratio ~ class_pb_ratio, main = "Philr vs DADA2 Classifier \nbacteroidaceae/prevotellaceae")
print(summary(lm(philr_pb_ratio ~ class_pb_ratio)))


plot(philr_pb_ratio ~ log(class_pb_ratio), main = "Philr vs DADA2 Classifier \nbacteroidaceae/prevotellaceae")
print(summary(lm(log(philr_pb_ratio) ~ log(class_pb_ratio))))

pdf(file = file.path(output_dir, "philr_vs_logClass.pdf"))
for(i in 1:ncol(metadata)){
  plot(philr_pb_ratio ~ log(class_pb_ratio), 
       main = paste("Philr vs DADA2 Classifier bacteroidaceae/prevotellaceae\n",
                    names(metadata)[i]),
       col = as.factor(metadata[,i]),
       )
  plot.new()
  legend('topleft',
         title = names(metadata)[i],
         legend = unique(as.factor(metadata[,i])), 
         col = unique(as.factor(metadata[,i])), 
         pch = 19,
         cex=0.5,
         
  )
}
dev.off()


#absolute abundance with counts - points colored
# pdf(file = file.path(output_dir, "philr_vs_Class.pdf"))
for(i in 1:ncol(metadata)){
  plot(philr_pb_ratio ~class_pb_ratio, 
       main = paste("Philr vs DADA2 Classifier bacteroidaceae/prevotellaceae\n",
                    names(metadata)[i]),
       col = as.factor(metadata[,i]),
  )
  plot.new()
  legend('topleft',
         title = names(metadata)[i],
         legend = unique(as.factor(metadata[,i])), 
         col = unique(as.factor(metadata[,i])), 
         pch = 19,
         cex=0.5,
  )
}
# dev.off()

# with(LifeCycleSavings[1:9,], text(sr~dpi, labels = row.names(LifeCycleSavings[1:9,]), pos = 4))
print(summary(lm(philr_pb_ratio ~ log(class_pb_ratio))))


#absolute abundance with counts -hidden points
rs = rowSums(ps@otu_table)

max_rs = max(rs)
min_rs = min(rs)

pdf(file = file.path(output_dir, "philr_vs_Class_abs_abund.pdf"))
for (i in seq(6000, min_rs, by = -500)){
  rm = rs > i
  
  plot(philr_pb_ratio ~ log(class_pb_ratio), 
       main = paste0("Philr vs DADA2 Classifier bacteroidaceae/prevotellaceae\n", 
                    "Hiding abs abund < ", i),
       col = rm,
  )
}
plot(philr_pb_ratio ~ log(class_pb_ratio), 
     main = paste("Philr vs DADA2 Classifier bacteroidaceae/prevotellaceae\nNothing hidden"),
)
dev.off()

#abundance of 

