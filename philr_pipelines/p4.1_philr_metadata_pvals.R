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
con <- gzfile(file.path( "philr_pipelines", "r_objects","philr_transform.rds"))
ps.philr = data.frame(readRDS(con))

con <- gzfile(file.path( "philr_pipelines", "r_objects","ps_philr_transform.rds"))
ps = readRDS(con)

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
le = ncol(ps.philr) * ncol(metadata)#number of rows in df
pValues <-vector(length = le)
philrColumns <- vector(length = le)
metaNames <- vector(length = le)
metaIndex <- vector(length = le)
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
  
  for( p in 1:ncol(metadata))
  {
    #print(m)
    #print(ps.philr[,m])
    aLm <- lm( ps.philr[,m] ~  metadata[,p])
    
    pValues[index]  <-1
    
    try( pValues[index] <- anova(aLm)$"Pr(>F)"[1])
    metaNames[index] <- names(metadata)[p]
    philrColumns[index] <- names(ps.philr)[m]
    metaIndex[index] <- p
    philrIndex[index] <- m
    philrAnnot[index] = name.balance(ps@phy_tree, ps@tax_table, names(ps.philr)[m])
  
    # genus[index] = paste0(names(which.max(votes$up.votes$Genus)), "/", names(which.max(votes$down.votes$Genus)))
    # family[index] = paste0(names(which.max(votes$up.votes$Family)), "/", names(which.max(votes$down.votes$Family)))
    # order[index] = paste0(names(which.max(votes$up.votes$Order)), "/", names(which.max(votes$down.votes$Order)))
    # class[index] = paste0(names(which.max(votes$up.votes$Class)), "/", names(which.max(votes$down.votes$Class)))
    # phylum[index] = paste0(names(which.max(votes$up.votes$Phylum)), "/", names(which.max(votes$down.votes$Phylum)))
    
    genus[index] = paste0(onlyOne(votes$up.votes$Genus), "/", onlyOne(votes$down.votes$Genus))
    family[index] = paste0(onlyOne(votes$up.votes$Family), "/", onlyOne(votes$down.votes$Family))
    order[index] = paste0(onlyOne(votes$up.votes$Order), "/", onlyOne(votes$down.votes$Order))
    class[index] = paste0(onlyOne(votes$up.votes$Class), "/", onlyOne(votes$down.votes$Class))
    phylum[index] = paste0(onlyOne(votes$up.votes$Phylum), "/", onlyOne(votes$down.votes$Phylum))
     
    index <- index + 1
  }
}

print(philrColumns[index])

pValAdj <- p.adjust( pValues, method = "BH" )
# dFrame <- data.frame(pValAdj,philrColumns,philrIndex,metaNames,metaIndex,philrIndex, philrAnnot, genus, family, order, class, phylum)
dFrame <- data.frame(pValAdj,philrIndex,metaNames, genus, family, order, class, phylum)
dFrame <- dFrame [order(dFrame$pValAdj),]
write.table(dFrame, file=file.path(output_dir, "philrVmetadata_only_one.tsv"), row.names=FALSE, sep="\t")

pdf(file.path(output_dir, "philrVmetadata.pdf"))
# par(mfrow=c(2,2))
for( i in 1:nrow(dFrame))
{
  #print(i)
  #print(paste(dFrame$philrIndex[index], " ",  names(ps.philr)[dFrame$metaIndex[i]], " ", dFrame$pValuesAdjusted[i]))
  aTitle <- paste(  dFrame$philrColumns[i], "vs",  dFrame$metaNames[i], "\nAdjusted Pvalue=",
                    dFrame$pValuesAdjusted[i])
  
  plot(ps.philr[,dFrame$philrIndex[i]] ~ as.factor(metadata[,dFrame$metaIndex[i]]), 
       las =2,
       main=aTitle, 
       xlab=paste(dFrame$philrIndex[i]), 
       ylab=names(ps.philr)[dFrame$metaIndex[i]])
}

dev.off()


