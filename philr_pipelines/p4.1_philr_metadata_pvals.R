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

##------------Import R objects and data preprocessing----------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects","philr_transform.rds"))
ps.philr = data.frame(readRDS(con))

con <- gzfile(file.path( "philr_pipelines", "r_objects","ps_philr_transform.rds"))
ps = readRDS(con)

metadata = data.frame(sample_data(ps))
metadata = metadata[, !apply(is.na(metadata), 2, all)]#remove columns where all are NA
metadata[c(1:3)] <- list(NULL)
drop = c("sample_accession","secondary_sample_accession","experiment_accession",             
        "tax_id","scientific_name","instrument_model","library_layout",
        "experiment_title","sample_title", "study_title", "run_alias",
        "fastq_ftp" ,"fastq_galaxy","submitted_ftp"                    
        ,"submitted_galaxy","sra_ftp", "sra_galaxy", "study_accession")
metadata = metadata[ , !(names(metadata) %in% drop)]
##----------------------t-test w metadata----------------------------##


pValues <-vector()
philrColumns <- vector()
metaNames <- vector()
metaIndex <- vector()
philrIndex <- vector()
index <- 1

for( m in 1:ncol(ps.philr))
{
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
    index <- index + 1
  }
}

dFrame <- data.frame(pValues,philrColumns,philrIndex,metaNames,metaIndex,philrIndex)

dFrame <- dFrame [order(dFrame$pValues),]
dFrame$pValuesAdjusted<- p.adjust( dFrame$pValues, method = "BH" )

write.table(dFrame, file=file.path(output_dir, "philrVmetadata.tsv"), row.names=FALSE, sep="\t")

pdf(file.path(output_dir, "philrVmetadata.pdf"))
# par(mfrow=c(2,2))

for( i in 1:nrow(dFrame))
{
  #print(i)
  #print(paste(dFrame$philrIndex[i], " ",  names(ps.philr)[dFrame$metaIndex[i]], " ", dFrame$pValuesAdjusted[i]))
  aTitle <- paste(  names(metadata)[as.numeric(dFrame$philrIndex[i])], "vs",  names(ps.philr)[dFrame$metaIndex[i]], "\nAdjusted Pvalue=",
                    dFrame$pValuesAdjusted[i])
  
  plot( ps.philr[,dFrame$philrIndex[i]] ~ as.factor(metadata[,dFrame$metaIndex[i]]),main=aTitle, xlab=paste(dFrame$philrIndex[i]), ylab=names(ps.philr)[dFrame$metaIndex[i]])
}

dev.off()


