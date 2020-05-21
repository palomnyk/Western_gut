# Author: Aaron Yerke
# This script is for exploring the B/V ratio and determining if it is the best ratio

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

##----------------------Make otu tables-----------------------------##

no = data.frame(ps@otu_table)
tax = data.frame(ps@tax_table)
otus = vector(length = ncol(no), mode = "character")
for (i in 1:ncol(no)){
  asv = names(no)[i]
  otus[i] = tolower(as.character(tax[asv, 6]))
}

##----------------------calculate ratios----------------------------##
bact = otus == "bacteroides"
bact = replace(bact, is.na(bact), F)
prev = otus == "prevotella"
prev = replace(prev, is.na(prev), F)

ratios = c(length = nrow(no))

for (i in 1:nrow(no)){
  ratios[i] = sum(no[i, bact])/ sum(no[i, prev])
}

pvals = vector(length = ncol(ps.philr), mode = "numeric")
rsqs = vector(length = ncol(ps.philr), mode = "numeric")
nod_nam = vector(length = ncol(ps.philr), mode = "integer")
index = vector(length = ncol(ps.philr), mode = "integer")

for (i in 1:ncol(ps.philr)){
  myLm = lm(ps.philr[,i] ~ ratios)
  pvals[i] = anova(myLm)$"Pr(>F)"[1]
  rsqs[i] = summary(myLm)$r.squared
  nod_nam[i] = names(ps.philr)[i]
  index[i] = i
}

hist(pvals, breaks = 50)

dFrame = data.frame(pvals, rsqs, nod_nam, index)
dFrame <- dFrame [order(dFrame$pvals),]
dFrame$pvalsAdj <- p.adjust( dFrame$pvals, method = "BH" )

write.table(dFrame, file=file.path(output_dir, "philrVbp.tsv"), row.names=FALSE, sep="\t")

pdf(file.path(output_dir, "philrVbp.pdf"))
# par(mfrow=c(2,2))
for( i in 1:nrow(dFrame))
{
  #print(paste(dFrame$philrIndex[i], " ",  names(ps.philr)[dFrame$metaIndex[i]], " ", dFrame$pValuesAdjusted[i]))
  aTitle <- paste(  dFrame$nod_nam[i], "Adjusted Pvalue=",
                    dFrame$pvalsAdj[i], "R^2:", dFrame$rsqs[i])
  plot(ps.philr[,dFrame$index[i]] ~ ratios, 
       las=2,
       main=aTitle, 
       xlab="bacteroides/prevolta ratios", 
       ylab="Balances"
  )
}

dev.off()

print(paste("number of significant hits:", sum(dFrame$pvalsAdj < 0.05), "out of:", length(dFrame$pvalsAdj) ))
