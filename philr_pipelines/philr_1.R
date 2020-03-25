# Author: Aaron Yerke
# This is a pipeline that was created to explore the philr transform.
# Many resources were used to create this pipeline:
#   Processing the sequences through DADA2, making trees with phangorn, and putting them into phyloseq objects:
#     This is by far the biggest source:
#     https://github.com/spholmes/F1000_workflow/blob/master/MicrobiomeWorkflow/MicrobiomeWorkflowII.Rmd
#     or (in article format)
#     https://f1000research.com/articles/5-1492/v2
#     
#   DADA2:
#     https://benjjneb.github.io/dada2/tutorial.html
#   
#   Using philr:
#     https://bioconductor.org/packages/release/bioc/html/philr.html
#     Code example:
#       https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.R

rm(list = ls()) #clear workspace

# ‘ape’, ‘dplyr’, ‘reshape2’, ‘plyr’
# .cran_packages <- c("ggplot2", "gridExtra")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2",type = "source", checkBuilt = TRUE)

library("dada2")

home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
project = "RDP Western Gut"
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
# list.files(f_path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq
fnFs <- sort(list.files(f_path, pattern="_1.fastq", full.names = TRUE))
sampleNames <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filt_path <- file.path(f_path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sampleNames, "_R1_filt.fastq"))

#Filter
out <- filterAndTrim(fnFs, filtFs,truncLen=240,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=FALSE, multithread=TRUE) 

dds <- vector("list", length(sampleNames))
names(dds) <- sampleNames

index <-1 

for (f in filtFs){
  errF <- learnErrors(f, multithread = FALSE)
  derepFs <- derepFastq(f, verbose=TRUE)
  dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
  dds[[index]]<-dadaFs
  index<-index+1
}

seqtab <- makeSequenceTable(dds)

#Removing chimeras
seqtab <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE)

# seqtabPath  <-file.path(output_dir,"dada2","ForwardReads_DADA2.rds")
# print(seqtabPath)
# fil <- tempfile(seqtabPath, fileext = ".rds")
setwd(file.path(home_dir, "philr_pipelines"))
saveRDS(seqtab, "ForwardReads_DADA2.rds")
write.table(seqtab,file.path(output_dir,"dada2","ForwardReads_DADA2.txt"),sep="\t")