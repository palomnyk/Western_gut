# Author: Aaron Yerke
# This is a pipeline that was created to process Western gut data through DADA2.
# Many resources were used to create this pipeline:
#   DADA2:
#     https://benjjneb.github.io/dada2/tutorial.html

rm(list = ls()) #clear workspace
library("dada2")

home_dir = file.path('~','git','Western_gut')
f_path <- "~/git/Western_gut/sequences" # CHANGE ME to the directory containing the fastq files after unzipping.
output_dir = file.path('~','git','Western_gut', 'output',"dada2")

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

saveRDS(seqtab, paste0(output_dir,"/ForwardsReads_DADA2.rds"))
write.table(seqtab,file=paste0(output_dir,"/ForwardReads_DADA2.txt"),sep="\t")



