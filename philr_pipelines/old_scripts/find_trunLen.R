# This script generates images that can be used to finetune the length of 
# the truncLen parameter of the filterAndTrim DADA2 function.

rm(list = ls()) #clear workspace

library(dada2); packageVersion("dada2")

# setwd(file.path('~','git','Western_gut'))
output_dir = file.path('~','git','Western_gut', 'output')
project = "RDP Western Gut"
#f_path = c('.', 'biolockj', 'rdp_2020Feb19', '04_BuildTaxaTables','output','rdp_2020Feb19_taxaCount_genus.tsv')
f_path <- "~/git/Western_gut/sequences" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(f_path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(f_path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(f_path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualF = plotQualityProfile(fnFs[1:2])

plotQualR = plotQualityProfile(fnRs[1:2])
dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir)

png(filename="plotQualF.png")
plot(plotQualF)
dev.off()

png(filename="plotQualR.png")
plot(plotQualR)
dev.off()