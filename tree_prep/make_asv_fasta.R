#creating fasta for blasting from DADA2 results
home_dir = file.path('~','git','Western_gut')

con <- gzfile(file.path( home_dir, "philr_pipelines", "r_objects", "ForwardReads_DADA2.rds"))
seqtab = readRDS(con)

setwd(file.path(home_dir, "tree_preprocessing_blast"))

for (i in 1:ncol(seqtab)){
  seq = colnames(seqtab)[i]
  cat(paste0(">", i, "\n", seq, "\n"), 
      file = "dada2seqs.fasta", 
      append=TRUE)
}