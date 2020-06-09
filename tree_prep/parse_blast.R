#!/usr/bin/env Rscript

# script for parsing

# Columns follow this pattern:
# qseqid sseqid pident length evalue bitscore score ppos
# 
# 1.	qseqid	 query (e.g., unknown gene) sequence id
# 2.	sseqid	 subject (e.g., reference genome) sequence id
# 3.	pident	 percentage of identical matches
# 4.	length	 alignment length (sequence overlap)
# 5.	evalue	 expect value
# 6.	bitscore	 bit score
# 7.  score     Raw score
# 8.  ppos      Percentage of positive-scoring matches
# Comparison stategy: compare ppos (if tie, go with biggest alignment length)
args = commandArgs(trailingOnly=TRUE)

home_dir = file.path('~','git','Western_gut')
setwd(file.path(home_dir, "tree_prep"))

if (exists(args[1])){
  input_file = args[1]
}else{
  input_file = "output.txt"
}

output_file = paste0("parsed_", input_file)

df <- data.frame(sseqid=character(),
                 ppos=double(),
                 len=integer(),
                 stringsAsFactors=FALSE)
tree_id = vector()
seqs = vector()

con = file(input_file, "r")
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  result = strsplit(line)
  print(paste(result[1],result[12],result[8],result[4]))
  qseq = result[1]
  sseq = result[2]
  ppos = as.numeric(result[8])
  len = as.integer(result[4])
  print(paste(sseq, ppos, len))
  if (! any(row.names(df) == qseq)){
    df[qseq] = list(sseq, ppos, len)
  }
  
  break
}
close(con)

