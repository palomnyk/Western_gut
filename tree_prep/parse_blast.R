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

df <- data.frame(#qseqid=character(),
                 sseqid=character(),
                 bitscore=integer(),
                 stringsAsFactors=FALSE)
con = file(input_file, "r")
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  result = unlist(strsplit(line, '\t'))
  qseq = result[1]
  sseq = unlist(strsplit(result[2],"\\|"))[2] #dbj|AB064923|
  bitsc = as.integer(result[6])
  # print(paste(sseq, bitsc))
  if (! qseq %in% row.names(df)){
    newRow <- data.frame(sseqid=sseq,
                         bitscore=bitsc) 
    row.names(newRow) = qseq
    df = rbind(df, newRow)
  }else{
    if (bitsc > df[qseq, "bitscore"]){
      df[qseq,] = list(sseq, bitsc)
    }
  }
  # break
}
close(con)

write.csv(df, file = "parsed_output.csv", sep = ",")
