#!/usr/bin/env Rscript

# makeblastdb was complaining about the length of the fasta headers
# writing this was quicker than fixing and rerunning download_seqs.R

# arg 1 is the file that we will need to fix

args = commandArgs(trailingOnly=TRUE)

home_dir = file.path('~','git','Western_gut')
setwd(file.path(home_dir, "tree_prep"))


if (exists(args[1])){
  oldFasta = args[1]
}else{
  oldFasta = "treeFasta.fasta"
}

fixedFasta = paste0("fixed_", oldFasta)

fastaHeaders = vector()

con = file(oldFasta, "r")
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  # print(line)
  if (grepl("_", line, fixed = TRUE)){
    id = strsplit(line, "_")[[1]][1]
    if (!id %in% fastaHeaders){
      fastaHeaders = c(id, fastaHeaders)
      cat(paste0(id, "\n"), file = fixedFasta, append=TRUE)
    }else{
      print(id)
    }
  } else{
    cat(paste0(line, "\n"), file = fixedFasta, append=TRUE)
  }
}
close(con)

