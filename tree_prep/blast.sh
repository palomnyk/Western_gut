module load blast/2.9.0+

cd ~/git/Western_gut/tree_prep

blastn -query dada2seqs.fasta \
  -db /users/amyerke/git/Western_gut/tree_prep/db/tree \
  -out output.txt \
  -outfmt "6 qseqid sseqid pident length evalue bitscore score ppos"

 # 1.	 qseqid	 query (e.g., unknown gene) sequence id
 # 2.	 sseqid	 subject (e.g., reference genome) sequence id
 # 3.	 pident	 percentage of identical matches
 # 4.	 length	 alignment length (sequence overlap)
 # 5.	 evalue	 expect value
 # 6.	 bitscore	 bit score
 # 7.  score     Raw score
 # 8.  ppos      Percentage of positive-scoring matches