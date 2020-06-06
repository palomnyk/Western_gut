# Author: Aaron Yerke
# Convert output from p1_dada2_rd1.R to silva based taxonomy. Loosely following :
#   http://benjjneb.github.io/dada2/training.html

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("dada2",type = "source", checkBuilt = TRUE)
  BiocManager::install("DECIPHER")
}

library("DECIPHER")
library("dada2")
library("phyloseq")

home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')

setwd(file.path(home_dir, "philr_pipelines", "taxonomy", "silva"))

# Download and extract the Mothur-formatted taxonomy files (see
#                https://mothur.org/wiki/silva_reference_files/),
# download.file(
#   "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138.tgz",
#   file.path( "silva.nr_v138.tgz")
# )
# command <- paste("tar -xvf", 
#                  file.path("silva.nr_v138.tgz"),
#                  "-C", 
#                  file.path(".")
# )
# system(command)
# list.files()
# 
# #Download the Silva v138 license,
# download.file(
#   "https://www.arb-silva.de/fileadmin/silva_databases/release_138/LICENSE.txt",
#   file.path("SILVA_LICENSE.txt")
# )
# 
# #get tree
# download.file(
#   "https://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_132/LTPs132_SSU_tree.newick",
#   file.path(home_dir,"philr_pipelines", "LTPs132_SSU_tree.newick")
# )
# 
# # Download the file `SILVA_138_SSURef_tax_silva.fasta.gz` from the Silva website,
# download.file(
#   "https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_tax_silva.fasta.gz",
#   file.path("SILVA_138_SSURef_tax_silva.fasta.gz")
# )
# Create the DADA2-formatted taxonomy database file,
dada2:::makeTaxonomyFasta_Silva(
  file.path( "silva.nr_v138.align"), 
  file.path( "silva.nr_v138.tax"), 
  file.path( "silva_nr_v138_train_set.fa.gz")
)

#Create the DADA2-formatted species database file,
dada2:::makeSpeciesFasta_Silva(
  file.path("SILVA_138_SSURef_tax_silva.fasta.gz"),
  file.path("silva_species_assignment_v138.fa.gz")
)

#
setwd(home_dir)
myMeta = read.table(file.path("fullMetadata.tsv"), 
                    sep="\t", 
                    header=TRUE, 
                    row.names = "run_accession", 
                    check.names = FALSE,
                    stringsAsFactors=FALSE)


con <- gzfile(file.path( home_dir, "philr_pipelines", "r_objects", "ForwardReads_DADA2.rds"))
seqtab = readRDS(con)

fastaRef <- file.path(home_dir, 'philr_pipelines', "taxonomy", "silva", "silva_nr_v138_train_set.fa.gz")
taxTab <- assignTaxonomy(seqtab, refFasta = fastaRef, multithread=TRUE)
species <- assignSpecies(seqtab, refFasta = file.path(home_dir, 'philr_pipelines', "taxonomy", "silva", "SILVA_138_SSURef_tax_silva.fasta.gz"))
unname(head(taxTab))

seqs <- getSequences(seqtab)

names(seqs) <- seqs # This propagates to the tip labels of the tree

# alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

# saveRDS(taxTab, file.path(home_dir, 'philr_pipelines', "r_objects", "silvaSpeciesForwardReads_DADA2_taxonomy.rds"))
# saveRDS(taxTab, file.path(home_dir, 'philr_pipelines', "r_objects", "silvaForwardReads_DADA2_taxonomy.rds"))
# saveRDS(alignment, file.path(home_dir, 'philr_pipelines', "r_objects","silvaForwardReads_DADA2_alignment.rds"))

species <- assignSpecies(seqtab, refFasta = file.path("philr_pipelines", "taxonomy", "silva", "silva_species_assignment_v138.fa.gz"))

print("Alignment completed")

#get the taxonmy to leave in tree
taxnames = vector()
for (i in 1:ncol(taxTab)){
  taxnames = c(taxnames, unique(taxTab[,i]))
}

for (i in 1:nrow(taxTab)){
  taxnames = c(taxnames, unique(taxTab[,i]))
}

tree = read_tree(file.path(home_dir,"philr_pipelines",  "taxonomy" , "silva","viFy10M5J2nvIBpCLM-QMQ_newick.txt"))
# plot_tree(tree, nodelabf=nodeplotblank, label.tips="taxa_names", ladderize="left")

tree = prune_taxa(taxnames, tree)

plot_tree(tree)

ps <- phyloseq(otu_table(seqTa, taxa_are_rows=FALSE),
               sample_data(myMeta),
               tax_table(taxTab),
               phy_tree(tree)
)

length(unique(row.names(taxTab)))
nrow(taxTab)

seqTa = seqtab[,row.names(taxTab)]

for (i in 1:nrow(taxTab)){
  seq = row.names(taxTab(i))
  cat(paste0(">", i, "\n", seq, "\n"), file = "dada2seqs.fasta", append=TRUE)
}

