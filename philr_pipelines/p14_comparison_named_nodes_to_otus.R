# Author: Aaron Yerke
# Compare the named nodes to named OTUs

##-Functions--------------------------------------------------------##
philr_node_annot_anova_pval <- function(philr_df, phylo_obj, meta_df, only_one = TRUE) {
  # gives philr nodes phylogenetic taxonomy
  # philr_df: a philr transformed df
  # phylo_obj: a phyloseq object from which we will get the taxonmy and the phylogenetic tree for name.balance
  # meta_df: metadata df that you want to iterate over for pvalues
  # INTERNAL FUNCTIONS:
  onlyOne = function(namedVector){
    # internal function for finding annotations with only one option (meaning perfect seperation)
    if(length(namedVector) == 1){
      return(names(namedVector)[1])
    }else {
      return("---")
    }
  }
  # MAIN METHOD OF FUNCTION
  le = ncol(philr_df) * ncol(metadata)#number of rows in df
  pValues <-vector(length = le)
  philrColumns <- vector(length = le)
  metaNames <- vector(length = le)
  metaIndex <- vector(length = le)
  philrIndex <- vector(length = le)
  
  #taxonomy
  philrAnnot = vector(length = le)
  phylum = vector(length = le)
  class = vector(length = le)
  order = vector(length = le)
  family = vector(length = le)
  genus = vector(length = le)
  
  index <- 1
  
  for( m in 1:ncol(philr_df))
  {
    votes <- name.balance(phylo_obj@phy_tree, phylo_obj@tax_table, names(philr_df)[m], return.votes = c('up', 'down'))
    
    for( p in 1:ncol(metadata))
    {
      aLm <- lm( philr_df[,m] ~  metadata[,p])
      pValues[index]  <- 1
      
      try( pValues[index] <- anova(aLm)$"Pr(>F)"[1])
      metaNames[index] <- names(metadata)[p]
      philrColumns[index] <- names(philr_df)[m]
      metaIndex[index] <- p
      philrIndex[index] <- m
      philrAnnot[index] = name.balance(phylo_obj@phy_tree, phylo_obj@tax_table, names(philr_df)[m])
      
      if (only_one == TRUE){
        genus[index] = paste0(onlyOne(votes$up.votes$Genus), "/", onlyOne(votes$down.votes$Genus))
        family[index] = paste0(onlyOne(votes$up.votes$Family), "/", onlyOne(votes$down.votes$Family))
        order[index] = paste0(onlyOne(votes$up.votes$Order), "/", onlyOne(votes$down.votes$Order))
        class[index] = paste0(onlyOne(votes$up.votes$Class), "/", onlyOne(votes$down.votes$Class))
        phylum[index] = paste0(onlyOne(votes$up.votes$Phylum), "/", onlyOne(votes$down.votes$Phylum))
      }else{
        genus[index] = paste0(names(which.max(votes$up.votes$Genus)), "/", names(which.max(votes$down.votes$Genus)))
        family[index] = paste0(names(which.max(votes$up.votes$Family)), "/", names(which.max(votes$down.votes$Family)))
        order[index] = paste0(names(which.max(votes$up.votes$Order)), "/", names(which.max(votes$down.votes$Order)))
        class[index] = paste0(names(which.max(votes$up.votes$Class)), "/", names(which.max(votes$down.votes$Class)))
        phylum[index] = paste0(names(which.max(votes$up.votes$Phylum)), "/", names(which.max(votes$down.votes$Phylum)))
      }
      index <- index + 1
    }
  }
  
  pValAdj <- p.adjust( pValues, method = "BH" )
  dFrame <- data.frame(pValAdj,philrColumns,philrIndex, philrAnnot,metaNames,metaIndex, genus, family, order, class, phylum)
  dFrame <- dFrame [order(dFrame$pValAdj),]
  return(dFrame)
}#end philr_node_annot_anova_pval

otu_anova_pval <- function(otu_df, meta_df){
  #gives pvalues and adjusted pvalues of otus and metadata
  le = ncol(otu_df) * ncol(meta_df)#number of rows in output df
  pValues <-vector(length = le)
  otuColumns <- vector(length = le)
  metaNames <- vector(length = le)
  metaIndex <- vector(length = le)
  otuIndex <- vector(length = le)
  
  index <- 1
  
  for( m in 1:ncol(otu_df))
  {
    for( p in 1:ncol(meta_df))
    {
      aLm <- lm( otu_df[,m] ~  meta_df[,p])
      pValues[index]  <- 1
      
      try( pValues[index] <- anova(aLm)$"Pr(>F)"[1])
      metaNames[index] <- names(meta_df)[p]
      otuColumns[index] <- names(otu_df)[m]
      metaIndex[index] <- p
      otuIndex[index] <- m
      index <- index + 1
    }
  }
  pValAdj <- p.adjust( pValues, method = "BH" )
  dFrame <- data.frame(pValAdj,otuColumns,otuIndex,metaNames,metaIndex)
  dFrame <- dFrame [order(dFrame$pValAdj),]
  return(dFrame)
}

##_Establish constants----------------------------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
# setwd(file.path(home_dir))

##-Load Dependencies------------------------------------------------##
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

##-Load data--------------------------------------------------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects","ps_philr_transform.rds"))
ps = readRDS(con)
close(con)

philr_trans = read.table(file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_philr_transform.csv"),
                         sep = ",",
                         header = TRUE)
philr_trans = philr_trans - min(philr_trans)

otu_tab = read.table(file.path(home_dir, "philr_pipelines", "tables", "otu_table.csv"),
                     sep = ",",
                     header = TRUE)

metadata = read.table(file.path("philr_pipelines", "tables", "ps_sample_data_reduced.csv"), 
                      sep=",", 
                      header=TRUE)
##-Run Comparison---------------------------------------------------##

philr_anova <- philr_node_annot_anova_pval(philr_df = philr_trans,
                                    phylo_obj = ps,
                                    meta_df = metadata)

otu_anova <- otu_anova_pval(otu_df = otu_tab,
                            meta_df = metadata)

count = 0

for (name1 in 1:nrow(otu_anova)){
  my_otu <- tolower(otu_anova$otuColumns[name1])
  for (name2 in 1:nrow(philr_anova)){
    my_nodes <- as.character(philr_anova$family[name2])
    my_node_split <- unlist(strsplit(my_nodes, "/"))
    for (nm in my_node_split){
      nm = tolower(nm)
      if (nm ==  my_otu &
          otu_anova$metaNames[name1] == philr_anova$metaNames[name2] &
          philr_anova$pValAdj[name2] < 0.05 &
          philr_anova$pValAdj[name2] < 0.05
          ){
        
        count = count + 1
        # print(paste("name1:", name1, "nm: ", nm, "nms:", my_nodes[1], my_nodes[2]))
        # for (met in 1:ncol(metadata)){
        # 
        # }
        # metadata columns where my_otu is significant in otu_anova
        # print(paste("my_otu:", my_otu, otu_anova$pValAdj[name1], otu_anova$metaNames[name1]))
        # print(paste("my_nodes:", my_nodes, philr_anova$pValAdj[name2], philr_anova$metaNames[name2]))
        
      }
    }
  }
}


#nested apply function
# apply(b, 2, function(x) apply(a, 2, function(y) summary(lm(x~y))$coefficients[2,4]))


