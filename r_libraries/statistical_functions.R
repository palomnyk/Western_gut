# Author: Aaron Yerke
# Library of statistical functions

##-Load Dependencies------------------------------------------------##

##_Establish constants----------------------------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
# setwd(file.path(home_dir))

##_Functions--------------------------------------------------------##
philr_node_annot_anova_pval <- function(philr_df, phylo_obj, meta_df, only_one = TRUE) {
  # gives philr nodes phylogenetic taxonomy
  # philr_df: a philr transformed df
  # phylo_obj: a phyloseq object from which we will get the taxonmy and the phylogenetic tree for name.balance
  # meta_df: metadata df that you want to iterate over for pvalues
  # Function depencencies:
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install("philr")
  }
  library(philr); packageVersion("philr")
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
  le = ncol(philr_df) * ncol(meta_df)#number of rows in df
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
    
    for( p in 1:ncol(meta_df))
    {
      aLm <- lm( philr_df[,m] ~  meta_df[,p])
      pValues[index]  <- 1
      
      try( pValues[index] <- anova(aLm)$"Pr(>F)"[1])
      metaNames[index] <- names(meta_df)[p]
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
  dFrame <- data.frame(pValues,pValAdj,philrColumns,philrIndex, philrAnnot,metaNames,metaIndex, genus, family, order, class, phylum)
  dFrame <- dFrame [order(dFrame$pValAdj),]
  
  detach(package:philr)
  
  return(dFrame)
}#end philr_node_annot_anova_pval

otu_anova_pval <- function(df1, df2){
  #gives pvalues and adjusted pvalues of otus and metadata
  le = ncol(df1) * ncol(df2)#number of rows in output df
  pValues <-vector(length = le)
  otuColumns <- vector(length = le)
  metaNames <- vector(length = le)
  metaIndex <- vector(length = le)
  otuIndex <- vector(length = le)
  
  index <- 1
  
  for( m in 1:ncol(df1))
  {
    for( p in 1:ncol(df2))
    {
      aLm <- lm( df1[,m] ~  df2[,p])
      pValues[index]  <- 1
      
      try( pValues[index] <- anova(aLm)$"Pr(>F)"[1])
      metaNames[index] <- names(df2)[p]
      otuColumns[index] <- names(df1)[m]
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

kendall_corr <- function(df1, df2){
  #gives pvalues and adjusted pvalues of otus and metadata
  le = ncol(df1) * ncol(df2)#number of rows in output df
  pValues <-vector(length = le)
  otuColumns <- vector(length = le)
  metaNames <- vector(length = le)
  metaIndex <- vector(length = le)
  otuIndex <- vector(length = le)
  
  index <- 1
  
  for( m in 1:ncol(df1))
  {
    for( p in 1:ncol(df2))
    {
      aLm <- lm( df1[,m] ~  df2[,p])
      pValues[index]  <- 1
      
      try( pValues[index] <- anova(aLm)$"Pr(>F)"[1])
      metaNames[index] <- names(df2)[p]
      otuColumns[index] <- names(df1)[m]
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

