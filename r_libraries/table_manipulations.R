# Author: Aaron Yerke
# Library of table manipulating functions

##-Load Dependencies------------------------------------------------##

##_Establish constants----------------------------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
# setwd(file.path(home_dir))

##_Functions--------------------------------------------------------##
equal_num_columns_top_abund <- function(dfs, percent_abund = 0.70) {
  # dfs: list of dfs
  # percent_abund: minimun percent of non-zero cells in column
  # INTERNAL FUNCTIONS:
  percent_abund_filt <- function(abund_percent, df) {
    # finds columns that are above X% populated by nonzeros
    # returns a bool vector same length as # cols
    or_abund = apply(df, 2, function(c){
      sum(c!=0) >= abund_percent*nrow(df)})
    return(df[,or_abund])
  }
  
  top_x_abund_filt <- function(top_x, df) {
    #returns a df with highest abundance top_x number of columns
    df_sums = apply(df, 2, sum)
    df_rank = rank(df_sums, ties.method = "random")
    return(df[,df_rank <= top_x])
  }
  # MAIN METHOD OF FUNCTION
  dfs_abund = lapply(dfs, percent_abund_filt, abund_percent = percent_abund)
  lowest_ncol = min(unlist(lapply(dfs, ncol)))
  dfs_top_x_abund = lapply(dfs_abund, top_x_abund_filt, top_x = lowest_ncol)

  return(dfs_top_x_abund)
}

philr_tutorial_normalization <- function(df) {
  # The philr tutorial https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.R uses this block for normalizations: 
  # ps <-  filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
  # ps <-  filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
  # ps <- transform_sample_counts(ps, function(x) x+1)
  df <- df[apply(df, 2, function(x) sum(x > 3) > (0.2*length(x)))]
  df <-  df[apply(df, 2, function(x) sd(x)/mean(x) > 3.0)]
  df <- df + 1
  return(df)
}
