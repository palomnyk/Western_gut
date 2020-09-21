# Author: Aaron Yerke
# Library of table manipulating functions

##-Load Dependencies------------------------------------------------##

##_Establish constants----------------------------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
# setwd(file.path(home_dir))

##_Functions--------------------------------------------------------##
equal_num_columns <- function(df1, df2, percent_abund = 0.70) {
  # INTERNAL FUNCTIONS:
  percent_abund_filt <- function(abund_percent, df) {
    # finds columns that are above X% populated by nonzeros
    # returns a bool vector same length as # cols
    or_abund = apply(df, 2, function(c){
      sum(c!=0) >= abund_percent*nrow(df)})
  }
  
  top_x_abund_filt <- function(top_x, df) {
    #returns a df with highest abundance top_x number of columns
    df_sums = apply(df, 2, sum)
    df_rank = rank(df_sums, ties.method = "random")
    return(df[,df_rank <= top_x])
  }
  # MAIN METHOD OF FUNCTION
  df1_abund <- df1[,percent_abund_filt(percent_abund, df1)]
  df2_abund <- df2[,percent_abund_filt(percent_abund,df2)]
  lowest_ncol <- min( c( ncol(df1_abund), ncol(df2_abund)))
  return(
    list(top_x_abund_filt(lowest_ncol, df2_abund),
         top_x_abund_filt(lowest_ncol, df1_abund))
  )
}


