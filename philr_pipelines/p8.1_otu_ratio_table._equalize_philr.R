# Author: Aaron Yerke
# Script for making ratio table of OTUS for ml purposes
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##
##---------------------Establish constants--------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
setwd(file.path(home_dir))

PERCENT_ABUND = 0.7
##---------------------------Functions------------------------------##

# finds columns that are above X% populated by nonzeros
# returns a bool vector same length as # cols
percent_abund_filt <- function(abund_percent, df) {
  or_abund = apply(df, 2, function(c){
    sum(c!=0) >= abund_percent*nrow(df)})
}

#returns a df with highest abundance top_x number of columns
top_x_abund_filt <- function(top_x, df) {
  df_sums = apply(df, 2, sum)
  df_rank = rank(df_sums, ties.method = "random")
  return(df[,df_rank <= top_x])
}
##------------------Import R objects and tables---------------------##
otu_ratios = read.table(file = file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table.csv"),
                        sep = ",")
# philr_trans = read.table(file = file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_philr_transform.csv"),
#                         sep = ",")
philr_trans = read.table(file = file.path(home_dir, "philr_pipelines", "tables", "denovo_ps_philr_transform.csv"),
                sep = ",")

##--------------------Create Ratio table----------------------------##

ps_abund <- philr_trans[,percent_abund_filt(PERCENT_ABUND, philr_trans)]
or_abund <- otu_ratios[,percent_abund_filt(PERCENT_ABUND,otu_ratios)]
lowest_ncol <- min( c( ncol(ps_abund), ncol(or_abund)))
or_top_x <- top_x_abund_filt(lowest_ncol, or_abund)
ps_top_x <- top_x_abund_filt(lowest_ncol, ps_abund)

write.table(or_top_x, 
            file = file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table_top_35.csv"),
            sep = ",")

# write.table(ps_top_x,
#             file = file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_top_35.csv"),
#             sep = ",")
write.table(ps_top_x, 
            file = file.path(home_dir, "philr_pipelines", "tables", "denovo_tree_ps_top_35.csv"),
            sep = ",")