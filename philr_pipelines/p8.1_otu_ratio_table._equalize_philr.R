# Author: Aaron Yerke
# Script for making ratio table of OTUS for ml purposes
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##
##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
setwd(file.path(home_dir))

PERCENT_ABUND = 0.7

##---------------------------Functions------------------------------##
#returns a bool vector same length as # cols
percent_abund_filt <- function(abund_percent, df) {
  or_abund = apply(df, 2, function(c){
    sum(c!=0) >= abund_percent*nrow(df)})
}
top_x_abund_filt <- function(top_x, df) {
  
}

##------------Import R objects and data preprocessing---------------##
otu_ratios = read.table(file = file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table.csv"),
                        sep = ",")
ps = read.table(file = file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_philr_transform.csv"),
                        sep = ",")
##--------------------Create Ratio table----------------------------##
print(dim(otu_ratios))
print(dim(ps))
ps_abund = ps[,percent_abund_filt(PERCENT_ABUND, ps)]
or_abund = otu_ratios[,percent_abund_filt(PERCENT_ABUND,otu_ratios)]

or_sums = apply(or_abund, 2, sum)
or_rank = rank(or_sums, ties.method = "random")

or_top_x = or_abund[,or_rank <= 35]

write.table(or_top_x, 
            file = file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table_top_35.csv"),
            sep = ",")


