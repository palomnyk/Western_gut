# Author: Aaron Yerke
# Script for making ratio table of OTUS for ml purposes
# This was helpful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace


##---------------------Establish constants--------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
setwd(file.path(home_dir))

PERCENT_ABUND = 0.7

##-Load Dependencies------------------------------------------------##
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

##------------------Import R objects and tables---------------------##
otu_ratios = read.table(file = file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table.csv"),
                        sep = ",")
# philr_trans = read.table(file = file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_philr_transform.csv"),
#                         sep = ",")
philr_trans = read.table(file = file.path(home_dir, "philr_pipelines", "tables", "denovo_ps_philr_transform.csv"),
                sep = ",")
##--------------------Create Ratio table----------------------------##
mylist = 
write.table(or_top_x,
            file = file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table_top_35.csv"),
            sep = ",")

# write.table(ps_top_x,
#             file = file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_top_35.csv"),
#             sep = ",")
write.table(ps_top_x,
            file = file.path(home_dir, "philr_pipelines", "tables", "denovo_tree_ps_top_35.csv"),
            sep = ",")



otu_ratio = read.table(file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table.csv"),
                       sep = ",",
                       header = TRUE)

philr_trans = read.table(file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_philr_transform.csv"),
                         sep = ",",
                         header = TRUE)

otu_tab = read.table(file.path(home_dir, "philr_pipelines", "tables", "otu_table.csv"),
                     sep = ",",
                     header = TRUE)

asv_table = read.table(file.path(home_dir, "philr_pipelines", "tables", "asv_table.csv"),
                       sep = ",",
                       header = TRUE)

my_datasets = list(asv_table, otu_tab, otu_ratio, philr_trans, cbind(otu_ratio, philr_trans), cbind(otu_tab, otu_ratio))
my_ds_names = c("ASV", "OTU", "OTU ratio", "philR", "OTU + philR", "OTU + OTU ratio")

test = equal_num_columns_top_abund(my_datasets)




