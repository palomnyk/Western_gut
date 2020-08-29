# Author: Aaron Yerke
# Script for ratio table of OTUS to 
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!requireNamespace("randomForest", quietly = TRUE)){
  install.packages("randomForest")
}
library(randomForest)
##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
setwd(file.path(home_dir))

##------------Import tables    and data preprocessing---------------##
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

metadata = read.table(file.path("philr_pipelines", "tables", "ps_sample_data.csv"), 
                    sep=",", 
                    header=TRUE)

##-----------------Create training/testing sets---------------------##
set.seed(36)
train_index = sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
test_index = c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]

resp_var_train = metadata$Sample.Group[train_index]
resp_var_test = metadata$Sample.Group[test_index]

# my_table = cbind(otu_tab, otu_ratio)
# my_table = asv_table
# my_table = otu_tab
my_table = cbind(otu_ratio, philr_trans)

my_table_train = my_table[train_index,]
my_table_test = my_table[test_index,]
##--------------------------Build model ----------------------------##
rf <- randomForest(
  my_table_train, resp_var_train
)
##-------------------------Test model------------------------------##
pred = predict(rf, newdata=my_table_test)
cm = table(resp_var_test, pred)
print(cm)
cm_percent = cm/sum(cm)

heatmap(cm, 
        symm = T,
        Rowv = FALSE,
        keep.dendro = FALSE,
        verbose = TRUE)
