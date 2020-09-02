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
if (!requireNamespace("pROC", quietly = TRUE)){
  install.packages("pROC")
}
library(pROC)
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
true_pos = c(0)
true_neg = c(0)
for (grp in unique(resp_var_test)){
  for (rw in 1:nrow(my_table_test)){
    pred = predict(rf, newdata=my_table_test[rw,])
    decision = pred == grp
    truth = resp_var_test[rw] == grp
    if (decision == truth & truth == T){
      true_pos = c(true_pos, tail(true_pos)+1)
      true_neg = c(true_neg, tail(true_neg))
    }
    else if (decision == truth & truth == F){
      true_neg = c(true_neg, tail(true_neg)+1)
      true_pos = c(true_pos, tail(true_pos))
    }
    else{
      true_neg = c(true_neg, tail(true_neg))
      true_pos = c(true_pos, tail(true_pos))
    }
    
  }#for (rw in 1:nrow(my_table_test)){
}

true_neg = true_neg/max(true_neg)
true_pos = true_pos/max(true_pos)

plot(true_pos ~ true_neg,
    type = "l")


