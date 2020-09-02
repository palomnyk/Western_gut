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

my_datasets = list(asv_table, otu_tab, philr_trans, cbind(otu_ratio, philr_trans), cbind(otu_tab, otu_ratio))

my_table_train = my_table[train_index,]
my_table_test = my_table[test_index,]
##--------------------------Build model ----------------------------##
rf <- randomForest(
  my_table_train, resp_var_train
)
##-------------------------Test model------------------------------##

roc_axes = function(groups = unique(resp_var_test),
                    test_data, 
                    true_resp = resp_var_test,
                    ml_model = rf){
  true_pos = c(0)
  false_pos = c(0)
  true_neg = c(0)
  for (grp in groups){
    for (rw in 1:nrow(test_data)){
      pred = predict(ml_model, newdata=test_data[rw,])
      decision = pred == grp
      truth = true_resp[rw] == grp
      if (decision == truth & truth == T){
        true_pos = c(true_pos, tail(true_pos)+1)
        true_neg = c(true_neg, tail(true_neg))
        false_pos = c(false_pos, tail(false_pos))
      }
      else if (decision == truth & truth == F){
        true_neg = c(true_neg, tail(true_neg)+1)
        true_pos = c(true_pos, tail(true_pos))
        false_pos = c(false_pos, tail(false_pos))
      }
      else{
        true_neg = c(true_neg, tail(true_neg))
        true_pos = c(true_pos, tail(true_pos))
        false_pos = c(false_pos, tail(false_pos)+1)
      }
    }#for (rw in 1:nrow(test_data)){
  }#for (grp in groups){
  true_neg = true_neg/max(true_neg)
  true_pos = true_pos/max(true_pos)
  false_pos = false_pos/max(false_pos)
  return(list(true_pos, false_pos))
}#end roc_data
my_t_roc = roc_axes(test_data = my_table_test,
                    true_resp = resp_var_test,
                    )


plot(my_t_roc[[1]] ~ my_t_roc[[2]],
    type = "l",
    xlab = "True positives",
    ylab = "False positives")


