# Author: Aaron Yerke
# Compare the named nodes to named OTUs

##-Functions--------------------------------------------------------##

##_Establish constants----------------------------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
# setwd(file.path(home_dir))

##-Load Dependencies------------------------------------------------##
source(file.path(home_dir, "r_libraries", "statistical_functions.R"))
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

##-Load data--------------------------------------------------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects","ps_philr_transform.rds"))
ps = readRDS(con)
close(con)

philr_trans = read.table(file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_philr_transform.csv"),
                         sep = ",",
                         header = TRUE)
# philr_trans = philr_trans - min(philr_trans)

otu_ratio = read.table(file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table.csv"),
                       sep = ",",
                       header = TRUE)

otu_tab = read.table(file.path(home_dir, "philr_pipelines", "tables", "otu_table.csv"),
                     sep = ",",
                     header = TRUE)

metadata = read.table(file.path("philr_pipelines", "tables", "ps_sample_data_reduced.csv"), 
                      sep=",", 
                      header=TRUE)
##-Run Comparison---------------------------------------------------##

philr_otu_corr <- kendall_corr(philr_trans, otu_ratio)
# philr_name <- philr_node_annot(philr_trans, ps, only_one = TRUE)

pdf(file = file.path(output_dir, "kendall_pval_histogram.pdf"))
hist(philr_otu_corr$pValues,
     breaks = 200,
     main = "Kendall Correlation Pvalue Histogram", 
     xlab = "Kendall pvalue")
dev.off()

pdf(file = file.path(output_dir, "kendall_stat_histogram.pdf"))
hist(philr_otu_corr$corr_stat,
     breaks = 120,
     main = "Kendall Correlation Statistic Histogram", 
     xlab = "Kendall Statistic")
dev.off()

pdf(file = file.path(output_dir, "kendall_est_histogram.pdf"))
hist(philr_otu_corr$est,
     breaks = 120,
     main = "Kendall Correlation Estimate Histogram", 
     xlab = "Kendall Estimate")
dev.off()

philr_otu_corr_extremes <- c(philr_otu_corr$est[philr_otu_corr$est < -0.4], philr_otu_corr$est[philr_otu_corr$est > 0.4])
pdf(file = file.path(output_dir, "kendall_est_extremes_histogram.pdf"))
hist(philr_otu_corr_extremes,
     breaks = 120,
     main = "Kendall Correlation Estimate Extremes Histogram", 
     xlab = "Kendall Estimate")
dev.off()

pdf(file = file.path(output_dir, "kendall_pval_histogram_log10.pdf"))
hist(-log10(philr_otu_corr$pValues),
     breaks = 200,
     main = "Kendall Correlation log10 Pvalue Histogram", 
     xlab = "Log10 Kendall pvalue")
dev.off()
# for 

# is this correlation table symmetrical?
length(philr_otu_corr$est[philr_otu_corr$est < -0.5]) - length(philr_otu_corr$est[philr_otu_corr$est > 0.5])

# determine how many correlations are significant at pval
l <- length(which(philr_otu_corr$pValAdj < 0.01))
n <- nrow(philr_otu_corr)

print(paste("percentage of significant pairs:", l/n))

# Determine if both ratios are inverses of each other test_df$test will be one if inverse is true
test_df = data.frame(otu_ratio$sutterellaceae.acidaminococcaceae, otu_ratio$acidaminococcaceae.sutterellaceae)

test_df$test = test_df$otu_ratio.sutterellaceae.acidaminococcaceae * test_df$otu_ratio.acidaminococcaceae.sutterellaceae

#look at the namespace of the highly correlated pairs


