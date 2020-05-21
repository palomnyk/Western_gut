# Author: Aaron Yerke
# This is a pipeline that was created to explore the philr transform.
# Many resources were used to create this pipeline:
#   Processing the sequences through DADA2, making trees with phangorn, and putting them into phyloseq objects:
#     This is by far the biggest source:
#     https://github.com/spholmes/F1000_workflow/blob/master/MicrobiomeWorkflow/MicrobiomeWorkflowII.Rmd
#     or (in article format)
#     https://f1000research.com/articles/5-1492/v2
#     
#   DADA2:
#     https://benjjneb.github.io/dada2/tutorial.html
#   
#   Using philr:
#     https://bioconductor.org/packages/release/bioc/html/philr.html
#     Code example:
#       https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.R

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("phyloseq")
  BiocManager::install("philr")
  BiocManager::install("ape")
}

library("phyloseq")
library("ape")
library(philr); packageVersion("philr")

##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
setwd(file.path(home_dir))

##------------------Import R objects and data-----------------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects", "phyloseq_obj.rds"))
ps = readRDS(con)

# con <- gzfile(file.path( "philr_pipelines","phyloseq_obj.rds"))
# ps = readRDS(con)

##------------Import R objects and data preprocessing----------------##
con <- gzfile(file.path( "philr_pipelines", "r_objects","philr_transform.rds"))
ps.philr = readRDS(con)

con <- gzfile(file.path( "philr_pipelines", "r_objects","ps_philr_transform.rds"))
ps = readRDS(con)


##---------------------PCOA transform and plots---------------------##
ps.dist <- dist(ps.philr, method="euclidean")
ps.pcoa <- ordinate(ps, 'PCoA', distance=ps.dist)

catagoriesString = "Recruitment.Date	Sample.Date	Recruitment.Location	Researcher	Sub.Study	Sample.Month	Birth.Year	Age	Public.Housing	Medical.Assistance	Children.Free.Lunch	Highest.Education	Ethnicity	Religion	Birth.Location	Type.Birth.Location	Arrival.in.US	Years.in.US	Location.before.US	Type.location.before.US	Years.lived.in.Location.before.US	Tobacco.Use	Alcohol.Use	Height	Weight	Waist	BMI	BMI.Class	Medications	Breastfed	Years.Breastfed	Notes.Participants	Age.at.Arrival	Weight.M6	BMI.M6	Waist.M6	Measurement.Date.M6	Sample.Order	Sample.Group	Waist.Height.Ratio"
catagoriesString = unlist(strsplit(catagoriesString, split = "\t"))

setwd(file.path(output_dir, "philr_pcoa"))

library("ggplot2")

theme_set(theme_classic())

#Shapes
mimicShapes = c(17, 21, 17, 16, 21, 16)

#Colors
mimicCols = c("gray", "maroon2", "palevioletred1", "red", "seagreen4", "seagreen4")
# Levels: Control Hmong1st Hmong2nd HmongThai Karen1st KarenThai
ellipsi_mimicCols = c("gray", "palevioletred1", "red", "seagreen4")
myTitle = paste("philr_PCOA_", "ethn", ".pdf", sep = "")

plt = plot_ordination(ps, ps.pcoa, color="Sample.Group", shape="Sample.Group") + 
  scale_colour_manual(values=mimicCols) + 
  scale_shape_manual(values = mimicShapes) +
  stat_ellipse(show.legend=F, type="t", level=.6) +
  scale_x_reverse() + 
  geom_point(size=3)
ggsave(plt, filename = myTitle)


#this code still works, just  commenting so it doesn't run every time I test this
# for (catag in catagoriesString){
#   print(catag)
#   catag = trimws(catag)
#   myTitle = paste("philr_PCOA_", catag, ".pdf", sep = "")
#   plt = plot_ordination(ps, ps.pcoa, color=catag)
#   ggsave(plt, filename = myTitle)
# }

print("done with PCOA")

library("vegan")

#create PCA
philr_prc = prcomp(ps.philr, 
                   center = TRUE,
                   scale = TRUE)

#extract PCA matrix and convert to dataframe
myPCA = data.frame(philr_prc$x)

myMeta = data.frame(sample_data(ps))
sa_gr = as.factor(myMeta[,"Sample.Group"])
# Levels: Control Hmong1st Hmong2nd HmongThai Karen1st KarenThai
#combine Hmong1st with Hmong2nd and Karen1st with KarenThai
comb_sa_gr = sa_gr
levels(comb_sa_gr) = c("Control", "Hmong1st", "Hmong1st", "HmongThai", "Karen1st", "Karen1st")
            

# plot_ordination

ggplot(myPCA) +
  scale_colour_manual(values=mimicCols) +
  scale_shape_manual(values = mimicShapes) +
  stat_ellipse(show.legend=F, type="t", 
               level=.6,
               inherit.aes = F,
               aes(PC1,
                   PC2,
                   color = comb_sa_gr,
                   # group = comb_sa_gr
                   )) +
  # scale_x_reverse() + 
  # scale_fill_discrete(name = "New Legend Title") +
  geom_point(data = myPCA, size=3, show.legend=T,
             # title = "Sample Group",
             aes(PC1,PC2,
                 col = sa_gr,
                 shape = sa_gr
                 ))
# ggsave(plt, filename = "sample_group_vegan.pdf")
# plt

# ggplot(data=myPCA, aes(PC1, PC2)) + geom_point(aes(colour=Sample.Group, shape=Sample.Group, size=Sample.Group, alpha=Sample.Group)) +
#   xlab(paste0("PC1 [",percent_var[1],"%]")) +
#   ylab(paste0("PC2 [",percent_var[2],"%]")) + 
#   scale_color_manual(name="Groups", values=cols) + #sets the color palette of the fill
#   scale_alpha_manual(name="Groups", values=alphas) +
#   scale_shape_manual(name="Groups", values=shapes) +
#   scale_size_manual(name="Groups", values=sizes) + 
#   stat_ellipse(data=d, aes(colour=Sample.Group), show.legend=F, type="t", level=.6)
