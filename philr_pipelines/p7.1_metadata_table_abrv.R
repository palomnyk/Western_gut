# Author: Aaron Yerke
# Script for making metadata table
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##

##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
setwd(file.path(home_dir))

##---------------------------Functions------------------------------##


##------------Import tables    and data preprocessing---------------#
metadata = read.table(file.path("philr_pipelines", "tables", "ps_sample_data.csv"), 
                      sep=",", 
                      header=TRUE)
drop = c("SampleID","Sample.Date", "Subject.ID","Old.Participant.ID","sample_accession",
         "secondary_sample_accession","experiment_accession",             
         "tax_id","scientific_name","instrument_model","library_layout",
         "experiment_title","sample_title", "study_title", "run_alias",
         "fastq_ftp" ,"fastq_galaxy","submitted_ftp","submitted_galaxy","sra_ftp", 
         "sra_galaxy", "study_accession",
         "Notes.Samples", "Sub.Study", "Notes.Participants", "Birth.Year")
metadata[metadata==""] <- NA
for (c in 1:ncol(metadata)){
  if (any(is.na(metadata[,c]))){
    drop = c(drop, colnames(metadata)[c])
  }
}
drop = unique(drop)
metadata = metadata[ , !(names(metadata) %in% drop)]
# lapply(metadata, class)
write.table(metadata, 
            file = file.path("philr_pipelines", "tables", "ps_sample_data_reduced.csv"), sep = ",")