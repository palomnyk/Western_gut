#file for combining patient metadata to sequencing metadata.

myMeta = read.table(file.path('metadata.txt'), 
                    sep="\t", 
                    header=TRUE,
                    # row.names = "run_accession",
                    check.names = FALSE,
                    stringsAsFactors=FALSE)
sampleData = read.table(file.path('PRJEB28687.txt'), 
                        sep="\t", 
                        header=TRUE, 
                        # row.names = "run_accession",
                        check.names = FALSE,
                        stringsAsFactors=FALSE)

combinedMetadata = merge(myMeta, sampleData, by = "SampleID")

write.table(combinedMetadata, 
            file = "fullMetadata.tsv", 
            sep = "\t",
            row.names = FALSE)