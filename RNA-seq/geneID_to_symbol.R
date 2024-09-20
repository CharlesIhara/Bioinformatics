# this script converts ENSEMBL gene ID's to gene symbols

# install / load packages
library(tibble)

# set variables
counts_filename <- "GSE171663_normalized_counts.txt"
extension_length <- 4 # length of extension including the "."
output_filename <- paste0(substr(counts_filename, 1, nchar(counts_filename)-extension_length), "_geneSymbols.txt") 

# load the counts and annotation as data frames
counts_df <- read.table(counts_filename)
class(counts_df)
annotation_file <- read.delim("Human.GRCh38.p13.annot.tsv", sep="\t")
class(annotation_file)

# convert the rownames to a column called "GeneID"
counts_df <- rownames_to_column(counts_df, "GeneID")

# subset the annotation file (only need the GeneID and Symbol columns)
annotation_file <- annotation_file[,c("GeneID", "Symbol")]

# merge the counts and annotation data frames using "GeneID" column
# ensure that all rows in x are kept, even if there's no match in y
counts_df <- merge(x=counts_df, y=annotation_file, by="GeneID", all.x=TRUE)

# replace the first column (GeneID) with the last column (Symbol)
counts_df$GeneID <- counts_df$Symbol

# remove last column (Symbol)
counts_df$Symbol <- NULL

# save the data table to a new counts file
write.table(counts_df, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE)
  
  
  