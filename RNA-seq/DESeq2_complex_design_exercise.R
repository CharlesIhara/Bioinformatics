# DESeq2 Exercise: design formula with complex variables 

# this was written following the tutorial from HBC: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
# a complex design formula can be used to study the effect of multiple known sources of variation on differential gene expression between samples

# load data
data <- "read.table() OR read.tsv(sep="\t") OR read.csv(sep=",")"

# STEP 1: design formula / metadata - this file tells DESeq2 the known sources of variation in the experiment
meta <- data.frame(
  sex=c("M", "M", "M", "M", "F", "F", "F", "F"),
  age=c(11, 13, 11, 13, 11, 13, 11, 13),
  litter=c(1, 2, 1, 1, 1, 1, 1, 2),
  treatment=c("ctrl", "ctrl", "treat", "treat", "ctrl", "ctrl", "treat", "treat")
)
row.names(meta)=c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7", "sample8")

# STEP 2: set "design" variable using "~ major sources of variation + variable of interest"
# variable of interest always comes last
# NOTE: any factors listed in the design formula must match a column name from the meta data frame
design <- ~ sex + age + treatment

# we could also explore the effect of one factor on another by creating a new factor in our meta data frame
# example: how does sex affect the treatment effect?

# create the sex_treatment column
sex_treatment <- c()
for (i in 1:nrow(meta)) {
  insertIndex <- length(sex_treatment)+1
  sex <- meta[["sex"]][i]
  treatment <- meta[["treatment"]][i]
  sex_treatment[insertIndex] <- paste0(sex, "_", treatment)
}
# add the sex_treatment column to meta data frame
meta <- cbind(meta, sex_treatment)
# create a new design formula (exclude variables that are confounded with the complex factor, in this case "sex" and "treatment")
design <- ~ age + sex_treatment

# STEP 3: create DESeqDataSet and run DESeq
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = design)
dds <- DEseq(dds)

