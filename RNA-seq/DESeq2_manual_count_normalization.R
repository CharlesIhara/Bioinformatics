# DESeq2


# Exercise 1: normalize the counts for a single gene, given the raw counts and size factors for each sample
# raw counts
PD1 <- c(21, 58, 17, 97, 83, 10)

# set the sample names for the vector values
names(PD1) <- c(paste0("sample", 1:6))

# convert into a data frame
PD1 <- data.frame(PD1)

# tranform/transpose the data frame so that genes = rows, samples = cols
PD1 <- t(PD1)

# size factors for each sample
size_factors <- c(1.32, 0.70, 1.04, 1.27, 1.11, 0.85)

# normalize the counts by dividing each sample's counts in the data frame by the corresponding size factor
PD1_normalized <- PD1/size_factors


# Exercise 2: manual calculation of size factors and normalization of counts using DESeq2's method with  dummy data

# STEP 1: create the data frame of counts
genes <- c("EF2A", "ABCD1", "MEFV1", "BAG1", "MOV10")
sample1 <- c(1489, 22, 793, 76, 521)
sample2 <- c(906, 13, 410, 42, 1196)
counts_df <- data.frame(sample1=sample1, sample2=sample2)
rownames(counts_df) <- genes

# STEP 2: create a pseudo-reference (row-wise mean)
# option A: use the logarithmic mean
counts_df$LogarithmicMean <- apply(counts_df, 1, function(x) exp(mean(log(x))))

# option B: use the geometric mean
geometricMean <- function(x) {
  n <- length(x) # length of the input row aka number of columns
  product <- prod(x) # product of all elements in the row vector
  nth_root <- product^(1/n) # taking the nth root = raising to the 1/nth power
}
counts_df$GeometricMean <- apply(counts_df, 1, geometricMean)

# STEP 3: calculate the ratio of each sample to the reference
samples <- c("sample1", "sample2")
for (sample in samples) {
  new_column_name <- paste0(sample, "_ratio")
  
  # option A: directly compute the new column by dividing columns
  counts_df[new_column_name] <- counts_df[sample]/counts_df$GeometricMean
  
  # option B: use the apply() function
  #counts_df[,new_column_name] <- apply(counts_df[,samples], 2, function(x) x/counts_df$GeometricMean)
}

# STEP 4: Calculate the normalization factor (median of ratios)




