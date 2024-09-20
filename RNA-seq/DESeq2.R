# DESeq2 DGE Analysis

# this analysis is performed using a dataset from NCBI GEO with the following accession ID: GSE171663
# link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171663
# this script was written using the following tutorials:
# https://hbctraining.github.io/DGE_workshop/schedule/1.5-day.html
# https://www.youtube.com/watch?v=Ul-9s8YOOSk

################################################################################
# PART 0: Setup

# install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("pheatmap")
BiocManager::install("apeglm")
install.packages("ashr")
install.packages("ggrepel")
BiocManager::install("DEGreport")

# load packages 
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(apeglm)
library(ashr)
library(dplyr)
library(tibble)
library(tidyverse)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)

# load data
data <- read.delim("GSE171663_normalized_counts_geneSymbols.txt", header=TRUE, row.names=1, sep="\t")
base <- "GSE171663"
colnames(data) <- c("Ctr_s2", "Dox1d_s4", "Dox21d_s5", "RS_s6", "Ctr_s7", "Dox1d_s10", "Dox21d_s11", "RS_s12", 
                    "Ctr_s13", "Dox1d_s14", "Dox21d_s15", "RS_s16", "Ctr_s1", "Dox1d_s3", "Dox21d_s8", "RS_s9")
# filter out genes with less than 50 total counts across all samples
data <- data[which(rowSums(data) > 50),]

# reorder the columns
data <- data[,c("Ctr_s1", "Ctr_s2", "Ctr_s7", "Ctr_s13", 
                "RS_s6", "RS_s9", "RS_s12", "RS_s16", 
                "Dox1d_s4", "Dox1d_s14", "Dox21d_s5", "Dox21d_s15")]

# create the meta data
condition <- factor(c("C", "C", "C", "C", "R", "R", "R", "R", "D", "D", "D", "D"))
meta <- data.frame(condition, row.names=colnames(data))

# check that the meta data rows match the columns of the data matrix
all(colnames(data) %in% rownames(meta))
all(colnames(data)==rownames(meta))

# variable for the number of samples (used later to gather all the samples/columns)
num_samples <- 12
################################################################################

################################################################################
# PART 1: Count Normalization

# create a DEseqDataSet object 
dds <- DESeqDataSetFromMatrix(countData=data, colData=meta, design=~condition)

# calculate size factors (median of ratios)
dds <- estimateSizeFactors(dds)

# scale counts using size factors
normalized_counts <- counts(dds, normalized=TRUE)

# save the normalized counts to a txt file
filename <- paste0(base, "_normalized_counts.txt")
write.table(normalized_counts, file=filename, sep="\t", quote=F, col.names=NA)
################################################################################

################################################################################
# PART 2: Quality Control

# transform the counts to log2 scale (reduces variation across samples for low count genes and normalizes for gene length, allowing for better clustering)
# blind=TRUE ensures that the samples are compared without considering experimental design or other information 
# rld is a DESeqTransform object, which not only contains the matrix of scaled count data, but also the size factors
rld <- rlog(dds, blind=TRUE)

# STEP 1: PCA plot (see how if samples cluster based on treatment, or another confounding variable)
# this function uses ggplot2 package under the hood
# the 'intgroup' specifies what variable from the meta data to use (in this case, sampletype could mean: treated vs untreated, mutant vs WT, 0h vs 12h vs 48h vs 72h)
# this returns the two largest principal components (pc1 and pc2)
factor_of_interest <- "condition"
plotPCA(rld, intgroup=factor_of_interest)
plot_name <- paste0(factor_of_interest, "_PC1_PC2_plot.png")
ggsave(filename=plot_name)

# to return the other principal components, we need to use additional functions
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
# create a data frame with meta data and PC's
df <- cbind(meta, pca$x)
# manually  plot PC3 and PC4
ggplot(data=df, mapping=aes(x=PC3, y=PC4, color=!!sym(factor_of_interest))) + geom_point()
plot_name <- paste0(factor_of_interest, "_PC3_PC4_plot.png")
ggsave(filename=plot_name)


# STEP 2: Hierarchical clustering heatmap
# extract the rlog matrix from the DESeqTransform object to use for heatmap generation
rld_mat <- assay(rld)
# compute the pairwise correlation values using cor()
rld_cor <- cor(rld_mat)
# plot the correlations in a heatmap
plot_name <- "QC_heatmap.png"
pheatmap(rld_cor, filename=plot_name, main="QC Heatmap")
################################################################################

################################################################################
# PART 3: (OPTIONAL) Re-create DESeqDataSet with New Design Formula if QC has Revealed Additional Sources of Variation

#dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sex + age + condition)
################################################################################

################################################################################
# PART 4: Running DESeq2 

dds <- DESeq(dds)

# Explanation of steps performed by this function:
# STEP 1: Estimate Size Factors (same as count normalization using median of ratios)
# look at size factors (samples with size factors > 1 will be scaled down; size factors < 1 will be scaled up)
sizeFactors(dds)
# total number of raw counts per sample
colSums(data)
# total number of normalized counts per sample (very similiar across samples)
colSums(normalized_counts)

# STEP 2: Estimate Dispersions
# dispersion is a measure of variability
# dispersion is directly related to variance
# dispersion is inversely related to mean counts (higher counts = less dispersion; lower counts = higher dispersion)
# dispersion estimates reflect the variance in gene expression for a given mean counts value
# a plot of variance vs mean counts shows that the correlation between these variables increases as mean counts increases
# this means that variance is well predicted by the mean counts for genes with higher average counts
# thus, the most likely estimate of dispersion for each gene is calculated based on its mean counts (across replicates)
# DESeq2 assumes that genes with similar mean counts have similar dispersion

# STEP 3: Fit a Curve to the Plot of Gene-wise Dispersion Estimates vs Average Counts
# the curve gives an expected dispersion value based on the average counts for a particular gene

# STEP 4: Shrink Gene-wise Dispersion Estimates to Values Predicted by the Curve
# genes that are further from the curve and have lower average counts will be shrunken more than points closer to the curve and with higher average counts
# genes with low dispersion estimates (below the curve) or slightly high dispersion estimates (slightly above the curve) are shrunken towards the curve
# genes with extremely high dispersion (very high above the curve) are not shrunken
# these genes likely have very high biological variation, so shrinking these to the curve could result in false positives
# plot of dispersions:
# in general, we expect higher dispersion for genes with lower average counts, and less dispersion for higher count genes
plot_name <- "dispersion_vs_mean_normalized_counts.png"
png(filename=plot_name)
plotDispEsts(dds)
dev.off()

# STEP 5: Fit the Counts Data to a Generalized Linear Model
# DESeq2 uses a negative binomial distribution to model counts data
# sequencing data exhibits over-dispersion, so the variance is not equal to the mean. This means we can't use Poisson
# negative binomial takes in a dispersion parameter, which is calculated in the previous step

# STEP 6: Shrink LFC Values
# LFC values are shrunken towards 0 if a gene has low counts or high dispersion
# NOTE: shrinking LFC values will not change whether genes are identified as differentially expressed
# instead, the shrinkage is used for downstream analysis, such as subsetting the differentially expressed genes based on shrunken LFC, or GSEA

# STEP 7: Hypothesis Testing Using Wald Test
# null hypothesis: there is no differential expression between groups X and Y
# Wald Test - commonly used for comparing two groups
# A Wald Test Statistic is computed along with a p-value (probability that the value of the statistic is observed due to chance)
# if the p-value is small, then it's unlikely that the observed difference is due to chance alone, so we reject the null hypothesis (there is differential expression between groups X and Y)
################################################################################

################################################################################
# PART 4: View Results

# set up variables to easily test different comparisons
factor_of_interest <- "condition"
levels_to_compare <- c("R", "C")
comparison <- paste0(levels_to_compare[1], "vs", levels_to_compare[2])
column_range <- 2:9 # variable for the range of columns being used in the comparison 

# STEP 1: use the results() function along with several parameters to build a results table
# the 'contrast' parameter is a vector used to specify a factor and the levels of that factor we would like to compare
# example: contrast <- c("condition", "WT", "CTR")
# NOTE: the last sample is always used as the "baseline" when reporting LFC values, so a negative LFC would mean lower expression in WT compared to CTR
# the 'alpha' parameter specifies the cutoff for the adjusted p-value
contrast <- c(factor_of_interest, levels_to_compare[1], levels_to_compare[2]) # compare Replicative Senesence vs Control cells
results_unshrunken <- results(dds, contrast=contrast, alpha=0.05)
# shrink LFC values towards zero if the mean count for that gene is low or dispersion is high
# NOTE: shrinking LFC values does not change the number of DE genes returned
results_shrunken <- lfcShrink(dds, contrast=contrast, res=results_unshrunken, type="ashr")

# STEP 2: MA plot
# this is a scatterplot that displays LFC values (y-axis) vs the mean counts (x-axis) for each gene
# genes identified as differentially expressed are colored as a separate series
# in general, we expect to see differentially expressed genes across all expression levels (high and low)
# MA plot with unshrunken LFC values:
plot_name <- paste0(factor_of_interest, "_", comparison, "_unshrunken_MAplot.png")
png(filename=plot_name)
plotMA(results_unshrunken, ylim=c(-2,2))
dev.off()
# MA plot with shrunken LFC values:
plot_name <- paste0(factor_of_interest, "_", comparison, "_shrunken_MAplot.png")
png(filename=plot_name)
plotMA(results_shrunken, ylim=c(-2,2))
dev.off()

# STEP 3: Multiple Test Correction
# with a large data set, even a fairly low false positive rate will result in a high number of false positives
# to address this issue, we filter the genes by their 'adjusted p-value' or 'false discovery rate' rather than basic 'p-value'

# STEP 4: Additional Filtering
# in addition to adjusted p-values, we may want to filter by LFC for additional stringency
padj.cutoff <- 0.05
lfc.cutoff <- 0.58 # this corresponds to a log2FC of 1.5
# convert the results table into a tibble
results_tibble <- results_shrunken %>% data.frame() %>% rownames_to_column("gene") %>% as_tibble()
# filter out genes with NA p-values
results_NA_filtered <- na.omit(results_tibble)
# subset the tibble by adjusted p-value and LFC cutoffs
results_filtered <- results_NA_filtered %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

# Alternative Filtering Approach
# we could also set a lfc threshold in the results() function directly, as we did for the adjusted p-value
# this performs a two-tailed statistical test against the absoulte value of the lfc threshold 
# this approach is more conservative and produces a much smaller set of differentially expressed genes (less false positives)
# the initial results table has the same number of DE genes (18,019), but with much higher p-values, so we will filter more of them out with the same p-value cutoffs
results_alt_filtered <- results(dds, contrast=contrast, alpha=0.05, lfcThreshold=0.58)
# now apply the same cutoffs as before
results_alt_filtered <- results_alt_filtered %>% 
                            data.frame() %>% 
                            rownames_to_column("gene") %>% 
                            as_tibble() %>% 
                            na.omit() %>% 
                            filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
# observe that this alternative filtering produces only 83 DE genes, compared to 697 found with the first method

# STEP 5: Summarize results
summary(results_shrunken)
################################################################################

################################################################################
# PART 5: Visualize Differentially Expressed Genes

# convert meta and normalized counts data into tibbles (a data structure that is used with tidyverse functions)
meta_tibble <- meta %>% rownames_to_column("samplename") %>% as_tibble()
normalized_counts_tibble <- normalized_counts %>% data.frame() %>% rownames_to_column("gene") %>% as_tibble()

# STEP 1: Scatter Plot
# plotting counts for one gene:
plotCounts(dds, gene="79915", intgroup="condition")

# we can also return the data frame from plotCounts() and use it to produce a more complex plot with ggplot2
plotCounts_df <- plotCounts(dds, gene="79915", intgroup="condition", returnData=TRUE)
ggplot(data=plotCounts_df, mapping=aes(x=condition, y=count, color=condition)) + 
  geom_point(position=position_jitter(width=0.1, height=0)) + 
  geom_text_repel(aes(label=rownames(plotCounts_df))) + 
  theme_bw() +
  ggtitle("79915") +
  theme(plot.title=element_text(hjust=0.5))
  
# we can also plot the normalized counts of the top 20 DE genes
# first identify the top 20 DE genes (20 lowest adjusted p-values) from the filtered results (LFC > 0.58, padj < 0.05)
top20_genes <- results_filtered %>% 
                  arrange(padj) %>% 
                  pull(gene) %>% 
                  head(n=20)
  
# get the normalized counts of these top 20 genes
# the expression will return true if a given row/gene is in the list of top 20 genes
# the filter() function will keep rows for which the expression is true
top20_normalized_counts_tibble <- normalized_counts_tibble %>% filter(gene %in% top20_genes)

# subset the relevant samples/columns (plus the gene column)
top20_gathered_tibble <- top20_normalized_counts_tibble[c(1,column_range)] 

# gather the samples into one column; specify a range of columns (could also be column names)
top20_gathered_tibble <- gather(data=top20_gathered_tibble, 2:ncol(top20_gathered_tibble), key="samplename", value="normalized_counts")

# merge the counts data for top 20 genes in all samples with the meta data 
# this adds a 'condition' column, allowing us to color data points based on condition
top20_gathered_tibble <- inner_join(top20_gathered_tibble, meta_tibble)

# plot using ggplot
# scale the counts using log to reduce the focus on outliers
title <- paste0("Counts of Top 20 DE Genes Between ", levels_to_compare[1]," and ", levels_to_compare[2])
ggplot(data=top20_gathered_tibble, mapping=aes(x=gene, y=normalized_counts, color=condition)) +
  geom_point() +
  scale_y_log10() +
  xlab("gene") + 
  ylab("log10(normalized counts)") +
  ggtitle(title) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(plot.title=element_text(hjust=0.5))
plot_name <- paste0(factor_of_interest, "_", comparison, "_top20_DEG_counts_scatterplot.png")
ggsave(filename=plot_name)

# STEP 2: Heat Map
# this is useful to visualize expression of many genes at once

# obtain a data frame of normalized counts for the significant genes (filtered by LFC and padj):
# 1. extract the relevant samples/columns from the normalized counts tibble (plus the genes column)
# 2. filter out rows/genes that are not in the filtered results table 
# 3. convert the tibble back into a data frame 
normalized_counts_sig <- normalized_counts_tibble[,c(1,column_range)] %>% 
                            filter(gene %in% results_filtered[["gene"]]) %>%
                            data.frame() %>%
                            column_to_rownames("gene")

# create data frame of annotation info
# NOTE: the "select" is only necessary if the meta data contains other factors besides condition (e.g. sex, age, sequencing group)
annotation <- meta_tibble %>% select(samplename, condition) %>% data.frame(row.names="samplename")

# set a color palette
heatmap_colors <- brewer.pal(n=6, name="YlOrRd")

# generate the heatmap
# scale=row computes a z-score for each count value (relative to the mean for that gene); this improves the color visualization
plot_name <- paste0(factor_of_interest, "_", comparison, "_heatmap.png")
title <- paste0("Significant (padj < 0.05, LFC > 0.58) DE Genes Between ", levels_to_compare[1], " and ", levels_to_compare[2])
pheatmap(mat=normalized_counts_sig,
         color=heatmap_colors,
         cluster_rows=T,
         cluster_cols=T,
         show_rownames=F,
         annotation=annotation,
         border_color=NA,
         fontsize=10,
         scale="row",
         fontsize_row = 10,
         height=20,
         main=title,
         filename=plot_name)

# STEP 3: Volcano Plot
# this is a scatter plot of log transformed padj values (y-axis) vs lFC values (x-axis)
# for this plot, we want to use the unfiltered results tibble because weant counts of both significant and insignificant DEGs
# we will color the significant DEGs as a different series to see how they compare to insignificant DEGs

# first, add a new column to our results tibble, which tells us if a gene is significant or not (does it meet the LFC and p-value cutoffs)
results_tibble_significant <- results_NA_filtered %>% mutate(significant=(padj<padj.cutoff & abs(log2FoldChange) > lfc.cutoff))
# alternatively, we could just check if that row is in the filtered results table
results_tibble_significant2 <- results_NA_filtered %>% mutate(significant=(gene %in% results_filtered$gene))
# we can verify that the two vectors are equivalent
identical(results_tibble_significant[["significant"]], results_tibble_significant2[["significant"]])

# create the volcano plot
plot_name <- paste0(factor_of_interest, "_", comparison, "_volcanoplot.png")
title <-paste0("Significant (padj < 0.05, LFC > 0.58) DE Genes Between ", levels_to_compare[1], " and ", levels_to_compare[2])
xlab <- paste0("log2(fold change) ", "(",levels_to_compare[1], "/", levels_to_compare[2], ")")
ggplot(data=results_tibble_significant, mapping=aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point() +
  ggtitle(title) +
  xlab(xlab) +
  ylab("log10(padj)") +
  theme(plot.title=element_text(size=rel(1.5), hjust=0.5),
        axis.title=element_text(size=rel(1.25)),
        plot.margin=margin(t=20, r=10, l=10, b=10, unit="pt"))
ggsave(filename=plot_name, scale=2.5)

# we can also label certain genes on the plot
# first order the results tibble of significant genes by padj value
results_tibble_significant_with_labels <- results_tibble_significant %>% arrange(padj) 
# then add a column for gene labels (initially all empty strings) 
results_tibble_significant_with_labels <- results_tibble_significant_with_labels %>% mutate(label="")
# populate gene label column for the top X genes only
x <- 10
results_tibble_significant_with_labels[1:x, "label"] <- results_tibble_significant_with_labels[1:x, "gene"]
# now create the plot
plot_name <- paste0(factor_of_interest, "_", comparison, "_volcanoplot_top10_labeled.png")
title <-paste0("Significant (padj < 0.05, LFC > 0.58) DE Genes Between ", levels_to_compare[1], " and ", levels_to_compare[2])
xlab <- paste0("log2(fold change) ", "(",levels_to_compare[1], "/", levels_to_compare[2], ")")
ggplot(data=results_tibble_significant_with_labels, mapping=aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(mapping=aes(color=significant)) +
  geom_text_repel(mapping=aes(label=label)) +
  ggtitle(title) +
  xlab(xlab) +
  ylab("log10(padj)") +
  theme(plot.title=element_text(size=rel(1.5), hjust=0.5),
        axis.title=element_text(size=rel(1.25)))
ggsave(filename=plot_name, scale=2.5)  
################################################################################

################################################################################
# PART 6 (OPTIONAL): Hypothesis Testing with Likelihood Ratio Test




