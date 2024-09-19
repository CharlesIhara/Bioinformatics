# Gene Set Enrichment Analysis (GSEA)
GSEA is a method to determine whether a pre-defined set of genes shows statistically significant differences between two biological states (i.e. whether a set of genes is enriched in one sample compared to another).

## Disclaimer
This is not intended to be a complete documentation. This document is for personal reference and records.

## Inputs
- Expression Data Set - normalized counts data obtained from running DESeq2 normalization on raw counts data
- Phenotype Labels - a text file mapping each sample to a particular phenotype/condition e.g. "CTRL" vs "MUT")
- Gene Sets - pre-defined sets of genes associated with particular biological pathways 

## Outputs
- List of gene sets that are enriched in one sample relative to another
- HTML report containing enrichment score plots and other metrics

## 1. Expression Dataset
An expression dataset is a tab-delimted text file that contains a row for each feature/gene and columns for each sample. 
The value of the cell at row i, column j is the "expression level of gene i in sample j".

### File types
There are several accepted file types for the expression dataset (.GCT, .RES, .PCL, and .TXT). We will use the .TXT (tab-delimited text) file type. 

To prepare the expression dataset file, open it in excel, modify to match the following format, and then save as a tab-delimited text file:

![alt text](https://software.broadinstitute.org/cancer/software/gsea/wiki/images/4/4c/Txt_format_snapshot.gif)

## 2. Phenotype Labels

The Phenotype Label file, also called a class file or template file, defines phenotype labels and assigns them to each sample in the expression dataset.
These labels can be either categorical or continuous. 
Categorical labels are used to define discrete phenotypes, and you need at least two. These enrichment of gene sets will be compared between phenotypes in a pairwise manner (ex: tumor cell expression vs healthy cell expression).
Continuous labels are used to define a phenotype profile for a time series experiment, or to find gene sets that correlate with a gene of interest. 

We will use categorical labels. Even though we have time-series data (0h, 12h, 48h, 72h, 96h, etc), we can perform pairwise comparisons between timepoints by treating each timepoint as a discrete category.
To prepare the categorical phenotype label file, open a new text file in textEdit, VS Code, or any other text editor, and format it according to the following template:

![alt text](https://software.broadinstitute.org/cancer/software/gsea/wiki/images/6/6d/Cls_format_snapshot.png)

Then, save this as a tab-delimted text file with the .cls extension.

## 3. Installing and Running GSEA
If needed, install the GSEA software from here: https://www.gsea-msigdb.org/gsea/index.jsp

### Load Data
- Once GSEA is installed, open the application, and click 'load data' icon.
- Find and upload your Expression Dataset and Phenotype Labels Files.

### Run GSEA
You can leave most of the settings set to their default values. The important ones to change have been outlined below:

#### Required Settings
- Select your Expression Dataset file from the dropdown menu (it should show your uploaded files)
- Select a gene set you would like to use from the 'Gene sets database' (e.g. Hallmark, KEGG)
- Number of permutations: 1000
- Select your Phenotype Labels file from the dropdown menu

#### Basic Fields
- Analysis Name: give the output a name based on the comparison you're running (e.g. 0_vs_12)
- Metric for ranking genes: Signal2Noise
- Save results in this folder: set the output folder path where you want to save the results

#### Advanced Fields
- Number of Markers: you may need to adjust this if you encounter an error (too many or too few)

After adjusting these settings:
- Click 'run' to initiate the analysis. The analysis name and status should appear in the bottom left window titled "GSEA reports Processes".
- After a few minuts, if/when the analysis is successfully completed, the status should change to 'success'.
- Click on 'success' to view the GSEA report in the browser (html file).
- You may need to rerun the analysis if it's unsuccessful. Click "error" to see what went wrong.

## 4. GSEA Algorithm
1. Rank all the genes in the current sample based on how strongly they correlate with the sample's phenotype
   - Example: samples 1-3 are of "WT" phenotype and samples 4-6 are of "MUT" phenotype
   - GSEA will use samples 1-3 to create a ranking of genes based on how strongly they correlate with the "WT" phenotype
   - GSEA will use samples 4-6 to create a ranking of genes based on how strongly they correlate with the "MUT" phenotype
   - A suitable metric such as "fold-change" is used for this ranking

2. Walk down the ranked gene list corresponding to the current sample's phenotype, starting with the most highly expressed gene
   - If the current gene is in the gene set, then the running enrichment score is increased.
   - If the current gene is not in the gene set, then the running enrichment score is decreased.

3. The Enrichment Score Plot shows how the Enrichment Score changes as we walk down the list of ranked genes

4. Calculate "Enrichment Score" as the maximum deviation from zero (absolute min or max)

5. Perform a test to determine whether the Enrichment Score is significant
   - The phenotype labels are randomly permuted among the samples and "random" Enrichment Scores are calculated
   - The real Enrichment Scores are compared to the random ones to generate a p-value

6. Normalize Enrichment Score to account for gene size

7. Adjust p-value is to account for multiple testing

## 5. GSEA Output Statistics and Interpretation

### Enrichment Score (ES)
- The enrichment score is the absolute maximum or minimum value of the running enrichment score obtained while walking down the list of ranked genes.
- For a gene set that is upregulated in the sample, genes from the gene set will likely fall among the top/highly ranked genes. This will lead to a positive enrichment score (peak) towards the "high rank" end the gene list.
- Example:
[insert image]

- For a gene set that is downregulated in the sample, genes from the gene set will likely fall among the bottom/lowly ranked genes. This will lead to a negative enrichment score (trough) towards the "low rank" end of the gene list.
- Example:
[insert image]

### Leading Edge Subset
- The leading edge subset is the subset of genes from the gene set that contribute most to the enrichment score.
- For a sample with a positive enrichment score, these will be highly expressed genes from the gene set which hike up the Enrichment Score, leading to a peak.
- Example:
[insert image]

- For a sample with a negative enrichment score, these will be the lowly expressed genes from the gene set that hike up the Enrichment Score after a "trough", forming the rising edge of a minimum.
- Example:
[insert image]

### Ranking Metric (bottom of Enrichment Plot)
- The Ranking Metric tells us the correlation between a gene from the ranked gene list and the phenotype of the sample.
- Highly expressed and therefore highly ranked genes towards the left end of the plot are positively correlated with the sample's phenotype.
- Lowly expressed and therefore lowly ranked genes towards the right end of the plot are negatively correlated with the sample's phenotype.
- Example: gene ABCD1 is highly expressed in a sample with a hairy phenotype. Thus, it will be highly ranked in the gene list and positively correlated with the hairy phenotype.

### Normalized Enrichment Score (NES)
- NES - Enrichment Score that has been normalized to account for gene set size
- Without normalization, larger gene sets would be more likely to have higher enrichment scores, because there are simply more genes that could be enriched and hike up the ES
- Normalization accounts for gene set size, allowing us to compare enrichment between gene sets
- NOTE: only the NES (not the ES) should be used for comparing gene sets

### p-value
- the p-value is the probability that the gene set was identified as "enriched" due to chance, rather than true overexpression
- The most significant gene sets are those with p-value < 0.05
- NOTE: if analyzing multiple gene sets, use the adjusted p-value or q-value instead 

### Adjusted p-value
- GSEA adjusts the p-value using the Benjamini-Hochberg (BH) Procedure
- Essentially, this inflates the p-values to account for multiple testing and gene set size
- Multiple testing is the idea that for a large number of gene sets, the number of falsely identified "enriched" gene sets will become high even if the false discovery rate itself is relatively small
- Thus, lowering the p-values corrects for underestimations in the likelihood of a false positive

### False Discovery Rate (FDR) aka q-value
- The q-value is the probability of a false positive (gene set being identified as "enriched" by chance, rather than true overexpression)
- This is another type of adjusted p-value and is less stringent than the (BH)-corrected p-value for filtering out insignificant gene sets
- The most significant gene sets are those with FDR (q-value) less than 0.25.

## Example Data

### Experimental Background
- The following analysis was performed on normalized counts data from an RNA-sequencing experiment involving a tumor organoid 
- The tumor organoid cells lack a specific protein called ARID1A, which is part of a chromatin remodeler that is important for tumor supression
- In the experiment, ARID1A expression was induced at 0h, and RNA expression levels were measured at subsequent time points (12h, 24h, 48h, 72h, 96h)
- Thus, the 12h data represents the RNA expression levels 12 hours after ARID1A addback in perviously ARID1A-deficient tumor cells
the 0h time point serves as a "control" or baseline for comparison, because this represents zero ARID1A expression (tumor cells just prior to ARID1A addback)

### GSEA Enrichment Score Plots


### Intrepretation



## References
1. https://www.gsea-msigdb.org/gsea/doc/subramanian_tamayo_gsea_pnas.pdf
2. https://www.youtube.com/watch?v=KY6SS4vRchY
3. https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
4. https://www.youtube.com/watch?v=Yi4d7JIlAsM
