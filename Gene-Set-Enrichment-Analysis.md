# Gene Set Enrichment Analysis (GSEA)
GSEA is a method to determine whether a pre-defined set of genes shows statistically significant differences between two biological states (i.e. whether a set of genes is enriched in one sample compared to another).

## Disclaimer
This is not intended to be a complete documentation. This document is for personal reference and records.

## Inputs
- Normalized counts data obtained from running DESeq2 (this is the Expression Dataset)
- Phenotype Labels

## Outputs
- list of gene sets that are enriched in one sample relative to another
- html report containing enrichment score plots

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

We will use categorical labels.
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

## 4. GSEA Statistics and Interpretation
### Enrichment Score (ES)

### Normalized Enrichment Score (NES)

### False Discovery Rate (FDR)

### Nominal P Value

## Example Data

### Experimental Background
this analysis used data that was collected from a tumor organoid; 
the organoid cells lack a specific protein called ARID1A, which is part of a chromatin remodeler that is important for tumor supression.
in the experiment, ARID1A expression was induced at 0h, and RNA expression levels were measured at subsequent time points (12h, 24h, 48h, 72h, 96h)
thus, the 12h data represents the RNA expression levels 12 hours after ARID1A addback in perviously ARID1A-deficient tumor cells
the 0h time point serves as a "control" or baseline for comparison, because this represents zero ARID1A expression (tumor cells just prior to ARID1A addback)

### GSEA Enrichment Score Plots


### Intrepretation



## References
1. https://www.youtube.com/watch?v=KY6SS4vRchY
2. https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
