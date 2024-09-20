#!/bin/bash

# DISCLAIMER: I did not write the code for this script. I just added comments.
# Alignment pipeline for RNA data. 
# Modified from /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/rna_alignment.sh
# Assumes that all fastq files are contained in a fastq folder, and have suffix _R1.fastq.gz. See areas marked "CHANGE FOR FILE EXTENSION"
# Requires the "alignment" and "R" environments, and the "qc" environment (for fastqc and multiqc)

# TODO: Ask/find out what the read length was
# TODO: Check to see if fastqc is recursive (Important!)

# 1. Setup
# generate required folders:
mkdir -p counts
mkdir -p clumped
mkdir -p trimmed
mkdir -p aligned
mkdir -p bigwig
mkdir -p results
mkdir -p PBS
mkdir -p log
# save location of current folder
folder=$(cd "$(dirname "$0")";pwd)  
# set the suffix name for the input fastq files
suffix1=_R1_001.fastq.gz

# 2. Find the number of fastq files with the specified suffix
# this command 'finds' files inside the fastq directory that have the specified suffix
# for each of these files, an 'x' is printed
# the resulting string of 'x's is piped into the wc (word count) command, which counts the number of characters when the -c flag is used
# the number of x's is equal to the number of fastq files
# This count is needed to know when to start qc. CHANGE FOR FILE EXTENSION
count=$(find ./fastq -mindepth 1 -type f -name "*${suffix1}" -printf x | wc -c)  
echo There are $count sets of files

# 3. Create a Meta file to know when qc will begin (empty file)
cat >${folder}/'alignMeta.txt' <<EOF
EOF
cd fastq

# 4. Loop over each file that has the specified suffix
for file in *${suffix1}; do
	base=$(basename "$file" "${suffix1}")
	smallBase=${base%_S*}
	echo ${smallBase}

	# for each file, create an R script in the counts folder
	# this R script produces a table of counts data using the rpkm packacge
	cat >${folder}/counts/${smallBase}_RPKM'.R' <<EOF
library(edgeR)
leng <- read.table("${smallBase}_featurecounts_Length.txt",header = TRUE,skip=1)
data <- read.table("${smallBase}_featurecounts_Count.txt",header = TRUE,skip=1)
geneid <- read.table("${smallBase}_featurecounts_Name.txt",header = TRUE,skip=1)

names(data)[names(data) == "aligned.${smallBase}_sorted.bam"] <- "Counts"
rpkm <- rpkm(data,leng)
names(rpkm)[names(rpkm) == "Length"] <- "RPKM"
final<- data.frame(geneid,data,rpkm)
write.table(final,file="${smallBase}_RPKM.csv",row.names = FALSE,quote = FALSE,sep = ",")
EOF

	# for each fastq file, create a PBS script that does the following:
		# Clumping, trimming, and alignment of reads, resulting in a BAM file
		# Generation of a BigWig track and a featureCounts file
		# Execution the R Script using the filtered featureCounts file
	# NOTE: DO NOT change cpus-per-task. Needs to be >= 8 for STAR alignment
	cat >${folder}/PBS/$smallBase'.sbatch' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${smallBase}_alignment # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=8

#Number of memory
#SBATCH --mem-per-cpu=16GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=10:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

################################
# Enter your code to run below #
################################

# Move into folder
cd ${folder}

# initialize and activate conda environment with necessary packages
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment
	
mkdir -p fastqc


# STEP 1: Clumpify 
# this sorts through the fastq files and finds reads that have similar portions and are thus likely overlapping 
# by grouping similar reads close to each other, compression into gz files will more efficient, and alignment will be easier
# CHANGE FOR FILE EXTENSION
# NOTE: backslashes allow long commands to be split up across multiple lines
clumpify.sh \
	in1=${folder}/fastq/${base}${suffix1} \
	in2=${folder}/fastq/${base}${suffix1/R1/R2} \
	out1=clumped/${smallBase}_R1_clumped.fastq.gz \
	out2=clumped/${smallBase}_R2_clumped.fastq.gz 

# decompress the gunzipped files containing clumped reads
gunzip clumped/${smallBase}*_clumped.fastq.gz

# STEP 2: BBDuk
# BBTools is a suite of bioinformatics tools used to analyze sequencing data
# BBDuk (Decontamination Using K-mers) - this combines quality trimming, adapter trimming, contaminant filtering via kmer matching
# a K-mer is a length k subsequence (ex: TAG is a 3-mer of TCAAG)
bbduk.sh \
	in1=clumped/${smallBase}_R1_clumped.fastq \
	in2=clumped/${smallBase}_R2_clumped.fastq \
	out1=trimmed/${smallBase}_R1_trimmed.fastq \
	out2=trimmed/${smallBase}_R2_trimmed.fastq \
	ref=/dartfs-hpc/rc/lab/W/WangX/Nicholas/bbmap/resources/adapters.fa \
	ktrim=r k=21 mink=11 hdist=1 tpe tbo

rm clumped/${smallBase}*

# STEP 3: STAR Alignment
# --genomeDir specifies the path to the directory containing the reference genome (in this case human genome 38)
# --runThreadN sets the number of threads to be used for genome generation (must match the number of cores allotted on server)
# --readFilesIn specifies the paths to files containing reads to be mapped (in this case the clumped and trimmed fastq files)
# --outFileNamePrefix specifies the prefix/name of the output file (suffix will be determined by output file type)
# --outSAMtype specifies the output file type (in this case BAM SortedByCoordinate with suffix '.Aligned.sortedByCoord.out.bam')
# --outSAMUnmapped Within outputs unmapped reads to the BAM file was well
# --outSAMattributes Standard outputs the standard SAM attributes (number of loci the read mapped to, multiple alignment index, local alignment score, number of mismatches)
STAR \
	--genomeDir /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/hg38_STAR_core \
	--runThreadN 8 \
	--readFilesIn trimmed/${smallBase}_R1_trimmed.fastq trimmed/${smallBase}_R2_trimmed.fastq \
	--outFileNamePrefix aligned/${smallBase} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard

rm -r aligned/${smallBase}_STARtmp
rm trimmed/${smallBase}*

# STEP 4: Sort and index the alignments using SamTools
# The 'sort -n' command-flag combination sorts aligned reads by read name (this ensures that paired end reads are adjacent to each other)
# The 'index' command creates indices for coordinate-sorted reads for fast access
# input: BAM file containing alignments "sorted by coordinate"
samtools sort \
	-n \
	-o aligned/${smallBase}_sorted.bam \
	aligned/${smallBase}Aligned.sortedByCoord.out.bam
samtools index aligned/${smallBase}Aligned.sortedByCoord.out.bam

# STEP 5: Generate BigWig file using the indexed BAM file
# BigWig files quantify the expression levels of certain regions in the genome based on the number of reads that aligned to those regions
# These files can be viewed in Genome Browsers
bamCoverage \
	-b aligned/${smallBase}Aligned.sortedByCoord.out.bam \
	-o bigwig/${smallBase}.bw

rm aligned/${smallBase}Aligned.sortedByCoord.out.bam*

# STEP 6: Generate featureCounts
# the 'featureCounts' command counts the number of reads that map to genomic features such as genes, exons, and promoters
# inputs: Gene Transfer File (annotation file with positions of genomic features), BAM file with alignments sorted by read name (from STEP 4)
# ouputs: a text file containing the number of reads that map to each genomic feature
# -T flag specifies the number of threads to use
# -s specifies whether reads are stranded (1), unstranded (0), or reverse stranded (2)
# -g specifies which attribute to use when grouping features (exons) into meta-features (genes)
# -a sets the filepath to the annotation file
# -o sets the path to the output counts file
# -p specifies that paired end reads are being used
featureCounts \
	-T 8 \
	-s 0 \
	-g gene_name \
	-a /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
	-o counts/${smallBase}_featurecounts.txt \
	-p \
	aligned/${smallBase}_sorted.bam

# STEP 7: Filter featureCounts
# grep command searches for text that matches a specified pattern
# -v flag inverts the match, selecting the text that does NOT contain the pattern 
# MT- is the pattern being searched for, which is a common prefix for mitchondrial DNA annotations
grep \
	-v \
	MT- \
	counts/${smallBase}_featurecounts.txt > counts/${smallBase}_featurecounts_MTfiltered.txt

# copy the 'featureCounts' text file into the 'results' folder
cp counts/${smallBase}_featurecounts.txt.summary results

# move into the counts directory
cd counts/
	
# STEP 8: Extract columns from 'featureCounts' text file
# the cut command is used to extract specific columns from a tab-delimited text file
# -f 7 extracts the 7th column, which contains the actual counts of each genomic feature
# -f 1 extracts the 1st column, which contains the names of each genomic feature
# -f 6 extracts the 6th column, which contains the length of each genomic feature (useful for normalizing counts by length of read)
cut -f 7 ${smallBase}_featurecounts_MTfiltered.txt > ${smallBase}_featurecounts_Count.txt
cut -f 1 ${smallBase}_featurecounts_MTfiltered.txt > ${smallBase}_featurecounts_Name.txt
cut -f 6 ${smallBase}_featurecounts_MTfiltered.txt > ${smallBase}_featurecounts_Length.txt

# STEP 9: Switch to a new conda environment that has the necessary packages for running the R script created earlier
source deactivate
source activate deseq
Rscript ${smallBase}_RPKM.R

echo "${smallBase} completed!" >> ${folder}/'alignMeta.txt'

# STEP 10: Check to see if all the files have been analyzed
# if so, will run fastqc and multiqc on fastq and bam files 
# (requires access to Nick's pipeline folder)

currLine=\$(wc -l < ${folder}/alignMeta.txt)
echo \${currLine}
if ((\$currLine == $count)); then
    cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/qc.sh ${folder}
    sh qc.sh
	cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/rna_multi.sh ${folder}
    sh rna_multi.sh
    rm ${folder}/alignMeta.txt
	rmdir clumped/
	rmdir trimmed/
fi
EOF

	cd ${folder}/log

	sbatch ${folder}/PBS/$smallBase'.sbatch'
	cd ${folder}
done