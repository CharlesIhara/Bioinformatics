#!/bin/bash

# RNA sequencing pipeline (reads to counts) to be run locally

# assumes the following directory system:
#   >folder
#       rna-seq_pipeline.sh (this script)
#       >fastq (inputs)
#           fastq1.fastq
#           fastq2.fastq
#       >clumped
#       >trimmed
#       >aligned
#           sam1.sam
#           sam2.sam
#       >fastqc
#           >raw_fastq
#           >trimmed_fastq
#           >bam
#       >multiqc 
#           >raw_fastq
#           >trimmed_fastq
#           >bam
#       >sorted
#           bam1_sorted.bam
#           bam2_sorted.bam
#       >HISAT2
#           >grch38
#       >featureCounts
#           >annotation_file.gtf
#           >results

# set up directories
mkdir -p fastqc/raw_fastq \
    multiqc/raw_fastq \
    clumped \
    trimmed \
    fastqc/trimmed_fastq \
    multiqc/trimmed_fastq \
    HISAT2 \
    aligned \
    sorted \
    fastqc/bam \
    multiqc/bam \
    HISAT2 \
    featureCounts/results

# get the path to parent folder
folder=$(cd "$(dirname "$0")"; pwd)
cd $folder

# get the basename for the fastq file
base='demo'

# activate an environment with the necessary packages: fastqc, multiqc, bbmap (bbduk.sh and clumpify.sh), hisat2, samtools, subread (featureCounts)
# conda init zsh
# conda activate working_env

# STEP 1: fastqc and multiqc on raw reads
fastqc -o fastqc/raw_fastq fastq/*.fastq
multiqc -f -o multiqc/raw_fastq -n raw_fastq fastqc/raw_fastq
echo "finished STEP 1: fastqc and multiqc on raw reads"

# STEP 2: clumpify and trim reads
# NOTE: the bbmap scripts have been installed in the following directory: /Users/charlesihara/bbmap
clumpify.sh in=${folder}/fastq/${base}.fastq out=${folder}/clumped/${base}_clumped.fastq.gz
gzip -d ${folder}/clumped/*.gz

# adapter and quality trimming
# trim bases from the right side with quality scores < 10 and only keep trimmed reads with atleast 50 bp
bbduk.sh \
    in=${folder}/clumped/${base}_clumped.fastq \
    out=${folder}/trimmed/${base}_trimmed.fastq \
    ref=adapters.fa \
    ktrim=r \
    k=21 mink=11 hdist=1 tpe tbo qtrim=r trimq=10 minlen=50
echo "finished STEP 2: clumping and trimming"

# STEP 3: fastqc and multiqc on trimmed reads
fastqc -o fastqc/trimmed_fastq trimmed/*.fastq
multiqc -f -o multiqc/trimmed_fastq -n trimmed_fastq fastqc/trimmed_fastq
echo "finished STEP 3: fastqc and multiqc on trimmed reads"

# STEP 4: alignment (because salmon, STAR, and HISAT2 require memory-intensive index generation, I will use a pre-made index from HISAT2)

# download the indices:
# mkdir -p ${folder}/HISAT2
# cd HISAT2
# wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# gunzip grch38_genome.tar.gz
# tar -xvf grch38_genome.tar
# cd ..

# -q flag specifies that the input file type is fastq
# -rna-strandedness R specifies that the fastq file consists of "reverse compliment" reads (as opposed to forward reads)
# -x specifies the directory path + basename for each genome index file (ex: ~/grch38/genome.1.ht2)
# -U argument specifies the filepath(s) to unpaired read files
# -S sets the output filepath for the SAM file 
hisat2 \
    -q \
    --rna-strandness R \
    -x HISAT2/grch38/genome \
    -U trimmed/${base}_trimmed.fastq \
    -S aligned/${base}_aligned.sam
echo "finished STEP 4: alignment with HISAT2"

# STEP 5: Use samtools to sort the reads by name and convert SAM to BAM
samtools sort \
    -o ${folder}/sorted/${base}_aligned_sorted.bam \
    ${folder}/aligned/${base}_aligned.sam
echo "finished STEP 5: samtools sort and BAM file generation"

# STEP 6: fastqc and multiqc on alignment data
fastqc -o fastqc/bam sorted/*.bam
multiqc -f -o multiqc/bam -n aligned_BAM fastqc/bam
echo "finished STEP 6: fastqc and multiqc on alignment data"

# STEP 7: featureCounts
# -a flag specifies the annotation file (GTF, GFF, or SAF file)
# -o flag specifies the path to the output counts file 
# -s 2 specifies "reverse-stranded" reads
featureCounts \
    -s 2 \
    -a /Users/charlesihara/Desktop/RNAseq/featureCounts/Homo_sapiens.GRCh38.106.gtf \
    -o featureCounts/results/${base}_featureCounts.txt \
    sorted/${base}_aligned_sorted.bam
echo "finished STEP 7: featureCounts"