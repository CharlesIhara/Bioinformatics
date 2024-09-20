#!/bin/bash

# this script creates a PBS script that runs salmon to align and quantify RNA sequencing data
# assumes that the following folder structure:
# >folder
#	>fastq
#		fastq1.fastq
#		fastq2.fastq
#	>PBS
#	>log
#	>salmon
#		>salmon_index
#		>results
#	this_script.sh

# set the samplename (common across fastq files)
# example filename: CR_H514_WT_Diff_A1A_R1_S54_R1_001.fastq
samplename='experiment1'

# get the path to the folder containing this script
folder=$(cd "$(dirname "$0")"; pwd)

# get the path to the reference transcriptome (fasta format)
reference_transcriptome='/path/to/reference/transcriptome.fasta'

# setup the folders
cd ${folder}
mkdir -p PBS
mkdir -p log
mkdir -p salmon/salmon_index salmon/results

# loop through each fastq file
for file in ${folder}/fastq/*; do

	# write the PBS script
	cat > ${folder}/PBS/salmon_alignment'.pbs' <<EOF
# Name of the job
#SBATCH --job-name=${samplename}_salmon_alignment

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=6

# Number of memory
#SBATCH --mem-per-cpu=8G

# Number of cores (number of tasks per node)
#SBATCH --ntasks-per-node=6

# Walltime (job duration)
#SBATCH -time=12:00:00

# Name of output files to be created
#SBATCH --output=${folder}/log/${samplename}_salmon.%j.out
#SBATCH --error=${folder}/log/${samplename}_salmon.%j.err

################################
# Enter your code to run below #
################################

# install salmon
# conda install bioconda::salmon

# activate an environment with the salmon package installed
conda activate working_env

# Salmon alignment

# get the file basename
base = $(basename "$file" '_R1_001.fastq')	

# quanitfy expression:
# -i flag specifies the path to the folder containing the indices
# -l A tells salmon that it should automatically detect the library type of the reads (stranded vs unstranded)
# -1 and -2 are needed to specify the two sets of reads for paired-end data
# -o flag sets the output file name
# --useEM tells the model to use the standard bayesian EM algorithm, rather than the default (variational bayesian EM), which is more accurate
# --seqBias enables the model to learn and correct for sequence-specific biases (e.g. random hexamer priming bias results in preferential sequencing of fragments that begin with certain nucleotide motifs)
# --validateMappings improves sensitivity and specificity of read mapping by removing loci that do not meet a threshold for base-to-base alignment
salmon quant \
	-i ${folder}/salmon_index \
	-l A \
	-1 $file
	-2 ${folder}/fastq/$base_R2_001.fastq
	-o ${folder}/salmon/results/base.salmon
	--seqBias
	--validateMappings
EOF
	# submit the PBS script
	sbatch ${folder}/PBS/salmon_alignment'.pbs'
done
