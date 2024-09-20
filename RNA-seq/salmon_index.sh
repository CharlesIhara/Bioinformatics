#!/bin/bash

# this script creates a PBS script that generates indices from the reference transcriptome for salmon

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
reference_transcriptome='path/to/reference/transcriptome.fasta'

# setup the folders
cd ${folder}
mkdir -p PBS
mkdir -p log
mkdir -p salmon/salmon_index salmon/results

# write the PBS script
cat > ${folder}/PBS/salmon_indexing'.pbs' <<EOF

# Name of the job
#SBATCH --job-name=${samplename}_salmon_indexing

# Number of compute nodes
#SBATCH --nodes=1

# Number of CPU cores
#SBATCH --cpus-per-task=6

# Number of memory
#SBATCH --mem-per-cpu=8G

# Number of cores (number of tasks per node)
#SBATCH --ntasks-per-node=6

# Walltime (job duration)
#SBATCH -time=12:00:00

# Name of output files to be created
#SBATCH --output=${folder}/log/${samplename}_salmon_index.%j.out
#SBATCH --error=${folder}/log/${samplename}_salmon_index.%j.err

################################
# Enter your code to run below #
################################

# install salmon
# conda install bioconda::salmon

# activate an environment with the salmon package
conda activate working_env

# STEP 1: Create an index
# -t flag sets the filepath for the reference transcriptome
# -i flag sets the output directory for the generated indices
# -k flag sets the k-mer size (default is 31, which is optimal for reads >= 75 bp; set to a smaller odd number for smaller reads)
salmon index \
	-t ${reference_transcriptome} \
	-i ${folder}/salmon_index \
	-k 31
EOF

sbatch ${folder}/PBS/salmon_index'.pbs'

