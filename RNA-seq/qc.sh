#!/bin/bash

# This script runs fastqc and multiqc for fastq and bam files. 
# modified from /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/qc.sh

# Assumes the following directory system:
#   >folder
#       qc.sh (this script)
#       >fastq (inputs)
#           fastq1.fastq
#           fastq2.fastq
#       >aligned (inputs)
#           bam1.bam
#           bam2.bam
#       >fastqc (outputs)
#           >fastq
#           >bam
#       >PBS
#       >logs

# 1. Run a script to initialize conda
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh

# 2. Get the name of the folder where this script is located
# "$0" is a special variable that contains the filename of this script itself
# dirname is a command used to extract the path to the enclosing directory from the full filepath
# thus, cd "$(dirname "$0")" moves the user into the directory containing this script
# pwd prints out the absolute path for this directory 
folder=$(cd "$(dirname "$0")";pwd)

# 3. Write the Portable Batch System (PBS) job script inside the PBS folder 
cat >${folder}/PBS/'fastqc.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=fastqc # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=8GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=24:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/fastqc.%j.out
#SBATCH --error=${folder}/log/fastqc.%j.err

################################
# Enter your code to run below #
################################

# run a script to initialize conda
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh

# activate a conda environment containing the 'fastqc' and 'multiqc' packages
source activate qc

# move into folder with script
cd ${folder}

# setup directories (-p flag allows the creation of nested directories)
mkdir -p fastqc
mkdir -p fastqc/fastq
mkdir -p fastqc/bam
mkdir -p multiqc

# fastq files
fastqc -o fastqc/fastq fastq/*
multiqc -f -o multiqc -n fastq fastqc/fastq/

# run fastqc - this generates html reports on the quality of fastq reads and alignment files
fastqc -o fastqc/fastq fastq/*.fastq
fastqc -o fastqc/bam aligned/*.bam

# run multiqc - this searches through the specified folders and generates a comprehensive quality report
# the -f flag overwrites old multiqc reports with the same name
# the -o flag specifies the output directory for the multiqc report
# the -n flag specifies the name for the report
# the last argument is the folders to search through
multiqc -f -o multiqc -n fastq fastqc/fastq
multiqc -f -o multiqc -n bam fastqc/bam/
EOF

# 4. Use the sbatch command to run the PBS script we just created
sbatch ${folder}/PBS/$base'fastqc.pbs'