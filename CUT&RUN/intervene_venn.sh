#!bin/bash/

# this script generates venn diagrams of intersection counts using Intervene 
# results output to a "Intervene_results" folder in the working directory

# initialize conda
eval "$(conda shell.bash hook)"

# install intervene or activate an environment that has it
conda activate working_env
# pip install intervene

# path to folder containing bed files
folder="$(pwd)/bed/merged/*.bed"

# run intervene
intervene venn -i $folder