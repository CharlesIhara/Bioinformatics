# this script generates heatmaps of protein binding signal at genomic regions defined by BED files

# assumes that you're comparing two proteins at a time
# assumes that original bed files had the "_merged.bed" suffix
#   example: IgG_1million_4M_merged_CTCF_1000_merged_unique.bed
# assumes that the bigwig files have the "_small.bw" suffix
#   example: CTCF_1000_small.bw
# assumes the following directory structure:
# >folder
#   >this_script.sh
#   >deeptools
#   >bed
#       >intersections
#           >a_merged_b_merged_intersection.bed
#       >unique
#           >a_merged_b_merged_unique.bed
#           >b_merged_a_merged_unique.bed
#   >bigwig
#       >a_small.bw
#       >b_small.bw

# initialize conda
eval "$(conda shell.bash hook)"

# install deeptools of activate an environment that already has it
source activate working_env
# pip install deeptools

# names of proteins to compare
names=("K27me3_200" "IgG_1million_4M")

# specify a "base" name used for file lookup and naming
base="${names[0]}_${names[1]}"

# create folder for the comparison data inside the deeptools directory
mkdir -p deeptools/${base}

# get the paths to the genomic regions (BED) files
# store as one string of space-separated filepaths
region_filepaths="bed/intersections/${names[0]}_merged_${names[1]}_merged_intersect.bed \
                bed/unique/${names[0]}_merged_${names[1]}_merged_unique.bed \
                bed/unique/${names[1]}_merged_${names[0]}_merged_unique.bed" 

# get the paths to the signal (BigWig) files
# store as one string of space-separated filepaths
signal_filepaths="bigwig/${names[0]}_small.bw bigwig/${names[1]}_small.bw" 

# compute the matrix for the comparison
# align each row/region by its center point
# include 3000 bp before and after each region
computeMatrix reference-point \
    --referencePoint center \
    -a 3000 \
    -b 3000 \
    -R $region_filepaths \
    -S $signal_filepaths \
    -o "$(pwd)/deeptools/${base}/${base}_matrix"

# create a heatmap from matrix
plotHeatmap \
    -m "deeptools/${base}/${base}_matrix" \
    -o "deeptools/${base}/${base}_heatmap.png" \
    --colorMap "magma" \
    --heatmapWidth 10 \
    --heatmapHeight 70 
