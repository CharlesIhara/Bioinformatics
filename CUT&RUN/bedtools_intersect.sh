# !bin/bash

# this script runs bedtools intersect to find overlapping genomic regions between two bed files
# this outputs all the regions in "b" that intersect with atleast 1 region in "a"
# NOTE: two "b"s could intersect with the same "a"

# assumes the following directory structure:
# >folder
#   >this_script.sh
#   >bed
#       >a_folder
#           >a1.bed
#           >a2.bed
#       >b_folder
#           >b1.bed
#           >b2.bed
#       >intersects
#           >a1_b1_intersects.bed
#       >unique
#           >a1_unique.bed
#           >a2_unique.bed
#           >b1_unique.bed
#           >b2_unique.bed

# install bedtools or activate an environment that has it
#conda install bioconda::bedtools

# initialize conda first so that this command will work
eval "$(conda shell.bash hook)"
source activate working_env

# array of bed filenames to use as "-a" (not full filepaths)
a_files=("IgG_1million_4M_merged.bed.txt" \
        "CTCF_1000_merged.bed.txt" \
        "K27me3_200_merged.bed.txt")

# array of corresponding filenames to use as "'b"
b_files=("IgG_1million_4M_merged.bed.txt" \
        "IgG_1million_4M_merged.bed.txt" \
        "IgG_1million_4M_merged.bed.txt")

# check that the number of a_files equals the number of b_files
a_files_length=${#a_files[@]}
b_files_length=${#b_files[@]}
if [ $a_files_length -ne $b_files_length ]; then
    echo "Error: the number of a_files does not match the number of b_files"
    exit 1
fi

# loop through and intersect corresponding a_files and b_files
for ((index=0; index<${a_files_length}; index++)); do

    # use %% operator to remove a pattern (in this case, everything after the first '.')
    a_filename=${a_files[$index]%%.*}
    b_filename=${b_files[$index]%%.*}

    # filepath for the output intersection file
    mkdir -p bed/intersections
    intersect_filepath="bed/intersections/${a_filename}_${b_filename}_intersect.bed"
    
    # filepath for the output unique files (regions in a or b that are not in the other)
    mkdir -p bed/unique
    unique_a_filepath="bed/unique/${a_filename}_${b_filename}_unique.bed"
    unique_b_filepath="bed/unique/${b_filename}_${a_filename}_unique.bed"

    # filepaths for input bed files
    a_filepath="$(pwd)/bed/a_files/${a_files[$index]}"
    b_filepath="$(pwd)/bed/b_files/${b_files[$index]}"
    
    # add permissions to read input files
    # chmod +r $a_filepath
    # chmod +r $b_filepath

    # perform intersection
    # -wa writes the original -a region that intersects with at least one region in -b
    # -v writes the regions that are unique to -a relative to -b
    bedtools intersect -a $a_filepath -b $b_filepath -wa > $intersect_filepath
    bedtools intersect -a $a_filepath -b $b_filepath -v > $unique_a_filepath
    bedtools intersect -a $b_filepath -b $a_filepath -v > $unique_b_filepath
    echo "complete!"
done

