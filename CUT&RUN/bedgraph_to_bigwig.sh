# this script converts a bedgraph files to a bigwig files

# initialize conda
eval "$(conda shell.bash hook)"

# activate environment
conda activate working_env

# install UCSC packages/scripts
#conda install -c bioconda ucsc-bedgraphtobigwig
#conda install bioconda::ucsc-fetchchromsizes

# generate chromosome size file (replace hg19 with your genome build)
#fetchChromSizes hg19 > "$(pwd)/bedgraph/hg19.chrom.sizes"

# path to chromsome size file
chrom_size_filepath="$(pwd)/bedgraph/hg19.chrom.sizes"

# folder of bedgraph files to convert
bedgraph_folder="$(pwd)/bedgraph/small"

# set up directory for output bigwig files
mkdir -p bigwig

# loop over bedgrpah files and convert to bigwigs
for bedgraph in $bedgraph_folder/*; do

    # create biwgwig filepath (%% operator removes everything after first '.')
    out_filepath="bigwig/$(basename "${bedgraph%%.*}").bw"
    echo $out_filepath

    # convert
    bedgraphToBigWig $bedgraph $chrom_size_filepath $out_filepath
done
