
library(FLAMES)
library(SingleCellExperiment)
library(optparse)

option_list = list(
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="Genome file", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help="GTF file", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Directory containing FASTQs", metavar="character"),
  make_option(c("-c", "--cell_count"), type="numeric", default=NULL,
              help="Expected number of cells per sample", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

genome_file <- opt$genome
annotation_file <- opt$annotation
cell_num <- opt$cell_count
output <- opt$output
fastqs_dir <- opt$input

fastqs_in <- list.files(fastqs_dir, full.names = TRUE)

expected_cells <- rep(cell_num, length(fastqs_in))

# command
# LD_PRELOAD="/home/josie/.cache/R/basilisk/1.18.0/FLAMES/2.0.1/flames_env/lib/libssl.so.3 /home/josie/.cache/R/basilisk/1.18.0/FLAMES/2.0.1/flames_env/lib/libnghttp2.so.14" Rscript flames-testing.R

if (length(fastqs_in) > 1) {
  message("Multiple fastq files detected")
  FLAMES::sc_long_multisample_pipeline(
    annotation_file,
    fastqs_in,
    genome_file,
    config_file = "/srv/shiny-server/data/config_file.json",
    barcodes_file = NULL,
    expect_cell_number = expected_cells,
    outdir = output
  )
} else {
  message("Single fastq file detected")
  FLAMES::sc_long_pipeline(
    annotation_file,
    fastqs_in,
    genome_fa = genome_file,
    config_file = "/srv/shiny-server/data/config_file.json",
    barcodes_file = NULL,
    expect_cell_number = expected_cells,
    outdir = output
  )
}





















