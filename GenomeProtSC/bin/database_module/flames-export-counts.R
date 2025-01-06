
library(Matrix)
library(dplyr)
library(data.table)
library(optparse)

option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL,
              help="FLAMES directory", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

flames_dir <- opt$directory
output <- opt$output

# function to read and create a sparse matrix from oarfish
read_oarfish <- function(mtx, feature_file, bc_file) {
  sparse_matrix <- t(readMM(mtx))
  rownames(sparse_matrix) <- read.table(feature_file)$V1
  colnames(sparse_matrix) <- read.table(bc_file)$V1
  return(sparse_matrix)
}

# function to process directory and combine sparse matrices
process_directory <- function(directory) {
  # directory <- "test"
  # find all barcode files
  barcode_files <- list.files(directory, pattern = "\\.barcodes\\.txt$", full.names = TRUE)
  
  # get sample IDs from file names
  sample_ids <- gsub("\\.barcodes\\.txt$", "", basename(barcode_files))
  
  # get paths for matrix and feature files
  sample_files <- lapply(sample_ids, function(sample_id) {
    list(
      mtx = file.path(directory, paste0(sample_id, ".count.mtx")),
      feature_file = file.path(directory, paste0(sample_id, ".features.txt")),
      bc_file = file.path(directory, paste0(sample_id, ".barcodes.txt"))
    )
  })
  names(sample_files) <- sample_ids
  
  # create an empty list to store all matrices
  matrices <- list()
  
  for (sample_id in names(sample_files)) {
    files <- sample_files[[sample_id]]
    
    # read sparse matrix
    sparse_matrix <- read_oarfish(mtx = files$mtx, feature_file = files$feature_file, bc_file = files$bc_file)
    
    # prefix cell barcodes with sample ID
    colnames(sparse_matrix) <- paste0(sample_id, "_", colnames(sparse_matrix))
    
    # add the sparse matrix to the list
    matrices[[sample_id]] <- sparse_matrix
  }
  
  # combine all sparse matrices
  combined_sparse <- do.call(cbind, matrices)
  
  # convert the combined sparse matrix to a dense matrix
  dense_matrix <- as.matrix(combined_sparse)
  
  return(dense_matrix)
  
}

# apply function
combined_matrix <- process_directory(flames_dir)
combined_df <- as.data.frame(combined_matrix)
combined_df <- tibble::rownames_to_column(combined_df, "transcript_id")

# export dense matrix
write.table(
  combined_df,
  file = paste0(output, "/transcript_counts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# function to process gene count files
process_gene_counts <- function(directory) {
  # find all gene count files
  gene_count_files <- list.files(directory, pattern = "_gene_count\\.csv$", full.names = TRUE)
  
  # get sample IDs from file names
  sample_ids <- gsub("_gene_count\\.csv$", "", basename(gene_count_files))
  
  # create a list to store data frames
  gene_counts_list <- list()
  
  # read and process each gene count file
  for (i in seq_along(gene_count_files)) {
    sample_id <- sample_ids[i]
    file_path <- gene_count_files[i]
    
    # read the gene count file
    gene_counts <- read.csv(file_path, row.names = 1)
    
    # prefix column names with sample ID
    colnames(gene_counts) <- paste0(sample_id, "_", colnames(gene_counts))
    
    # add the processed data frame to the list
    gene_counts_list[[sample_id]] <- gene_counts
  }
  
  # get all unique row names across all dfs
  all_genes <- unique(unlist(lapply(gene_counts_list, rownames)))
  
  # function to add missing rows with zeros
  fill_missing_rows <- function(df, all_genes) {
    missing_genes <- setdiff(all_genes, rownames(df))
    if (length(missing_genes) > 0) {
      # create new rows with zeros
      new_rows <- matrix(0, 
                         nrow = length(missing_genes), 
                         ncol = ncol(df),
                         dimnames = list(missing_genes, colnames(df)))
      # bind the new rows to the original df
      df_filled <- rbind(df, new_rows)
      # sort rows to maintain consistent order
      df_filled[all_genes,]
    } else {
      df[all_genes,]
    }
  }
  
  # apply the function to each df in the list
  filled_list <- lapply(gene_counts_list, fill_missing_rows, all_genes)
  names(filled_list) <- NULL
  
  # combine all dfs
  combined_gene_counts <- do.call(cbind, filled_list)
  
  return(combined_gene_counts)
  
}

# apply function
combined_gene_count_matrix <- process_gene_counts(flames_dir)
combined_gene_df <- as.data.frame(combined_gene_count_matrix)
combined_gene_df <- tibble::rownames_to_column(combined_gene_df, "gene_id")

# export gene counts
write.table(
  combined_gene_df,
  file = paste0(output, "/gene_counts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


