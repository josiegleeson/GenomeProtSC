
library(Seurat)
library(Matrix)
library(dplyr)
library(data.table)
library(tibble)
library(optparse)

option_list = list(
  make_option(c("-g", "--gene_counts"), type="character", default=NULL,
              help="Gene counts file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

gene_counts_file <- opt$gene_counts
output <- opt$output

gene_counts <- fread(gene_counts_file)
gene_counts <- gene_counts %>% remove_rownames %>% column_to_rownames("gene_id")

combined_seurat <- CreateSeuratObject(gene_counts, project = "SeuratProject")

# standard workflow
combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)
combined_seurat <- ScaleData(combined_seurat)
combined_seurat <- RunPCA(combined_seurat)

# neighbors
# calculate cumulative variance explained in PCA
stdev <- combined_seurat[["pca"]]@stdev
# get vector of cumulative variance
cumvar <- cumsum(stdev^2 / sum(stdev^2))
# choose dims based on cumulative variance >90%
dims_to_use <- which(cumvar >= 0.9)[1]
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:dims_to_use)

# clusters
# set resolution based on number of cells
num_cells <- ncol(combined_seurat)
# under 1k, use 0.4. under 5k, use 0.8, over 5k use 1.2
dynamic_resolution <- ifelse(num_cells < 1000, 0.4, ifelse(num_cells > 5000, 1.2, 0.8))
combined_seurat <- FindClusters(combined_seurat, resolution = dynamic_resolution)

# UMAP
combined_seurat <- RunUMAP(combined_seurat, dims = 1:dims_to_use)
# plot
pdf(file = paste0(output, "/UMAP.pdf"), width=8, height=6)
DimPlot(combined_seurat, reduction = "umap", label = TRUE)
dev.off()

# extract cell clusters
metadata_seurat <- combined_seurat@meta.data
metadata_seurat$sample_cellbarcode <- row.names(metadata_seurat)
metadata_seurat$sample <- metadata_seurat$orig.ident 

metadata_seurat <- metadata_seurat %>% dplyr::select(sample_cellbarcode, sample, seurat_clusters)
row.names(metadata_seurat) <- NULL

# export 
write.table(
  metadata_seurat,
  file = paste0(output, "/sample_cellbarcode_cellcluster.txt"),
  sep = "\t",
  quote = FALSE)

saveRDS(combined_seurat, file = paste0(output, "/seurat_object.rds"))


