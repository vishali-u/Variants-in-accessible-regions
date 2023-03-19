# Find the clusters in the multiome atac-seq data

library(Seurat)
library(data.table)
library(tidyverse)
library(Signac)

# 1. Load data

file_path <- '/nethome/kcni/vumaiyalan/BCB330/chromatin_accessibility/trevino_2021/'

# The name of the atac-seq counts file (either multiome_atac_counts.tsv.gz or )
counts_file_name <- 'multiome_atac_counts.tsv.gz'

# Either multiome_cell_metadata.txt.gz or 
metadata_name <- 'multioime_cell_metadata.txt.gz'

# Either multiome_atac_consensus_peaks.txt or
consensus_peaks_name <- 'multiome_atac_consensus_peaks.txt'


multiome_atac_counts <- fread(file = paste(file_path, 
                                           counts_file_name,
                                           sep = ''))
meta_data <- fread(file = paste(file_path, 
                                metadata_name,
                                sep = ''), 
                   data.table = FALSE)

# Set column names in counts file to match row names in metadata file
rownames(meta_data) <- meta_data$Cell.ID
colnames(multiome_atac_counts) <- meta_data$Cell.ID 

consensus_peaks <- fread(file = paste(file_path, 
                                      consensus_peaks_name,
                                      sep = ''))


atac_seurat <- CreateSeuratObject(counts = multiome_atac_counts, 
                                  assay = 'peaks', 
                                  meta.data = meta_data,
                                  row.names = consensus_peaks$name)


# 2. Data normalization and linear dimensional reduction

atac_seurat <- RunTFIDF(atac_seurat)
# Warning message:
# In RunTFIDF.default(object = GetAssayData(object = object, slot = "counts"),
# Some features contain 0 total counts

atac_seurat <- FindTopFeatures(atac_seurat, min.cutoff = 'q0')
atac_seurat <- RunSVD(atac_seurat)


# 3. Non-linear dimension reduction and clustering
atac_seurat <- RunUMAP(object = atac_seurat, 
                       reduction = 'lsi',
                       dims = 1:10,
                       n.neighbors = 50,
                       min.dist = 0.6, 
                       metric = 'cosine')

atac_seurat <- FindNeighbors(object = atac_seurat,
                             reduction = 'lsi', 
                             dims = 1:10)

atac_seurat <- FindClusters(object = atac_seurat, 
                            verbose = FALSE,
                            algorithm = 1, 
                            resolution = 0.6)


DimPlot(object = atac_seurat, label = TRUE) + NoLegend() 
