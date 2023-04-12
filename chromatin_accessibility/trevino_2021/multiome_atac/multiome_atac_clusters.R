# Find the clusters in the multiome atac-seq data

library(Seurat)
library(data.table)
library(tidyverse)
library(Signac)
library(ggplot2)

# 1. Load data

file_path = '/nethome/kcni/vumaiyalan/BCB330/chromatin_accessibility/trevino_2021/'

multiome_atac_counts <- fread(file = paste(file_path, 
                                           'multiome_atac_counts.tsv.gz',
                                           sep = ''))
meta_data <- fread(file = paste(file_path, 
                                'multiome_cell_metadata.txt.gz',
                                sep = ''), 
                   data.table = FALSE)

# Set column names in counts file to match row names in metadata file
rownames(meta_data) <- meta_data$Cell.ID
colnames(multiome_atac_counts) <- meta_data$Cell.ID 


atac_seurat <- CreateSeuratObject(counts = multiome_atac_counts, 
                                  assay = 'peaks', 
                                  meta.data = meta_data)


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

DimPlot(object = atac_seurat, label = TRUE)


# 4. Compare generated clusters with clusters from paper

# Add the original clusters into the seurat object
atac_seurat@meta.data$original_clusters = meta_data$seurat_clusters

DimPlot(object = atac_seurat, label = TRUE, group.by = 'original_clusters') +
  ggtitle(NULL)







