# Use Seurat to find differentially accessible regions in multiome atac_counts

library(Seurat)
library(data.table)
library(tidyverse)
library(Matrix)
library(Signac)

# 1. Load data and create Seurat object

file_path = '/nethome/kcni/vumaiyalan/BCB330/chromatin_accessibility/trevino_2021/'

multiome_atac_counts <- fread(file = paste(file_path, 
                                           'multiome_atac_counts.tsv.gz',
                                           sep = ''))
meta_data <- fread(file = paste(file_path, 
                                'multiome_cell_metadata.txt.gz',
                                sep = ''), 
                   data.table = FALSE)

consensus_peaks <- fread(file = paste(file_path, 
                                      'multiome_atac_consensus_peaks.txt', 
                                      sep = ''))

# Set column names in counts file to match row names in metadata file
rownames(meta_data) <- meta_data$Cell.ID
colnames(multiome_atac_counts) <- meta_data$Cell.ID 

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

# 3. Set identities to clusters 
Idents(atac_seurat) <- "seurat_clusters"


# 4. Find differentially accessible peaks between clusters

# Rename seurat clusters
atac_seurat <- subset(atac_seurat, idents = 16, invert = TRUE)
atac_seurat <- RenameIdents(object = atac_seurat,
                            'c0' = 'CluN3',
                            'c1' = 'GluN4',
                            'c2' = 'IN1',
                            'c3' = 'GluN2',
                            'c4' = 'IN2', 
                            'c5' = 'GluN6', 
                            'c6' = 'GluN5', 
                            'c7' = 'RG', 
                            'c8' = 'nIPC/GluN1', 
                            'c9' = 'GluN1', 
                            'c10' = 'mGPC/OPC', 
                            'c11' = 'IN3',
                            'c12' = 'IN4',
                            'c13' = 'SP',
                            'c14' = 'GluN7',
                            'c15' = 'MG/EC/Peric.')

#gluN4_markers <- FindMarkers(object = atac_seurat, 
                            # ident.1 = "GluN4", 
                            # min.pct = 0.01)

# Find differentially accessible peaks
da_peaks <- FindAllMarkers(object = atac_seurat)

# 5. Merge da_peaks with consensus_peaks to get locations

da_peaks$gene <- gsub('-', '_', da_peaks$gene)

da_peaks <- merge(consensus_peaks, da_peaks, by.x = 'name', by.y = 'gene')

saveRDS(object = da_peaks, file = paste(file_path, 
                                        'multiome_dar.rds', 
                                        sep = ''))

