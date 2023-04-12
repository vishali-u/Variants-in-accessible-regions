# Use Seurat to find differentially accessible regions in singleome atac_counts

if (! requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}

if (! requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

# Remove the current BiocManager, BiocGenerics, and GenomeInfoDb libraries
remove.packages('BiocManager')
remove.packages('BiocGenerics')
remove.packages('GenomeInfoDb')

# Install version 1.30.29 of BiocManager
package_url <- 'https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.18.tar.gz'
install.packages(package_url, repos=NULL, type="source")

# Install version 0.45.0 of BiocGenerics (needed for GenomeInfoDb)
bioc_generics_url <- 'https://bioconductor.org/packages/devel/bioc/src/contrib/BiocGenerics_0.45.3.tar.gz'
install.packages(bioc_generics_url, repos=NULL, type="source")

# Install version 1.34.9 of GenomeInfoDb
genome_db_url <- 'https://bioconductor.org/packages/release/bioc/src/contrib/GenomeInfoDb_1.34.9.tar.gz'
install.packages(genome_db_url, repos=NULL, type="source")

setRepositories(ind=1:3)
if (! requireNamespace("Signac", quietly = TRUE)) {
  install.packages("Signac")
}

library(Seurat)
library(data.table)
library(BiocManager)
library(BiocGenerics)
library(GenomeInfoDb)
library(Signac)
library(future)

plan(strategy=multicore, workers=12)

# 1. Load data and create Seurat object

file_path <- '/nethome/kcni/vumaiyalan/BCB330/chromatin_accessibility/trevino_2021/single_atac/'
counts_name <- paste(file_path, 'atac_counts.tsv.gz', sep = '')
metadata_name <- paste(file_path, 'atac_cell_metadata.txt', sep = '')
consensus_peaks_name <- paste(file_path, 'atac_consensus_peaks.bed.gz', sep = '')

atac_counts <- fread(file = counts_name)

meta_data <- fread(file = metadata_name, data.table = FALSE)

consensus_peaks <- fread(file = consensus_peaks_name, data.table = FALSE)

peak_names <- c()

# Add a column to consensus peaks to create a unique name for each peak
for (i in 1:nrow(consensus_peaks)) {
  peak_names[i] <- paste('peak_', i, sep = '')
}
consensus_peaks <- cbind(consensus_peaks, peak_names)

# Set column names in counts file to match row names in metadata file
rownames(meta_data) <- meta_data$Cell.ID
colnames(atac_counts) <- meta_data$Cell.ID 

atac_seurat <- CreateSeuratObject(counts = atac_counts, 
                                  assay = 'peaks', 
                                  meta.data = meta_data,
                                  row.names = consensus_peaks$peak_names)


# 2. Data normalization and linear dimensional reduction

atac_seurat <- RunTFIDF(atac_seurat)
atac_seurat <- FindTopFeatures(atac_seurat, min.cutoff = 'q0')
atac_seurat <- RunSVD(atac_seurat)

# 3. Set identities to seurat_clusters 
Idents(atac_seurat) <- "Iterative.LSI.Clusters"


# 4. Find differentially accessible peaks between clusters

# Rename seurat clusters
atac_seurat <- RenameIdents(object = atac_seurat,
                            'c0' = 'GluN6',
                            'c1' = 'GluN4',
                            'c2' = 'GluN7',
                            'c3' = 'IN2',
                            'c4' = 'IN4',
                            'c5' = 'GluN9',
                            'c6' = 'GluN2',
                            'c7' = 'GluN3',
                            'c8' = 'nIPC',
                            'c9' = 'late RG',
                            'c10' = 'oIPC',
                            'c11' = 'early RG',
                            'c12' = 'IN3',
                            'c13' = 'GluN5',
                            'c14' = 'GluN8',
                            'c15' = 'OPC/Oligo',
                            'c16' = 'IN5',
                            'c17' = 'Peric',
                            'c18' = 'GluN1',
                            'c19' = 'MG',
                            'c20' = 'IN1',
                            'c21' = 'EC')

# gluN4_markers <- FindMarkers(object = atac_seurat, 
#                              ident.1 = "GluN4", 
#                              min.pct = 0.01)

# cluster_names <- c('GluN6', 'GluN4', 'GluN7', 'IN2', 'IN4', 'GluN9', 'GluN2',
#                    'GluN3', 'nIPC', 'late RG', 'oIPC', 'early RG', 'IN3',
#                    'GluN5', 'GluN8', 'OPC/Oligo', 'IN5', 'Peric', 'GluN1',
#                    'MG', 'IN1', 'EC')
# 
# for (i in length(cluster_names)) {
#   markers <- findMarkers(object = atac_seurat, ident.1 = cluster_names[i],
#                          thresh.use = 0.1, min.pct = 0)
#   
#   saveRDS(markers, file = paste(file_path, cluser_names[i], sep = ''))
#   
# }


# Find differentially accessible peaks
da_peaks <- FindAllMarkers(object = atac_seurat,
                           logfc.threshold = 0.1,
                           min.pct = 0)

saveRDS(object = da_peaks, file = paste(file_path, 
                                        'no_location_dar.rds', 
                                        sep = ''))

# Merge da_peaks with consensus_peaks to get locations
da_peaks$gene <- gsub('-', '_', da_peaks$gene)

da_peaks <- merge(consensus_peaks,
                  da_peaks,
                  by.x = 'peak_names',
                  by.y = 'gene')

# Save the file of differentially accessible peaks
saveRDS(object = da_peaks, file = paste(file_path,
                                        'single_dar.rds',
                                        sep = ''))



