# Use Seurat to find differentially accessible regions 
library(Seurat)
library(data.table)
library(tidyverse)
library(Matrix)

file_path = '~/Documents/3rd_Year/BCB330/bcb330_code/chromatin_accessibility/trevino_2021/'

# Read the counts file
multiome_atac_counts <- fread(file = paste(file_path, 
                                           'multiome_atac_counts.tsv',
                                           sep = ''))

metadata <- fread(file = paste(f))

# Create a matrix object from the atac_counts data
count_matrix <- Matrix(as.matrix(multiome_atac_counts), sparse = TRUE)
