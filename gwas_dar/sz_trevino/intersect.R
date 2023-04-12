# Intersect differentially accessible peaks with the SZ GWAS (data from
# Trevino et al. paper)

library(tidyverse)
library(data.table)

# 1. Load the data
gwas_file_path <- '/nethome/kcni/vumaiyalan/BCB330/gwas/sz_gwas/sz_gwas_hg38.rds'
dar_file_path <- '/nethome/kcni/vumaiyalan/BCB330/chromatin_accessibility/trevino_2021/'

# The name of the dar file (either multiome_dar.rds or singleome_dar.rds)
dar_file_name <- 'multiome_dar.rds'

gwas <- readRDS(file = gwas_file_path)

multiome_dar <- readRDS(file = paste(dar_file_path, 
                                     dar_file_name, 
                                     sep = ''))

# Add a SNP to the GWAS for testing
gwas[1,] <- list('fake_entry', 0, 0, 'fake_entry', 
                      2, 184692037, "B", "C")

# Read the file containing the cluster names and store these names
clusters_file <- fread(paste(dar_file_path, 
                             'multiome_cluster_names.txt', 
                             sep = ''))

cluster_names <- clusters_file$Cluster.Name[15:30]

# Access argument passed in from shell script
args <- as.numeric(commandArgs(trailingOnly = TRUE))
ending <- args + 4
intersection <- data.frame()

# Iterate over the five rows determined by args
for (i in args:ending) {
  if (i > nrow(gwas)) {
    break
  }
  
  if (i == args) {
    
    # Get the gwas data and store the chromosome and location 
    row <- as.data.frame(gwas[i,])
    snp_chr <- paste('chr', row[['chromosome']], sep = '')
    snp_location <- row[['position']]
    
    # Initialize the accessible region and cell type columns
    row[['accessible']] <- 'unavailable'
    for (p in 1:length(cluster_names)) {
      row[[cluster_names[p]]] <- 0.00
    }
    
    # Loop over the accessible peaks and search for overlaps
    for (j in 1:nrow(multiome_dar)) {
      peak_chr <- multiome_dar$seqnames[j]
      peak_start <- multiome_dar$start[j]
      peak_end <- multiome_dar$end[j]
      peak_cluster <- as.character(multiome_dar$cluster[j])
      peak_logFC <- multiome_dar$avg_log2FC[j]
      
      if (peak_chr == snp_chr & peak_start <= snp_location & 
          peak_end >= snp_location) {
        row[['accessible']] <- paste(peak_chr, ':', peak_start, '-', 
                                     peak_end, sep = '')
        row[[peak_cluster]] <- peak_logFC
      }
    }
    intersection <- rbind(intersection, row)
  }
  
  else {
    row <- as.data.frame(gwas[i,])
    
    snp_chr <- paste('chr', row[['chromosome']], sep = '')
    snp_location <- row[['position']]
    
    row[['accessible']] <- 'unavailable'
    for (p in 1:length(cluster_names)) {
      row[[cluster_names[p]]] <- 0.00
    }
    
    for (j in 1:nrow(multiome_dar)) {
      peak_chr <- multiome_dar$seqnames[j]
      peak_start <- multiome_dar$start[j]
      peak_end <- multiome_dar$end[j]
      peak_cluster <- as.character(multiome_dar$cluster[j])
      peak_logFC <- multiome_dar$avg_log2FC[j]    
   
      if (peak_chr == snp_chr & peak_start <= snp_location & 
          peak_end >= snp_location) {
        row[['accessible']] <- paste(peak_chr, ':', peak_start, '-', 
                                     peak_end, sep = '')
        row[[peak_cluster]] <- peak_logFC
      }
    }
    intersection <- rbind(intersection, row)
  }
}

# Save the data into a file
saveRDS(intersection, file = 
          paste('/nethome/kcni/vumaiyalan/BCB330/gwas_dar/sz_trevino/intersected_files/intersection', 
                as.character(args),'.rds', sep=''))





