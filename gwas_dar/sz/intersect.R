# Overlap the differentially accessible chromatin peaks with the SNPs from the
# GWAS. 

# Access argument passed in from shell script
args <- as.numeric(commandArgs(trailingOnly = TRUE))
ending <- args + 4

# Import chromatin accessibility & gwas data
dev_dar <- readRDS('~/BCB330/chromatin_accessibility/ziffra_2021/developmental_dar.rds')
gwas_data <- readRDS('~/BCB330/gwas/sz_gwas/sz_gwas_hg38.rds')

# Add fake entry to GWAS list to verify that code is working
gwas_data[1,] <- list('fake_entry', 0, 0, 'fake_entry', 
                      1, 629199, "B", "C")
gwas_data[args : ending, "accessible_region"] <- "unavailable"
intersection <- data.frame()
all_cell_names <- colnames(dev_dar)[2:13]

# Loop over the arguments passed in from shell script
for (i in args:ending) {
  
  # Stop looping if beyond the number of GWAS
  if (i > nrow(gwas_data)){
    break
  }
  
  if (i == args) {
    
    # Store all the information stored in one row (from the GWAS)
    row <- as.list(gwas_data[i,])
    
    # Store the location of each SNP
    SNP_chr <- row[['chromosome']] 
    SNP_location <- row[['position']]
    
    # Add columns for accessible region and all cell types
    row[['accessible_region']] <- 'unavailable'
    for (p in 1:length(all_cell_names)) {
      row[[all_cell_names[p]]] <- 0
    }
    
    # Loop over each peak (i.e. row) in the ATAC-seq data
    for (j in 1:nrow(dev_dar)) {
      
      # Store the location of each peak
      peak_chr <- dev_dar$chr[j]
      peak_start <- dev_dar$start[j]
      peak_end <- dev_dar$end[j]
      
      # Check if SNP overlaps with ATAC-seq peak
      if (peak_chr == SNP_chr & peak_start <= SNP_location &
          peak_end >= SNP_location) {
        
        # Iterate through the cell types in dev_dar
        for (k in 2:(length(all_cell_names) + 1)) {
          
          # If the column contains a 1, update the cell type in row to also
          # contain a 1
          if (dev_dar[j,][k] == 1) {
            row[[colnames(dev_dar)[k]]] <- 1
            
            # Add the peak location to accessible_region
            row['accessible_region'] <- paste('chr', peak_chr,':', peak_start, 
                                              '-', peak_end, sep = '')
          }
        }
      }
      next
    }
    
    # Insert the row (containing SNP information and the cell type in which
    # the SNP is an accessible peak)
    intersection <- rbind(intersection, row)
    
  }
  
  else {
    
    # Store all the information stored in one row (from the GWAS data)
    row <- as.list(gwas_data[i,])
    
    # Store the location of each SNP
    SNP_chr <- row[['chromosome']] 
    SNP_location <- row[['position']]
    
    row[['accessible_region']] <- 'unavailable'
    for (p in 1:length(all_cell_names)) {
      row[[all_cell_names[p]]] <- 0
    }
    
    
    # Loop over each peak in the ATAC-seq data
    for (j in 1:nrow(dev_dar)) {
      
      peak_chr <- dev_dar$chr[j]
      peak_start <- dev_dar$start[j]
      peak_end <- dev_dar$end[j]
      
      # Check if SNP overlaps with ATAC-seq peak
      if (peak_chr == SNP_chr & peak_start <= SNP_location &
          peak_end >= SNP_location) {
        
        # Iterate through the cell types
        for (k in 2:(length(all_cell_names) + 1)) {
          
          # If the column contains a 1, update the cell type in row to also
          # contain a 1
          if (dev_dar[j,][k] == 1) {
            row[[colnames(dev_dar)[k]]] <- 1
            row['accessible_region'] <- paste('chr', peak_chr,':', peak_start, '-', 
                                              peak_end, sep = '')
          }
        }
      }
      next
    }
    
    # Insert the row (containing SNP information and the cell type in which
    # the SNP is in an accessible peak)
    intersection <- rbind(intersection, row)
  }
}


# Save the data into a file
saveRDS(intersection, file = 
          paste("~/BCB330/gwas_dar/sz/intersected_files/intersection", 
                                   as.character(args),".rds", sep=""))





#qbatch -w 00:10:00 --ppj 1 ~/SZ_parallel/script_list.txt
#ppj is the number of cores
#showq
#scancel jobnumber