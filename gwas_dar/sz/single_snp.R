# Check if a single SNP is located in an accessible region.
# Location of SNP in SZ_finemapped.csv
# Another SNP: rs758749, 19, 56678350

dev_dar <- read_rds(url('https://github.com/jess-xia/chromatin-accessibility/blob/main/Developmental/Developmental_DAR.rds?raw=true'))

# Store the information about the SNP
SNP_name <- "rs758749"
SNP_chr <- 19
SNP_location <- 56678350
cell_name <- numeric(12)
all_cell_names <- colnames(dev_dar)[2:13]

# Iterate over the rows in the CA data
for (i in 1:nrow(dev_dar)) {
  
  # Check if peak and SNP are on the same chromosome
  if (dev_dar$chr[i] == SNP_chr) {
    
    # Check if peak start is before SNP
    if (dev_dar$start[i] < SNP_location) {
      
      # Check if peak end is after SNP
      if (dev_dar$end[i] > SNP_location) {
        
        # Iterate through the row and select the first column containing a 1
        for (j in 2:(length(all_cell_names) + 1)) {
          if (dev_dar[i,][j] == 1) {
            cell_name[j-1] <- colnames(dev_dar)[j]
          }
        }
      }
    }
  }
}
cell_name[cell_name != 0]

