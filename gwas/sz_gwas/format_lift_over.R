# Return the SNPs in the needed format for the lift over in a txt file. 
# The format is: chr#:position-position

library(readxl)
library(dplyr)
library(stringr)

# Read the SNP data contained in the excel file on the specified sheet
# and insert a new column to with the formatted positions
SNP_data <- read_excel('~/BCB330/gwas/sz_gwas/gwas_data_hg37.xlsx', 
                       sheet = 'ST11b 95% Credible Sets k≤3.5',
                       guess_max = 20800) %>%
  mutate(liftover = paste("chr", as.character(chromosome), 
                          ":", as.character(position), "-", 
                          as.character(position), sep=""))

# Export the formatted locations into a txt file for the lift over
write.table(SNP_data$liftover, file = 'snp_locations_for_lift_over.txt', 
            sep = '\n', row.names = FALSE, col.names = FALSE, quote = FALSE)

# The SNPs that could not be lifted over
failures <- c('chr1:150007105-150007105', 'chr1:150147150-150147150', 
              'chr3:16862873-16862873', 'chr6:165183911-165183911', 
              'chr6:165183961-165183961', 'chr6:165209316-165209316')

# Read the new positions for the SNPs that were given by the lift over
new_positions <- read.table("hg38_liftover.bed", header=FALSE)

# Choose rows not in failures
corrected_SNP_data <- filter(SNP_data, !(liftover %in% failures)) %>%
  
  # Add the new positions to the data frame (with column header 'V1')
  cbind(., new_positions) %>% 
  
  # Change the chromosome and position columns to match the lift over results
  mutate(chromosome = str_match(V1, "chr([A-Z0-9]{1,2}):")[,2], 
         position = as.numeric(str_match(V1, ":([0-9]+)-([0-9]+)")[,2]))

# Rename the column named 'V1'
return_data <- subset(corrected_SNP_data, select = -c(V1, liftover))

# Save the data into a file with the given name
saveRDS(return_data, file = '~/BCB330/gwas/sz_gwas/SZ_finemapped_hg38.rds')
