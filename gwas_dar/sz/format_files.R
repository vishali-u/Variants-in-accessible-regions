# Combine all gwas_dar_intersection.rds files. 
# Change the format of the created table (change from wide to long) and select 
# the specific cell types in which a SNP is accessible

library(tidyverse)

# Combine the gwas_dar_intersection files into one data fram
combined_rds <- list.files(path = '~/BCB330/gwas_dar/sz/intersection_files', 
                           pattern = '*.rds', full.names = TRUE) %>%
  map_dfr(read_rds)


# Get all the cell names (these are the columns that need to be pivoted to 
# longer format)
cell_names = colnames(combined_rds)[25:35]

# Change table from wide to long
formatted <- pivot_longer(combined_rds, cols = cell_names, 
                          names_to = 'cell_type', values_to = 'accessible')

# Only keep the cell types in which the SNP is accessible
filtered <- filter(formatted, accessible == 1)

saveRDS(filtered, file = '~/BCB330/gwas_dar/sz/intersection.rds')