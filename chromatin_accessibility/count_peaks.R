# Create a bar graph plotting the number of cell-type specific 
# peaks each cell has from the chromatin accessibility data

library(ggplot2)
library(tidyverse)

dev_dar = readRDS(url('https://github.com/jess-xia/chromatin-accessibility/blob/main/Developmental/Developmental_DAR.rds?raw=true'))

# Store the names of included cell types
cellNames <- colnames(dev_dar)[2:13]

# Initialize an array to store the number of peaks in each cell type 
peak_counts <- numeric(length(cellNames))

# Iterate through the columns in the data file and count up the number
# of peaks in each cell type
for (i in 1:length(cellNames)) {
  # select the column for the cell type
  cell_df <- select(dev_dar, cellNames[i]) 
  
  cell_name_as_symbol <- as.symbol(cellNames[i])
  
  # select rows that have a 1 in the column for that cell type
  peak_cols <- filter(cell_df, !!cell_name_as_symbol == 1) 
  
  cell_peak_count <- nrow(peak_cols) 
  peak_counts[i] <- cell_peak_count 
}

cells <- 1:length(cell_names)

data <- data.frame(peak_counts = peak_counts)
data$cell_names = cellNames
ggplot(data, aes(x = cellNames, y = peak_counts)) + geom_bar(stat = "identity") +
  labs(title = "Number of cell-type specific peaks per cell type") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

