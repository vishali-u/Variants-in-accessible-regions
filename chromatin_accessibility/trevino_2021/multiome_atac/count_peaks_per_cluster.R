# Find the number of DA peaks per cell type in the multiome ATAC-seq data

library(plyr)
library(ggplot2)

dar_path <- "~/BCB330/chromatin_accessibility/trevino_2021/multiome_atac/multiome_dar.rds"

multiome_dar <- readRDS(file = dar_path)

frequencies <- count(multiome_dar, "cluster")

ggplot(frequencies, aes(x = cluster, y = freq)) + geom_bar(stat = "identity") +
  labs(title = "Number of DA peaks per cluster") + ylab("peak_count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



