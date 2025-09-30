### Normalizing the read count using DESEq2 package
library(DESeq2)
count_data <- read.csv("your_count_data.csv", row.names = 1)  # Replace with your file path
col_data <- read.csv("your_sample_info.csv", row.names = 1)  # Replace with your file path
all(colnames(count_data) %in% rownames(col_data))
all(colnames(count_data) == rownames(col_data))
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design    = ~ Age + Sex + RIN + Group) 
dds <- estimateSizeFactors(dds)  
normalized_counts <- counts(dds, normalized = TRUE)
head(normalized_counts)
write.csv(normalized_counts, "normalized_counts.csv")