library(readxl)
library(DESeq2)
counts_file   <- "/path/to/read_count_data.xlsx"
metadata_file <- "/path/to/read_count_data_metadata_file.xlsx"
counts <- as.data.frame(read_excel(counts_file, sheet = "Sheet1"))
rownames(counts) <- counts[[1]]
counts <- counts[, -1, drop = FALSE]
metadata <- as.data.frame(read_excel(metadata_file))
colnames(counts)        <- trimws(colnames(counts))
metadata$sample_ID      <- trimws(as.character(metadata$sample_ID))
metadata                <- metadata[match(colnames(counts), metadata$sample_ID), , drop = FALSE]
stopifnot(all(colnames(counts) == metadata$sample_ID))
rownames(metadata)      <- metadata$sample_ID
metadata$Age <- as.numeric(metadata$Age)
metadata$RIN <- as.numeric(metadata$RIN)
metadata$Age_scaled <- as.numeric(scale(metadata$Age))
metadata$RIN_scaled <- as.numeric(scale(metadata$RIN))

# ensure Sex is a categorical factor (not numeric 1/2) ---
if (is.numeric(metadata$Sex)) {
  metadata$Sex <- factor(metadata$Sex, levels = c(1, 2), labels = c("Male", "Female"))
} else {
  metadata$Sex <- factor(trimws(as.character(metadata$Sex)))
}

metadata$Group <- factor(trimws(as.character(metadata$Group)))
metadata$Group <- relevel(metadata$Group, ref = "CTRL")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ Age_scaled + Sex + RIN_scaled + Group
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("Group", "SPMI", "CTRL"), independentFiltering = FALSE)
out_csv <- "_DESeq2_results.csv"
write.csv(as.data.frame(res), file = out_csv, row.names = TRUE)
cat("DESeq2 results written to:\n", out_csv, "\n")
