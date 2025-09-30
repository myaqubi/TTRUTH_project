library(clusterProfiler)
counts_data <- read.csv("DESEq2_output_file.csv", check.names = FALSE)
row.names(counts_data) = counts_data[[1]] 
counts_data <- counts_data[, -1]

### GO for the up reg genes
sigs <- counts_data[counts_data$log2FoldChange > 1,]
up_reg_genes <- row.names(sigs)
Go_results_treat_vs_control_up_reg <- enrichGO(gene = up_reg_genes, OrgDb = 'org.Hs.eg.db',
                                          keyType = "SYMBOL", ont = "BP")
tiff('file_name.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
barplot(Go_results_treat_vs_control_up_reg, showCategory = 20)
dev.off()

### GO for the down reg genes
sigs <- counts_data[counts_data$log2FoldChange < 1,]
down_reg_genes <- row.names(sigs)
Go_results_treat_vs_control_down_reg <- enrichGO(gene = down_reg_genes, OrgDb = 'org.Hs.eg.db',
                                                  keyType = "SYMBOL", ont = "BP")
tiff('file_name.tiff', units="in", width=10, height=2, res=300, compression = 'lzw')
barplot(Go_results_treat_vs_control_down_reg, showCategory = 20)
dev.off()