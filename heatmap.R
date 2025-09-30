normalized_counts <- read.csv("normalized_counts_matrix.csv", header = T)
heatmap.2(normalized_counts,
          scale = "row",                   
          dendrogram = "both",             
          trace = "none",                  
          col = bluered(100),              
          margins = c(8, 8),               
          cexRow = 0.5,                   
          cexCol = 0.7,                    
          key = TRUE,                      
          key.title = "Expression",        
          key.xlab = "Z-score",            
          main = "Heatmap of Normalized Counts")  