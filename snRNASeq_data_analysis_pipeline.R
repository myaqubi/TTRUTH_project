## ── 1) load packages ──────────────────────────────
library(scFlow)
library(SingleCellExperiment)
library(scater)
library(DropletUtils)
library(Seurat)
library(DoubletFinder)
library(tidyverse) 
library(ggplot2)

## ── 2) Paths ────────────────────────────────────────────────
base_dir <- "C:/snRNAseq_project"
raw_dir  <- file.path(base_dir, "raw_matrices")
outdir   <- file.path(base_dir, "results_scflow_TTRUTH")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

sample_dirs <- list(
  "sample1" = file.path(raw_dir, "sample1"),
  "sample2" = file.path(raw_dir, "sample2"),
  "sample3" = file.path(raw_dir, "sample3"),
  "sample4" = file.path(raw_dir, "sample4"),
  "sample5" = file.path(raw_dir, "sample5")
)

## ── 3) Load data ─────────────────────────────────────────────
sce_list <- lapply(names(sample_dirs), function(snm){
  sce <- load_sce(sample_dirs[[snm]])
  colnames(sce) <- paste0(snm, "_", colnames(sce))
  sce$manifest <- snm
  sce
})
sce <- merge_sces(sce_list)

## ── 4) Ambient RNA removal ───────────────────────────────────
set.seed(1)
ed <- emptyDrops(counts(sce), lower = 100, retain = "auto", alpha = 0.001, niters = 10000)
sce <- sce[, ed$FDR <= 0.001 & !is.na(ed$FDR)]

## ── 5) QC filtering ─────────────────────────────────────────
sce <- addPerCellQC(sce)
qc_filter <- sce$sum >= 400 & sce$sum <= 16000 &
  sce$detected >= 200 & sce$detected <= 8000 &
  sce$subsets_Mito_percent <= 20
sce <- sce[, qc_filter]
sce <- sce[!grepl("^MT-", rowData(sce)$Symbol), ]   # remove mitochondrial genes

## ── 6) Normalisation & PCA ───────────────────────────────────
sce <- normalise_sce(sce, method = "log_scran")
sce <- select_features(sce, nfeatures = 2000)
sce <- run_dimred(sce, method = "PCA", ncomponents = 50)

## ── 7) DoubletFinder ─────────────────────────────────────────
seu <- as.Seurat(sce, counts = "counts", data = "logcounts")
sweep.res   <- paramSweep_v3(seu, PCs = 1:10, sct = FALSE, pN = 0.008, pK = NULL)
sweep.stats <- summarizeSweep(sweep.res)
pk.tab      <- find.pK(sweep.stats)
best.pK     <- as.numeric(as.character(pk.tab$pK[which.max(pk.tab$BCmvn)]))
cat("Optimal pK =", best.pK, "\n")
seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.008, pK = best.pK,
                        nExp = round(0.008 * ncol(seu)),
                        reuse.pANN = FALSE, sct = FALSE)
dbl.col <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)
sce$DoubletFinder <- seu@meta.data[[dbl.col]]
sce <- sce[, sce$DoubletFinder == "Singlet"]

## ── 8) UMAP & clustering ─────────────────────────────────────
sce <- run_dimred(sce, method = "UMAP")
sce <- cluster_sce(sce, method = "leiden", resolution = 0.6)

## ── 9) Cell-type annotation ──────────────────────────────────
sce <- annotate_celltypes(sce, reference = "Azimuth", species = "human")
sce$Cell_Type <- sce$celltype_annotation
Idents(sce)   <- "Cell_Type"

## ── 10) TTRUTH gene-profile placeholder ──────────────────────
TTRUTH_profile <- c()    # TTRUTH gene symbols 

## ── 11) Module score & FeaturePlot ───────────────────────────
sce <- AddModuleScore(object = sce, features = list(TTRUTH_profile), name = "TTRUTH_Score")

FeaturePlot(sce, features = "TTRUTH_Score1", reduction = "umap") +
  ggtitle("Average TTRUTH-module expression on UMAP")
ggsave(file.path(outdir, "TTRUTH_Module_FeaturePlot_UMAP.png"),
       width = 7, height = 6, dpi = 300)

## ── 12) PMI-based simple DotPlots ────────────────────────────
Glutamatergic   <- subset(sce, idents = "Glutamatergic")
Oligodendrocyte <- subset(sce, subset = Cell_Type == "Oligodendrocyte")

manifest_to_PMI <- c("sample1"="0h","sample2"="0h","sample3"="6h","sample4"="24h","sample5"="24h")
Glutamatergic@meta.data$PMI   <- manifest_to_PMI[Glutamatergic@meta.data$manifest]
Oligodendrocyte@meta.data$PMI <- manifest_to_PMI[Oligodendrocyte@meta.data$manifest]

p_glut <- DotPlot(Glutamatergic, features = TTRUTH_profile,
                  group.by = "PMI", cols = c("blue","red","green"),
                  scale = TRUE, dot.scale = 16)
ggsave(file.path(outdir, "Glutamatergic_TTRUTH_DotPlot.tiff"), plot = p_glut,
       device = "tiff", width = 18, height = 4, units = "in", dpi = 300, compression = "lzw")

p_olig <- DotPlot(Oligodendrocyte, features = TTRUTH_profile,
                  group.by = "PMI", cols = c("blue","red","green"),
                  scale = TRUE, dot.scale = 16)
ggsave(file.path(outdir, "Oligodendrocyte_TTRUTH_DotPlot.tiff"), plot = p_olig,
       device = "tiff", width = 18, height = 4, units = "in", dpi = 300, compression = "lzw")

## ── 13) Save object & session info ───────────────────────────
saveRDS(sce, file.path(outdir, "scflow_TTRUTH_processed_sce.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

message("✅ Full pipeline complete – outputs saved in: ", outdir)
