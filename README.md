RNA-seq & snRNA-seq • Machine-Learning Analysis

Code for “Characterising Processing Conditions that Artifactually Bias Human Brain Tissue Transcriptomes”

This repository contains all scripts used to produce the analyses and figures for <Yaqubi et al., 2025. Nat. Comm.>.



## Repository Structure

- `rnaseq/` – bulk RNA-seq pipeline 
- `snRNAseq/` – single-nucleus RNA-seq data analysis
- `ML_analysis/` – machine-learning gene selection, SHAP feature importance, RF+NN models
- `utils/` – helper functions


## Requirements

- R >=4.3 with packages:

1- DESeq2 – for normalization and differential expression
2- clusterProfiler – for GO enrichment analysis
3- org.Hs.eg.db – human gene annotation database used by clusterProfiler
4- gplots – (for heatmap.2 function) – heatmap plotting
5- scFlow – single-nucleus RNA-seq workflow
6- SingleCellExperiment – single-cell data container
7- scater – QC and visualization for single-cell data
8- DropletUtils – ambient RNA removal (emptyDrops)
9- Seurat – single-cell data integration, clustering, visualization
10- DoubletFinder – detection of doublets in single-cell data
11-tidyverse – for data wrangling and plotting (includes dplyr, tibble, etc.)
12- ggplot2 – plotting 


- Python >=3.9 with packages:

1-pandas – data manipulation and file I/O
2-numpy – numerical computing
3- matplotlib.pyplot – visualization (plots, volcano plot, etc.)
4- seaborn – statistical data visualization
5- os – file system operations (built-in Python module)
6- scikit-learn (sklearn) – machine learning:
    a. sklearn.decomposition.PCA – principal component analysis
    b. sklearn.cluster.KMeans – clustering
    c. sklearn.preprocessing.StandardScaler – feature scaling
    d. sklearn.ensemble.RandomForestClassifier – feature importance
    e. sklearn.model_selection.StratifiedKFold – cross-validation
    f. sklearn.metrics.confusion_matrix – evaluation metric

7- tensorflow / tensorflow.keras – deep learning:
    a. tensorflow.keras.layers.Dense, Dropout, Input, Model – neural network layers and model definition

8- matplotlib.lines.Line2D – for custom volcano plot legend



HPC & Workflow Management

Uses Nextflow + Singularity to run nf-core/rnaseq on HPC clusters.




Contact

Maintainer: Moein Yaqubi
📧 moein.yaqubi@mail.mcgill.ca