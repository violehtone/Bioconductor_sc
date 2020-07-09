###################################
# Chapter 5 - Overview of a scRNA-seq workflow
###################################

# Load data
library(scRNAseq)
sce <- MacoskoRetinaData() # Macosko mouse retina data is a scRNA-seq data from Macosko et al. (2016)

# Quality control
library(scater)
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets = list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets = "subsets_Mito_percent")
sce <- sce[, !filtered$discard]

# Normalization
sce <- logNormCounts(sce)

# Feature selection
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)

# Dimensionality reduction
set.seed(1234)
sce <- runPCA(sce, ncomponents = 25, subset_row=hvg)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors = TRUE)

# Clustering
g <- buildSNNGraph(sce, use.dimred = 'PCA')
colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)

# Visualization
plotUMAP(sce, colour_by="label")