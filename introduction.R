#################################
# Coding parts of chapters 1-4
#################################

library(BiocManager)
library(SingleCellExperiment)
library(scater)
library(AnnotationHub)
library(ensembldb)
library(uwot)

# Installing packages with bioconductor:
# i.e. BiocManager::install("AnnotationHub")

# Define a gene expression profile
counts_matrix <- data.frame(cell_1 = rpois(10,10),
                            cell_2 = rpois(10,10),
                            cell_3 = rpois(10,30))

rownames(counts_matrix) <- paste0("gene", 1:10)
counts_matrix <- as.matrix(counts_matrix)

# Create a single cell experiment object from the data
sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

# Get gene expressions from sce
assay(sce, "counts")
counts(sce)

# Compute normalized & log-transformed rep. of data
sce <- scater::logNormCounts(sce)

# Inspect the dataset
logcounts(sce) # show the logcounts assay
assays(sce) # show all assays

# Create a new assay
counts_100 <- counts(sce) + 100
assay(sce, "counts_100") <- counts_100 # assign a new entry to assays slot
assays(sce) # new assay is now visible

# Create metadata
cell_metadata <- data.frame(batch = c(1,1,2))
rownames(cell_metadata) <- paste0("cell_", 1:3)

# Define sce again with the metadata included
sce <- SingleCellExperiment(assays = list(counts = counts_matrix),
                            colData = cell_metadata)

# Inspect coldata
colData(sce)
sce$batch

# Append quality control metrics
sce <- scater::addPerCellQC(sce)
colData(sce)[, 1:5]

# Add more fields to metadata
sce$more_stuff <- runif(ncol(sce))
colnames(colData(sce))

# Using colData values for subsetting
# i.e. choose only cells from batch 1
sce[, sce$batch == 1]

# Inspecting rows (rowChanges, rowData)
rowRanges(sce) # <- empty slot

# insert values to rowData slot of sce object
sce <- scater::addPerFeatureQC(sce)
rowData(sce)

# Pull down an Ensembl annotation object
edb <- AnnotationHub()[["AH73881"]] # human, ensembl v97
genes(edb)[,2]

sce[c("gene_1", "gene_4"), ]
sce[c(1,4), ]

# Other metadata
my_genes <- c("gene_1", "gene_5")
metadata(sce) <- list(favorite_genes = my_genes)
metadata(sce)

your_genes <- c("gene_4", "gene_8")
metadata(sce)$your_genes <- your_genes
metadata(sce)

# Single-cell -specific fields
# PCA
sce <- scater::logNormCounts(sce)
sce <- scater::runPCA(sce)
reducedDim(sce, "PCA")

# t-SNE
sce <- scater::runTSNE(sce, perplexity = 0.1)
reducedDim(sce, "TSNE")
reducedDims(sce)

# umap
u <- uwot::umap(t(logcounts(sce)), n_neighbors = 2)
reducedDim(sce, "UMAP_uwot") <- u
reducedDims(sce)
reducedDim(sce, "UMAP_uwot")

# Alternative experiments
spike_counts <- cbind(cell_1 = rpois(5, 10),
                      cell_2 = rpois(5, 10),
                      cell_3 = rpois(5, 30))
rownames(spike_counts) <- paste0("spike_", 1:5)
spike_se <- SummarizedExperiment(list(counts = spike_counts))
spike_se

# Store summarized experiment into sce
altExp(sce, "spike") <- spike_se
altExps(sce)

sub <- sce[, 1:2]
altExp(sub, "spike")

# Size factors
sce <- scran::computeSumFactors(sce)
sizeFactors(sce)

# Manually add size factors
sizeFactors(sce) <- scater::librarySizeFactors(sce)
sizeFactors(sce)

# Column labels

