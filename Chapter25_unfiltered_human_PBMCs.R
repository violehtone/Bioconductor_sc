##########################
# 25.1 - Introduction
##########################
# - Analysis of the peripheral blood mononuclear cell (PBMC) dataset from 10X genomics

##########################
# 25.2 - Data loading
##########################
# BFC is used to Create a on-disk cache of files
library(BiocFileCache)

# Create a cache object specifying a location. 
bfc <- BiocFileCache("raw_data", ask = FALSE)

# Define the path of the data set file
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
# Extract files from a tar archive
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

# DropletUtils provides a number of utility functions for handling data from doplet-baesd seq. tech
library(DropletUtils)

# Create a file path
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")

# Create a SingleCellExperiment object from the raw data
# - rows = features (genes, transcripts)
# - columns = cells
sce.pbmc <- read10xCounts(fname, col.names = TRUE)

# Scater is a collection of tools for performing various scRNA-seq analysis with focus on QC and visualization
library(scater)

# Transform the rownames into unique names
# - If name is unique, use it. If not, append the ID to any non-unique value
# - Missing names will be replaced by ID
rownames(sce.pbmc) <- uniquifyFeatureNames(rowData(sce.pbmc)$ID,
                                           rowData(sce.pbmc)$Symbol)

# Load the annotation database
library(EnsDb.Hsapiens.v86)

# use mapIds to establish a mapping between ids and values
location <- mapIds(x = EnsDb.Hsapiens.v86,
                   keys = rowData(sce.pbmc)$ID,
                   column = "SEQNAME", # column to search on (for mapIds)
                   keytype = "GENEID") # keytype that matches the keys used. 


##########################
# 25.3 - Quality control
##########################
# Distinguish between droplets containing cells and ambient RNA
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]

unfiltered <- sce.pbmc

# Calculate the cell quality metrics with perCellQCMetrics
# - results are stored in the 'subsets' field of the stats object
stats <- perCellQCMetrics(x = sce.pbmc, # sce data
                          subsets = list(Mito = which(location == "MT"))) # list used to identify interesting subsets (i.e. mito genes)

# Detect outliers based on mitochondrial transcription
# high.mito is a logical vector (TRUE / FALSE)
high.mito <- isOutlier(stats$subsets_Mito_percent,
                       type = "higher") # higher = find outliers only from the upper tail

# Remove cells with large mitochondrial proportions using it as a proxy for cell damage
sce.pbmc <- sce.pbmc[, !high.mito]

# Inspect the cells with high mito proportions
summary(high.mito) # -> 311 (TRUE) cells were removed

# Combine the perCellQCMetrics (stats) to the pbmc data
colData(unfiltered) <- cbind(colData(unfiltered), stats)
# Define a 'discard' field with the cells that should be removed
unfiltered$discard <- high.mito

# Plot QC metrics on a cell-level
gridExtra::grid.arrange(
  # sum = total number of counts for the cell (= library size)
  plotColData(unfiltered,
              y = "sum",
              colour_by = "discard") +
    scale_y_log10() +
    ggtitle("Total count"),
  # detected = number of features for the cell that have counts above the detection limit (default zero)
  plotColData(unfiltered,
              y = "detected",
              colour_by = "discard") +
    scale_y_log10() +
    ggtitle("Detected features"),
  # discard = cells with high mitochondrial content
  plotColData(unfiltered,
              y = "subsets_Mito_percent",
              colour_by = "discard") +
    scale_y_log10() +
    ggtitle("Mito percent"),
  ncol = 2
)

# -> Each point in the plot is a cell
# -> Colour indicates whether the cell is removed or not
# -> In other words, following cells were removed:
#     - Cells with low libary sizes
#     - Cells with low amount of detected gene transcripts
#     - Cells with high level of mitochondrial transcription

# Plot the relationship between the library size and mitochondrial content
plotColData(unfiltered,
            x = "sum",
            y = "subsets_Mito_percent",
            colour_by = "discard") +
  scale_x_log10()

# -> When library size is low, mitochondrial content is high
# -> Cells are coloured based on whether they were discarded or not


##########################
# 25.4 - Normalization
##########################
library(scran)
set.seed(1000)

# Cluster similar cells based on their expression profiles
clusters <- quickCluster(sce.pbmc)

# Perform scaling-normalization
# = divide all counts for each cell by a cell-specific "size factor"
# Size factor = estimate of the relative bias in that cell
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)

# Compute log-transformed normalized expression values
sce.pbmc <- logNormCounts(sce.pbmc)

# Inspect the size factors
summary(sizeFactors(sce.pbmc))

# Plot the relationship of library size factors vs. size factors
plot(x = librarySizeFactors(sce.pbmc),
     y = sizeFactors(sce.pbmc),
     pch = 16,
     xlab = "Library size factors",
     ylab = "Deconvolution factors",
     log = "xy")


##########################
# 25.5 - Variance modelling
##########################
set.seed(1001)
# Model the variance of log-expression profile for each gene
# Decompose the variation into technical and biological components
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)

# Define a set of highly variable genes (HVGs)
top.pbmc <-getTopHVGs(stats = dec.pbmc, 
                      prop = 0.1)

# plot the mean against total
# - mean: mean normalized log-expression per gene
# - total: variance of the normalized log-exp. per gene
# - bio: biological component of variance
# - tech: technical component of variance

plot(x = dec.pbmc$mean,
     y = dec.pbmc$total,
     pch = 16,
     cex = 0.5,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression")
curfit <- metadata(dec.pbmc)
curve(curfit$trend(x),
      col = 'dodgerblue',
      add = TRUE,
      lwd = 2)


##########################
# 25.6 - Dimensionality reduction
##########################
set.seed(10000)

# Remove principal components corresponding to technical noise
sce.pbmc <- denoisePCA(x = sce.pbmc, # data set
                       subset.row = top.pbmc, #HVGs
                       technical = dec.pbmc) # object containign the technical variation ('tech' field)
                       
set.seed(100000)

# Perform t-SNE for the cells
# t-SNE is used for visualization purposes
# - dimred = specify existing dim.red. results to use ('PCA' field in the sce.pbmc data)
sce.pbmc <- runTSNE(sce.pbmc, dimred = 'PCA')

# UMAP is an alternative for tSNE
set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred = "PCA")

# Inspect the number of PCs retained
ncol(reducedDim(sce.pbmc, 'PCA'))



##########################
# 25.7 - Clustering
##########################
# build a k-nearest neighbor graph of cells based on similarities in expression profiles
g <- buildSNNGraph(x = sce.pbmc, # data
                   k = 10, #number of nearest neighbors
                   use.dimred = 'PCA') # use existing PCA values

# find connected subgraphs (communities) by random walks
# -> Output is a vector of cluster labels
clust <- igraph::cluster_walktrap(g)$membership

# Set the cluster labels for each cell
colLabels(sce.pbmc) <- factor(clust)

# Inspect the amount of cells in each cluster
table(colLabels(sce.pbmc))

# make a t-SNE plot of the data and colour by clustering label
plotTSNE(sce.pbmc, colour_by = "label")


##########################
# 25.8 - Interpretation
##########################

# Find candidate marker genes for groups of cells (clusters) by testing DE between pairs of groups
markers <- findMarkers(sce.pbmc,
                       pval.type = "some", #how p-values are combined
                       direction = "up") #up-regulated genes

# Examiner markers for cluster 7 in more detail
marker.set <- markers[["7"]]
as.data.frame(marker.set[1:30, 1:3])

# Plot the expression of selected genes
plotExpression(sce.pbmc, 
               features=c("CD14", "CD68", "MNDA", "FCGR3A"), 
               x="label",
               colour_by="label")
