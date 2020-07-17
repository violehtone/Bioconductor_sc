library(scRNAseq)
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)

#################
# 7.1 Motivation
#################
# list available data sets
browseVignettes("scRNAseq")

# Load the Zeisel et al. (2015) mouse brain scRNA-seq data set
sce.zeisel <- ZeiselBrainData()
sce.zeisel

#################
# 7.2 Library size normalization
#################
# Library size = sum(counts for all genes of each cell)
# Library size factor = counts of a cell / library size
# Sum(library size factors of all cells) = 1

# Calculate the library size factors
lib.sf.zeisel <- librarySizeFactors(sce.zeisel)
summary(lib.sf.zeisel)

# Plot the library size factors
hist(x = log10(lib.sf.zeisel),
     xlab = "Log10[size factor]",
     col = 'grey80')


#################
# 7.3 - Normalization by deconvolution
#################
# Pool-based size factors are deconvolved into cell-based factors for normalization
# of each cell's expression profile.

set.seed(100)
clust.zeisel <- quickCluster(sce.zeisel) # performs clustering of cells
table(clust.zeisel) # inspect the size of clusters (14 in total)

deconv.sf.zeisel <- calculateSumFactors(x = sce.zeisel,
                                        cluster = clust.zeisel) # specify which cells belong to which cluster

summary(deconv.sf.zeisel)

# Plot library size factors vs. deconvolution size factor
plot(x = lib.sf.zeisel,
     y = deconv.sf.zeisel,
     xlab = "Library size factor",
     ylab = "Deconvolution size factor",
     log = 'xy',
     pch = 16,
     col = as.integer(factor(sce.zeisel$level1class)))

abline(a=0, b=1, col='red')


#################
# 7.4 - Normalization by spike-ins
#################
# Spike-in normalization is based on the assumption that the same amount of
# spike-in RNA was added to each cell. Systematic differences in the coverage of
# spike-in transcripts can only be due to cell-specific biases, i.e. capture efficiency
# or sequencing depth

# Practically, spike-in normalization should be used if differences in the 
# total RNA content of individual cells are of interest and must be preserved 
#in downstream analyses.

# Spike-in normalization with the Richard et al. (2018) dataset
sce.richard <- RichardTCellData()
sce.richard <- sce.richard[, sce.richard$`single cell quality` == "OK"]
sce.richard

# Use computeSpikeFactors() to estimate spike-in size factors for all cells
sce.richard <- computeSpikeFactors(sce.richard, "ERCC")
summary(sizeFactors(sce.richard))

# Inspect the relationship between spike-in size factors and deconvolution size factors
to.plot <- data.frame(
  DeconvFactor = calculateSumFactors(sce.richard),
  SpikeFactor = sizeFactors(sce.richard),
  Stimulus = sce.richard$stimulus,
  Time = sce.richard$time
)

# Plot the spike-in size factors vs. deconvolution size factors
ggplot(data = to.plot,
       aes(x = DeconvFactor,
           y = SpikeFactor,
           color=Time)) +
  geom_point() +
  facet_wrap(~Stimulus) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0,
              slope = 1,
              color = "red")


# Inspect the effect of the two normalization strategies
sce.richard.deconv <- logNormCounts(sce.richard, size_factors = to.plot$DeconvFactor)
sce.richard.spike <- logNormCounts(sce.richard, size_factors = to.plot$SpikeFactor)

# Plot the two normalization method results for one gene
gridExtra::grid.arrange(
  plotExpression(sce.richard.deconv,
                 x = "stimulus",
                 colour_by = "time",
                 features = "ENSMUSG00000092341") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("After deconvolution"),
  plotExpression(sce.richard.spike,
                 x = "stimulus",
                 colour_by = "time",
                 features = "ENSMUSG00000092341") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("After spike-in normalization"),
  ncol = 2
)


#################
# 7.5 - Applying the size factors
#################
# 7.5.1 - Scaling and log-transforming

# Perform clustering of cells
set.seed(100)
clust.zeisel <- quickCluster(sce.zeisel)

# Perform scaling normalization by deconvolving size factors from cell pools
sce.zeisel <- computeSumFactors(sce.zeisel, cluster = clust.zeisel, min.mean = 0.1)

# Compute log-transformed normalized expression values from a count matrix
sce.zeisel <- logNormCounts(sce.zeisel)

# Get the names of assay elements
assayNames(sce.zeisel)


# 7.5.2 - Downsampling and log-transforming
# In rare cases, direct scaling of counts is not appropriate
# i.e. mean of log-norm counts != log-transformed mean of normalized counts

# Load a dataset to demonstrate this
library(BiocFileCache)
bfc <- BiocFileCache(ask = FALSE)
qcdata <- bfcrpath(bfc, "https://github.com/LuyiTian/CellBench_data/blob/master/data/mRNAmix_qc.RData?raw=true")

env <- new.env()
load(qcdata, envir=env)
sce.8qc <- env$sce8_qc

# Perform library size normalization and log-transformation
sce.8qc <- logNormCounts(sce.8qc) # Compute log-transformed normalized exp. values
sce.8qc <- runPCA(sce.8qc) # perform dimension reduction with PCA

# Plot the data
gridExtra::grid.arrange(
  plotPCA(sce.8qc,
          colour_by = I(factor(sce.8qc$mix))),
  plotPCA(sce.8qc,
          colour_by = I(librarySizeFactors(sce.8qc))),
  ncol = 2
)

# Perform downsampling of counts of the high-coverage cells to match those of low-coverage cells
sce.8qc2 <- logNormCounts(sce.8qc, downsample=TRUE)
sce.8qc2 <- runPCA(sce.8qc2)

gridExtra::grid.arrange(
  plotPCA(sce.8qc2, 
          colour_by=I(factor(sce.8qc2$mix))),
  plotPCA(sce.8qc2, 
          colour_by=I(librarySizeFactors(sce.8qc2))),
  ncol=2
)

# NOTE! downsampling should only be used after first trying scaled counts and if
# this reveals suspicious trajectories that are strongly correlated with the size
# factors, then downsampling should be tried.

