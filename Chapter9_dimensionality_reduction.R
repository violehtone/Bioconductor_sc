#################
# 9.1 - Overview
#################
# Load libraries
library(scRNAseq)
library(scater)
library(org.Mm.eg.db)
library(scran)

# list available data sets
#browseVignettes("scRNAseq")

# Upload Zeisel et al. (2015) mouse brain data set
sce.zeisel <- ZeiselBrainData()

sce.zeisel <- aggregateAcrossFeatures(sce.zeisel, 
                                      id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))

# Gene annotation
rowData(sce.zeisel)$Ensembl <- mapIds(x = org.Mm.eg.db,
                                      keys = rownames(sce.zeisel),
                                      keytype = "SYMBOL",
                                      column = "ENSEMBL")

# Quality control
stats <- perCellQCMetrics(x = sce.zeisel,
                          subsets = list(Mt = rowData(sce.zeisel)$featureType == "mito"))

qc <- quickPerCellQC(df = stats,
                     percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent"))

sce.zeisel <- sce.zeisel[, !qc$discard]


# Normalization
set.seed(1000)
clusters <- quickCluster(sce.zeisel) #perform clustering of cells
sce.zeisel <- computeSumFactors(sce.zeisel, clusters = clusters) #scaling normalization
sce.zeisel <- logNormCounts(sce.zeisel) #log-transformation

# Variance modeling (technical & biological variation based on spike-ins)
dec.zeisel <- modelGeneVarWithSpikes(sce.zeisel, "ERCC")
top.hvgs <- getTopHVGs(dec.zeisel, prop = 0.1) #get a set of highly variable genes

# Inspcect the dataset
sce.zeisel


#################
# 9.2 - Principal component analysis (PCA)
#################
# Perform PCA on the log-normalized expression values using runPCA()
# runPCA captures top 50 PCs by default and stores them in the reducedDims() of the SingleCellExperiment object
# We use only top 2000 genes for the PCA

top.zeisel <- getTopHVGs(dec.zeisel, n = 2000)
set.seed(100)
sce.zeisel <- runPCA(sce.zeisel, subset_row = top.zeisel)
reducedDimNames(sce.zeisel)

dim(reducedDim(sce.zeisel, "PCA")) # 2816 cells x 50 PCs

# Alternative: Use approximate SVD algorithms that only compute the top PCs
# Useful for large datasets
library(BiocSingular)
set.seed(1000)
sce.zeisel <- runPCA(sce.zeisel,
                     subset_row = top.zeisel,
                     BSPARAM = RandomParam(), name = "IRLBA")

reducedDimNames(sce.zeisel)


#################
# 9.3 - Choosing the number of PCs
#################
# How many PCs to use?
# - Many PCs: more noise but avoids discarding biological signal in later PCs
# - Few PCs: less noise but risk of discarding biological signal
# Usually #PCs (=d) is set between 10 to 50.

#9.3.2 - Using the elbow point: identifying the 'elbow point' in the % of variance explained by PCs
library(PCAtools)
percent.var <- attr(reducedDim(sce.zeisel), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var) #<- 7 PCs selected

plot(percent.var, xlab = "PC", ylab = "Variance explained (%)") # plot the elbow
abline(v=chosen.elbow, col='red')


#9.3.3 - Using technical noise to select d
# retain all PCs until the percentage of total var. explained reaches some threshold T
# -> Calculate the proportion of biological variation in the data 

#library(scran)
#set.seed(111001001)
#denoised.pbmc <- denoisePCA(sce.pbmc, technical=dec.pbmc, subset.row=top.pbmc)
#ncol(reducedDim(denoised.pbmc))

# 9.3.4 - d based on population structure
# uses information about the number of subpopulations in the data
# i.e.  if we have 10 different cell types (clusters), we would set d = 10

pcs <- reducedDim(sce.zeisel)
choices <- getClusteredPCs(pcs)
metadata(choices)$chosen # <- 17 PCs selected

# plot the number of clusters as a function of number of PCs
plot(choices$n.pcs, choices$n.clusters,
     xlab="Number of PCs", ylab="Number of clusters")
abline(a=1, b=1, col="red")
abline(v=metadata(choices)$chosen, col="grey80", lty=2)


#9.3.5 - Putting it together
# when d is selected, next thing to do is to subset the SingleCellExperiment object
# i.e. select top 20 PCs
reducedDim(sce.zeisel, "PCA") <- reducedDim(sce.zeisel, "PCA")[, 1:20]

# alternatively, full set of PCs can be retained and the top set can be assigned another name
reducedDim(sce.zeisel, "PCA_20") <- reducedDim(sce.zeisel, "PCA")[,1:20]
reducedDimNames(sce.zeisel)


#################
# 9.4 - Non-negative matrix factorization (NMF)
#################
# NMF is an alternative for PCA for performing dimensionality reduction
# NMF can be performed with the runNMF() function

set.seed(101001)
nmf.zeisel <- runNMF(sce.zeisel,
                     ncomponents = 10,
                     subset_row = top.zeisel)

#---#

# Extracting the basis matrix of per-gene contributions to each factor.
nmf.out <- reducedDim(nmf.zeisel, "NMF")
nmf.basis <- attr(nmf.out, "basis")
colnames(nmf.out) <- colnames(nmf.basis) <- 1:10

# Creating a heatmap where each row is a cell:
per.cell <- pheatmap::pheatmap(nmf.out, silent=TRUE, 
                               main="By cell", show_rownames=FALSE,
                               color=rev(viridis::magma(100)), cluster_cols=FALSE) 

# Creating a heatmap where each row is a gene:
per.gene <- pheatmap::pheatmap(nmf.basis, silent=TRUE, 
                               main="By gene", cluster_cols=FALSE, show_rownames=FALSE,
                               color=rev(viridis::magma(100)))

gridExtra::grid.arrange(per.cell[[4]], per.gene[[4]], ncol=2)

by.factor <- list()
for (x in colnames(nmf.basis)) {
  by.factor[[x]] <- sort(nmf.basis[,x], decreasing=TRUE)
}
lapply(by.factor, head, n=10)


#################
# 9.5 - Dimensionality reduction for visualization
#################
# In addition to preparing the data for downstream analysis, DR can also be used for visualization

# 9.5.2 - Visualizing with pCA
# Plot top 2 PCs
plotReducedDim(sce.zeisel,
               dimred = "PCA",
               colour_by = "level1class")

# Plot several PCs against each other in pairwise plots
plotReducedDim(sce.zeisel,
               dimred = "PCA",
               ncomponents = 4,
               colour_by = "level1class")


# 9.5.3 - tsNE for visualization (instead of PCA)
# tSNE is 'de facto standard' for visualization of scRNA-seq data

set.seed(00101001101)
sce.zeisel <- runTSNE(sce.zeisel, dimred = "PCA")
plotReducedDim(sce.zeisel,
               dimred = "TSNE",
               colour_by = "level1class")

# tSNE involves random initialization so we need to repeat the plotting several times to ensure the results are correct

set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA", perplexity=5)
out5 <- plotReducedDim(sce.zeisel, dimred="TSNE",
                       colour_by="level1class") + ggtitle("perplexity = 5")

set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA", perplexity=20)
out20 <- plotReducedDim(sce.zeisel, dimred="TSNE",
                        colour_by="level1class") + ggtitle("perplexity = 20")

set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA", perplexity=80)
out80 <- plotReducedDim(sce.zeisel, dimred="TSNE", 
                        colour_by="level1class") + ggtitle("perplexity = 80")

multiplot(out5, out20, out80, cols=3)



# 9.5.4 - Uniform manifold approximation and projection (UMAP)
# UMAP is alternative for tSNE

set.seed(1100101001)
sce.zeisel <- runUMAP(sce.zeisel, dimred="PCA")
plotReducedDim(sce.zeisel, dimred="UMAP", colour_by="level1class")

# compared to tSNE ,UMAP is unarguably much faster, and for that reason
# alone, it is increasingly displacing t-SNE as the method of choice
# for visualizing large scRNA-seq data sets.

# NOTES:
# - Clustering shouldn't be performed on t-SNE coordinates
# - Instead, we should perform clustering on the first 10-50 PCs and then
#   visualize the clusters on the tSNE plot

