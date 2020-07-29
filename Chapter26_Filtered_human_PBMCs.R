########################
# 26.1 - Introduction
########################
# - Analysis of the PBMC ID dataset
# - Analysis starts from the filtered count matrix

########################
# 26.2 - Data loading
########################
library(TENxPBMCData)

# load 3 different pbmc data sets into a list
all.sce <- list(
  pbmc3k = TENxPBMCData('pbmc3k'),
  pbmc4k = TENxPBMCData('pbmc4k'),
  pbmc8k = TENxPBMCData('pbmc8k')
)

########################
# 26.3 - Quality control
########################
unfiltered <- all.sce

# Filter on mitochondrial proportion
library(scater)
stats <- high.mito <- list()

# For all three datasets, remove the low quality cells based on mitochondiral content
for (n in names(all.sce)) { #loop through all 3 datasets
  current <- all.sce[[n]] # current data set (i.e. pbmc3k)
  is.mito <- grep("MT", rowData(current)$Symbol_TENx)
  stats[[n]] <- perCellQCMetrics(current, subsets = list(Mito = is.mito))
  high.mito[[n]] <- isOutlier(stats[[n]]$subsets_Mito_percent,
                              type = "higher")
  all.sce[[n]] <- current[, !high.mito[[n]]]
}

# for all datasets, add the stats data and add a 'discard' field
# also plot the data
qcplots <- list()
for (n in names(all.sce)) {
  current <- unfiltered[[n]]
  colData(current) <- cbind(colData(current), stats[[n]])
  current$discard <- high.mito[[n]]
  qcplots[[n]] <- plotColData(current, x="sum", y="subsets_Mito_percent",
                              colour_by="discard") + scale_x_log10()
}

#do.call constructs and executes a function call from a name or a function
do.call(gridExtra::grid.arrange, c(qcplots, ncol=3))

# lapply() applies a funcion over a list / vector
lapply(X = high.mito, # list
       FUN = summary) # function to apply for each element


########################
# 26.4 - Normalization
########################
# Perform library size normalization by log-transform
# Perform log-transformation for counts
all.sce <- lapply(X = all.sce,
                  FUN = logNormCounts)

# Inspect the size factors for each data set (= relative bias)
lapply(X = all.sce,
       FUN = function(x) summary(sizeFactors(x)))


########################
# 26.5 - Variance modelling
########################
library(scran)

# model the variance for each gene decomposint it to bio/tech components
all.dec <- lapply(X = all.sce,
                  FUN = modelGeneVar)

# Identify HVGs
all.hvgs <- lapply(X = all.dec,
                   FUN = getTopHVGs,
                   prop = 0.1) # proportion of genes to report as HVGs (10%)

# Plot the mean (avg. expression) and total (variance) per gene
par(mfrow=c(1,3))
for (n in names(all.dec)) {
  curdec <- all.dec[[n]]
  plot(curdec$mean, curdec$total, pch=16, cex=0.5, main=n,
       xlab="Mean of log-expression", ylab="Variance of log-expression")
  curfit <- metadata(curdec)
  curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
}

########################
# 26.6 - Dimensionality reduction
########################
library(BiocSingular)
set.seed(10000)

#mapply applies FUN to all the elements of a list
all.sce <- mapply(x = all.sce, # data
                  FUN = runPCA, # perform PCA
                  subset_row = all.hvgs,
                  MoreArgs = list(ncomponents = 25,
                                  BSPARAM = RandomParam()), #Randomized SVD
                  SIMPLIFY = FALSE) # don't simplify to a vector / matrix

# Perform t-SNE
set.seed(100000)
all.sce <- lapply(X = all.sce,
                  FUN = runTSNE,
                  dimred = "PCA")


# perform UMAP
set.seed(1000000)
all.sce <- lapply(X = all.sce,
                  FUN = runUMAP,
                  dimred = "PCA")



########################
# 26.7 - Clustering
########################
# Perform clustering for each dataset
for (n in names(all.sce)) {
  g <- buildSNNGraph(all.sce[[n]], k=10, use.dimred='PCA')
  clust <- igraph::cluster_walktrap(g)$membership
  colLabels(all.sce[[n]])  <- factor(clust)
}

# inspect the amount of cells in each cluster in each data set
lapply(all.sce, function(x) table(colLabels(x)))

# make t-SNE plots for each dataset
# - colour by clusters (label)
all.tsne <- list()
for (n in names(all.sce)) {
  all.tsne[[n]] <- plotTSNE(all.sce[[n]], colour_by="label") + ggtitle(n)
}
do.call(gridExtra::grid.arrange, c(all.tsne, list(ncol=2)))


########################
# 26.8 - Data integration
########################
# Untill now, all processing has been performed per-dataset
# We will repeat the processes after merging the 3 batches together

# Intersectig the common genes
universe <- Reduce(f = intersect, # intersect function
                   x = lapply(all.sce, rownames))

# FUN = '[' makes the object a nested list
all.sce2 <- lapply(X = all.sce,
                   FUN = "[",
                   i = universe)

all.dec2 <- lapply(X = all.dec,
                   FUN = "[",
                   i=universe,)

# Renormalizing to adjust for differences in depth
library(batchelor)
# Perform scaling normalization within each batch
normed.sce <- do.call(multiBatchNorm, all.sce2)

# Identifying a set of HVGs using stats from all batches
combined.dec <- do.call(what = combineVar, # combine the results of multiple variance decompositions 
                        args = all.dec2)

combined.hvg <- getTopHVGs(combined.dec, n = 5000)

# Merge the data
set.seed(1000101)
merged.pbmc <- do.call(what = fastMNN, # correct for batch effects by using MNN method
                       args = c(normed.sce,
                                list(subset.row = combined.hvg,
                                     BSPARAM = RandomParam())))

# Inspect the % of lost variance
metadata(merged.pbmc)$merge.info$lost.var


# Perform clustering
g <- buildSNNGraph(merged.pbmc, use.dimred="corrected")
colLabels(merged.pbmc) <- factor(igraph::cluster_louvain(g)$membership)

# Inspect clusters by batches (datasets)
table(colLabels(merged.pbmc), merged.pbmc$batch)

# use a TSNE plot to visualize the clusters
set.seed(10101010)
merged.pbmc <- runTSNE(merged.pbmc, dimred="corrected")
gridExtra::grid.arrange(
  plotTSNE(merged.pbmc, colour_by="label", text_by="label", text_colour="red"),
  plotTSNE(merged.pbmc, colour_by="batch")
)




