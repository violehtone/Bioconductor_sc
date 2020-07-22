###############################
# 13.1 - Motivation
###############################
# Large single-cell RNA sequencing (scRNA-seq) projects usually need to generate data across
# multiple batches due to logistical constraints.

# > There are systematic differences in the observed expression in cells from different batches
# > BATCH EFFECT!
# > Handling and eliminating batch effect is critical for downstream analyses!!!

###############################
# 13.2 - Setting up the data
###############################
# Use two separata 10X genomics PBMC datasets generated in 2 different batches

#--- loading ---#
library(TENxPBMCData)
all.sce <- list(
  pbmc3k=TENxPBMCData('pbmc3k'),
  pbmc4k=TENxPBMCData('pbmc4k'),
  pbmc8k=TENxPBMCData('pbmc8k')
)

#--- quality-control ---#
library(scater)
stats <- high.mito <- list()
for (n in names(all.sce)) {
  current <- all.sce[[n]]
  is.mito <- grep("MT", rowData(current)$Symbol_TENx)
  stats[[n]] <- perCellQCMetrics(current, subsets=list(Mito=is.mito))
  high.mito[[n]] <- isOutlier(stats[[n]]$subsets_Mito_percent, type="higher")
  all.sce[[n]] <- current[,!high.mito[[n]]]
}

#--- normalization ---#
all.sce <- lapply(all.sce, logNormCounts)

#--- variance-modelling ---#
library(scran)
all.dec <- lapply(all.sce, modelGeneVar)
all.hvgs <- lapply(all.dec, getTopHVGs, prop=0.1)

#--- dimensionality-reduction ---#
library(BiocSingular)
set.seed(10000)
all.sce <- mapply(FUN=runPCA, x=all.sce, subset_row=all.hvgs, 
                  MoreArgs=list(ncomponents=25, BSPARAM=RandomParam()), 
                  SIMPLIFY=FALSE)

set.seed(100000)
all.sce <- lapply(all.sce, runTSNE, dimred="PCA")

set.seed(1000000)
all.sce <- lapply(all.sce, runUMAP, dimred="PCA")

#--- clustering ---#
for (n in names(all.sce)) {
  g <- buildSNNGraph(all.sce[[n]], k=10, use.dimred='PCA')
  clust <- igraph::cluster_walktrap(g)$membership
  colLabels(all.sce[[n]])  <- factor(clust)
}

# --------- End loading data sets ------------- #
# PBMC3K dataset (and variance modeling results)
pbmc3k <- all.sce$pbmc3k
dec3k <- all.dec$pbmc3k

# PBMC4K dataset (and variance modeling results)
pbmc4k <- all.sce$pbmc4k
dec4k <- all.dec$pbmc4k


# Batch correction preparations:

# Step 1: Subset all batches to the common "universe" of features
universe <- intersect(rownames(pbmc3k), rownames(pbmc4k))
length(universe)

# Subsetting the datasets
pbmc3k <- pbmc3k[universe,]
pbmc4k <- pbmc4k[universe,]

# Subsetting the variance modelling results
dec3k <- dec3k[universe,]
dec4k <- dec4k[universe,]


# Step 2: Rescale each batch to adjust for differences in seq. depth between batches
library(batchelor)
rescaled <- multiBatchNorm(pbmc3k, pbmc4k) # Perform scaling normalization within each batch

pbmc3k <- rescaled[[1]]
pbmc4k <- rescaled[[2]]

# Step 3: Feature selection by averaging the variance components across all batches
library(scran)
combined.dec <- combineVar(dec3k, dec4k)
chosen.hvgs <- combined.dec$bio > 0
sum(chosen.hvgs)



###############################
# 13.3 - Diagnosing batch effects
###############################
# Before correcting batch effects, we should first look at the data and see
# whether there even are any batch effects
rowData(pbmc3k) <- rowData(pbmc4k)
pbmc3k$batch <- "3k"
pbmc4k$batch <- "4k"
uncorrected <- cbind(pbmc3k, pbmc4k)

set.seed(0010101010)
uncorrected <- runPCA(uncorrected,
                      subset_row = chosen.hvgs,
                      BSPARAM = BiocSingular::RandomParam())

# > Ideally each cluster should contain cells from both batches
# > If not, then there is a batch effect present in the data

snn.gr <- buildSNNGraph(uncorrected, use.dimred = "PCA")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
tab <- table(Clusters = clusters,
             Batch = uncorrected$batch)

# inspect table and see if cells are equally distributed between the batches
tab

# plot using tSNE
set.seed(1111001)
uncorrected <- runTSNE(uncorrected, dimred="PCA")
plotTSNE(uncorrected, colour_by="batch")


###############################
# 13.4 - Linear regression
###############################
# Batch effects in bulk RNA sequencing studies are commonly removed with linear regression.
# Linear modelling is the basis of the removeBatchEffect() function from the limma package
# and for comBat() function from the sva package

# We use the rescaleBatches() function from the batchelor package to remove the batch effect.
rescaled <- rescaleBatches(pbmc3k, pbmc4k)

# Inspect new clusters after removing the batch effect
set.seed(1010101010)
rescaled <- runPCA(rescaled, subset_row=chosen.hvgs, exprs_values="corrected")

snn.gr <- buildSNNGraph(rescaled, use.dimred="PCA")
clusters.resc <- igraph::cluster_walktrap(snn.gr)$membership
tab.resc <- table(Cluster=clusters.resc, Batch=rescaled$batch)
tab.resc

rescaled <- runTSNE(rescaled, dimred="PCA")
rescaled$batch <- factor(rescaled$batch)
plotTSNE(rescaled, colour_by="batch")

###############################
# 13.5 - Performing MNN correction
###############################
# Mutual nearest neighbors are pairs of cells from different batches that belong in each other’s set of nearest neighbors
# MNN pairs represent cells from the same biological state prior to the application of a batch effect
# the difference between cells in MNN pairs can be used as an estimate of the batch effec

# 13.5.2 - Application of MNN to the PBMC data
set.seed(1000101001)
mnn.out <- fastMNN(pbmc3k,
                   pbmc4k,
                   d = 50,
                   k = 20,
                   subset.row = chosen.hvgs,
                   BSPARAM = BiocSingular::RandomParam(deferred = TRUE))

# Corrected matrix contains the low-dim corrected coordinates for all cells
dim(reducedDim(mnn.out, "corrected"))

# Reconstructed matrix contains the corrected expression values for each gene in each cell
assay(mnn.out, "reconstructed")


# 13.5.3 - Correction diagnostics
snn.gr <- buildSNNGraph(mnn.out, use.dimred = "corrected")
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership

# Inspect from table
tab.mnn <- table(Cluster = clusters.mnn,
                 Batch = mnn.out$batch)

# Produce a plot to see batch effect
set.seed(0010101010)
mnn.out <- runTSNE(mnn.out, dimred="corrected")

mnn.out$batch <- factor(mnn.out$batch)
plotTSNE(mnn.out, colour_by="batch")

# Inspect the proportion of variance within each batch that is lost during MNN correction
metadata(mnn.out)$merge.info$lost.var

# > Large proportions of lost variance (>10%) suggest that correction is removing genuine biological heterogeneity


###############################
# 13.6 - Preserving bioogical heterogeneity
###############################
# 13.6.1 - Comparison to within-batch clusters
# > Another useful diagnostic check is to compare the clustering within each batch to the clustering of the merged data. 

# For the first batch (adding +10 for a smoother color transition
# from zero to non-zero counts for any given matrix entry).
tab <- table(paste("after", clusters.mnn[rescaled$batch==1]),
             paste("before", colLabels(pbmc3k)))
heat3k <- pheatmap(log10(tab+10), cluster_row=FALSE, cluster_col=FALSE,
                   main="PBMC 3K comparison", silent=TRUE)

# For the second batch.
tab <- table(paste("after", clusters.mnn[rescaled$batch==2]),
             paste("before", colLabels(pbmc4k)))
heat4k <- pheatmap(log10(tab+10), cluster_row=FALSE, cluster_col=FALSE,
                   main="PBMC 4K comparison", silent=TRUE)

gridExtra::grid.arrange(heat3k[[4]], heat4k[[4]])


# Another evaluation approach is to compute the coassignment probabilities
# = the probability that cells from two within-batch clusters are clustered together in the across-batch clustering

# For the first batch.
tab <- coassignProb(colLabels(pbmc3k), clusters.mnn[rescaled$batch==1])

library(pheatmap)
heat3k <- pheatmap(tab, 
                   cluster_row=FALSE, 
                   cluster_col=FALSE,
                   col=rev(viridis::magma(100)), 
                   main="PBMC 3K probabilities",
                   silent=TRUE)

# For the second batch.
tab <- coassignProb(colLabels(pbmc4k),clusters.mnn[rescaled$batch==2])

heat4k <- pheatmap(tab,
                   cluster_row=FALSE, 
                   cluster_col=FALSE,
                   col=rev(viridis::magma(100)), 
                   main="PBMC 4K probabilities", 
                   silent=TRUE)

gridExtra::grid.arrange(heat3k[[4]], heat4k[[4]])

#  summarize the agreement between clusterings by computing the Rand index
# = a simple metric that we can use to assess the preservation of variation by different correction methods.
# > Larger rand indices (i.e., closer to 1) are more desirable

library(fossil)
ri3k <- rand.index(as.integer(clusters.mnn[rescaled$batch==1]),
                   as.integer(colLabels(pbmc3k)))

ri4k <- rand.index(as.integer(clusters.mnn[rescaled$batch==2]),
                   as.integer(colLabels(pbmc4k)))

ri3k
ri4k


# 13.6.2 - Encouraging consistency with marker genes
# > using the marker genes within each dataset as our selected feature set for fastMNN()
# > This represents a semi-supervised approach

# identify the top marker genes from pairwise Wilcoxon ranked sum tests between every pair of clusters within each batch
stats3 <- pairwiseWilcox(pbmc3k, direction="up")
markers3 <- getTopMarkers(stats3[[1]], stats3[[2]], n=10)

stats4 <- pairwiseWilcox(pbmc4k, direction="up")
markers4 <- getTopMarkers(stats4[[1]], stats4[[2]], n=10)

marker.set <- unique(unlist(c(unlist(markers3), unlist(markers4))))
length(marker.set) # getting the total number of genes selected in this manner.

# Use these selected genes for MNN
set.seed(1000110)
mnn.out2 <- fastMNN(pbmc3k, 
                    pbmc4k, 
                    subset.row=marker.set,
                    BSPARAM=BiocSingular::RandomParam(deferred=TRUE))


# Make a tSNE plot of the merged PBMC datasets
mnn.out2 <- runTSNE(mnn.out2, dimred="corrected")

gridExtra::grid.arrange(
  plotTSNE(mnn.out2[,mnn.out2$batch==1], colour_by=I(colLabels(pbmc3k))),
  plotTSNE(mnn.out2[,mnn.out2$batch==2], colour_by=I(colLabels(pbmc4k))),
  ncol=2
)

###############################
# 13.7 - Using the corrected values
###############################
#  it is preferable to perform DE analyses using the uncorrected expression values with blocking on the batch

m.out <- findMarkers(uncorrected,
                     clusters.mnn,
                     block = uncorrected$batch,
                     direction = "up",
                     lfc = 1,
                     row.data = rowData(uncorrected)[,3, drop = FALSE])

demo <- m.out[["10"]]
as.data.frame(demo[1:20,c("Symbol", "Top", "p.value", "FDR")]) 

plotExpression(uncorrected,
               x=I(factor(clusters.mnn)), 
               features="ENSG00000177954",
               colour_by="batch") 
+ facet_wrap(~colour_by)
