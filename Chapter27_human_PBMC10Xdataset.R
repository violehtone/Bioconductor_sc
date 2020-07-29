######################
# 27.1 - Introduction
######################
# - analysis of PBMC dataset


######################
# 27.2.1 - Data loading
######################
library(BiocFileCache)
bfc <- BiocFileCache(ask = FALSE)
exprs.data <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples/cell-vdj/3.1.0",
  "vdj_v1_hs_pbmc3",
  "vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz"))
untar(exprs.data, exdir=tempdir())

library(DropletUtils)

# Create a singleCellExperiment object from reads
sce.pbmc <- read10xCounts(file.path(tempdir(), "filtered_feature_bc_matrix"))
# split the object based on the feature type, creating alternative experiments to hold features that are not in the majority set
sce.pbmc <- splitAltExps(x = sce.pbmc, #sce object
                         f = rowData(sce.pbmc)$Type) # factor (gene exp., or antibody capture)



######################
# 27.2.2 - Quality control
######################
unfiltered <- sce.pbmc

# Discard cells with high mitochondrial proportion and low ADT counts
library(scater)
is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=is.mito))

high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
low.adt <- stats$`altexps_Antibody Capture_detected` < nrow(altExp(sce.pbmc))/2

discard <- high.mito | low.adt
sce.pbmc <- sce.pbmc[,!discard]

# Examine the distribution of each QC metric:
#  - total counts (sum),
#  - detected features (detected)
#  - Mito% (subsets_Mito_percent)
#  - ADT detected (altexps_Antibody Capture_detected)

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- discard

gridExtra::grid.arrange(
  plotColData(unfiltered, y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, y="subsets_Mito_percent",
              colour_by="discard") + ggtitle("Mito percent"),
  plotColData(unfiltered, y="altexps_Antibody Capture_detected",
              colour_by="discard") + ggtitle("ADT detected"),
  ncol=2
)

# Plot mito% against total count for each cell
plotColData(object = unfiltered,
            x="sum",
            y="subsets_Mito_percent",
            colour_by="discard")
+ scale_x_log10()


######################
# 27.2.3 - Normalization
######################
# Compute size factors for gene exp. and ADT counts

library(scran)

set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)

# altExp() - store alternative experiments inside the sce object
# computeMedianFactors() - Define per-cell size factors by taking the median of ratios to a reference exp. profile
altExp(sce.pbmc) <- computeMedianFactors(altExp(sce.pbmc))
sce.pbmc <- logNormCounts(sce.pbmc, use_altexps=TRUE)

# Look a the distribution of size factors vs. library size for each set of features
par(mfrow=c(1,2))
plot(librarySizeFactors(sce.pbmc), sizeFactors(sce.pbmc), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", 
     main="Gene expression", log="xy")
plot(librarySizeFactors(altExp(sce.pbmc)), sizeFactors(altExp(sce.pbmc)), pch=16,
     xlab="Library size factors", ylab="Median-based factors", 
     main="Antibody capture", log="xy")


######################
# 27.2.4 - Dimensionality reduction
######################
# Skip PCA and just compute tSNE and UMAP for visualization purposes
set.seed(100000)
altExp(sce.pbmc) <- runTSNE(altExp(sce.pbmc))

set.seed(1000000)
altExp(sce.pbmc) <- runUMAP(altExp(sce.pbmc))


######################
# 27.2.5 - Clustering
######################
# Perform graph-based clustering on ADT data
g.adt <- buildSNNGraph(altExp(sce.pbmc),
                       k = 10,
                       d = NA)

clust.adt <- igraph::cluster_walktrap(g.adt)$membership
colLabels(altExp(sce.pbmc)) <- factor(clust.adt)

# Inspect amount of cells in each cluster in the ADT data
table(colLabels(altExp(sce.pbmc)))

# Calculate modularity of each cluster from a graph
#  - modularity is a measure of the structure of networks / graphs
#  - it measures the strength of division into clusters
#  - clusters with high modularity have dense connections between the nodes within modules but sparese connections between nodes in different modules
#  ==> High modularity means _good clustering_

mod <- clusterModularity(graph = g.adt,
                         clusters = clust.adt,
                         as.ratio = TRUE)

# Make a heatmap of the modularity
library(pheatmap)
pheatmap::pheatmap(mat = log10(mod + 10),
                   # when both = FALSE ->
                   # hclust object is clustered
                   cluster_row = FALSE,
                   cluster_col = FALSE,
                   color = c("white", "firebrick"))

# Make a t-SNE plot based on ADT exp. values
plotTSNE(altExp(sce.pbmc),
         colour_by="label",
         text_by="label",
         text_col="red")


# Perform additional subclustering to mimic in-silico FACS experiment

set.seed(1010010)

# quickSubCluster(): performs subclustering for all cells within each group
subclusters <- quickSubCluster(sce.pbmc,
                               clust.adt,
                               prepFUN = function(x) {
                                 # detect bio/tech variance
                                 dec <- modelGeneVarByPoisson(x)
                                 # find highly variable genes
                                 top <- getTopHVGs(dec, prop=0.1)
                                 # perform pCA
                                 denoisePCA(x, dec, subset.row=top)
                               },
                               clusterFUN = function(x) {
                                 # clustering
                                 g.gene <- buildSNNGraph(x, k=10, use.dimred = 'PCA')
                                 igraph::cluster_walktrap(g.gene)$membership
                                 })

# Inspect the subclusters in each cluster
data.frame(
  Cluster=names(subclusters),
  Ncells=vapply(subclusters, ncol, 0L),
  Nsub=vapply(subclusters, function(x) length(unique(x$subcluster)), 0L)
)
