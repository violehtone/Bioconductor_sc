#################
# 10.1 - Motivation
#################

# After annotation based on marker genes, the clusters can be
# treated as proxies for more abstract biological concepts 
# such as cell types or states.

# Load libraries
library(scRNAseq)
library(scater)
library(org.Mm.eg.db)
library(scran)
library(DropletUtils)
library(BiocFileCache)
library(EnsDb.Hsapiens.v86)


# Load and read the data set
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)

# Gene annotation
rownames(sce.pbmc) <- uniquifyFeatureNames(rowData(sce.pbmc)$ID,
                                           rowData(sce.pbmc)$Symbol)

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID")

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID")

# Cell detection (distinguish between true cells and ambient RNA)
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]

# Quality control
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.pbmc <- sce.pbmc[,!high.mito]

# Normalization
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

# Variance modelling
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)

# Dimensionality reduction
set.seed(100000)
sce.pbmc <- denoisePCA(sce.pbmc, subset.row = top.pbmc, technical = dec.pbmc)

set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred="PCA")

set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred="PCA")

# Inspect data
sce.pbmc


#################
# 10.3 - Graph-based clustering
#################
# Graph where each node = cell, edges = similarities between cells

# Use 10 nearest neighbours of each cell to construct a graph
# Walktrap method is used to identify communities

g <- buildSNNGraph(x = sce.pbmc, # data
                   k = 10, # number of nearest neighbors
                   use.dimred = "PCA") #use values in reducedDims(sce.pbmc)

clust <- igraph::cluster_walktrap(g)$membership

table(clust)

# Visualize with t-SNE
colLabels(sce.pbmc) <- factor(clust)
plotReducedDim(sce.pbmc, "TSNE", colour_by = "label")

# Use different values of k
g.5 <- buildSNNGraph(sce.pbmc, k=5, use.dimred = 'PCA')
clust.5 <- igraph::cluster_walktrap(g.5)$membership
table(clust.5)

g.50 <- buildSNNGraph(sce.pbmc, k=50, use.dimred = 'PCA')
clust.50 <- igraph::cluster_walktrap(g.50)$membership
table(clust.50)

# visualize with force-directed layout
set.seed(2000)
reducedDim(sce.pbmc, "force") <- igraph::layout_with_fr(g)
plotReducedDim(sce.pbmc, colour_by = "label", dimred = "force")


#10.3.3 - Other parameters (further tweaking of parameters)
# type = "number" -> weight edges based on the # of nearest neighbors
# type = "jaccard" -> weight edges based on jaccard index

# Build graph objects based on different measures
g.num <- buildSNNGraph(sce.pbmc,
                       use.dimred = "PCA",
                       type = "number")

g.jaccard <- buildSNNGraph(sce.pbmc,
                           use.dimred = "PCA",
                           type = "jaccard")

g.none <- buildKNNGraph(sce.pbmc, use.dimred = "PCA")

# Use different community detection algorithms to build clusters
clust.louvain <- igraph::cluster_louvain(g)$membership
clust.infomap <- igraph::cluster_infomap(g)$membership
clust.fast <- igraph::cluster_fast_greedy(g)$membership
clust.labprop <- igraph::cluster_label_prop(g)$membership
clust.eigen <- igraph::cluster_leading_eigen(g)$membership

# Compare how different clustering strategies differ
library(pheatmap)

tab <- table(paste("Infomap", clust.infomap),
             paste("Walktrap", clust))

ivw <- pheatmap(log10(tab+10),
                main = "Infomap vs. Walktrap",
                color = viridis::viridis(100),
                silent = TRUE)

tab <- table(paste("Fast", clust.fast),
             paste("Walktrap", clust))

fvw <- pheatmap(log10(tab+10),
                main = "Fast-greedy vs. Walktrap",
                color = viridis::viridis(100),
                silent = TRUE)

gridExtra::grid.arrange(ivw[[4]], fvw[[4]])


# Manually tuning the clustering resolution by generating nested clustering
# using the cut_at() function
community.walktrap <- igraph::cluster_walktrap(g)
table(igraph::cut_at(community.walktrap, n = 5))

table(igraph::cut_at(community.walktrap, n = 20))


# 10.3.4 - Assessing cluster separation
# clusterModularity() with as.ratio=TRUE returns
# the ratio of the observed to expected sum of weights between each pair of clusters. 

ratio <- clusterModularity(graph = g,
                           clusters = clust,
                           as.ratio = TRUE)

# In this matrix, each row / column corresponds to a cluster and each
# entry contains the ratio of the observed to total weight of edges
# between cells in the respective clusters.
dim(ratio)

pheatmap(log2(ratio+1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "blue"))(100))


# Form a graph where nodes = clusters to explore the relationships between clusters
cluster.gr <- igraph::graph_from_adjacency_matrix(log2(ratio+1),
                                                  mode = "upper",
                                                  weighted = TRUE,
                                                  diag=FALSE)
set.seed(11001010)
plot(x = cluster.gr,
     edge.width = igraph::E(cluster.gr)$weight*5,
     layout = igraph::layout_with_lgl)



#################
# 10.4 - K-means clustering
#################
# Call kmeans() with the top PCs and set k=10
set.seed(100)
clust.kmeans <- kmeans(reducedDim(sce.pbmc, "PCA"),
                                  centers=10)
table(clust.kmeans$cluster)

colLabels(sce.pbmc) <- factor(clust.kmeans$cluster)
plotReducedDim(sce.pbmc, "TSNE", colour_by = "label")

# Find suitable k value
# -> calculate the gap statistic (goodness of clustering measure) for each 
#    number fo k to find the optimal value
library(cluster)
gaps <- clusGap(reducedDim(sce.pbmc, "PCA"),
                kmeans,
                K.max = 20) 

# Calculate the gap statistic
best.k <- maxSE(gaps$Tab[, "gap"],
                gaps$Tab[, "SE.sim"])

best.k # best value for k = 10!

# Plot to verify
plot(gaps$Tab[,"gap"], xlab="Number of clusters", ylab="Gap statistic")
abline(v=best.k, col="red")


# Set k to a large value to achieve overclustering
set.seed(100)
clust.kmeans2 <- kmeans(reducedDim(sce.pbmc, "PCA"),
                        centers = 20)

table(clust.kmeans2$cluster)

# Plot the clusters
colLabels(sce.pbmc) <- factor(clust.kmeans2$cluster)
plotTSNE(sce.pbmc, colour_by="label", text_by="label")


# 10.4.3 - Assessing cluster separation
#  within-cluster sum of squares (WCSS) for each cluster is the most relevant diagnostic 
# for k-means. WCSS can be used to calculate the RMSD (root-mean squared deviation).

ncells <- tabulate(clust.kmeans2$cluster)

tab <- data.frame(wcss = clust.kmeans2$withinss,
                  ncells = ncells)

tab$rms <- sqrt(tab$wcss / tab$ncells)

tab

# To explore relationships between K-means clusters, calculate distanes between centroids
cent.tree <- hclust(dist(clust.kmeans2$centers),
                    "ward.D2")

plot(cent.tree)


# 10.4.4 - In two-step procedures
# clusterSNNGraph() function from scran can use k-means as an initial step to obtain
# representative centroids that are then subjected to graph-based clustering
set.seed(0101010)
kgraph.clusters <- clusterSNNGraph(sce.pbmc,
                                   use.dimred = "PCA",
                                   use.kmeans = TRUE,
                                   kmeans.centers = 1000,
                                   k = 5)

table(kgraph.clusters)
plotTSNE(sce.pbmc, colour_by = I(kgraph.clusters))



#################
# 10.5 - Hierarchical clustering
#################
# Hclust produces a dendrogram -> this can be used to describe the relationships
# between cells and subpopulations

# However, Hclust. is too slow to be used for larger scRNA-seq datasets!
# i.e. pbmc dataset is too large, but we can use the sce.416b dataset for this:

#################
# PREPARING THE 416B DATASET
#################
#--- loading ---#
library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b") 
sce.416b$block <- factor(sce.416b$block)

#--- gene-annotation ---#
library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.416b)$ENSEMBL <- rownames(sce.416b)
rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
                                   keytype="GENEID", column="SYMBOL")
rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
                                    keytype="GENEID", column="SEQNAME")

library(scater)
rownames(sce.416b) <- uniquifyFeatureNames(rowData(sce.416b)$ENSEMBL, 
                                           rowData(sce.416b)$SYMBOL)

#--- quality-control ---#
mito <- which(rowData(sce.416b)$SEQNAME=="MT")
stats <- perCellQCMetrics(sce.416b, subsets=list(Mt=mito))
qc <- quickPerCellQC(stats, percent_subsets=c("subsets_Mt_percent",
                                              "altexps_ERCC_percent"), batch=sce.416b$block)
sce.416b <- sce.416b[,!qc$discard]

#--- normalization ---#
library(scran)
sce.416b <- computeSumFactors(sce.416b)
sce.416b <- logNormCounts(sce.416b)

#--- variance-modelling ---#
dec.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC", block=sce.416b$block)
chosen.hvgs <- getTopHVGs(dec.416b, prop=0.1)

#--- batch-correction ---#
library(limma)
assay(sce.416b, "corrected") <- removeBatchEffect(logcounts(sce.416b), 
                                                  design=model.matrix(~sce.416b$phenotype), batch=sce.416b$block)

#--- dimensionality-reduction ---#
sce.416b <- runPCA(sce.416b, ncomponents=10, subset_row=chosen.hvgs,
                   exprs_values="corrected", BSPARAM=BiocSingular::ExactParam())

set.seed(1010)
sce.416b <- runTSNE(sce.416b, dimred="PCA", perplexity=10)

#################
# ...
#################
# Compute cell-cell distance matrix using the top PCs and apply Hclust with Ward's method
dist.416b <- dist(reducedDim(sce.416b, "PCA"))
tree.416b <- hclust(dist.416b, "ward.D2")

# Making a prettier dendrogram.
library(dendextend)
tree.416b$labels <- seq_along(tree.416b$labels)

dend <- as.dendrogram(tree.416b, 
                      hang=0.1)

combined.fac <- paste0(sce.416b$block, ".", 
                       sub(" .*", "", sce.416b$phenotype))

labels_colors(dend) <- c(
  `20160113.wild`="blue",
  `20160113.induced`="red",
  `20160325.wild`="dodgerblue",
  `20160325.induced`="salmon"
)[combined.fac][order.dendrogram(dend)]

plot(dend)

# To obtain explicit clusters, we "cut" the tree by removing internal branches
# Performed with cutree() function or dynamicTreeCut -package

library(dynamicTreeCut)
clust.416b <- cutreeDynamic(tree.416b, 
                            distM=as.matrix(dist.416b),
                            minClusterSize=10, 
                            deepSplit=1)

table(clust.416b)
labels_colors(dend) <- clust.416b[order.dendrogram(dend)]
plot(dend)

# t-SNE plot for the 416b dataset
colLabels(sce.416b) <- factor(clust.416b)
plotReducedDim(sce.416b, "TSNE", colour_by="label")


# 10.5.3 - Assessing cluster separation
# We check the separation of the clusters using the silhouette width.
# For each cell, we compute the average distance to cells in each other cluster.

# cells with large positive silhouette widths are closer to other cells in the same
# cluster than to cells in different clusters.

sil <- silhouette(clust.416b, dist = dist.416b)
plot(sil)

# Identify the closest neighboring cluster for cells with neg. widths
neg.widths <- sil[,3] < 0
table(Cluster = sil[neg.widths, 1], Neighbor = sil[neg.widths, 2])


#################
# 10.6 - Evaluating cluster stability
#################
# Cluster stability = how prone the clustering is to small changes (perturbations) in the data
myClusterFUN <- function(x) {
  g <- buildSNNGraph(x, use.dimred="PCA", type="jaccard")
  igraph::cluster_louvain(g)$membership
}

originals <- myClusterFUN(sce.pbmc)

set.seed(0010010100)
coassign <- bootstrapCluster(sce.pbmc, FUN=myClusterFUN, clusters=originals)
dim(coassign)

pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE,
         color=rev(viridis::magma(100)))


#################
# 10.7 - Subclustering
#################
# Another simple approach to improving resolution is to repeat the
# feature selection and clustering within a single cluster

# plot the expression values for several T cell markers
g.full <- buildSNNGraph(sce.pbmc, use.dimred = "PCA")
clust.full <- igraph::cluster_walktrap(g.full)$membership
plotExpression(sce.pbmc,
               features = c("CD3E", "CCR7", "CD69", "CD44"),
               x = I(factor(clust.full)),
               colour_by = I(factor(clust.full)))

# Repeating modelling and PCA on the subset
memory <- 6L
sce.memory <- sce.pbmc[, clust.full == memory]
dec.memory <- modelGeneVar(sce.memory)
sce.memory <- denoisePCA(sce.memory,
                         technical = dec.memory,
                         subset.row = getTopHVGs(dec.memory,
                                                 prop = 0.1))

# Apply graph-based clustering on memory subset to obtain subclusters
g.memory <- buildSNNGraph(sce.memory, use.dimred = "PCA")
clust.memory <- igraph::cluster_walktrap(g.memory)$membership

plotExpression(sce.memory,
               features=c("CD8A", "CD4"),
               x = I(factor(clust.memory)))

# It is useful to define a custom function that calls desired algorithms to obtain subclustering
# i.e. quickSubCluster()

library(BiocSingular)
set.seed(1000010)
subcluster.out <- quickSubCluster(sce.pbmc, 
                                  groups=clust.full,
                                  prepFUN = function(x) {
                                    dec <- modelGeneVar(x)
                                    input <- denoisePCA(x, 
                                                        technical = dec,
                                                        subset.row = getTopHVGs(dec,
                                                                                prop = 0.1),
                                                        BSPARAM = BiocSingular::IrlbaParam())
                                  },
                                  clusterFUN = function(x) {
                                    g <- buildSNNGraph(x,
                                                       use.dimred = "PCA",
                                                       k = 20)
                                    igraph::cluster_walktrap(g)$membership
                                  })

# Inspect the subclustering
names(subcluster.out)
table(subcluster.out[[1]]$subcluster)








