#############################
# 17.1 - Overview
#############################
# Small changes in the cells' transcriptome happen during biological processes
# These processes from sc expression data can be characterized by identifying a "trajectory"
# trajectory = path through the high-dim. expression space that traverses the various cellular states
#  - simple case: trajectory is a path from A to B
#  - complex case: trajectory with multiple end points

# pseudotime = positioning of cells along the trajectory that quantifies the relative activity
#             of the underlying biological process
# i.e. pseudotime might be used to represent the degree of differentiation from a pluripootent cell to a terminal state



# load data set
#--- data-loading ---#
library(scRNAseq)
sce.nest <- NestorowaHSCData()

#--- gene-annotation ---#
library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
anno <- select(ens.mm.v97, keys=rownames(sce.nest), 
               keytype="GENEID", columns=c("SYMBOL", "SEQNAME"))
rowData(sce.nest) <- anno[match(rownames(sce.nest), anno$GENEID),]

#--- quality-control-grun ---#
library(scater)
stats <- perCellQCMetrics(sce.nest)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent")
sce.nest <- sce.nest[,!qc$discard]

#--- normalization ---#
library(scran)
set.seed(101000110)
clusters <- quickCluster(sce.nest)
sce.nest <- computeSumFactors(sce.nest, clusters=clusters)
sce.nest <- logNormCounts(sce.nest)

#--- variance-modelling ---#
set.seed(00010101)
dec.nest <- modelGeneVarWithSpikes(sce.nest, "ERCC")
top.nest <- getTopHVGs(dec.nest, prop=0.1)

#--- dimensionality-reduction ---#
set.seed(101010011)
sce.nest <- denoisePCA(sce.nest, technical=dec.nest, subset.row=top.nest)
sce.nest <- runTSNE(sce.nest, dimred="PCA")

#--- clustering ---#
snn.gr <- buildSNNGraph(sce.nest, use.dimred="PCA")
colLabels(sce.nest) <- factor(igraph::cluster_walktrap(snn.gr)$membership)


#############################
# 17.2 - Obtaining pseudo-times
#############################
# TSCAN package can be used for trajectory inference
#  1) clusters cells and compute cluster centroids
#  2) builds a minimum spanning tree (MST) across centroids

# Perform clustering and calculate centroids
by.cluster <- aggregateAcrossCells(sce.nest,
                                   ids = colLabels(sce.nest))

centroids <- reducedDim(by.cluster, "PCA")

# Calculate distance matrix
dmat <- dist(centroids)
dmat <- as.matrix(dmat)

# Build a minimum spanning tree
g <- igraph::graph.adjacency(dmat,
                             mode = "undirected",
                             weighted = TRUE)
mst <- igraph::minimum.spanning.tree(g)

set.seed(1000)
plot(mst)

# We can build the same lines in a t-SNE plot
pairs <- Matrix::which(mst[] > 0, arr.ind = TRUE)
coords <- reducedDim(by.cluster, "TSNE")
group <- rep(seq_len(nrow(pairs)), 2)
stuff <- data.frame(rbind(coords[pairs[, 1], ],
                          coords[pairs[, 2], ]),
                    group)

plotTSNE(sce.nest,
         colour_by = "label") +
  geom_line(data = stuff,
            mapping = aes(x = X1,
                          y = X2,
                          group = group))


# pseudotime ordering is obtained by projecting cells onto the MST
# > calculated as the distance along the MST from a position X to a "root node"
# > let's choose one of the endpoint nodes as the root

.map2edges <- function(points, center, edge.ends, previous) {
  all.distances <- list()
  all.pseudo <- list()
  edge.len <- list()
  
  # Computing distance of each point from each edge.
  # Edges defined from 'center' to 'edge.ends'.
  for (i in rownames(edge.ends)) {
    edge.end <- edge.ends[i,]
    delta <- center - edge.end
    max.d <- sqrt(sum(delta^2))
    delta <- delta/max.d
    
    centered <- t(t(points) - center)
    proj <- as.numeric(centered %*% delta)
    proj <- pmax(0, pmin(proj, max.d))
    mapped <- outer(proj, delta)
    
    dist <- sqrt(rowSums((centered - mapped)^2))
    all.distances[[i]] <- dist
    all.pseudo[[i]] <- proj
    edge.len[[i]] <- max.d
  }
  
  all.distances <- do.call(cbind, all.distances)
  all.pseudo <- do.call(cbind, all.pseudo)
  chosen <- colnames(all.distances)[max.col(-all.distances)]
  
  # Flipping the distance of points to the previous node,
  # in order to enforce a directional pseudo-time.
  dist.previous <- 0
  if (!is.na(previous)) {
    on.previous <- chosen==previous
    dist.previous <- edge.len[[previous]]
    previous.proj <- dist.previous - all.pseudo[on.previous,previous,drop=FALSE]
    
    if (all(on.previous)) {
      return(list(dist=dist.previous, pseudo=list(previous.proj)))
    }
  }
  
  # Filling out the branches, where points are NA for a branch's
  # pseudo-time if they were assigned to another branch.
  output <- list()
  for (leftover in setdiff(rownames(edge.ends), previous)) {
    empty <- rep(NA_real_, nrow(points))
    if (!is.na(previous)) {
      empty[on.previous] <- previous.proj
    }
    current <- chosen==leftover
    empty[current] <- all.pseudo[current,leftover]
    output[[leftover]] <- empty
  }
  
  list(dist=dist.previous, pseudo=output)
}

originals <- reducedDim(sce.nest, "PCA")
cluster <- colLabels(sce.nest)
starting.cluster <- names(igraph::V(mst)[igraph::degree(mst)==1])[1]
collated <- list()

latest <- starting.cluster
parents <- NA_character_ 
progress <- list(rep(NA_real_, length(cluster)))
cumulative <- 0

while (length(latest)) {
  new.latest <- new.parents <- character(0)
  new.progress <- list()
  new.cumulative <- numeric(0)
  
  for (i in seq_along(latest)) {
    curnode <- latest[i]
    all.neighbors <- names(igraph::adjacent_vertices(mst, curnode, mode="all")[[1]])
    in.cluster <- cluster==curnode 
    
    mapped <- .map2edges(originals[in.cluster,,drop=FALSE], center=centroids[curnode,], 
                         edge.ends=centroids[all.neighbors,,drop=FALSE], previous=parents[i])
    edge.len <- mapped$dist
    pseudo <- mapped$pseudo
    
    collected.progress <- list()
    for (j in seq_along(pseudo)) {
      sofar <- progress[[i]] # yes, using 'i' here.
      sofar[in.cluster] <- pseudo[[j]] + cumulative[i]
      collected.progress[[j]] <- sofar
    }
    
    all.children <- setdiff(all.neighbors, parents[i])
    if (length(all.children)==0) {
      collated[[curnode]] <- collected.progress[[1]]
    } else {
      new.latest <- c(new.latest, all.children)
      new.parents <- c(new.parents, rep(curnode, length(all.children)))
      new.progress <- c(new.progress, collected.progress)
      new.cumulative <- c(new.cumulative, rep(cumulative[i] + edge.len, length(all.children)))
    }
  }
  
  latest <- new.latest
  parents <- new.parents
  progress <- new.progress
  cumulative <- new.cumulative
}
tscan.pseudo <- do.call(cbind, collated)

plotTSNE(sce.nest, colour_by=I(rowMeans(tscan.pseudo, na.rm=TRUE)), text_by="label") +
  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))


#############################
# 17.2.2 - Principal curves
#############################
# “fitting” a one-dimensional curve so that it passes through the cloud of cells in the high-dimensional expression space
# Principal curves is a non-linear generalization of PCA where the axes of most variation are allowed to bend
# Slingshot package can be used to fit a Principal Curve, which yields a pseudotime ordering of cells

library(slingshot)
sce.sling <- slingshot(sce.nest,
                       reducedDim = 'PCA')

head(sce.sling$slingPseudotime_1)

# Visualize the path taken by the fitted curve
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce.sling$slingPseudotime_1, breaks=100)]

# Creating a PCA plot.
plot(reducedDim(sce.sling, "PCA"), col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce.sling), lwd=2, col='black')

# UMAP plot of the data set, where each cell is colored by the slingshot pseudotime ordering
sce.sling <- runUMAP(sce.sling, dimred="PCA")
sce.sling$cell.type <- sce.sling$FACS <- NULL

library(viridis)
ggcells(sce.sling, mapping=aes(x=UMAP.1, 
                               y=UMAP.2, col=slingPseudotime_1)) +
  geom_point() + scale_color_viridis()


# Multiple trajectories (instead of a single one-dimensional trajectory)
sce.sling2 <- slingshot(sce.nest,
                        cluster = colLabels(sce.nest),
                        reducedDim = 'PCA')

plot(reducedDim(sce.sling2, "PCA"), col="grey80", pch=16, asp = 1)
lines(SlingshotDataSet(sce.sling2), lwd=2, col='black')


# We can use slingshotBranchID() to determine whether a particular cell is shared across multiple curves or is unique to a subset of curves (i.e., is located “after” branching).
curve.assignments <- slingBranchID(sce.sling2)
table(curve.assignments)

# For larger datasets, we can speed up the algorithm by approximating each principal curve with a fixed number of points
sce.sling3 <- slingshot(sce.nest,
                        cluster = colLabels(sce.nest),
                        reducedDim = 'PCA',
                        approx_points = 100)

plot(reducedDim(sce.sling3, "PCA"), col="grey80", pch=16, asp = 1)
lines(SlingshotDataSet(sce.sling3), lwd=2, col='black')