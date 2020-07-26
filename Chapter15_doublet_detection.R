############################
# 15.1. - Overview
############################
# doublets are artifactual libraries generated from two cells
# it is desirable to identify and remove doublet libraries so that they do not compromise interpretation of the results.

# Strategies for doublet removal:

# 1) Exploit natural genetic variation when pooling cells from multiple donor individuals.
#   Doublets can be identified as libraries with allele combinations that do not exist in any single donor

# 2) Mark a subset of cells (e.g., all cells from one sample) with an antibody conjugated to a different oligonucleotide
#   Upon pooling, libraries that are observed to have different oligonucleotides are considered to be doublets and removed.

# 3) Infer doublets from the expression profiles alone.

##### LOAD DATA ########
library(scRNAseq)
sce.mam <- BachMammaryData(samples="G_1")

#--- gene-annotation ---#
library(scater)
rownames(sce.mam) <- uniquifyFeatureNames(
  rowData(sce.mam)$Ensembl, rowData(sce.mam)$Symbol)

library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.mam)$SEQNAME <- mapIds(ens.mm.v97, keys=rowData(sce.mam)$Ensembl,
                                   keytype="GENEID", column="SEQNAME")

#--- quality-control ---#
is.mito <- rowData(sce.mam)$SEQNAME == "MT"
stats <- perCellQCMetrics(sce.mam, subsets=list(Mito=which(is.mito)))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent")
sce.mam <- sce.mam[,!qc$discard]

#--- normalization ---#
library(scran)
set.seed(101000110)
clusters <- quickCluster(sce.mam)
sce.mam <- computeSumFactors(sce.mam, clusters=clusters)
sce.mam <- logNormCounts(sce.mam)

#--- variance-modelling ---#
set.seed(00010101)
dec.mam <- modelGeneVarByPoisson(sce.mam)
top.mam <- getTopHVGs(dec.mam, prop=0.1)

#--- dimensionality-reduction ---#
library(BiocSingular)
set.seed(101010011)
sce.mam <- denoisePCA(sce.mam, technical=dec.mam, subset.row=top.mam)
sce.mam <- runTSNE(sce.mam, dimred="PCA")

#--- clustering ---#
snn.gr <- buildSNNGraph(sce.mam, use.dimred="PCA", k=25)
colLabels(sce.mam) <- factor(igraph::cluster_walktrap(snn.gr)$membership)



############################
# 15.2 - Doublet detection with clusters
############################
# doubletCluster() function identifies clusters with expression profiles lying between two other clusters 
library(scran)
dbl.out <- doubletCluster(sce.mam)

# Identify clusters that have unusually low N
library(scater)
chosen.doublet <- rownames(dbl.out)[isOutlier(dbl.out$N,
                                    type = "lower",
                                    log = TRUE)]
# > Cluster 6 has the fewest unique genes and library sizes that are comparable to or greater than its sources

markers <- findMarkers(sce.mam, direction = "up")
dbl.markers <- markers[[chosen.doublet]]

chosen <- rownames(dbl.markers)[dbl.markers$Top <= 10]

plotHeatmap(sce.mam,
            order_columns_by = "label",
            features = chosen,
            center = TRUE,
            symmetric = TRUE,
            zlim=c(-5,5))

plotExpression(sce.mam,
               features = c("Acta2", "Csn2"),
               x = "label",
               colour_by = "label")


############################
# 15.3 - Doublet detection by simulation
############################
# doubletCells() function will
# 1. Simulate thousands of doublets by adding together two randomly chosen single-cell profiles.
# 2. For each original cell, compute the density of simulated doublets in the surrounding neighborhood.
# 3. For each original cell, compute the density of other observed cells in the neighborhood.
# 4. Return the ratio between the two densities as a “doublet score” for each cell.
set.seed(100)
dbl.dens <- doubletCells(sce.mam,
                         subset.row = top.mam,
                         d = ncol(reducedDim(sce.mam)))

sce.mam$DoubletScore <- log10(dbl.dens+1)

plotTSNE(sce.mam, colour_by = "DoubletScore")
plotColData(sce.mam,
            x = "label",
            y = "DoubletScore",
            colour_by = "label")

# Simply removing cells with high doublet scores will not be sufficient to eliminate real doublets from the data set.
# Simulation approach is more robust than doubletClusters() to the quality of the clustering as the scores are computed on a per-cell basis.

############################
# 15.4 Doublet detection in multiplexed experiments
############################
# For multiplexed samples, we can identify doublet cells based on the cells that have multiple labels. 

library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
hash.tar <- bfcrpath(bfc, "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108313&format=file")

fname <- "GSM2895283_Hashtag-HTO-count.csv.gz"
untar(hash.tar, files=fname, exdir=tempdir())
hto.counts <- read.csv(file.path(tempdir(), fname), row.names=1)
hto.counts <- as.matrix(hto.counts[1:8,])

dim(hto.counts)

# 15.4.2. - Identifying inter-sample doublets
# Before we proceed to doublet detection, we simplify the problem by first identifying the barcodes that contain cells.
library(DropletUtils)
set.seed(101)
hash.calls <- emptyDrops(hto.counts, lower = 200)
is.cell <- which(hash.calls$FDR <= 0.001)
length(is.cell)

# Visualize
par(mfrow=c(1,2))
r <- rank(-hash.calls$Total)
plot(r,
     hash.calls$Total,
     log="xy",
     xlab = "Rank",
     ylab = "Total HTO count",
     main = "")

hist(log10(hash.calls$Total[is.cell]),
     xlab = "Log[10] HTO count",
     main = "")


# Run hashedDrops() on the subset of cell barcode libraries that actually contain cells.
# This returns the likely sample of origin for each barcode library based on its most abundant HTO
hash.stats <- hashedDrops(hto.counts[, is.cell],
                          ambient = metadata(hash.calls)$ambient)

par(mfrow=c(1,1))
hist(hash.stats$LogFC, xlab="Log fold-change from best to second HTO", main="")

# Raw assignments
table(hash.stats$Best)
# Confident assignments
table(hash.stats$Best[hash.stats$Confident])


# Hashing information can be used to detect doublets by reporting the log-fold change between the count for the second HTO and the estimated contribution from ambient contamination
# > a large log-fold change indicates that the second HTO still has an above-expected abundance
# > use outlier detection to explicitly identify putative doublets as those barcode libraries that have large log-fold changes

summary(hash.stats$Doublet)

# Plot the log-fold change of the second most abundant HTO over ambient contamination
# compared to the log-fold change of the first HTO over the second HTO
# > each point = cell
# > potential doublets = red
# > confidently assigned = black
colors <- rep("grey", nrow(hash.stats))
colors[hash.stats$Doublet] <- "red"
colors[hash.stats$Confident] <- "black"

plot(hash.stats$LogFC, hash.stats$LogFC2,
     xlab="Log fold-change from best to second HTO",
     ylab="Log fold-change of second HTO over ambient",
     col=colors)



# 15.4.3 - Guilt by association for unmarked doublets
# = recover the remaining intra-sample doublets based on their similarity with known doublets in gene expression space (hence, “guilt by association”)

# Load data set
gname <- "GSM2895282_Hashtag-RNA.umi.txt.gz"
untar(hash.tar, files=gname, exdir=tempdir())

# Reading it in as a sparse matrix in a SingleCellExperiment.
library(scater)
gene.counts <- readSparseCounts(file.path(tempdir(), gname))
sce.hash <- SingleCellExperiment(list(counts=gene.counts))

# Subsetting to all barcodes detected as cells. Requires an intersection,
# because `hto.counts` and `gene.counts` are not the same dimensions! 
common <- intersect(colnames(sce.hash), rownames(hash.stats))
sce.hash <- sce.hash[,common]
colData(sce.hash) <- hash.stats[common,]

sce.hash

# For each cell, we calculate the proportion of its nearest neighbors that are known doublets
# > Intra-sample doublets should have high proportions

# Performing a quick-and-dirty analysis to get some PCs to use
# for nearest neighbor detection inside doubletRecovery().
library(scran)
sce.hash <- logNormCounts(sce.hash)
dec.hash <- modelGeneVar(sce.hash)
top.hash <- getTopHVGs(dec.hash, n = 1000)

set.seed(1011110)
sce.hash <- runPCA(sce.hash,
                   subset_row = top.hash,
                   ncomponents = 20)

# Recovering the intra-sample doublets
hashed.doublets <- doubletRecovery(sce.hash,
                                   use.dimred = "PCA",
                                   doublets = sce.hash$Doublet,
                                   samples = table(sce.hash$Best))

set.seed(1000101001)
sce.hash <- runTSNE(sce.hash, dimred="PCA")
sce.hash$proportion <- hashed.doublets$proportion
sce.hash$predicted <- hashed.doublets$predicted

gridExtra::grid.arrange(
  plotTSNE(sce.hash, colour_by="proportion") + ggtitle("Doublet proportions"),
  plotTSNE(sce.hash, colour_by="Doublet") + ggtitle("Known doublets"),
  ggcells(sce.hash) +
    geom_point(aes(x=TSNE.1, y=TSNE.2), color="grey") +
    geom_point(aes(x=TSNE.1, y=TSNE.2), color="red", 
               data=function(x) x[x$predicted,]) +
    ggtitle("Predicted intra-sample doublets"),
  ncol=2        
)

# worth noting that even known doublets may not necessarily have high doublet neighbor proportions
state <- ifelse(hashed.doublets$predicted, "predicted",
                ifelse(hashed.doublets$known, "known", "singlet"))
ggplot(as.data.frame(hashed.doublets)) + 
  geom_violin(aes(x=state, y=proportion)) 

# Notes:
# Doublet detection procedures should only be applied to libraries generated in the same experimental batch.
# It is obviously impossible for doublets to form between two cells that were captured separately.









