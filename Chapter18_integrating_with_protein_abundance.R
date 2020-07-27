########################################
# 18.1 - Motivation
########################################
# = study of proteome with transcriptome data

# Cellular indexing of transcriptomes and epitopes by sequencing (CITE-seq)
# is a technique that quantifies both gene expression and the abundance of selected
# surface proteins in each cell simultaneously

# 1) cells are labelled with antibodies (RNA tags)
# 2) cells are separated into their own reaction chambers using droplet-based microfluidics
# 3) abundance of each protein is quantified by sequencing each set of features
# --> provides a powerful tool for interrogating aspects of the proteome

# * ADT = antibody derived tag
# here ==> strategies for integrated analysis of ADT and transcript data in CITE-seq experiments!

# Load data
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
stuff <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com",
                                 "samples/cell-exp/3.0.0/pbmc_10k_protein_v3",
                                 "pbmc_10k_protein_v3_filtered_feature_bc_matrix.tar.gz"))
untar(stuff, exdir=tempdir())

# Loading it in as a SingleCellExperiment object.
library(DropletUtils)
sce <- read10xCounts(file.path(tempdir(), "filtered_feature_bc_matrix"))
sce


########################################
# 18.2 - Preprocessing
########################################
# 18.2.1 - Setting up the data
# “alternative Experiment” can be used to store data for different sets of features but the same cells

sce <- splitAltExps(sce,
                    rowData(sce)$Type)
altExpNames(sce)
altExp(sce)

# coerce the sparse matrix for ADTs into a dense matrix
counts(altExp(sce)) <- as.matrix(counts(altExp(sce)))
counts(altExp(sce))[, 1:10]


# 18.2.2 - Quality control
# remove empty droplets and low-quality cells
#  filter on the mitochondrial proportions to remove putative low-quality cells

mito <- grep("^MT-", rowData(sce)$Symbol)
df <- perCellQCMetrics(sce, subsets = list(Mito = mito))
mito.discard <- isOutlier(df$subsets_Mito_percent, type = "higher")
summary(mito.discard)

# only retain cells that actually have ADT counts
# ->  remove cells that have unusually low numbers of detected ADTs
ab.discard <- isOutlier(df$`altexps_Antibody Capture_detected`,
                        log=TRUE, 
                        type="lower", 
                        min_diff=1)

summary(ab.discard)

# make a histogram of the # of detected ADTs across all cells
hist(df$`altexps_Antibody Capture_detected`, col='grey', 
     main="", xlab="Number of detected ADTs")
abline(v=attr(ab.discard, "thresholds")["lower"], col="red", lty=2)

# to remove the low-quality cells, we subset the SingleCellExperiment
discard <- ab.discard | mito.discard
sce <- sce[,!discard]



# 18.2.3 - Normalization
# - Capture efficiency varies from cell to cell though the differences in biophysical properties between endogenous transcripts and the (much shorter) ADTs
# - Composition biases are also much more pronounced in ADT data
# ==> normalize on the total ADT counts, effectively the library size for the ADTs
sf.lib <- librarySizeFactors(altExp(sce))
summary(sf.lib)

# Alternatively: taking the geometric mean of all counts as the size factor for each cell
sf.geo <- librarySizeFactors(altExp(sce), geometric=TRUE)
summary(sf.geo)

#  obtain an estimate of the ambient profile from the barcodes that were identified as empty droplets
ambient <- rowMeans(counts(altExp(sce)))
sf.amb <- medianSizeFactors(altExp(sce), reference=ambient)
summary(sf.amb)

# Plot DESeq-like size factors for each cell vs. ADT library size factors
tagdata <- logNormCounts(altExp(sce)) # library size factors by default.
g <- buildSNNGraph(tagdata, k=20, d=NA) # no need for PCA, see below.
clusters <- igraph::cluster_walktrap(g)$membership

plot(sf.lib, sf.amb, log="xy", col=clusters, 
     xlab="Library size factors (tag)",
     ylab="DESeq-like size factors (tag)")
abline(0, 1, col="grey", lty=2)

# computing size factors from the immunoglobulin (IgG) controls
controls <- grep("^Ig", rownames(altExp(sce)))
sf.control <- librarySizeFactors(altExp(sce), subset_row=controls) 
summary(sf.control)

plot(sf.amb, sf.control, log="xy", 
     xlab="DESeq-like size factors (tag)",
     ylab="Control size factors (tag)")
abline(0, 1, col="grey", lty=2)

#  scaling normalization and log-transformation
sizeFactors(altExp(sce)) <- sf.amb
sce <- logNormCounts(sce, use_altexps=TRUE)


########################################
# 18.3 - Clustering and interpretation
########################################
# feature selection is largely unnecessary for analyzing ADT data
# -> we should directly apply downstream procedures like clustering and visualization
g.adt <- buildSNNGraph(altExp(sce),
                       d = NA)

clusters.adt <- igraph::cluster_walktrap(g.adt)$membership

# generate a t-SNE plot
set.seed(1010010)
altExp(sce) <- runTSNE(altExp(sce))
colLabels(altExp(sce)) <- factor(clusters.adt)
plotTSNE(altExp(sce), colour_by="label", text_by="label", text_col="red")

# With only a few ADTs, characterization of each cluster is most efficiently achieved by
# creating a heatmap of the average log-abundance of each tag 
se.averaged <- sumCountsAcrossCells(altExp(sce), clusters.adt,
                                    exprs_values="logcounts", average=TRUE)

library(pheatmap)
averaged <- assay(se.averaged)
pheatmap(averaged - rowMeans(averaged),
         breaks=seq(-3, 3, length.out=101))




########################################
# 18.4 - Integrating with gene expression data
########################################
# we take cells in each of the ADT-derived clusters and perform subclustering using the transcript data
# -> using quickSubCluster() to loop over all of the ADT-derived clusters and subcluster on gene expression
set.seed(101010)
all.sce <- quickSubCluster(sce, clusters.adt,
                           prepFUN=function(x) {
                             dec <- modelGeneVar(x)
                             top <- getTopHVGs(dec, prop=0.1)
                             x <- runPCA(x, subset_row=top, ncomponents=25)
                           },
                           clusterFUN=function(x) {
                             g.trans <- buildSNNGraph(x, use.dimred="PCA")
                             igraph::cluster_walktrap(g.trans)$membership
                           }
)

ncells <- vapply(all.sce, ncol, 0L)
nsubclusters <- vapply(all.sce, FUN=function(x) length(unique(x$subcluster)), 0L)
plot(ncells, nsubclusters, xlab="Number of cells", type="n",
     ylab="Number of subclusters", log="xy")
text(ncells, nsubclusters, names(all.sce))

#  identify internal subclusters based on granzyme expression
of.interest <- "12"
plotExpression(all.sce[[of.interest]], x="subcluster",
               features=c("ENSG00000100450", "ENSG00000113088"))


# perform some additional checks to ensure that each subcluster has similar protein abundances
sce.cd8 <- all.sce[[of.interest]]
plotExpression(altExp(sce.cd8), x=I(sce.cd8$subcluster),
               features=c("CD3", "CD8a"))


# 18.4.2 - By combined clustering
# > Alternatively, we can combine the information from both sets of features into a single matrix for use in downstream analyses

# perform some standard steps on the transcript count matrix
sce.main <- logNormCounts(sce)
dec.main <- modelGeneVar(sce.main)
top.main <- getTopHVGs(dec.main, prop=0.1)
sce.main <- runPCA(sce.main, subset_row=top.main, ncomponents=25)


# combining the log-normalized abundance matrix for the ADTs
# with the log-expression matrix to obtain a single matrix for use in downstream analyses
library(DelayedMatrixStats)
transcript.data <- logcounts(sce.main)[top.main,,drop=FALSE]
transcript.var <- sum(rowVars(DelayedArray(transcript.data)))
tag.data <- logcounts(altExp(sce.main))
tag.var <- sum(rowVars(DelayedArray(tag.data)))

reweight <- sqrt(transcript.var/tag.var)
combined <- rbind(transcript.data, tag.data*reweight)

# build SNN graph
set.seed(100010)
g.com <- buildSNNGraph(combined, d=50) 
clusters.com <- igraph::cluster_walktrap(g.com)$membership
table(clusters.com)


# Alternatively: use UMAP to integrate information from two sets of features
set.seed(1001010)
combined2 <- runMultiUMAP(
  list(reducedDim(sce.main, "PCA"), 
       t(logcounts(altExp(sce.main)))),
  n_components=20, n_neighbors=30, min_dist=0
)

g.com2 <- buildSNNGraph(combined2, d=NA, transposed=TRUE)
clusters.com2 <- igraph::cluster_walktrap(g.com2)$membership
table(clusters.com2)

# visualize results
set.seed(0101110)
reducedDim(sce.main, "combinedUMAP") <- runMultiUMAP(
  list(reducedDim(sce.main, "PCA"), 
       t(logcounts(altExp(sce.main))))
)
colLabels(sce.main) <- clusters.com2
plotReducedDim(sce.main, "combinedUMAP", 
               colour_by="label", text_by="label")


# 18.4.3 - By differential testing
# > protein targets are chosen that reflect some functional activity rather than cell type
# > use the transcript data for clustering and perform differential testing between clusters or conditions for the relevant ADTs.

# Performing a quick analysis of the gene expression data.
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
top <- getTopHVGs(dec, prop=0.1)

set.seed(1001010)
sce <- runPCA(sce, subset_row=top, ncomponents=25)

g <- buildSNNGraph(sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(g)$membership
colLabels(sce) <- factor(clusters)

set.seed(1000010)
sce <- runTSNE(sce, dimred="PCA")
plotTSNE(sce, colour_by="label", text_by="label")

# test for differences in tag abundance between clusters
markers <- findMarkers(altExp(sce), colLabels(sce))
of.interest <- markers[[16]]
pheatmap(getMarkerEffects(of.interest), breaks=seq(-3, 3, length.out=101))
