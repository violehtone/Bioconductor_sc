#################
# 11.1 Motivation
#################

# Marker genes = genes that drive separation between clusters (~cell types)
# genes that are more strongly DE are more likely to have caused separate clustering of cells

# Load data set
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir=file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)

# Gene annotation
library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)

library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID")

# Cell detection
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

# Normalization
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

# Variance modeling
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)

# dimensionality reduction
set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc, subset.row=top.pbmc, technical=dec.pbmc)

set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred="PCA")

set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred="PCA")

# Clustering
g <- buildSNNGraph(sce.pbmc, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.pbmc) <- factor(clust)

# Inspect dataset
sce.pbmc

#################
# 11.2 - Using pairwise t-tests
#################

# 11.2.1 - Standard application

# findMarkers() performs pairwise comparison between clusters for each gene,
# which returns a list of dataframes containing ranked candidate markers for each cluster
library(scran)
markers.pbmc <- findMarkers(sce.pbmc)

# Alternative way of producing the same result as markers.pbmc
same.markers <- findMarkers(sce.pbmc, groups = colLabels(sce.pbmc))

# Inspect markers in cluster 9
chosen <- "9"
interesting <- markers.pbmc[[chosen]]

colnames(interesting)

# Top values of 1 contains the gene with the lowest p-value from each comparison
interesting[1:10, 1:4]

# Top field can be used to identify a set of genes that distinguis cluster 9 from other clusters
# Here we examine top 6 marker genes
best.set <- interesting[interesting$Top <= 6, ]
logFCs <- getMarkerEffects(best.set)

library(pheatmap)
pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))


# 11.2.2 - Using the log-fold change

# downregulated genes are less appealing as markers as it
# is more difficult to interpret and experimentally validate an absence of expression
# -> To focus on up-regulated markers, we can instead perform a one-sided t-test 

# direction = "up" is used to find up-regulated genes
markers.pbmc.up <- findMarkers(sce.pbmc, direction = "up")
interesting.up <- markers.pbmc.up[[chosen]]
interesting.up[1:10,1:4]

# lfc = 1 and direction = "up" -> find genes with log-fold changes that are significantly greater than 1
markers.pbmc.up2 <- findMarkers(sce.pbmc,
                                direction = "up",
                                lfc = 1)

interesting.up2 <- markers.pbmc.up2[[chosen]]
interesting.up2[1:10,1:4]

best.set <- interesting.up2[interesting.up2$Top <= 5,]
logFCs <- getMarkerEffects(best.set)
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))


# 11.2.3 - Finding cluster-specific markers
# consider genes that are differentially expressed in all pariwsie comparisons
# -> pval.type = "all"
# -> reports only genes that are highly specific to the cluster of interest
markers.pbmc.up3 <- findMarkers(sce.pbmc,
                                pval.type="all",
                                direction="up")

interesting.up3 <- markers.pbmc.up3[[chosen]]
interesting.up3[1:10,1:3]



# pval.type="some" -> at least 50% of the individual pairwise comparisons exhibit no DE
markers.pbmc.up4 <- findMarkers(sce.pbmc, pval.type="some", direction="up")
interesting.up4 <- markers.pbmc.up4[[chosen]]
interesting.up4[1:10,1:3]


#################
# 11.3 - Alternative testing regimes
#################

# 11.3.1 - Using the Wilcoxon rank sum test (= Wilcoxon-Mann-Whitney test)
# Alternative for Welch t-test pairwise comparison between groups of observations

# The WMW test statistic is proportional to the area-under-the-curve (AUC),
# In a pairwise comparison, AUCs of 1 or 0 indicate that the two clusters have
# perfectly separated expression distributions.

markers.pbmc.wmw <- findMarkers(sce.pbmc, 
                                test="wilcox", 
                                direction="up")
names(markers.pbmc.wmw)
interesting.wmw <- markers.pbmc.wmw[[chosen]]
interesting.wmw[1:10,1:4]

# AUC value greater than 0.5 indicates that the gene is upregulated in the current
# cluster compared to the other cluster, while values less than 0.5 correspond to
# downregulation

best.set <- interesting.wmw[interesting.wmw$Top <= 5,]
AUCs <- getMarkerEffects(best.set, prefix="AUC")

pheatmap(AUCs, 
         breaks=seq(0, 1, length.out=21),
         color=viridis::viridis(21))

# The main disadvantage of the WMW test is that the AUCs are much slower 
# to compute compared to t-statistics.


# 11.3.2 - Using a binomial test
# The binomial test identifies genes that differ in the proportion of expressing cells between clusters

markers.pbmc.binom <- findMarkers(sce.pbmc, test="binom", direction="up")
names(markers.pbmc.binom)

interesting.binom <- markers.pbmc.binom[[chosen]]
colnames(interesting.binom)

top.genes <- head(rownames(interesting.binom))
plotExpression(sce.pbmc, x="label", features=top.genes)


# 11.3.3 - Using custom DE methods
# i.e. voom() from Limma

library(limma)
design <- model.matrix(~0 + label, data=colData(sce.pbmc))
colnames(design)

# Removing very low-abundance genes.
keep <- calculateAverage(sce.pbmc) > 0.1 
summary(keep)

y <- convertTo(sce.pbmc, subset.row=keep)
# Voom() is used to prepare data for linear modeling
v <- voom(y, design) # transform count data to log2 counts per million (logCPM)
# fit a linear model
fit <- lmFit(v, design)


# 11.3.4 - Combining multiple marker statistics
# i.e. merge results from t-test, wilcoxon rank-test, and binomial test
# This allows us to easily inspect multiple statistics at once to verify that a particular gene is a strong candidate marker

combined <- multiMarkerStats(
  t=findMarkers(sce.pbmc, direction="up"),
  wilcox=findMarkers(sce.pbmc, test="wilcox", direction="up"),
  binom=findMarkers(sce.pbmc, test="binom", direction="up")
)

# Interleaved marker statistics from both tests for each cluster.
colnames(combined[["1"]])


#################
# 11.4 - Handling blocking factors 
#################
# Using the block -argument
# > block can be used to ignore  factors of variation that are known and not interesting
# i.e. m.out <- findMarkers(sce.416b, block=sce.416b$block, direction="up") 

# Using the design -argument
# situations where multiple batches contain unique clusters, 
#as comparisons can be implicitly performed via shared cell types in each batch.

# Setting up the design matrix (we remove intercept for full rank
# in the final design matrix with the cluster-specific terms).
#i.e.:
# design <- model.matrix(~sce.416b$block)
# design <- design[,-1,drop=FALSE]
# 
# m.alt <- findMarkers(sce.416b, design=design, direction="up")
# demo <- m.alt[["1"]]
# demo[demo$Top <= 5,1:4]


#################
# 11.4 - Unvalidity of p-values
#################
# The DE analysis is performed on the same data used to obtain the 
# clusters, which represents “data dredging” (also known as fishing or data snooping)
