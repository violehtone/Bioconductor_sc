################################
# 36.1 - Introduction
################################
# - Analysis of the Bach et al. (2017) 10X genomics dataset
# - single sample of epithelial cells from the mouse mammary gland

################################
# 36.2 - Data loading
################################
library(scRNAseq)
library(scater)
library(AnnotationHub)

# Obtain the mouse mammary gland scRNA-seq data
sce.mam <- BachMammaryData(samples="G_1")

# Combine the feature name with a standard identifier
rownames(sce.mam) <- uniquifyFeatureNames(
  rowData(sce.mam)$Ensembl, rowData(sce.mam)$Symbol)

# Use bioconductor's annotation service to annotate the data
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.mam)$SEQNAME <- mapIds(ens.mm.v97, 
                                   keys=rowData(sce.mam)$Ensembl,
                                   keytype="GENEID", 
                                   column="SEQNAME")


################################
# 36.3 - Quality control
################################
unfiltered <- sce.mam

# Remove where mitochondrial content is high
is.mito <- rowData(sce.mam)$SEQNAME == "MT"
stats <- perCellQCMetrics(sce.mam, subsets=list(Mito=which(is.mito)))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent")
sce.mam <- sce.mam[,!qc$discard]

# Combine unfiltered data with perCellQCMetrics
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

# plot sum, detected, and mito%
gridExtra::grid.arrange(
  plotColData(unfiltered, y="sum", colour_by="discard") + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, y="detected", colour_by="discard") + 
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, y="subsets_Mito_percent", 
              colour_by="discard") + ggtitle("Mito percent"),
  ncol=2
)

# plot mito% against total count
plotColData(unfiltered, x="sum", y="subsets_Mito_percent", 
            colour_by="discard") + scale_x_log10()

# Inspect discarded cells
colSums(as.matrix(qc))


################################
# 36.4 - Normalization
################################
library(scran)
set.seed(101000110)

# clustering, scaling, and log normalization
clusters <- quickCluster(sce.mam)
sce.mam <- computeSumFactors(sce.mam, clusters=clusters)
sce.mam <- logNormCounts(sce.mam)

# plot library size factors vs. size factors
plot(librarySizeFactors(sce.mam), sizeFactors(sce.mam), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy")


################################
# 36.5 - Variance modelling
################################
set.seed(00010101)
dec.mam <- modelGeneVarByPoisson(sce.mam)
top.mam <- getTopHVGs(dec.mam, prop=0.1)

################################
# 36.6 - Dimensionality reduction
################################
library(BiocSingular)
set.seed(101010011)
sce.mam <- denoisePCA(sce.mam, technical=dec.mam, subset.row=top.mam)
sce.mam <- runTSNE(sce.mam, dimred="PCA")


################################
# 36.7 - Clustering
################################
snn.gr <- buildSNNGraph(sce.mam, use.dimred="PCA", k=25)
colLabels(sce.mam) <- factor(igraph::cluster_walktrap(snn.gr)$membership)

# Inspect clustering
table(colLabels(sce.mam))
plotTSNE(sce.mam, colour_by="label")

