##############################
# 21.1 - Motivation
##############################
# - Number of cells that can be assayed in a routine experiment can be huge
# - how to tune our analysis pipeline for greater speed and efficiency?

##############################
# 21.2 - Fast approximations
##############################
# Load dataset
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir=file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)

#--- gene-annotation ---#
library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)

library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID")

#--- cell-detection ---#
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

#--- quality-control ---#
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.pbmc <- sce.pbmc[,!high.mito]

#--- normalization ---#
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

#--- variance-modelling ---#
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)

#--- dimensionality-reduction ---#
set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc, subset.row=top.pbmc, technical=dec.pbmc)

set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred="PCA")

set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred="PCA")

#--- clustering ---#
g <- buildSNNGraph(sce.pbmc, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.pbmc) <- factor(clust)


# Clustering using an approximate search (Annoy algorithm)
library(BiocNeighbors)
snn.gr <- buildSNNGraph(sce.pbmc,
                        BNPARAM = AnnoyParam(),
                        use.dimred = 'PCA')

clusters <- igraph::cluster_walktrap(snn.gr)
table(Exact=colLabels(sce.pbmc), Approx=clusters$membership)

set.seed(1000)
y1 <- matrix(rnorm(50000), nrow=1000)
y2 <- matrix(rnorm(50000), nrow=1000)
Y <- rbind(y1, y2)
exact <- findKNN(Y, k=20)
approx <- findKNN(Y, k=20, BNPARAM=AnnoyParam())
mean(exact$index!=approx$index)


# 21.2.2 - Singular value decomposition (SVD)
# SVD is a factorization of a matrix. PCA uses SVD in its calculation
# - base::svd() function performs an exact SVD
# - irlba and rsvd packages offer faster approximations of SVD

set.seed(101000)
r.out <- runPCA(sce.pbmc,
                ncomponents = 20,
                BSPARAM = RandomParam())

set.seed(101001)
i.out <- runPCA(sce.pbmc, 
                ncomponents=20, 
                BSPARAM=IrlbaParam())



##############################
# 21.3 - Parallelization
##############################
# Parallelization can be used to speed up scRNA-seq calculations
# We can pick from a diverse range of parallelization backends depending on the available hardware and operating system
library(BiocParallel)

# Use forking across 2 cores to parallelize the variance calculations
dec.pbmc.mc <- modelGeneVar(sce.pbmc,
                            BPPARAM = MulticoreParam(2))


# Alternative: distribute jobs across a network of computers
dec.pbmc.snow <- modelGeneVar(sce.pbmc,
                              BPPARAM = SnowParam(5))


# Alternative: Distribute jobs via the job scheduler
bpp <- BatchtoolsParam(10, cluster = "slurm",
                       resources = list(walltime = 7200,
                                        memory = 8000,
                                        ncpus = 1))

# NOTE: Parallelization is best suited for CPU-intensive calculations where the division of labor results in a concomitant reduction in compute time


##############################
# 21.4 - Out of memory representations
##############################
# - in memory count matrix may not be feasible for very large datasets
# - i.e. the 1.3 million brain cell data set from would require over 100 GB of RAM to hold as a matrix
# Solution? -> use a file-backed matrix representation
#   * data are held on disk
#   * subsets are retrieved into memory as requested

# Load test data
library(TENxBrainData)
sce.brain <- TENxBrainData20k() 
class(sce.brain)

tmp <- counts(sce.brain)
tmp <- log2(tmp + 1)
class(tmp)

# Compute QC metrics
is.mito <- grepl("^mt-", rowData(sce.brain)$Symbol)
qcstats <- perCellQCMetrics(sce.brain, subsets = list(Mt=is.mito))













