#########################
# 16.1 - Motivation
#########################

# Load data set
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

#--- clustering ---#
my.dist <- dist(reducedDim(sce.416b, "PCA"))
my.tree <- hclust(my.dist, method="ward.D2")

library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist),
                                    minClusterSize=10, verbose=0))
colLabels(sce.416b) <- factor(my.clusters)

#########################
# 16.2 - Using the cyclins
#########################
# * Cyclin is a family of proteins that controls the progression of a cell through the cell cycle by activating cyclin-dependent kinase (CDK) enzymes or group of enzymes required for synthesis of cell cycle.
# The cyclins control progression through the cell cycle and have well-characterized patterns of expression across cell cycle phases
# - Cyclin D = peaks at G1
# - Cyclin E = G1/S transition
# - Cyclin A = S and G2
# - Cyclin B = late G2 and mitosis
# -->  Inspection of the relative expression of cyclins across the population can often be sufficient to determine the relative cell cycle activity in each cluster

cyclin.genes <- grep("^Ccn[abde][0-9]$", rowData(sce.416b)$SYMBOL)
cyclin.genes <- rownames(sce.416b)[cyclin.genes]

plotHeatmap(sce.416b,
            order_columns_by = "label",
            cluster_rows = FALSE,
            features = sort(cyclin.genes))


# we can use standard DE methods to look for upregulation of each cyclin
markers <- findMarkers(sce.416b,
                       subset.row = cyclin.genes,
                       test.type = "wilcox",
                       direction = "up")

markers[[4]]


#########################
# 16.3 - Using reference profiles
#########################
# Cell cycle assignment can be considered a specialized case of cell annotation
# given a reference dataset  with known cell cycle phases, we can determine the phase of each cell in a test dataset.

# Load data set
sce.ref <- BuettnerESCData()

# Find genes that are present in both datasets and are cell cycle-related.
library(org.Mm.eg.db)
cycle.anno <- select(org.Mm.eg.db,
                     keytype = "GOALL",
                     keys = "GO:0007049",
                     columns = "ENSEMBL")[, "ENSEMBL"]

candidates <- Reduce(intersect,
                     list(rownames(sce.ref),
                          rowData(sce.416b)$ENSEMBL,
                          cycle.anno))

str(candidates)

# Identifying markers between cell cycle phases.
sce.ref <- logNormCounts(sce.ref)
phase.stats <- pairwiseWilcox(logcounts(sce.ref),
                              sce.ref$phase,
                              direction = "up",
                              subset.row = candidates)

cycle.markers <- getTopMarkers(phase.stats[[1]],
                               phase.stats[[2]])

# Switching row names back to Ensembl to match the reference.
test.data <- logcounts(sce.416b)
rownames(test.data) <- rowData(sce.416b)$ENSEMBL

library(SingleR)
assignments <- SingleR(test.data, ref=sce.ref,
                       label=sce.ref$phase, genes=cycle.markers)
tab <- table(assignments$labels, colLabels(sce.416b))
tab

# Lef1 is detected as one of the top markers to distinguish between G1 from G2/M in the reference but has no detectable expression in the 416B dataset
gridExtra::grid.arrange(
  plotExpression(sce.ref, features="ENSMUSG00000027985", x="phase"),
  plotExpression(sce.416b, features="Lef1", x="label"),
  ncol=2)

# Interpret frequencies in a relative sense: chi-squared test
chisq.test(tab[, 1:2])



#########################
# 16.4 - Using the cyclone() classifier
#########################
# prediction method: 
# - Using a reference dataset, we first compute the sign of the difference in expression between each pair of genes
# - Pairs with changes in the sign across cell cycle phases are chosen as markers
# - Cells in a test dataset can then be classified into the appropriate phase, based on whether the observed sign for each marker pair is consistent with one phase or another.

set.seed(100)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds",
                                package = "scran"))

# Classify cells into their cell cycle phases
assignments <- cyclone(sce.416b,
                       mm.pairs,
                       gene.names = rowData(sce.416b)$ENSEMBL)


# Make a plot of the cell cycle phase G1 and G2/M scores
# > For each cell, a higher score for a phase corresponds to a higher probability that the cell is in that phase
plot(assignments$score$G1, assignments$score$G2M,
     xlab="G1 score", ylab="G2/M score", pch=16)

# .. results compared to ones from SingleR
table(assignments$phases, colLabels(sce.416b))


#########################
# 16.5 - Regressing out cell cycle phase
#########################
# aim is to remove uninteresting variation due to cell cycle, thus improving resolution of other biological processes of interest.
# The most common approach is to use a linear model to simply regress out the phase effect, e.g., via regressBatches()

library(batchelor)
sce.nocycle <- regressBatches(sce.416b,
                              batch = assignments$phases)

dec.nocycle <- modelGeneVarWithSpikes(sce.416b, "ERCC", 
                                      block=assignments$phases)
marker.nocycle <- findMarkers(sce.416b, block=assignments$phases)


# NOTES:
# * we do not consider adjusting for cell cycle to be a necessary step in routine scRNA-seq analyses
# * the cell cycle is a minor factor of variation
# * We suggest only performing cell cycle adjustment on an as-needed basis in populations with clear cell cycle effects.