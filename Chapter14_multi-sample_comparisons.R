###############################
# 14.1 - Motivation
###############################
# Multi-sample = many samples of cells from different conditions
# Differential analyses of multi-condition scRNA-seq experiments can be
# * differential expression (DE): cells of the same type in diff. conditions
# * differential abundance (DA) analyses: changes in cell types (or states) between conditions

###############################
# 14.2 - Setting up the data
###############################
#--- loading ---#
library(MouseGastrulationData)
sce.chimera <- WTChimeraData(samples=5:10)
sce.chimera

#--- feature-annotation ---#
library(scater)
rownames(sce.chimera) <- uniquifyFeatureNames(
  rowData(sce.chimera)$ENSEMBL, rowData(sce.chimera)$SYMBOL)

#--- quality-control ---#
drop <- sce.chimera$celltype.mapped %in% c("stripped", "Doublet")
sce.chimera <- sce.chimera[,!drop]

#--- normalization ---#
sce.chimera <- logNormCounts(sce.chimera)

#--- variance-modelling ---#
library(scran)
dec.chimera <- modelGeneVar(sce.chimera, block=sce.chimera$sample)
chosen.hvgs <- dec.chimera$bio > 0

#--- merging ---#
library(batchelor)
set.seed(01001001)
merged <- correctExperiments(sce.chimera, 
                             batch=sce.chimera$sample, 
                             subset.row=chosen.hvgs,
                             PARAM=FastMnnParam(
                               merge.order=list(
                                 list(1,3,5), # WT (3 replicates)
                                 list(2,4,6)  # td-Tomato (3 replicates)
                               )
                             )
)

#--- clustering ---#
g <- buildSNNGraph(merged, use.dimred="corrected")
clusters <- igraph::cluster_louvain(g)
colLabels(merged) <- factor(clusters$membership)

#--- dimensionality-reduction ---#
merged <- runTSNE(merged, dimred="corrected", external_neighbors=TRUE)
merged <- runUMAP(merged, dimred="corrected", external_neighbors=TRUE)

# Inspect batch effects
table(colLabels(merged), merged$tomato)
table(colLabels(merged), merged$pool)

gridExtra::grid.arrange(
  plotTSNE(merged,
           colour_by = "tomato",
           text_by = "label"),
  plotTSNE(merged,
           colour_by = data.frame(pool = factor(merged$pool))),
  ncol = 2
)

# Assign cell type labels
by.label <- table(colLabels(merged), merged$celltype.mapped)
pheatmap::pheatmap(log2(by.label + 1),
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   color = viridis::viridis(101))


###############################
# 14.3 - Differential expression between conditions
###############################
# 14.3.1 - Creating pseudo-bulk samples
# Perform DE for each label and sample separately
summed <- aggregateAcrossCells(merged,
                               id = colData(merged)[, c("celltype.mapped", "sample")])


# 14.3.2 - Performing DE analysis
# 14.3.2.1 - Introduction
# DE analysis will be performed using quasi-likelihood (QL) methods

# Pick randomly one label
label <- "Mesenchyme"
current <- summed[, label == summed$celltype.mapped]

# Create a DGElist object
library(edgeR)
y <- DGEList(counts(current), samples = colData(current))


#14.3.2.2 - Pre-processing
# removing label-sample combinations that have very few or lowly-sequenced cells
# i.e. fewer than 20 cells
discarded <- current$ncells < 20
y <- y[, !discarded]
summary(discarded)

#remove genes that are lowly expressed
keep <- filterByExpr(y, group = current$tomato)
y <- y[keep, ]
summary(keep)

# Correct for composition biases -> scale normalization
y <- calcNormFactors(y) #calculate normalization factors
y$samples


#14.3.2.3 - Statistical modelling
#  set up the design matrix to block on the batch-to-batch differences across different embryo pools
design <- model.matrix(~factor(pool) + factor(tomato),
                       y$samples)

# Estimate negative binomial (NB) dispersions
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
plotBCV(y)

# estimate the quasi-likelihood dispersions
fit <- glmQLFit(y, design, robust = TRUE)
summary(fit$var.prior)
plotQLDisp(fit)

# test for differences in expression due to injection
res <- glmQLFTest(fit,
                  coef = ncol(design))
summary(decideTests(res))

topTags(res)


#14.3.3 - Putting it all together
# repeat this process for each of the labels to identify injection-induced DE in each cell type

# Filter out all sample-label combinations with insufficient cells
summed.filt <- summed[,summed$ncells >= 20]

# Construct a common design matrix that can be used for each label
targets <- colData(merged)[!duplicated(merged$sample),]
design <-  model.matrix(~factor(pool) + factor(tomato), data=targets)
rownames(design) <- targets$sample

# obtain a list of injection-induced DE genes for each label
de.results <- pseudoBulkDGE(summed.filt,
                            sample = summed.filt$sample,
                            label = summed.filt$celltype.mapped,
                            design = design,
                            coef = ncol(design),
                            condition = targets$tomato)

# examine the numbers of DEGs at a FDR of 5% for each label
is.de <- decideTestsPerLabel(de.results, threshold = 0.05)
summarizeTestsPerLabel(is.de)

# Examine upregulated genes across all cells
up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing = TRUE), 10)

# Examine downregulated genes across all cells
down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 10)

#  identify label-specific DE genes that are significant in our label of interest yet not DE in any other label.
remotely.de <- decideTestsPerLabel(de.results, threshold = 0.5)
not.de <- remotely.de==0 | is.na(remotely.de)

other.labels <- setdiff(colnames(not.de), "Allantois")
unique.degs <- is.de[,"Allantois"]!=0 & rowMeans(not.de[,other.labels])==1
unique.degs <- names(which(unique.degs))
unique.degs

# Choosing the top ranked gene for inspection
de.allantois <- de.results$Allantois
de.allantois <- de.allantois[order(de.allantois$PValue),]
de.allantois <- de.allantois[rownames(de.allantois) %in% unique.degs,]

sizeFactors(summed.filt) <- NULL
plotExpression(logNormCounts(summed.filt), 
               features=rownames(de.allantois)[1],
               x="tomato", colour_by="tomato", 
               other_fields="celltype.mapped") + 
  facet_wrap(~celltype.mapped)

# list the labels that were skipped due to the absence of replicates or contrasts
metadata(de.results)$failed



###############################
# 14.4 - Differential abundance between conditions
###############################
# 14.4.1 - Overview
# Quantify the # of cells assigned to each label
abundances <- table(merged$celltype.mapped, merged$sample)
abundances <- unclass(abundances)

# 14.4.2 - Performing the DA analysis
# Counts are here cells per label (not reads per gene)
extra.info <- colData(merged)[match(colnames(abundances),
                                    merged$sample),]

y.ab <- DGEList(abundances, samples = extra.info)

# Filter out low-abundance labels
keep <- filterByExpr(y.ab, group = y.ab$samples$tomato)
y.ab <- y.ab[keep, ]
summary(keep)

# Normalization based on the library size (= total # of cells in each sample)
design <- model.matrix(~factor(pool) + factor(tomato), y.ab$samples)

# Estimate the NB dispersion for each cluster
y.ab <- estimateDisp(y.ab, design, trend = "none")

summary(y.ab$common.dispersion)
plotBCV(y.ab, cex = 1)

# QL dispersion
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
summary(fit.ab$df.prior)
plotQLDisp(fit.ab, cex=1)

# Test for differences in abundance between td-Tomato-positive and negative samples
res <- glmQLFTest(fit.ab, coef = ncol(design))
summary(decideTests(res))
topTags(res)


# 14.4.3 - Handling composition effects
# 1. Assuming most labels (= cell types) do not change
y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples$norm.factors
y.ab2 <- estimateDisp(y.ab2, design, trend="none")
fit.ab2 <- glmQLFit(y.ab2, design, robust=TRUE, abundance.trend=FALSE)
res2 <- glmQLFTest(fit.ab2, coef=ncol(design))
topTags(res2, n=10)

#2. Removing the offending labels
# Another approach is to repeat the analysis after removing DA clusters containing many cells
offenders <- "ExE ectoderm"
y.ab3 <- y.ab[setdiff(rownames(y.ab), offenders),,
              keep.lib.sizes = FALSE]
y.ab3$samples

y.ab3 <- estimateDisp(y.ab3, design, trend="none")
fit.ab3 <- glmQLFit(y.ab3, design, robust=TRUE, abundance.trend=FALSE)
res3 <- glmQLFTest(fit.ab3, coef=ncol(design))
topTags(res3, n=10)

#3. Testing against a log-fold chagne threshold
res.lfc <- glmTreat(fit.ab, coef=ncol(design), lfc=1)
summary(decideTests(res.lfc))

topTags(res.lfc)


###############################
# 14.5 6 - Avoiding problems with ambient RNA
###############################
# Ambient contamination =
# extracellular RNA (most commonly released upon cell lysis) is captured along with each cell
# in its reaction chamber, contributing counts to genes that are not otherwise expressed in that cell 

# load data
library(MouseGastrulationData)
sce.tal1 <- Tal1ChimeraData()

rownames(sce.tal1) <- uniquifyFeatureNames(
  rowData(sce.tal1)$ENSEMBL, 
  rowData(sce.tal1)$SYMBOL
)
sce.tal1

# perform DE analysis with cells labeled as "neural crest"
summed.tal1 <- aggregateAcrossCells(sce.tal1, 
                                    ids=DataFrame(sample=sce.tal1$sample,
                                                  label=sce.tal1$celltype.mapped)
)
summed.neural <- summed.tal1[,summed.tal1$label=="Neural crest"]
summed.neural

# Standard edgeR analysis, as described above.
y.neural <- DGEList(counts(summed.neural), samples=colData(summed.neural))
keep.neural <- filterByExpr(y.neural, group=y.neural$samples$tomato)
y.neural <- y.neural[keep.neural,]
y.neural <- calcNormFactors(y.neural)

block <- y.neural$samples$sample %% 2 == 0
design <- model.matrix(~factor(block) + factor(tomato), y.neural$samples)
y.neural <- estimateDisp(y.neural, design)
fit.neural <- glmQLFit(y.neural, design, robust=TRUE)
res.neural <- glmQLFTest(fit.neural, coef=ncol(design))
summary(decideTests(res.neural))

topTags(res.neural, n=10)

# 14.6.2 - Discarding ambient DEGs
library(DropletUtils)
ambient <- vector("list", ncol(summed.neural))

# Looping over all raw (unfiltered) count matrices and
# computing the ambient profile based on its low-count barcodes.
# Turning off rounding, as we know this is count data.
for (s in seq_along(ambient)) {
  raw.tal1 <- Tal1ChimeraData(type="raw", samples=s)[[1]]
  ambient[[s]] <- estimateAmbience(counts(raw.tal1), 
                                   good.turing=FALSE, round=FALSE)
}

# Cleaning up the output for pretty printing.
ambient <- do.call(cbind, ambient)
colnames(ambient) <- seq_len(ncol(ambient))
rownames(ambient) <- uniquifyFeatureNames(
  rowData(raw.tal1)$ENSEMBL, 
  rowData(raw.tal1)$SYMBOL
)
head(ambient)

# Looping over all samples and computing the maximum proportion of 
# counts explained by ambient contamination in each sample.
max.ambient <- list()
for (i in seq_len(ncol(ambient))) {
  max.ambient[[i]] <- maximumAmbience(counts(summed.neural)[,i], 
                                      ambient[,i], mode="proportion")
}

max.ambient <- do.call(cbind, max.ambient)
dimnames(max.ambient) <- dimnames(ambient)
head(max.ambient)

non.ambient <- rowMeans(max.ambient, na.rm=TRUE) <= 0.1
summary(non.ambient)

okay.genes <- names(non.ambient)[which(non.ambient)]
res.neural2 <- res.neural[rownames(res.neural) %in% okay.genes,]
summary(decideTests(res.neural2))

topTags(res.neural2)

# Use of prior knowledge of mutually exlcusive gene expression profiles
is.hbb <- grep("^Hb[ab]-", rownames(summed.neural))
neural.hb <- colSums(counts(summed.neural)[is.hbb,])
ambient.hb <- colSums(ambient[is.hbb,])
scaled.ambient <- t(t(ambient) * neural.hb/ambient.hb)
head(scaled.ambient)

alt.prop <- scaled.ambient/counts(summed.neural)
alt.prop[!is.finite(alt.prop)] <- NA
alt.non.ambient <- rowMeans(alt.prop, na.rm=TRUE) <= 0.1
summary(alt.non.ambient)


# 14.6.3 Subtracting ambient counts
subtracted <- counts(summed.neural) - scaled.ambient
subtracted <- round(subtracted)
subtracted[subtracted < 0] <- 0
subtracted[is.hbb,]

# Another tempting approach is to use interaction models to implicitly subtract the ambient effect during GLM fitting. 
# Re-using keep.neural to simplify comparison.
y.ambient <- DGEList(ambient)
y.ambient <- y.ambient[keep.neural,]
y.ambient <- calcNormFactors(y.ambient)
y.ambient <- estimateDisp(y.ambient, design)
fit.ambient <- glmQLFit(y.ambient, design, robust=TRUE)
res.ambient <- glmQLFTest(fit.ambient, coef=ncol(design))
summary(decideTests(res.ambient))

topTags(res.ambient, n=10)
