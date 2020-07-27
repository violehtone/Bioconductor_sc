############################
# 19.1 - Motivation
############################
# An organism’s immune repertoire is defined as:
# "the set of T and B cell subtypes that contain genetic diversity in the T cell receptor (TCR) components or immunoglobin chains, respectively"

# > We can profile the immune repertoire by simply sequencing the relevant transcripts
# > This data can then be used to characterize an individual’s immune response based on the expansion of T or B cell clones

# single-cell repertoire sequencing data can be readily analyzed using tools from the ImmCantation suite

# Load data
#--- loading ---#
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
exprs.data <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples/cell-vdj/3.1.0",
  "vdj_v1_hs_pbmc3",
  "vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz"))
untar(exprs.data, exdir=tempdir())

library(DropletUtils)
sce.pbmc <- read10xCounts(file.path(tempdir(), "filtered_feature_bc_matrix"))
sce.pbmc <- splitAltExps(sce.pbmc, rowData(sce.pbmc)$Type)

#--- quality-control ---#
library(scater)
is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=is.mito))

high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
low.adt <- stats$`altexps_Antibody Capture_detected` < nrow(altExp(sce.pbmc))/2

discard <- high.mito | low.adt
sce.pbmc <- sce.pbmc[,!discard]

#--- normalization ---#
library(scran)

set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
altExp(sce.pbmc) <- computeMedianFactors(altExp(sce.pbmc))
sce.pbmc <- logNormCounts(sce.pbmc, use_altexps=TRUE)

#--- dimensionality-reduction ---#
set.seed(100000)
altExp(sce.pbmc) <- runTSNE(altExp(sce.pbmc))

set.seed(1000000)
altExp(sce.pbmc) <- runUMAP(altExp(sce.pbmc))

#--- clustering ---#
g.adt <- buildSNNGraph(altExp(sce.pbmc), k=10, d=NA)
clust.adt <- igraph::cluster_walktrap(g.adt)$membership
colLabels(altExp(sce.pbmc)) <- factor(clust.adt)

colLabels(sce.pbmc) <- colLabels(altExp(sce.pbmc))


#-> goal: define a single data structure that captures both the expression profile and repertoire state for each cell.

############################
# 19.2 - Analyzing the T cell receptor repertoire
############################
# 19.2.1 - Data processing
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
tcr.data <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples/cell-vdj/3.1.0",
  "vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_t_filtered_contig_annotations.csv"))
tcr <- read.csv(tcr.data, stringsAsFactors=FALSE)

# store the repertoire data as a SplitDataFrameList object 
tra <- tcr[tcr$chain=="TRA",]
sce.pbmc$TRA <- split(DataFrame(tra), factor(tra$barcode, sce.pbmc$Barcode))
length(sce.pbmc$TRA) # Now the same as the number of cells.


trb <- tcr[tcr$chain=="TRB",]
sce.pbmc$TRB <- split(DataFrame(trb), factor(trb$barcode, sce.pbmc$Barcode))
length(sce.pbmc$TRB)


# 19.2.2 - Basic diagnostics
# determine the proportion of cells that have at least one sequence of a TCR component
ncells <- table(colLabels(sce.pbmc))

at.least.one.A <- lengths(sce.pbmc$TRA) > 0
tra.counts.any <- table(colLabels(sce.pbmc)[at.least.one.A])

at.least.one.B <- lengths(sce.pbmc$TRB) > 0
trb.counts.any <- table(colLabels(sce.pbmc)[at.least.one.B])

# plot proportion of cells in each cluster that express at least one seq of the TCR a or b chains
barplot(rbind(TRA=tra.counts.any/ncells,
              TRB=trb.counts.any/ncells),
        beside=TRUE)

# only consider the productive sequences, i.e., contigs that are likely to produce a functional protein
is.prod.A <- sce.pbmc$TRA[,"productive"]=="True" 
is.prod.A

has.prod.A <- any(is.prod.A)

# count the number of cells in each cluster.
tra.counts.prod <- table(colLabels(sce.pbmc)[has.prod.A])

is.prod.B <- sce.pbmc$TRB[,"productive"]=="True"
has.prod.B <- any(is.prod.B)
trb.counts.prod <- table(colLabels(sce.pbmc)[has.prod.B])

barplot(rbind(TRA=tra.counts.prod/ncells,
              TRB=trb.counts.prod/ncells),
        legend=TRUE,
        beside=TRUE)

# count the number of cells in each cluster that have multiple sequences for a component
tra.counts.multi <- table(colLabels(sce.pbmc)[lengths(sce.pbmc$TRA) > 1])
trb.counts.multi <- table(colLabels(sce.pbmc)[lengths(sce.pbmc$TRB) > 1])
barplot(rbind(TRA=tra.counts.multi/ncells,
              TRB=trb.counts.multi/ncells), 
        beside=TRUE)


# 19.2.3 - complex diagnostics
# use the α-chain data to extract some complex features
# -> identification of seq. that have UMI counts >= 50% of the largest UMI count for the same cell

tra <- sce.pbmc$TRA
max.umi <- max(tra[, "umis"])
keep <- tra[, "umis"] >= max.umi / 2

# -> identify sequences that are full length, productive and have the largest UMI count in the cell
keep <- tra[,"full_length"]=="True" &
  tra[,"productive"]=="True" &
  tra[,"umis"] == max(tra[,"umis"])

tra.sub <- tra[keep]

# Quantify all combinations of V and J genes
combined <- paste(tra[, "v_gene"],
                  tra[, "j_gene"])

combo.freq <- table(unlist(combined))

# Reover the original seq-level data frame
tra.seq <- unlist(tra)

# Add extra annotations
extra.anno <- DataFrame(anno = sample(LETTERS,
                                      nrow(tra.seq),
                                      replace = TRUE))

tra.seq <- cbind(tra.seq, extra.anno)

# Regenerate the SplitDataFrameList form the df
tra2 <- relist(tra.seq, tra)


# 19.2.4 - Quantifying clonal expression
# we can gain some insights into the immune activity of each T cell cluster by counting the number of expanded clonotypes in each cluster
# * clonotype = The phenotype of a clone of a cell.

clone.id.A <- unlist(unique(sce.pbmc$TRA[,"raw_clonotype_id"]))
expanded.id.A <- setdiff(clone.id.A[duplicated(clone.id.A)], "None")

clone.id.B <- unlist(unique(sce.pbmc$TRB[,"raw_clonotype_id"]))
expanded.id.B <- setdiff(clone.id.B[duplicated(clone.id.B)], "None")

is.clone.A <- any(sce.pbmc$TRA[,"raw_clonotype_id"] %in% expanded.id.A)
tra.counts.clonal <- table(colLabels(sce.pbmc)[is.clone.A])
is.clone.B <- any(sce.pbmc$TRB[,"raw_clonotype_id"] %in% expanded.id.B)
trb.counts.clonal <- table(colLabels(sce.pbmc)[is.clone.B])

# Plot the proportion of cels in each cluster that have multiple clonotypes
barplot(rbind(TRA=tra.counts.clonal/ncells, TRB=trb.counts.clonal/ncells), 
        legend=TRUE, beside=TRUE)


# determine whether a particular T cell cluster is enriched for expanding clonotypes
# -> use Fisher’s exact test to identify a significant increase in the proportion of expanded clonotypes
# -> This provides some relative measure of the average immune activity of each cluster
tclust.1 <- "2"
tclust.2 <- "6"

mat <- cbind(Expanded=tra.counts.clonal,
             Unexpanded=tra.counts.any - tra.counts.clonal)[c(tclust.1, tclust.2),]

stats <- fisher.test(mat)

# Check identities of the relevant clusters
of.interest <- colLabels(sce.pbmc) %in% c(tclust.1, tclust.2)
plotExpression(altExp(sce.pbmc)[,of.interest], 
               features=rownames(altExp(sce.pbmc)),
               other_fields="label") + facet_wrap(~label, ncol=1)



# 19.2.5 - Quantifying gene expression and properties
# Inspection of expression of specific TCR genes, which can provide some insight into the type of antigens being targeted
gene.id.A <- sce.pbmc$TRA[,"v_gene"]
expanded.cluster <- rep(colLabels(sce.pbmc), lengths(gene.id.A))
gene.tab.A <- table(unlist(gene.id.A), expanded.cluster)

# Testing for differences between our clusters.
collected <- list()
totals <- colSums(gene.tab.A)
for (i in rownames(gene.tab.A)) {
  mat <- cbind(
    Gene=gene.tab.A[i,],
    Other=totals - gene.tab.A[i,])[c(tclust.1, tclust.2),]
  stats <- fisher.test(mat)
  collected[[i]] <- DataFrame(OR=stats$estimate, p.value=stats$p.value,
                              row.names=i)
}

collected <- do.call(rbind, collected)
collected$FDR <- p.adjust(collected$p.value, method="BH")
collected <- collected[order(collected$p.value),]

############################
# 19.4 - Multi-sample analyses
############################
# Multimple samples and treatment conditions
# -> aim is to determine which clusters contain T cell clonotypes that expand in response to treatment
pretend.samples <- sample(letters[1:4], ncol(sce.pbmc), replace=TRUE)

# Creating a count matrix.
clone.counts <- any.counts <- list()
for (i in sort(unique(pretend.samples))) {
  current.sample <- sce.pbmc[,i==pretend.samples]
  clone.id.A <- unlist(unique(current.sample$TRA[,"raw_clonotype_id"]))
  expanded.id.A <- setdiff(clone.id.A[duplicated(clone.id.A)], "None")
  is.clone.A <- any(current.sample$TRA[,"raw_clonotype_id"] %in% expanded.id.A)
  clone.counts[[i]] <- table(colLabels(current.sample)[is.clone.A])
  any.counts[[i]] <- table(colLabels(current.sample)[lengths(current.sample$TRA) > 0])
}

clone.counts <- do.call(cbind, clone.counts)
any.counts <- do.call(cbind, any.counts)

# test for condition-specific differences in the proportion of clonotypes that are expanded
prop <- clone.counts/any.counts
wilcox.test(prop[tclust.1,1:2], prop[tclust.1,3:4])
wilcox.test(prop[tclust.2,1:2], prop[tclust.2,3:4])
