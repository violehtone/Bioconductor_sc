###############################
# 12.1 - Motivation
###############################

# Interpreting results (i.e. what biological state each cluster represents)
# is the hardest task in scRNA-seq analysis

# we can use various computational approaches that exploit prior information to
# assign meaning to an uncharacterized scRNA-seq dataset.
#  > Gene ontology (GO)
#  > Kyoto Encyclopedia of Genes and Genomes (KEGG)
#  > Directly compare expression profiles to published reference datasets
#    , where each cell has already been annotated

# Load pbmc dataset
#----------------------------------------------------------------
#--- loading ---#
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

sce.pbmc
#----------------------------------------------------------------

###############################
# 12.2 - Assigning cell labels from reference data
###############################

# SingleR is a tool for performing cell type annotation
# We label the data with singleR() function
library(SingleR)
ref <- BlueprintEncodeData()
pred <- SingleR(test = sce.pbmc,
                ref = ref,
                labels = ref$label.main)

table(pred$labels)
plotScoreHeatmap(pred)

sum(is.na(pred$pruned.labels))
plotScoreDistribution(pred)

# Inspect the amount of different cell types in each cluster
tab <- table(Assigned = pred$pruned.labels,
             Cluster = colLabels(sce.pbmc))

# Make a heatmap (add pseudocount of +10 to avoid color jumps with just 1 cell)
pheatmap(log2(tab+10),
         color = colorRampPalette(c("white", "blue"))(101))


# 13.2.3 - Using custom references
# > Use an existing data set by the user

## -- Load muraro data set -- ##
#--- loading ---#
library(scRNAseq)
sce.muraro <- MuraroPancreasData()

#--- gene-annotation ---#
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
gene.symb <- sub("__chr.*$", "", rownames(sce.muraro))
gene.ids <- mapIds(edb, keys=gene.symb, 
                   keytype="SYMBOL", column="GENEID")

# Removing duplicated genes or genes without Ensembl IDs.
keep <- !is.na(gene.ids) & !duplicated(gene.ids)
sce.muraro <- sce.muraro[keep,]
rownames(sce.muraro) <- gene.ids[keep]

#--- quality-control ---#
library(scater)
stats <- perCellQCMetrics(sce.muraro)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent",
                     batch=sce.muraro$donor, subset=sce.muraro$donor!="D28")
sce.muraro <- sce.muraro[,!qc$discard]

#--- normalization ---#
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.muraro)
sce.muraro <- computeSumFactors(sce.muraro, clusters=clusters)
sce.muraro <- logNormCounts(sce.muraro)
#### -----------  #####

# Inspect cell types in muraro data
sce.muraro <- sce.muraro[,!is.na(sce.muraro$label) & 
                           sce.muraro$label!="unclear"]
table(sce.muraro$label)

# Use the muraro dataset to set cell type labels to the Seger dataset
# ---- Load seger dataset ------ #
#--- loading ---#
library(scRNAseq)
sce.seger <- SegerstolpePancreasData()

#--- gene-annotation ---#
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
symbols <- rowData(sce.seger)$symbol
ens.id <- mapIds(edb, keys=symbols, keytype="SYMBOL", column="GENEID")
ens.id <- ifelse(is.na(ens.id), symbols, ens.id)

# Removing duplicated rows.
keep <- !duplicated(ens.id)
sce.seger <- sce.seger[keep,]
rownames(sce.seger) <- ens.id[keep]

#--- sample-annotation ---#
emtab.meta <- colData(sce.seger)[,c("cell type", 
                                    "individual", "single cell well quality")]
colnames(emtab.meta) <- c("CellType", "Donor", "Quality")
colData(sce.seger) <- emtab.meta

sce.seger$CellType <- gsub(" cell", "", sce.seger$CellType)
sce.seger$CellType <- paste0(
  toupper(substr(sce.seger$CellType, 1, 1)),
  substring(sce.seger$CellType, 2))

#--- quality-control ---#
low.qual <- sce.seger$Quality == "low quality cell"

library(scater)
stats <- perCellQCMetrics(sce.seger)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent",
                     batch=sce.seger$Donor,
                     subset=!sce.seger$Donor %in% c("HP1504901", "HP1509101"))

sce.seger <- sce.seger[,!(qc$discard | low.qual)]

#--- normalization ---#
library(scran)
clusters <- quickCluster(sce.seger)
sce.seger <- computeSumFactors(sce.seger, clusters=clusters)
sce.seger <- logNormCounts(sce.seger) 

pred.seger <- SingleR(test = sce.seger,
                      ref = sce.muraro,
                      labels = sce.muraro$label,
                      de.method = "wilcox")

# Inspect cell types
table(pred.seger$labels)

# Compare the predicted labels to the independently defined labels
tab <- table(pred.seger$pruned.labels, sce.seger$CellType)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))


###############################
# 12.3 - Assigning cell labels from gene sets
###############################
# Load the neuronal cell type markers from Zeisel et al. study
#--- loading ---#
library(scRNAseq)
sce.zeisel <- ZeiselBrainData()

library(scater)
sce.zeisel <- aggregateAcrossFeatures(sce.zeisel, 
                                      id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))

#--- gene-annotation ---#
library(org.Mm.eg.db)
rowData(sce.zeisel)$Ensembl <- mapIds(org.Mm.eg.db, 
                                      keys=rownames(sce.zeisel), keytype="SYMBOL", column="ENSEMBL")

#--- quality-control ---#
stats <- perCellQCMetrics(sce.zeisel, subsets=list(
  Mt=rowData(sce.zeisel)$featureType=="mito"))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", 
                                              "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[,!qc$discard]

#--- normalization ---#
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.zeisel)
sce.zeisel <- computeSumFactors(sce.zeisel, cluster=clusters) 
sce.zeisel <- logNormCounts(sce.zeisel)


## -------------- ##

# Identify sets of marker genes that are highly expressed in each individual cell
wilcox.z <- pairwiseWilcox(sce.zeisel,
                           sce.zeisel$level1class,
                           lfc=1,
                           direction = "up")

markers.z <- getTopMarkers(wilcox.z$statistics,
                           wilcox.z$pairs,
                           pairwise = FALSE,
                           n = 50)

lengths(markers.z)

# Load tasic brain dataset (test set)
sce.tasic <- TasicBrainData()

# Use AUCell package to identify marker sets
library(GSEABase)
all.sets <- lapply(names(markers.z),
                   function(x) {
                     GeneSet(markers.z[[x]],
                             setName = x)
                   })

all.sets <- GeneSetCollection(all.sets)

#BiocManager::install("AUCell")
library(AUCell)
rankings <- AUCell_buildRankings(counts(sce.tasic),
                                 plotStats = FALSE,
                                 verbose = FALSE)

cell.aucs <- AUCell_calcAUC(all.sets, rankings)
results <- t(assay(cell.aucs))
head(results)

new.labels <- colnames(results)[max.col(results)]
tab <- table(new.labels, sce.tasic$broad_type)
tab

#BiocManager::install("DelayedMatrixStats")
library(DelayedMatrixStats)

deltas <- rowMaxs(results) - rowMedians(results)
discard <- isOutlier(deltas, type="lower", batch=new.labels)
table(new.labels[discard])

# Make a boxplot of the distribution of differences between the max and median AUCs
# for each cell
par(mar=c(10,4,1,1))
boxplot(split(deltas, new.labels), las=2)
points(attr(discard, "thresholds")[1,], col="red", pch=4, cex=2)


# Use single-cell signatures defined from MSigDB
# Downloading the signatures and caching them locally.
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
scsig.path <- bfcrpath(bfc, file.path("http://software.broadinstitute.org",
                                      "gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt"))
scsigs <- getGmt(scsig.path)

muraro.mat <- counts(sce.muraro)
rownames(muraro.mat) <- rowData(sce.muraro)$symbol
muraro.rankings <- AUCell_buildRankings(muraro.mat,
                                        plotStats=FALSE, verbose=FALSE)

# Applying MsigDB to the Muraro dataset, because it's human:
scsig.aucs <- AUCell_calcAUC(scsigs, muraro.rankings)
scsig.results <- t(assay(scsig.aucs))
full.labels <- colnames(scsig.results)[max.col(scsig.results)]
tab <- table(full.labels, sce.muraro$label)
fullheat <- pheatmap(log10(tab+10), color=viridis::viridis(100), silent=TRUE)

# Restricting to the subset of Muraro-derived gene sets:
scsigs.sub <- scsigs[grep("Pancreas", names(scsigs))]
sub.aucs <- AUCell_calcAUC(scsigs.sub, muraro.rankings)
sub.results <- t(assay(sub.aucs))
sub.labels <- colnames(sub.results)[max.col(sub.results)]
tab <- table(sub.labels, sce.muraro$label)
subheat <- pheatmap(log10(tab+10), color=viridis::viridis(100), silent=TRUE)

gridExtra::grid.arrange(fullheat[[4]], subheat[[4]])



###############################
# 12.4 - Assigning cluster labels from markers
###############################
# another strategy for annotation is to perform a gene set enrichment analysis on 
# the marker genes defining each cluster.
# > This identifies the pathways and processes that are (relatively) active
#   in each cluster based on upregulation of the associated genes compared to other clusters.

# Load Bach et al. (2017) data set
## ------------------------ ##
#--- loading ---#
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

## ------------------------ ##

# Find up-regulated marker genes
markers.mam <- findMarkers(sce.mam, 
                           direction="up",
                           lfc=1)

# annotations for the marker genes that define cluster 2
# > use gene sets defined by the Gene Ontology (GO) project
# > relevant marker genes have FDR <= 5%
# > goana() identifies GO terms that are overrepresented in our marker subset

chosen <- "2"
cur.markers <- markers.mam[[chosen]]
is.de <- cur.markers$FDR <= 0.05 
summary(is.de)

library(org.Mm.eg.db)
entrez.ids <- mapIds(org.Mm.eg.db, keys=rownames(cur.markers), 
                     column="ENTREZID", keytype="SYMBOL")

library(limma)
go.out <- goana(unique(entrez.ids[is.de]), species="Mm", 
                universe=unique(entrez.ids))

# Only keeping biological process terms that are not overly general.
go.out <- go.out[order(go.out$P.DE),]
go.useful <- go.out[go.out$Ont=="BP" & go.out$N <= 200,]
head(go.useful, 20)

# Plot genes in cell types
plotExpression(sce.mam, features=c("Csn2", "Csn3"), 
               x="label", colour_by="label")

# Extract the relevant genes
# Extract symbols for each GO term; done once.
tab <- select(org.Mm.eg.db, keytype="SYMBOL", 
              keys=rownames(sce.mam), columns="GOALL")
by.go <- split(tab[,1], tab[,2])

# Identify genes associated with an interesting term.
adhesion <- unique(by.go[["GO:0022408"]])
head(cur.markers[rownames(cur.markers) %in% adhesion,1:4], 10)

# Gene set testing of marker lists is a reliable approach for determining 
# if pathways are up- or down-regulated between clusters


###############################
# 12.5 - Computing gene set activities
###############################
# we can also quantify gene set activity on a per-cell level and test for differences in activity. 
#  we simply compute the average of the log-expression values across all genes in the set for each cell.

aggregated <- sumCountsAcrossFeatures(sce.mam, by.go,
                                      exprs_values="logcounts", 
                                      average=TRUE)

dim(aggregated) # rows are gene sets, columns are cells
aggregated[1:10,1:5]

# We can then identify “differential gene set activity” between clusters
# by looking for significant differences in the per-set averages of the relevant cells
# i.e. cluster 2 has the highest average expression for the triacylglycerol biosynthesis GO term
plotColData(sce.mam, y=I(aggregated["GO:0019432",]), x="label")

# Choose the top-ranking gene in GO:0019432.
plotExpression(sce.mam, "Thrsp", x="label")

