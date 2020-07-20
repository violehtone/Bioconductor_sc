library(scRNAseq)

#################
# 8.1 - Motivation
#################
# list available data sets
#browseVignettes("scRNAseq")

# Load data sets:
# Dataset of Lun et al., (2017)
sce.416b <- LunSpikeInData('416b')

# Peripheral blood mononuclear cell (PBMC) dataset from 10X Genomics (Zheng et al., 2017)
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

# View datasets
sce.pbmc
sce.416b


#################
# 8.2 - Quantifying per-gene variation
#################
#8.2.1 Variance of the log-counts
# Use modelGeneVar() to fit a trend to the variance with respect to abundance across all genes
library(scran)
dec.pbmc <- modelGeneVar(sce.pbmc)

# Visualizing the fit
fit.pbmc <- metadata(dec.pbmc)
plot(x = fit.pbmc$mean,
     y = fit.pbmc$var,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression")
curve(fit.pbmc$trend(x),
      col = "dodgerblue",
      add = TRUE,
      lwd = 2)


# Ordering by most interesting genes for inspection
dec.pbmc[order(dec.pbmc$bio, decreasing = TRUE), ]

#8.2.2 - Coefficient of variation
#CV^2 = squared coefficient of variation -> metric for describing variation in non-neg data

# Compute CV^2 for each gene in pbmc dataset using modelGeneCV2() function
dec.cv2.pbmc <- modelGeneCV2(sce.pbmc)

fit.cv2.pbmc <- metadata(dec.cv2.pbmc)
plot(x = fit.cv2.pbmc$mean,
     y = fit.cv2.pbmc$cv2,
     log="xy")
curve(fit.cv2.pbmc$trend(x),
      col = "dodgerblue",
      add=TRUE,
      lwd=2)

# Quantify the deviation from the trend for each gene
dec.cv2.pbmc[order(dec.cv2.pbmc$ratio, decreasing = TRUE), ]


#8.2.3 - Quantifying technical noise
dec.spike.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC")
dec.spike.416b[order(dec.spike.416b$bio, decreasing = TRUE), ]


# plot variance vs. mean
plot(x = dec.spike.416b$mean, 
     y = dec.spike.416b$total, 
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression")

fit.spike.416b <- metadata(dec.spike.416b)

points(fit.spike.416b$mean, fit.spike.416b$var, col="red", pch=16)
curve(fit.spike.416b$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# Create a mean-variance trend
set.seed(0010101)
dec.pois.pbmc <- modelGeneVarByPoisson(sce.pbmc)
dec.pois.pbmc <- dec.pois.pbmc[order(dec.pois.pbmc$bio, decreasing = TRUE), ]
head(dec.pois.pbmc)

plot(dec.pois.pbmc$mean,
     dec.pois.pbmc$total,
     pch=16,
     xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(metadata(dec.pois.pbmc)$trend(x), col = "dodgerblue", add = TRUE)

# 8.2.4 - Accounting for blocking factors
# 8.2.4.1 Fitting block-specific trends

# Find genes that are highly variable (HVGs) within each batch
dec.block.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC", block = sce.416b$block)
head(dec.block.416b[order(dec.block.416b$bio, decreasing = TRUE), 1:6])

# Plot the variance of 416B data set as a function of the mean after blocking on the plate of origin
par(mfrow = c(1,2))

blocked.stats <- dec.block.416b$per.block

for(i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  plot(current$mean,
       current$total,
       main = i,
       pch = 16,
       cex = 0.5,
       xlab = "Mean of log-expression",
       ylab = "Variance of log-expression")
  
  curfit <- metadata(current)
  points(curfit$mean, curfit$var, col="red", pch=16)
  curve(curfit$trend(x), col = 'dodgerblue', add = TRUE, lwd = 2)
}


# 8.2.4.2 - Using a design matrix
# > when there is a large number of blocking factors
# > design argument specifies a design matrix with uninteresting factors of variation

design <- model.matrix(~factor(block) + phenotype,
                       colData(sce.416b))

dec.design.416b <- modelGeneVarWithSpikes(sce.416b,
                                          "ERCC",
                                          design = design)

dec.design.416b[order(dec.design.416b$bio, decreasing = TRUE), ]


#################
# 8.3 - Selecting highly variable genes
#################
# Once the variation per gene is quantified, next we want to select
# the HVGs for downstream analyses

# 8.3.2 - Simple: Feature selection based on largest metrics
# i.e. take top 1000 genes
hvg.pbmc.var <- getTopHVGs(dec.pbmc, n = 1000)
str(hvg.pbmc.var)

hvg.pbmc.cv2 <- getTopHVGs(dec.cv2.pbmc, var.field = "ratio", n = 1000)
str(hvg.pbmc.cv2)

# -> not a suitable approach for very heterogenous populations
# -> choice of n is random are isn't based on any statistical test



# 8.3.3 - Feature selection based on significance (p-value)
hvg.pbmc.var.2 <- getTopHVGs(dec.pbmc, fdr.threshold = 0.05)
length(hvg.pbmc.var.2)

# > There is no reason to think that .05 threshold on FDR gives
# better results as compared to selecting top n=1000 features based on variance



# 8.3.4 - Feature selection by keeping all genes above the trend
# > aim is to only remove the obviously uninteresting genes with variances below the trend
hvg.pbmc.var.3 <- getTopHVGs(dec.pbmc, var.threshold = 0)
length(hvg.pbmc.var.3)

hvg.pbmc.cv2.3 <- getTopHVGs(dec.cv2.pbmc,
                             var.field = "ratio",
                             var.threshold = 1)
length(hvg.pbmc.cv2.3)

# > Most useful for rare subpopulations and heterogenous populations


#################
# 8.4 - Selecting a priori genes of interst
#################
# Use a pre-defined set of interesting genes
# i.e. using only genes associated with a specific metabolic pathway

library(msigdbr)
c7.sets <- msigdbr(species = "Homo sapiens",
                   category = "C7")
head(unique(c7.sets$gs_name))

# Use the Goldrath sets to distinguish CD8 subtypes
cd8.sets <- c7.sets[grep("GOLDRATH", c7.sets$gs_name), ]
cd8.genes <- rowData(sce.pbmc)$Symbol %in% cd8.sets$human_gene_symbol
summary(cd8.genes)

# Using GSE11924 to distinguish between T helper subtypes
th.sets <- c7.sets[grep("GSE11924", c7.sets$gs_name), ]
th.genes <- rowData(sce.pbmc)$Symbol %in% th.sets$human_gene_symbol
summary(th.genes)

# Using GSE11961 to distinguish between B cell subtypes
b.sets <- c7.sets[grep("GSE11961", c7.sets$gs_name),]
b.genes <- rowData(sce.pbmc)$Symbol %in% b.sets$human_gene_symbol
summary(b.genes)

# Identify ribosomal protein genes
ribo.discard <- grepl("^RP[SL]\\d+", rownames(sce.pbmc))
sum(ribo.discard)

# More curated approach for identifying ribosomal protein genes
c2.sets <- msigdbr(species = "Homo sapiens", category = "C2")
ribo.set <- c2.sets[c2.sets$gs_name=="KEGG_RIBOSOME",]$human_gene_symbol
ribo.discard <- rownames(sce.pbmc) %in% ribo.set
sum(ribo.discard)

# Removing immunoglobulin variable chains
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
anno <- select(x = edb,
               keys = rowData(sce.pbmc)$ID,
               keytype = "GENEID",
               columns = "TXBIOTYPE")

igv.set <- anno$GENEID[anno$TXBIOTYPE %in% c("IG_V_gene", "IG_V_pseudogene")]
igv.discard <- rowData(sce.pbmc)$ID %in% igv.set
sum(igv.discard)

# Removing TCR variable chains
tcr.set <- anno$GENEID[anno$TXBIOTYPE %in% c("TR_V_gene", "TR_V_pseudogene")]
tcr.discard <- rowData(sce.pbmc)$ID %in% tcr.set
sum(tcr.discard)



#################
# 8.5 - Putting it all together
#################
# Select top 10% of genes with the highest biological components
dec.pbmc <- modelGeneVar(sce.pbmc)
chosen <- getTopHVGs(dec.pbmc, prop=0.1)
str(chosen)

# Option 1 - retain only our selection of HVGs.
sce.pbmc.hvg <- sce.pbmc[chosen, ]
dim(sce.pbmc.hvg)

# Option 2 - keep the original object and specify the genes to use
# in downstream analyses
library(scater)
sce.pbmc <- runPCA(sce.pbmc, subset_row = chosen)
reducedDimNames(sce.pbmc)

# Option 3 - retain full dataset while using HVGs for downstream analyses
altExp(sce.pbmc.hvg, "original") <- sce.pbmc
altExpNames(sce.pbmc.hvg)

sce.pbmc.hvg <- runPCA(sce.pbmc.hvg)
sce.pbmc.original <- altExp(sce.pbmc.hvg, "original", withColData = TRUE)
