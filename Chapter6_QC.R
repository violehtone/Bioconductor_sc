# BiocManager::install("package_name")

library(scRNAseq)
library(SingleCellExperiment)
library(scater)

#################
# Introduction
#################

# List the available datasets in the scRNA-seq package
# > browseVignettes("scRNAseq")
# -> HTML
# -> 2 - Available data sets

# Spike-in scRNA-seq dataset from Lun et al. (2017)
sce.416b <- LunSpikeInData('416b')

#################
# 6.2 - Choice of QC metrics
#################
# Retrieve the mitochondrial transcripts
location <- rowRanges(sce.416b)
is.mito <- any(seqnames(location) == "MT")
df <- perCellQCMetrics(sce.416b, subsets = list(Mito = is.mito))
head(df)

# Alternative way: use addPerCellQC()
sce.416b <- addPerCellQC(sce.416b,
                         subsets = list(Mito = is.mito))
colnames(colData(sce.416b))

#################
# 6.3 - Identifying low-quality cells
#################

# Identifying and discard low quality cells
qc.lib <- df$sum < 1e5
qc.nexprs <- df$detected < 5e3
qc.spike <- df$altexps_ERCC_percent > 10
qc.mito <- df$subsets_Mito_percent > 10

discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

# Summarize the number of cells removed for each reason
DataFrame(LibSize = sum(qc.lib), Nexprs = sum(qc.nexprs),
          SpikeProp = sum(qc.spike), MitoProp = sum(qc.mito),
          Total = sum(discard))

# Identifying outliers
qc.lib2 <- isOutlier(df$sum, log = TRUE, type = "lower")
qc.nexprs2 <- isOutlier(df$detected, log = TRUE, type = "lower")
qc.spike2 <- isOutlier(df$altexps_ERCC_percent, type = "higher")
qc.mito2 <- isOutlier(df$subsets_Mito_percent, type = "higher")

# Check the exact filter tresholds
attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")
attr(qc.spike2, "thresholds")
attr(qc.mito2,  "thresholds")

# A cell is considered outlier if it is considered to be low quality in any of the above metrics
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib2), NExprs=sum(qc.nexprs2),
          SpikeProp=sum(qc.spike2), MitoProp=sum(qc.mito2), Total=sum(discard2))

# An alternative way to do the detection of outliers is using the quickPerCellQC()
reasons <- quickPerCellQC(df = df,
                          percent_subsets = c("subsets_Mito_percent",
                                              "altexps_ERCC_percent"))


colSums(as.matrix(reasons))

# Considering experimental factors
batch <- paste0(sce.416b$phenotype, "-", sce.416b$Plate)
batch.reasons <- quickPerCellQC(df,  
                                percent_subsets = c("subsets_Mito_percent", "altexps_ERCC_percent"),
                                batch = batch)

colSums(as.matrix(batch.reasons))


# Human pancreas scRNA-seq data from Grun et al. (2016)
sce.grun <- GrunPancreasData()
sce.grun <- addPerCellQC(sce.grun)

# First attempt with batch-specific thresholds
discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
                          type="higher", batch=sce.grun$donor)


discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
                          type = "higher",
                          batch = sce.grun$donor)

with.blocking <- plotColData(object = sce.grun,
                             x = "donor",
                             y = "altexps_ERCC_percent",
                             colour_by = I(discard.ercc))

# Second attempt, sharing information across batches
discard.ercc2 <- isOutlier(sce.grun$altexps_ERCC_percent,
                          type = "higher",
                          batch = sce.grun$donor,
                          subset = sce.grun$donor %in% c("D17", "D2", "D7"))

without.blocking <- plotColData(object = sce.grun,
                                x = "donor",
                                y = "altexps_ERCC_percent",
                                colour_by = I(discard.ercc2))

# Plot and compare the two results
gridExtra::grid.arrange(with.blocking, without.blocking, ncol = 2)

#  Finding batches with QC thresholds that are themselves outliers compared to the thresholds of other batches
ercc.thresholds <- attr(discard.ercc, "thresholds")["higher", ]
ercc.tresholds
names(ercc.thresholds)[isOutlier(ercc.thresholds, type="higher")]

# Other approaches for outlier detection
stats <- cbind(log10(df$sum), log10(df$detected),
               df$subsets_Mito_percent, df$altexps_ERCC_percent)

library(robustbase)
outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")
summary(multi.outlier)


#################
# 6.4 - Checking diagnostic plots
#################

colData(sce.416b) <- cbind(colData(sce.416b), df)
sce.416b$block <- factor(sce.416b$block)
sce.416b$phenotype <- ifelse(grepl("induced", sce.416b$phenotype),
                             "induced", "wild type")
sce.416b$discard <- reasons$discard

gridExtra::grid.arrange(
  plotColData(sce.416b, 
              x="block", 
              y="sum", 
              colour_by="discard",
              other_fields="phenotype") 
  + facet_wrap(~phenotype) 
  + scale_y_log10() 
  + ggtitle("Total count"),
  plotColData(sce.416b, 
              x="block", 
              y="detected", 
              colour_by="discard", 
              other_fields="phenotype") 
  + facet_wrap(~phenotype) 
  + scale_y_log10() 
  + ggtitle("Detected features"),
  plotColData(sce.416b, 
              x="block", 
              y="subsets_Mito_percent", 
              colour_by="discard",
              other_fields="phenotype") 
  + facet_wrap(~phenotype) 
  + ggtitle("Mito percent"),
  plotColData(sce.416b, 
              x="block", 
              y="altexps_ERCC_percent", 
              colour_by="discard", 
              other_fields="phenotype") 
  + facet_wrap(~phenotype) 
  + ggtitle("ERCC percent"),
  ncol=1
)

# Plotting the proportion of mitochondrial counts against some of the other QC metrics
# Mouse brain data set (Zeisel et al., 2015)

sce.zeisel <- ZeiselBrainData()

sce.zeisel <- addPerCellQC(sce.zeisel,
                           subsets = list(Mt = rowData(sce.zeisel)$featureType == "mito"))

qc <- quickPerCellQC(colData(sce.zeisel),
                     percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent"))

sce.zeisel$discard <- qc$discard

# Plot the data
# Mt percent (% of UMIs assigned to mito transcripts vs. total # of UMIs)
plotColData(sce.zeisel, x="sum", y="subsets_Mt_percent", colour_by="discard")
# ERCC percent (% of UMIs assigned to mito transcripts vs. % UMIs assigned to spike-in transcripts)
plotColData(sce.zeisel, x="altexps_ERCC_percent", y="subsets_Mt_percent", colour_by="discard")



#################
# 6.5 - Cell calling for droplet data
#################

# Obtain unfiltered count matrix for the PBMC dataset from 10X Genomics
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc,
                     file.path("http://cf.10xgenomics.com/samples",
                               "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
library(Matrix)

fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = T)
sce.pbmc

# Compute barcode rank statistics from count data (# reads)
bcrank <- barcodeRanks(counts(sce.pbmc))

uniq <- !duplicated(bcrank$rank)

# Make a plot of total UMI counts vs. Rank
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)

abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)


# Testing for empty droplets
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
summary(e.out$FDR <= 0.001)

table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)

set.seed(100)
limit <- 100
all.out <- emptyDrops(counts(sce.pbmc),
                      lower = limit,
                      test.ambient = TRUE)

hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
     xlab="P-value",
     main="",
     col="grey80") 

# Subset the sc-experiment object to retain only the detected cells
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]

# Other QC metrics (filtering out mitochondrial proportion)
is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)
pbmc.qc <- perCellQCMetrics(sce.pbmc, subsets=list(MT=is.mito))
discard.mito <- isOutlier(pbmc.qc$subsets_MT_percent, type="higher")
summary(discard.mito)

# Plot the mito % vs. total counts
plot(pbmc.qc$sum,
     pbmc.qc$subsets_MT_percent,
     log="x",
     xlab="Total count", 
     ylab='Mitochondrial %')

abline(h=attr(discard.mito, "thresholds")["higher"], col="red")


#################
# 6.6 - Removing low-quality cells
#################

# Keep columns (cells) we want to keep
filtered <- sce.416b[, !reasons$discard]

# compute avg. count across the discarded and retained pools
lost <- calculateAverage(counts(sce.416b)[,!discard])
kept <- calculateAverage(counts(sce.416b)[,discard])

library(edgeR)
logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

# Plot the fold change of lost/kept vs. avg. count
# Each point represents a gene with mitochondrial transcripts in blue
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[is.mito], logFC[is.mito], col="dodgerblue", pch=16)

# If a fixed treshold on the library size to filter cells was applied:
alt.discard <- colSums(counts(sce.pbmc)) < 500
lost <- calculateAverage(counts(sce.pbmc)[,alt.discard])
kept <- calculateAverage(counts(sce.pbmc)[,!alt.discard])

logged <- edgeR::cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

# Reproduce plot. Platelet-related genes are highlighted in orange
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
platelet <- c("PF4", "PPBP", "CAVIN2")

library(org.Hs.eg.db)
ids <- mapIds(org.Hs.eg.db, keys=platelet, column="ENSEMBL", keytype="SYMBOL")
points(abundance[ids], logFC[ids], col="orange", pch=16)


#################
# 6.7 - Marking low-quality cells
#################
# Another option is to just mark (and not discard) the low quality cells and thus retain them for the downstream analysis
marked <- sce.416b
marked$discard <- batch.reasons$discard


