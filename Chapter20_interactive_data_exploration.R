###########################
# 20.1 - Motivation
###########################

# - Many GUI tools exist for exploratory data analysis (EDA)
# - Here: iSEE (bioconductor package)

# Load data
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

###########################
# 20.2 - Quick start
###########################
library(iSEE)
app <- iSEE(sce.pbmc)

###########################
# 20.3 - Usage examples
###########################
# 20.3.1 - Quality control
# iSEE app can be configured to focus on quality control metrics.
#  - Library size of each cell
#  - Dimensionality reduction result where cells are coloured by the log-library size

copy.pbmc <- sce.pbmc

# computing QC metrics
copy.pbmc <- addPerCellQC(copy.pbmc, exprs_values="counts")
copy.pbmc$log10_total_counts <- log10(copy.pbmc$total)
copy.pbmc$total_counts_rank <- rank(-copy.pbmc$total)

initial.state <- list(
  # Configure a "Column data plot" panel
  ColumnDataPlot(YAxis="log10_total_counts",
                 XAxis="Column data",
                 XAxisColumnData="total_counts_rank",
                 DataBoxOpen=TRUE,
                 PanelId=1L),
  
  # Configure a "Reduced dimension plot " panel
  ReducedDimensionPlot(
    Type="TSNE",
    VisualBoxOpen=TRUE,
    DataBoxOpen=TRUE,
    ColorBy="Column data",
    ColorByColumnData="log10_total_counts",
    SelectionBoxOpen=TRUE,
    ColumnSelectionSource="ColumnDataPlot1")
)

app <- iSEE(copy.pbmc,
            initial = initial.state)


# 20.3.2 - Annotation of cell populations
# use iSEE to interactively examine the marker genes to conveniently determine cell identities

copy.pbmc <- sce.pbmc
# identify upregulated markers in each cluster and collect log-p-value for each gene in each cluster
markers.pbmc.up <- findMarkers(copy.pbmc,
                               direction = "up",
                               log.p = TRUE,
                               sorted = FALSE)

# Collate the log-p-value for each marker in a single table
all.p <- lapply(markers.pbmc.up, FUN = "[[", i="log.p.value")
all.p <- DataFrame(all.p, check.names=FALSE)
colnames(all.p) <- paste0("cluster", colnames(all.p))

# Store the table of results as row metadata
rowData(copy.pbmc) <- cbind(rowData(copy.pbmc), all.p)

# Create an app that contains:
# - table of feature statistics
# - plot showing the distribution of expression values for a chosne gene in each cluster
# - plot showing the result of UMAP dim.red.

initial.state <- list(
  RowDataTable(PanelId=1L),
  
  # Configure a "Feature assay plot" panel
  FeatureAssayPlot(
    YAxisFeatureSource="RowDataTable1",
    XAxis="Column data",
    XAxisColumnData="label",
    Assay="logcounts",
    DataBoxOpen=TRUE
  ),
  
  # Configure a "Reduced dimension plot" panel
  ReducedDimensionPlot(
    Type="UMAP",
    ColorBy="Feature name",
    ColorByFeatureSource="RowDataTable1",
    ColorByFeatureNameAssay="logcounts"
  )
)

app <- iSEE(copy.pbmc,
            initial = initial.state)



# 20.3.3 - Querying features of interest
# > gene-centric exploratory analyses

# add variance modelling statistics to the rowData() of our SingleCellExperiment object.
copy.pbmc <- sce.pbmc
dec <- modelGeneVarByPoisson(copy.pbmc)
rowData(copy.pbmc) <- cbind(rowData(copy.pbmc), dec)

# Create an app that contains
# - mean-var trend where each point is a cell
# - feature statistics
# - heatmap for the genes of first plot
initial.state <- list(
  # Configure a "Feature assay plot" panel
  RowDataPlot(
    YAxis="total",
    XAxis="Row data",
    XAxisRowData="mean",
    PanelId=1L
  ),
  
  RowDataTable(
    RowSelectionSource="RowDataPlot1"
  ),
  
  # Configure a "ComplexHeatmap" panel
  ComplexHeatmapPlot(
    RowSelectionSource="RowDataPlot1",
    CustomRows=FALSE,
    ColumnData="label",
    Assay="logcounts",
    ClusterRows=TRUE,
    PanelHeight=800L,
    AssayCenterRows=TRUE
  )
)

app <- iSEE(copy.pbmc,
            initial = initial.state)


# set up an app where selecting a single HVG in the mean-variance plot causes the neighboring  
# t-SNE to be colored by the expression of the selected gene 
initial.state <- list(
  # Configure a "Feature assay plot" panel
  RowDataPlot(
    YAxis="total",
    XAxis="Row data",
    XAxisRowData="mean",
    PanelId=1L
  ),
  
  # Configure a "Reduced dimension plot" panel
  ReducedDimensionPlot(
    Type="TSNE",
    ColorBy="Feature name",
    ColorByFeatureSource="RowDataPlot1",
    ColorByFeatureNameAssay="logcounts"
  )
)

app <- iSEE(copy.pbmc, initial=initial.state)


###########################
# 20.5 - Dissemination of analysis results
###########################
# writing a step-by-step walkthrough of the different panels with explanations to facilitate their interpretation
# -> add a data frame with two columns named “element” and “intro” ot the iSEE object
# -> provide this df with the 'tour' - argument
tour <- data.frame(
  element = c(
    "#Welcome",
    "#ReducedDimensionPlot1",
    "#ColumnDataPlot1",
    "#ColumnDataPlot1_DataBoxOpen",
    "#Conclusion"),
  intro = c(
    "Welcome to this tour!",
    "This is a <i>Reduced dimension plot.</i>",
    "And this is a <i>Column data plot.</i>",
    "<b>Action:</b> Click on this collapsible box to open and close it.",
    "Thank you for taking this tour!"),
  stringsAsFactors = FALSE)

initial.state <- list(
  ReducedDimensionPlot(PanelWidth=6L), 
  ColumnDataPlot(PanelWidth=6L)
)

app <- iSEE(sce.pbmc,
            initial = initial.state,
            tour = tour)
