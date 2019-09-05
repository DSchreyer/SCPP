#!/usr/bin/env R
library(dplyr)
library(ggplot2)
library(tidyr)
library(scater)
library(scran)
library(Seurat)


args <- commandArgs(TRUE)

library(org.Mm.eg.db)

all.ensembl <- unique(toTable(org.Mm.egENSEMBL)$ensembl_id)

count.file <- as.character(args[1])
metadata.file <- as.character(args[2])

data <- read.table(count.file, sep = "\t", header = T, stringsAsFactors = F)


### Test --- Excel mistake of authors ###
data <- data %>% distinct(X, .keep_all = T)
###
genes <- data[,1]

data <- data[-1]
rownames(data) <- genes
n.genes <- nrow(data)

sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(data))
)
n.cells.original <- ncol(sce)

# remove genes with zero counts in each cell
keep.feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep.feature, ]
n.genes.filt <- nrow(counts(sce))

print("Filter genes with zero counts across all cells:")
cat("Before:\t", n.genes)
cat("After:\t", n.genes.filt)

# generate log normalized counts and counts per million (cpm)
# generate size.factor
size.factor <- librarySizeFactors(sce)

sce <- normalize(sce, exprs_values = "counts", return_log = F)
sce <- normalize(sce, exprs_values = "counts", return_log = T)
cpm(sce) <- calculateCPM(sce)


# Compute quality control (QC) metrics for each feature and cell - Controls: Spike ins or mt genes
isSpike <- grepl("^Ercc", rownames(sce))
is.mito <- grepl("^mt-", rownames(sce))

cat("Number of mt genes:", sum(is.mito))
cat("Number of ERCC Spike ins:", sum(isSpike))
qc <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls = list(ERCC=isSpike,
                                                                               MT=is.mito))

# Drop cells with log-library sizes with more that 3 median absloute deviations (MADs) below median log-library sizes 
libsize.drop <- isOutlier(qc$total_counts, nmads = 5, type = "lower", log = T)
cat("Drop", sum(libsize.drop), "cells, because of its low library size.")
# Drop cells with 3 MADs below the median of log-transformed number of expressed genes
feature.drop <- isOutlier(qc$total_features_by_counts, nmads = 5, type = "lower", log = T)
cat("Drop", sum(feature.drop), "cells, because of its low number of expressed genes.")
# Drop cells with high mt gene / Spike in counts - 3 MADs above Median
mito.drop <- isOutlier(qc$pct_counts_MT, nmads=5, type="higher")
cat("Drop", sum(mito.drop), "cells, because of high mt gene counts.")
spike.drop <- isOutlier(qc$pct_counts_ERCC, nmads=5, type="higher")
cat("Drop", sum(spike.drop), "cells, because of high spike in counts.")


sce <- sce[ ,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
df <- data.frame(Original = n.cells.original, LibSize = sum(libsize.drop), Features = sum(feature.drop),
           Mt=sum(mito.drop), SpikeIns = sum(spike.drop), Remaining = ncol(sce))

print.data.frame(df)


# plot highest expressed genes
plotHighestExprs(qc)

# plot expression values for each gene for different features --> violin plot, specify marker genes
plotExpression(qc, features = c("Gm26870"))


plotScater(sce)
plotPCA(sce)
plotTSNE(sce)

print("DONE: Quality Control with Scater")

# Normalization and Cell Cycle Asignment with Scran -- Normalization possible with Deconvolution method or with spike-in counts
# Default: Deconvolution method
print("START: Normalization with Scran")

print("Cell cycle phase assignment")
# using mouse cycle markers file which is implemented in the scran package
# genes need ensemble id instead of normal gene names
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assigned <- cyclone(as.matrix(data), pairs=mm.pairs)
head(assigned)
table(assigned$phases)

# after that we can filter out low abundance genes 
# Cells should be in G1 phase --> only use G1 cells
# filter out cells with other cell cycle phases


# Normalization with large data sets use "quickCluster" before Normalization
# clusters <- quickCluster(sce, min.size = 100)

# Normalization with computeSumFactors (Scran)
sce <- normalize(sce)


# Batch correction possible

# Seurat 
# Create Seurat object
genes <- rownames(sce)
cell <- colnames(sce)

# extract normalized counts from sce object
table <- unlist(sce@assays$data[2])
colnames(table) <- cell
rownames(table) <- genes

# Initialize Seurat object with normalized data
SO <- CreateSeuratObject(counts = table, project = "Test",
                         min.cells = 3, min.features = 200)

# Find top variable features. Selection.method: vst, mean.var.plot, dispersion
# vst default 
SO <- FindVariableFeatures(object = SO, selection.method = "vst",
                           nfeatures = 2000)

# identify the 20 most highly variable genes
top10 <- head(VariableFeatures(SO), 10)

# plot vairable features with and without labels
plot1 <- VariableFeaturePlot(SO)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1
plot2
CombinePlots(plots = list(plot1, plot2))



# scaling the data --> linear transfromation 
# preprocessing step prior to dimensional reduction techniques
SO <- ScaleData(SO, features = genes)

# Perform linear dimensional reduction
# PCA
SO <- RunPCA(SO, features = VariableFeatures(SO))

# examine and visualize PCA results a few different ways
print(SO[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(SO, dims = 1:2, reduction = "pca")

DimPlot(SO, reduction = "pca")

DimHeatmap(SO, dims = 1, cells = 234, balanced = TRUE)

DimHeatmap(SO, dims = 1:15, cells = 234, balanced = TRUE)

# Determine the dimensionality of the dataset
# randomly permute a subset of the data (1% default) and rerun PCA
# constructiong a "null distribution" of feature scores

SO <- JackStraw(SO, num.replicate = 100)
SO <- ScoreJackStraw(SO, dims = 1:20)

JackStrawPlot(SO, dims = 1:15)

# Elbow plot
ElbowPlot(SO)

# Cluster the cells
SO <- FindNeighbors(SO, dims = 1:10)
SO <- FindClusters(SO, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(x = Idents(SO), 5)

# Run non-linear dimensional reduction (tSNE/UMAP)
SO <- RunUMAP(SO, dims = 1:10)
DimPlot(SO, reduction = "umap")

# find differently expressed features in all clusters
SO.markers <- FindAllMarkers(SO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SO.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val)
