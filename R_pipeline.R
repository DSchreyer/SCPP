#!/usr/bin/env R
library(dplyr)
library(ggplot2)
library(tidyr)
library(scater)

args <- commandArgs(TRUE)

## Testing ##
args <- c("/home/daniel/master_thesis/single_cell_data/GSE80232_vsx2.RSEM.genes.counts.matrix.txt", "")
####

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
is.spike <- grepl("^Ercc", rownames(sce))
is.mito <- grepl("^mt-", rownames(sce))

cat("Number of mt genes:", sum(is.mito))
cat("Number of ERCC Spike ins:", sum(is.spike))
qc <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls = list(ERCC=is.spike,
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
