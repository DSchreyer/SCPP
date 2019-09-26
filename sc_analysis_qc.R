library(scater)
library(Seurat)
library(Matrix)
library(dplyr)
file = "/home/daniel/master_thesis/bassoon_data/test.counts.tsv"
# Remove cells with less than [int] UMI counts and less than [int] expressed genes
filterCountTable <- function(
  file,
  umi.count = 100,
  expressed.genes = 0,
  umitools = FALSE
){
  if (umitools){
    data <- read.table(file, sep = "\t", stringsAsFactors = F, header = T)
    barcodes <- unique(data$cell)
    genes <- unique(data$gene)
    data$gene_index <- match(data$gene, genes)
    data$cell_index <- match(data$cell, barcodes)
    data <- Matrix::sparseMatrix(
      i = data$gene_index,
      j = data$cell_index,
      x = data$count)
    rownames(data) <- genes
    colnames(data) <- barcodes
  }else{
    data <- Read10X(data.dir = file, gene.column = 2, unique.features = TRUE)
  }
  filter.umi.count <- Matrix::colSums(data) >= umi.count
  data <- data[ , filter.umi.count]
  filter.expressed.genes <- Matrix::colSums(data > 0) >= expressed.genes
  data <- data[ , filter.expressed.genes]
  return(data)
}

# Remove outlier cells based on MAD
# Outliers are determined on library size, MT genes expressed, and expressed genes
qcControl <- function(
  count.table,
  MAD = 3, 
  abundant.mt = 1,
  generate.info = FALSE
){
  # Create SingleCellExperiment object with sparse matrix
  sce <- SingleCellExperiment(assays = list(counts = count.table))
  
  # Quality Control with scater
  is.mito <- grepl("^mt-", rownames(sce))
  # quality control with caluclateQCMetrics
  sce <- calculateQCMetrics(sce, exprs_values = "counts",
                            feature_controls = list(MT = is.mito))
  n.cells.pre <- ncol(sce)
  # Drop cells above [n] median absolute deviations: libsize, mitochondrial gene count, expressed genes
  libsize.drop <- isOutlier(sce$total_counts, nmads = MAD, type = "higher", log = T)
  feature.drop <- isOutlier(sce$total_features_by_counts, nmads = MAD, type = "higher", log = T)
  mito.drop <- isOutlier(sce$pct_counts_MT, nmads = MAD, type = "higher")
  sce <- sce[, !(libsize.drop | feature.drop | mito.drop)]
  
  # remove cells with percentage of mt genes expressed to total number of expressed gebes
  pct.mt.outliers <- sce$pct_counts_MT > abundant.mt*100
  sce <- sce[ , !pct.mt.outliers]
  n.cells.post <- ncol(sce)
  
  # generate an information table with additional information
  if (generate.info){
    sample <-  sample
    n.libsize.drop <- sum(libsize.drop)
    n.feature.drop <- sum(feature.drop)
    n.mito.drop <- sum(mito.drop)
    before.qc <- n.cells.pre
    after.qc <- n.cells.post
    median.umi <- median(sce$total_counts)
    median.gene <- median(sce$total_features_by_counts)
    mean.umi <- mean(sce$total_counts)
    mean.gene <- mean(sce$total_features_by_counts)
    n.pct.mt.outlier <- sum(pct.mt.outliers)
    pct.mt.outlier.median <- median(sce$pct_counts_MT)
  }
  i <- i + 1
  if (generate.info){
    qc.info <- cbind(sample, before.qc, after.qc, n.libsize.drop, n.feature.drop, n.mito.drop, median.umi,
                     median.gene, mean.umi, mean.gene, pct.mt.outlier.median, n.pct.mt.outlier)
    qc.info <- as.data.frame(qc.info)
    write.table(qc.info, "samples_qc_info.csv", sep = ",", col.names = F)
  }
  return(sce)
}

# Remove genes, which are expressed in no more than [fraction] number of cells
geneFiltering <- function(
  object,
  gene.expr,
  sample.names = names(object)
){
  filter.gene.expr.cells <- ncol(object)*gene.expr
  keep.features <- Matrix::rowSums(counts(object) > 0) >= filter.gene.expr.cells
  object <- object[keep.features, ]
  return(object)
}

writeCountMatrix <- function(
  sce,
  log.normalize = "yes",
  filter.genes = TRUE,
  output.dir=getwd(),
  gene.expressed.cells = 0.001
){
  setwd(output.dir)
  if (log.normalize == "yes"){
    seurat.object <- as.Seurat(sce, data = "counts")
    seurat.object <- NormalizeData(seurat.object, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 1000, 
                                   verbose = FALSE)
    log.counts <- GetAssayData(seurat.object)
    if (filter.genes){
      log.counts <- SingleCellExperiment(assays = list(counts = log.counts))
      log.counts <- geneFiltering(log.counts, gene.expr = gene.expressed.cells)
      log.counts <- counts(log.counts)
    }
    writeMM(log.counts, file = "log.normalized.matrix.mtx")
  }
  if (filter.genes){
    sce <- geneFiltering(object = sce, gene.expr = gene.expressed.cells)
  }
  features <- rownames(sce)
  barcodes <- colnames(sce)
  write(features, file = "features.tsv")
  write(barcodes, file = "barcodes.tsv")
  writeMM(counts(sce), file = "raw.matrix.mtx")
}

options <- commandArgs(trailingOnly = TRUE)

branch <- as.character(options[1])
counts <- as.character(options[2])
n.genes <- as.numeric(options[3])
n.umis <-as.numeric(options[4])
MAD <- as.numeric(options[5])
mt.threshold <- as.numeric(options[6])
normalize <- as.character(options[7])
filter.genes <- as.numeric(options[8])
output.dir <- options[9]

if (branch == "UMItools"){
  table <- read.table(counts, sep = "\t", header=TRUE, stringsAsFactors = F)
  rownames(table) <- table[, 1]
  table[, 1] <- NULL
  count.table <- filterCountTable(file = counts, 
                                  umi.count = n.umis, 
                                  expressed.genes = n.genes,
                                  matrix = TRUE)
  
}else{
  count.table <- filterCountTable(file = counts, 
                                  umi.count = n.umis, 
                                  expressed.genes = n.genes)
}

sce <- qcControl(count.table = count.table,
                      MAD = MAD,
                      generate.info = TRUE)

writeCountMatrix(sce = sce, 
                 log.normalize = normalize, 
                 filter.genes = TRUE,
                 gene.expressed.cells = filter.genes,
                 output.dir = output.dir)



