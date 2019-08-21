library(scater)
library(scran)
library(Seurat)
library(Matrix)

# laptop
# path.to.cellranger.outs <- paste("/home/daniel/master_thesis/bassoon_data/Cellranger output/", sample, sep = "")
# cluster
# path.to.cellranger.outs <- paste("/media/data2/Daniel/cellranger_output/", sample, "/outs/", sep = "/")

# Read in sparse matrix with cellranger output
# data.dir = paste(path.to.cellranger.outs, "/raw_feature_bc_matrix/", sep = "/")

# cluster
# data.dir <- "/media/data2/Daniel/cellranger_aggr/AGG_Bsn/outs/raw_feature_bc_matrix/"
# output.dir <- "media/data2/Daniel/cellranger_aggr/R_output/post_qc/"

# laptop
data.dir <- "/home/daniel/master_thesis/bassoon_data/Cellranger output/Aggr/raw_feature_bc_matrix/"
output.dir <- "/home/daniel/master_thesis/bassoon_data/Output/post_qc/"

# Create Output Directory
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
setwd(output.dir)

sample.names <- paste0("Bsn.", seq(221931,221940))

generateCountTables <- function(
  aggr.dir,
  umi.count = 100,
  expressed.genes = 0,
  sample.names = NULL
){
  data <- Read10X(data.dir = aggr.dir, gene.column = 2, unique.features = TRUE)
  
  if (is.null(sample.names)){
    m <- regexpr("-.*", colnames(data), perl = TRUE)
    n.samples <- length(unique(regmatches(x = colnames(data), m)))
    sample.names <- paste0("sample", as.character(seq(1, n.samples)))
  } else {
    n.samples <- length(sample.names)
  }
  
  filter.umi.count <- Matrix::colSums(data) >= umi.count
  data <- data[ , filter.umi.count]
  
  filter.expressed.genes <- Matrix::colSums(data > 0) >= expressed.genes
  data <- data[ , filter.expressed.genes]
  
  all.count.tables <- lapply(c(as.character(seq(1, n.samples))),
                             function(x) data[, grep(x, colnames(data))])
  names(all.count.tables) <- sample.names
  return(all.count.tables)
}

count.tables <- generateCountTables(aggr.dir = data.dir, umi.count = 120, expressed.genes = 100, sample.names = sample.names)


qcControl <- function(
  count.tables,
  MAD = 3, 
  sample.names = names(count.tables)
){
  if (is.null(sample.names)){
    sample.names <- paste0("sample", as.character(seq(1, length(count.tables))))
  }
  i <- 1
  sce.list <- list()
  for (count.table in count.tables){
    sample <- sample.names[i]
    i <- i + 1
    # Create SingleCellExperiment object with sparse matrix
    sce <- SingleCellExperiment(assays = list(counts = count.table))
    
    # Quality Control with scater
    is.mito <- grepl("^mt-", rownames(sce))
    # quality control with caluclateQCMetrics
    sce <- calculateQCMetrics(sce, exprs_values = "counts",
                              feature_controls = list(MT = is.mito))
    
    # Drop below and above 3 median absolute deviations: libsize, mitochondrial gene count, expressed genes
    libsize.drop <- isOutlier(sce$total_counts, nmads = MAD, type = "both", log = T)
    
    feature.drop <- isOutlier(sce$total_features_by_counts, nmads = MAD, type = "both", log = T)
    
    mito.drop <- isOutlier(sce$pct_counts_MT, nmads = MAD, type = "higher")
    
    sce <- sce[, !(libsize.drop | feature.drop | mito.drop)]

    sce.list[[sample]] <- sce
    print(paste(sample, "Done"))
  }
  return(sce.list)
}

sce.list <- qcControl(count.tables = count.tables, MAD = 3, sample.names = names(count.tables))


geneFiltering <- function(
  object,
  gene.expr,
  sample.names = names(object)
){
  if (is.list(object)){
    obj.list <- list()
    i <- 1
    for (matrix in object){
      filter.gene.expr.cells <- ncol(sce)*gene.expr
      keep.features <- Matrix::rowSums(counts(matrix) > 0) >= filter.gene.expr.cells
      matrix <- matrix[keep.features, ]
      
      sample <- sample.names[i]
      obj.list[[sample]] <- matrix
      i <- i + 1
      print(paste(sample, "Done"))
    }
    return(obj.list)
  } else {
    matrix <- object
    filter.gene.expr.cells <- ncol(matrix)*gene.expr
    keep.features <- Matrix::rowSums(counts(matrix) > 0) >= filter.gene.expr.cells
    matrix <- matrix[keep.features, ]
    return(matrix)
  }
}

WriteCountMetrices <- function(
  sce.list,
  log.normalize = TRUE,
  filter.genes = TRUE,
  sample.names = names(sce.list)
){
  i <- 1
  for (sce in sce.list){
    sample <- sample.names[i]
    i <- i + 1
    if (log.normalize){
      seurat.object <- as.Seurat(sce, data = "counts")
      seurat.object <- NormalizeData(seurat.object, 
                                     normalization.method = "LogNormalize", 
                                     scale.factor = 1000, 
                                     verbose = FALSE)
      log.counts <- GetAssayData(seurat.object)
      if (filter.genes){
        log.counts <- SingleCellExperiment(assays = list(counts = log.counts))
        log.counts <- geneFiltering(log.counts, gene.expr = 0.001)
        log.counts <- counts(log.counts)
      }
    }
    if (filter.genes){
      sce <- geneFiltering(object = sce, gene.expr = 0.001)
    }
    features <- rownames(sce)
    barcodes <- colnames(sce)
    write(features, file = paste0(sample, ".features.txt"))
    write(barcodes, file = paste0(sample, ".barcodes.txt"))
    writeMM(counts(sce), file = paste0(sample, ".raw.counts.mtx"))
    writeMM(log.counts, file = paste0(sample, ".log.counts.mtx"))
    print(paste(sample, "Done"))
  }
}

WriteCountMetrices(sce.list = sce.list, log.normalize = TRUE, filter.genes = TRUE)
