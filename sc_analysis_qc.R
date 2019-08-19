library(scater)
library(scran)
library(Seurat)
library(Matrix)

# laptop
# path.to.cellranger.outs <- paste("/home/daniel/master_thesis/bassoon_data/Cellranger output/", sample, sep = "")
# cluster
# path.to.cellranger.outs <- paste("/media/data2/Daniel/cellranger_output/", sample, "/outs/", sep = "/")

# Read in sparse matrix with cellranger output
data.dir = paste(path.to.cellranger.outs, "/raw_feature_bc_matrix/", sep = "/")

# cluster
# data.dir <- "/media/data2/Daniel/cellranger_aggr/AGG_Bsn/outs/raw_feature_bc_matrix/"
# output.dir <- "media/data2/Daniel/cellranger_aggr/R_output/post_qc/"

# laptop
data.dir <- "/home/daniel/master_thesis/bassoon_data/Cellranger output/Aggr/raw_feature_bc_matrix/"
output.dir <- "/home/daniel/master_thesis/bassoon_data/Cellranger output/Aggr/R_output"

# Create Output Directory
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
setwd(output.dir)

n.samples <- 10
sample.names <- paste0("Bsn.", seq(221931,221940))


generateCountTables <- function(
  aggr.dir,
  umi.count = 100,
  sample.names = seq(1,n.samples)
){
  data <- Read10X(data.dir = aggr.dir, gene.column = 2, unique.features = TRUE)
  
  umi.count <- umi.count
  filter.umi.count <- Matrix::colSums(data) >= umi.count
  data <- data[ , filter.umi.count]
  
  all.count.tables <- lapply(c(as.character(seq(1,n.samples))),
                             function(x) data[, grep(x, colnames(data))])
  names(all.count.tables) <- sample.names
  rm(data)
  return(all.count.tables)
}

count.tables <- GenerateCountTables(aggr.dir = data.dir, umi.count = 120, sample.names = sample.names)


qcControl <- function(
  count.tables,
  MAD = 3,
  feature.expr = 0.01
){
  i <- 1
  sce.list <- list()
  for (count.table in count.tables){
    sample <- names(count.tables)[i]
    i <- i + 1
    filter.feature.expr.cells <- ncol(count.table)*feature.expr
    keep.features <- Matrix::rowSums(count.table > 0) >= filter.feature.expr.cells
    count.table <- count.table[keep.features, ]
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
  
    # Normalize counts with Deconvolution method
    sce <- computeSumFactors(sce, min.mean = 0.1)
    
    # Remove cells with negative size factors 
    neg.size.factors <- sizeFactors(sce) < 0 
    sce <- sce[, !neg.size.factors]
    
    sce <- normalize(sce, return_log = FALSE)
    
    # Transform counts to natural log counts
    logcounts(sce) <- as(log(normcounts(sce) + 1), "sparseMatrix")
    
    features <- rownames(sce)
    # write count matrix and feature list to file
    sce.list[[sample]] <- sce
  }
  rm(features, libsize.drop, feature.drop, mito.drop, is.mito, keep.features, sample,
     filter.feature.expr.cells, i, count.table)
  return(sce.list)
}

sce <- qcControl(count.tables = count.tables, MAD = 3, feature.expr = 0.001 )


# write files
write(features, file = paste0(sample, ".features.txt"))
writeMM(normcounts(sce), file = paste0(sample, ".log.counts.mtx"))
print(paste(sample, "Done"))


sce
# deconvolution method normalization
# larger dataset --> quick cluster
sce <- normalize(sce, return_log = FALSE)
sce <- as(log(normcounts(sce) + 1), "sparseMatrix")
