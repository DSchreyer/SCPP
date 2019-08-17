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
data.dir = "/home/daniel/master_thesis/bassoon_data/Cellranger output/Aggr/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data.dir, gene.column = 2, unique.features = TRUE)


umi.count <- 200
filter.umi.count <- Matrix::colSums(data) >= umi.count
data <- data[ , filter.umi.count]

# gene.count <- 80
# filter.gene.count <- Matrix::colSums(data > 0) >= gene.count
# data <- data[ , filter.gene.count]

all.count.tables <- lapply(c(as.character(seq(1,10))),
                           function(x) data[, grep(x, colnames(data))])
names(all.count.tables) <- paste0("Bsn.", seq(221931,221940))

i <- 1
for (count.table in all.count.tables){
  sample <- names(all.count.tables)[i]
  i <- i + 1
  filter.feature.expr.cells <- ncol(count.table)*0.001
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
  libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "both", log = T)
  
  feature.drop <- isOutlier(sce$total_features_by_counts, nmads = 3, type = "both", log = T)
  
  mito.drop <- isOutlier(sce$pct_counts_MT, nmads = 3, type = "higher")
  
  sce <- sce[, !(libsize.drop | feature.drop | mito.drop)]
  print(dim(sce))
  features <- rownames(sce)
  write(features, file = paste0("/home/daniel/master_thesis/bassoon_data/Output/post_qc/", sample, ".features.txt"))
  writeMM(count.table, file = paste0("/home/daniel/master_thesis/bassoon_data/Output/post_qc/", sample, ".raw.counts.mtx"))
}

so <- as.Seurat(sce, assay = "logcounts")

sce
# deconvolution method normalization
# larger dataset --> quick cluster
clusters <- quickCluster(sce, min.size = 100, use.ranks = FALSE)
sce <- computeSumFactors(sce, min.mean = 0.1, clusters = clusters)
