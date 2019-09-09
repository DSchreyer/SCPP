library(Seurat)
library(Matrix)

input.dir <- "/home/daniel/master_thesis/bassoon_data/Output/test_qc_120//"
output.dir <- "/home/daniel/master_thesis/bassoon_data/Output/test_ct_ident/"
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
setwd(input.dir)

samples <- paste0("Bsn.", seq(221931, 221940))

ReadCountMetrices <- function(
  dir,
  log.counts = TRUE,
  samples,
  counts = "log.counts.mtx",
  feature = "features.txt",
  barcode = "barcodes.txt"
){
  seurat.objects <- list()
  ls <- list.files(dir)
  for (sample in samples){
    ind <- grep(pattern = sample, ls)
    count.file <- grep(counts, ls[ind], value = TRUE)
    feature.file <- grep(feature, ls[ind], value = TRUE)
    barcode.file <- grep(barcode, ls[ind], value = TRUE)
    count.table <- readMM(count.file)
    features <- read.csv(feature.file, stringsAsFactors = F, header = F)
    barcodes <- read.csv(barcode.file, stringsAsFactors = F, header = F)
    rownames(count.table) <- features$V1
    colnames(count.table) <- barcodes$V1
    seurat.object <- CreateSeuratObject(project = sample, counts = count.table)
    if (!log.counts){
      seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 1000)
    }
    seurat.objects[[sample]] <- seurat.object
  }
  return(seurat.objects)
}

seurat.object.list <- ReadCountMetrices(dir = input.dir, log.counts = TRUE, samples = samples)

# ### TRANSFER CELL LABELS TO OTHER OBJECTS
ref <- seurat.object.list$Bsn.221932
ref <- FindVariableFeatures(ref, nfeatures = 1000)
ref <- ScaleData(ref)
ref <- RunPCA(ref)
ref <- FindNeighbors(ref, dims = 1:10)
ref <- FindClusters(ref, resolution = 0.03)
markers <- FindAllMarkers(ref, only.pos = T)
anchors <- FindTransferAnchors(ref, seurat.object.list$Bsn.221931)
predictions <- TransferData(anchorset = anchors, refdata = ref$seurat_clusters, 
                            dims = 1:10)
queries <- AddMetaData(seurat.object.list$Bsn.221931, metadata = predictions)
Idents(queries) <- queries$predicted.id
FindMarkers(queries, ident.1 = "0", ident.2 = "2")

new.so.list <- list()
for (object in seurat.object.list){
  name <- object@project.name
  if (ncol(object) < 500){
    next
  }
  anchors <- FindTransferAnchors(ref, object)
  predictions <- TransferData(anchorset = anchors, refdata = ref$seurat_clusters, dims = 1:10)
  object <- AddMetaData(object, metadata = predictions)
  Idents(object) <- object$predicted.id
  new.so.list[[name]] <- object
}

for (so in new.so.list){
  # print(FindMarkers(so, ident.1 = "0", ident.2 = "2"))
  print(table(Idents(so)))
}


## INTEGRATE DATA
features <- SelectIntegrationFeatures(object.list = seurat.object.list, nfeatures = 1000)
reference <- which(names(seurat.object.list) == "Bsn.221932")
subset <- seurat.object.list
subset$Bsn.221935 <- NULL 


so.anchors <- FindIntegrationAnchors(object.list = subset, anchor.features = features )
so.integrated <- IntegrateData(anchorset = so.anchors)
so.integrated <- ScaleData(so.integrated)
so.integrated <- RunPCA(so.integrated)
so.integrated <- RunUMAP(so.integrated, dims = 1:10)
so.integrated <- FindNeighbors(so.integrated, dims = 1:10)
# resolution 0.1 works
so.integrated <- FindClusters(so.integrated, resolution = 0.15)
FindAllMarkers(so.integrated, only.pos = T)


DimPlot(so.integrated, group.by = "celltype", reduction = "umap")
save(so.integrated, file = "/home/daniel/master_thesis/bassoon_data/Output/test_qc_120/post_integration_so.Robj")

  
Rbad <- seurat.object.list$Bsn.221938
anchors <- FindTransferAnchors(ref, bad)
predictions <- TransferData(anchorset = anchors, refdata = ref$seurat_clusters, dims = 1:10)
bad <- AddMetaData(bad, metadata = predictions)
Idents(bad) <- bad$predicted.id
# m <- FindAllMarkers(bad, only.pos = T)
FindMarkers(bad, ident.1 = "0", ident.2 = "2")

# CellTypeScoring function is based on ScoreCellCycle function from Seurat
# Input are marker genes of 2 different cell types and it scores each cell based on the marker gene expression level
# With these scores it defines a cell type for each cell --> Cell Type 1, Cell Type 2 or Unknown
# Unknown cells have no definite marker gene expression or no marker is expressed in these cells
CellTypeScoring <- function(
  object,
  type.1.features,
  type.2.features,
  type.1,
  type.2,
  name,
  set.ident = TRUE
){
  features <- list("Type.1" = type.1.features, "Type.2" = type.2.features)
  object.ct <- AddModuleScore(
    object = object, 
    features = features,
    name = name,
    ctrl = min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  )
  ct.columns <- grep(pattern = name, x = colnames(x = object.ct[[]]), value = TRUE)
  ct.scores <- object.ct[[ct.columns]]
  rm(object.ct)
  
  assignments <- apply(
    X = ct.scores,
    MARGIN = 1,
    FUN = function(scores, first = type.1, second = type.2, null = 'Unknown') {
      if (all(scores < 0)) {
        return(null)
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return('Unknown')
        } else {
          return(c(first, second)[which(x = scores == max(scores))])
        }
      }
    }
  )
  ct.scores <- merge(x = ct.scores, y = data.frame(assignments), by = 0)
  type.1.score <- paste(type.1, "Score", sep = ".")
  type.2.score <- paste(type.2, "Score", sep = ".")
  colnames(x = ct.scores) <- c("rownames",
                               type.1.score,
                               type.2.score,
                               "CellType")
  rownames(x = ct.scores) <- ct.scores$rownames
  ct.scores <- ct.scores[, c(type.1.score, type.2.score, "CellType")]
  object[[colnames(x = ct.scores)]] <- ct.scores
  if (set.ident) {
    object[['old.ident']] <- Idents(object = object)
    Idents(object = object) <- "CellType"
  }
  return(object)
}

one <- FindVariableFeatures(seurat.object.list$Bsn.221932, nfeatures = 1000)
one <- ScaleData(one)
one <- RunPCA(one)
one <- FindNeighbors(one, dims = 1:6)
one <- FindClusters(one, resolution = 0.05)
FindAllMarkers(one, only.pos = T) %>% filter(p_val_adj < 0.05)
table(Idents(one))
one <- RunTSNE(one)
DimPlot(one, reduction = "tsne")

anchors <- FindTransferAnchors(one, seurat.object.list$Bsn.221931)
predictions <- TransferData(anchorset = anchors, refdata = one$seurat_clusters, 
                            dims = 1:6)
queries <- AddMetaData(seurat.object.list$Bsn.221931, metadata = predictions)

# Alex approved Markers
rod.markers <- c("Rho", "Nt5e", "Nr2e3", "Gnat1", "Cngb1", "Crx", "Nrl", "Pde6a")
cone.markers <- c("Opn1sw", "Opn1mw", "Arr3", "Gnat2", "Pde6h", "Gngt2", "Gnb3", "Rora")


for (seurat.object in seurat.object.list){
  seurat.object <- CellTypeScoring(object = seurat.object, type.1 = "Cone", type.2 = "Rod", type.1.features = cone.markers,
                         type.2.features = rod.markers, name = "Cell Type")
  celltype <- as.vector(seurat.object@meta.data$CellType)
  count.table <- GetAssayData(seurat.object)
  barcodes <- colnames(count.table)
  features <- rownames(count.table)
  sample <- seurat.object@project.name
  write(features, file =  paste0(output.dir, sample, ".features.txt"))
  write(barcodes, file =  paste0(output.dir, sample, ".barcodes.txt"))
  write(celltype, file =  paste0(output.dir, sample, ".celltype.txt"))
  writeMM(count.table, file = paste0(output.dir, sample, ".logcounts.mtx"))
  print(paste(sample, "Done"))
}

Idents(seurat.object) <- celltype
FindAllMarkers(seurat.object)

rm(seurat.object, celltype, count.table, barcodes, features, sample)
# write feature list, barcode list, and raw count table to a file
