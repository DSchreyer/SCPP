library(Seurat)
library(Matrix)

input.dir <- "/home/daniel/master_thesis/bassoon_data/Output/post_qc/"
output.dir <- "/home/daniel/master_thesis/bassoon_data/Output/new_ct_ident/"
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
  first <- TRUE
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
    if (first){
      seurat.objects <- seurat.object
      first = FALSE
    } else {
      seurat.objects <- merge(seurat.objects, seurat.object)
    }
  }
  return(seurat.objects)
}

seurat.objects <- ReadCountMetrices(dir = input.dir, log.counts = TRUE, samples = samples)

# ### TRANSFER CELL LABELS TO OTHER OBJECTS
ref <- subset(seurat.objects, subset = orig.ident == "Bsn.221932")
ref <- FindVariableFeatures(ref, nfeatures = 1500)
ref <- ScaleData(ref)
ref <- RunPCA(ref)
ElbowPlot(ref)
ref <- FindNeighbors(ref, dims = 1:7)
# 0.1 works
ref <- FindClusters(ref, resolution = 0.1)
m <- FindAllMarkers(ref, only.pos = T)
# 2 is cones 0 is rods
ref <- RunTSNE(ref, dims = 1:7)
DimPlot(ref, reduction = "tsne")

seurat.object.list <- SplitObject(seurat.objects, split.by = "orig.ident")
seurat.object.list$Bsn.221932 <- NULL
seurat.object.list$Bsn.221935 <- NULL
new.so.list <- list()
for (object in seurat.object.list){
  name <- unique(object$orig.ident)
  anchors <- FindTransferAnchors(ref, object)
  predictions <- TransferData(anchorset = anchors, 
                              refdata = ref$seurat_clusters)
  object <- AddMetaData(object, metadata = predictions)
  Idents(object) <- object$predicted.id
  new.so.list[[name]] <- object
}
new.so.list$Bsn.221932 <- ref
save(new.so.list, file = "/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/so.transfer.Robj")
load("/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/so.transfer.Robj")

for (so in new.so.list){
  print(unique(so$orig.ident))
  print(table(Idents(so)))
}

celltypes <- c("rod", "bipolar_cell", "cone", "Müller_glia", "RPE")
names(celltypes) <- c("0", "1", "2", "3", "4")

for (seurat.object in new.so.list){
  # remove cells with a prediction score below 70%
  # try(seurat.object <- seurat.object[, seurat.object$prediction.score.max >= 0.7])
  
  seurat.object <- RenameIdents(seurat.object, celltypes)
  seurat.object$celltype <- Idents(seurat.object)
  celltype <- as.vector(seurat.object$celltype)
  count.table <- GetAssayData(seurat.object)
  barcodes <- colnames(count.table)
  features <- rownames(count.table)
  sample <- unique(seurat.object$orig.ident)
  write(features, file =  paste0(output.dir, sample, ".features.txt"))
  write(barcodes, file =  paste0(output.dir, sample, ".barcodes.txt"))
  write(celltype, file =  paste0(output.dir, sample, ".celltype.txt"))
  writeMM(count.table, file = paste0(output.dir, sample, ".logcounts.mtx"))
  print(paste(sample, "Done"))
}


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
