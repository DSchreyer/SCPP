library(Seurat)
library(Matrix)
library(ggplot2)
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
                               "ct_score_celltype")
  rownames(x = ct.scores) <- ct.scores$rownames
  ct.scores <- ct.scores[, c(type.1.score, type.2.score, "ct_score_celltype")]
  object[[colnames(x = ct.scores)]] <- ct.scores
  if (set.ident) {
    object[['old.ident']] <- Idents(object = object)
    Idents(object = object) <- "ct_score_celltype"
  }
  return(object)
}

input.dir <- "/home/daniel/master_thesis/bassoon_data/Output/final_qc/"
output.dir <- "/home/daniel/master_thesis/bassoon_data/Output/final_ct_class/"
analysis.out <- "/home/daniel/master_thesis/bassoon_data/Output/final_downstream_analysis/"
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

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
  ls <- list.files(dir, full.names = T)
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

# load("/home/daniel/master_thesis/bassoon_data/Output/final_ct_class//final.new.so.list.Robj")

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
celltypes <- c("rod", "bipolar cell", "cone", "Müller glia", "RPE")
names(celltypes) <- c("0", "1", "2", "3", "4")
ref$celltypes <- Idents(ref)
ref <- RenameIdents(ref, celltypes)
ref <- RunTSNE(ref, dims = 1:7)
# pdf("reference_plots.pdf")
tsne <- DimPlot(ref, reduction = "tsne")

ggsave(plot = tsne, filename = "/home/daniel/master_thesis/bassoon_data/Output/Tables_Graphs/tsne.png")

# DimPlot(ref, reduction = "pca")

seurat.object.list <- SplitObject(seurat.objects, split.by = "orig.ident")
seurat.object.list$Bsn.221932 <- NULL
seurat.object.list$Bsn.221935 <- NULL
new.so.list <- list()
for (object in seurat.object.list){
  name <- unique(object$orig.ident)
  anchors <- FindTransferAnchors(ref, object, dims = 1:7)
  predictions <- TransferData(anchorset = anchors,
                              refdata = ref$seurat_clusters, dims = 1:7)
  object <- AddMetaData(object, metadata = predictions)
  Idents(object) <- object$predicted.id
  new.so.list[[name]] <- object
}
new.so.list$Bsn.221932 <- ref
# save(new.so.list, file = "/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/final.new.so.list.Robj")
# load("/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/so.transfer.Robj")

first = TRUE
for (so in new.so.list){
  print(unique(so$orig.ident))
  print(table(Idents(so)))
  if (first){
    seurat.objects <- so
    first = FALSE
  } else {
    seurat.objects <- merge(seurat.objects, so)
  }
}
'sos <- sos[, sos$prediction.score.max >= 0.50 | is.na(sos$prediction.score.max)]
celltypes <- c("rod", "bipolar_cell", "cone", "Müller_glia", "RPE")
names(celltypes) <- c("0", "1", "2", "3", "4")
sos <- RenameIdents(sos, celltypes)
sos$celltype <- Idents(sos)
'
# seurat.objects <- sos


seurat.objects <- seurat.objects[, seurat.objects$prediction.score.max >= 0.95 | 
                                   is.na(seurat.objects$prediction.score.max)]

celltypes <- c("rod", "bipolar cell", "cone", "Müller glia", "RPE")
names(celltypes) <- c("0", "1", "2", "3", "4")
seurat.objects <- RenameIdents(seurat.objects, celltypes)
seurat.objects$celltype <- Idents(seurat.objects)

seurat.objects.keep <- seurat.objects
seurat.objects <- subset(seurat.objects.keep, subset = celltype == "rod" | celltype == "cone")

################

# generate table with number of cells per cell type for each sample
ct_class_info <- table(seurat.objects$orig.ident, seurat.objects$celltype)
write.table(file = paste0(analysis.out, "/ct_class_info.csv"), x = ct_class_info, sep = ",")

seurat.object.list <- SplitObject(seurat.objects, split.by = "orig.ident")
for (seurat.object in seurat.object.list){
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







# seurat.objects <- CellTypeScoring(object = seurat.objects, type.1 = "cone", type.2 = "rod", type.1.features = cone.markers,
#                                 type.2.features = rod.markers, name = "Cell Type")
# cones <- subset(seurat.objects, subset = celltype == "cone" & ct_score_celltype == "cone")
# rods <- subset(seurat.objects, subset = celltype == "rod" & ct_score_celltype == "rod")

# seurat.objects <- merge(cones, rods)




seurat.object <- CellTypeScoring(object = seurat.objects, type.1 = "Cone", type.2 = "Rod", type.1.features = cone.markers,
                                 type.2.features = rod.markers, name = "Cell Type")




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

ref <- new.so.list$Bsn.221932
ref <- RenameIdents(ref, celltypes)


