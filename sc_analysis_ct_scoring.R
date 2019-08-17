library(Seurat)
library(scater)
library(scran)
library(dplyr)
library(tidyr)
library(Matrix)


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

# normalize with scater function normalize using library sizes as size factors
# sce <- normalize(sce)
# plotHighestExprs(sce)
# plotExpression(sce, features = c("Rho", "Opn1sw", "Opn1mw", "Arr3", "Nrl", "Gnat2"))


# Load into Seurat objet
Bsn <- as.Seurat(sce, data = "counts", assay = "RNA", counts = "counts", project = sample)

raw.counts <- Bsn[["RNA"]]@counts

# normalize the data
Bsn <- NormalizeData(Bsn, normalization.method = "LogNormalize", scale.factor = 1000)

# Alex approved Markers
rod.markers <- c("Rho", "Nt5e", "Nr2e3", "Gnat1", "Cngb1")
cone.markers <- c("Opn1sw", "Opn1mw", "Arr3", "Gnat2", "Pde6h", "Gngt2", "Gnb3")


Bsn <- CellTypeScoring(object = Bsn, type.1 = "Cone", type.2 = "Rod", type.1.features = cone.markers,
                   type.2.features = rod.markers, name = "Cell Type")

CellType <- as.vector(Bsn@meta.data$CellType)
barcodes <- colnames(raw.counts)
features <- rownames(raw.counts)

cell.type <- cbind(barcodes, CellType)

# write feature list, barcode list, and raw count table to a file
write(features, file = paste0("/home/daniel/master_thesis/bassoon_data/Cellranger output/Samples_ct_scored/", sample, ".features.txt"))
write.table(cell.type, file = paste0("/home/daniel/master_thesis/bassoon_data/Cellranger output/Samples_ct_scored/", sample, ".barcode.celltype.csv"), sep = ",", quote = F, row.names = F)
writeMM(raw.counts, file = paste0("/home/daniel/master_thesis/bassoon_data/Cellranger output/Samples_ct_scored/", sample, ".logcounts.mtx"))
}






