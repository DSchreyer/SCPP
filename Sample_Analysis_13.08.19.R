library(Seurat)
library(scater)
library(scran)
library(dplyr)
library(ggpubr)
library(UpSetR)
library(ggsci)
library(tidyr)
library(cluster)
library(fastcluster)

# Project name: Don't use special symbols for project name: /?^%$ --> Is used for file names

# all.samples <- seq(221931, 221940, 1)
# for (sample in all.samples){
sample <- 221940
# laptop
path.to.cellranger.outs <- paste("/home/daniel/master_thesis/bassoon_data/Cellranger output/", sample, sep = "")
# cluster
# path.to.cellranger.outs <- paste("/media/data2/Daniel/cellranger_output/", sample, "/outs/", sep = "/")
project <- as.character(sample)

f <- paste0(project, ".markers.supervised.txt")
write(paste("Project:", project, sep = " "), file = f)

# Read in sparse matrix with cellranger output
# data.dir = path to cellranger output matrix.mtx, features.tsv, barcodes.tsv
data.dir = paste(path.to.cellranger.outs, "/raw_feature_bc_matrix/", sep = "/")
data <- Read10X(data.dir = data.dir, gene.column = 2, unique.features = TRUE)


filter.umi.count <- 120
keep.cells <- Matrix::colSums(data) >= filter.umi.count
filter.feature.expr.cells <- sum(keep.cells)*0.005
keep.features <- Matrix::rowSums(data > 0) >= filter.feature.expr.cells


write(paste("Data:", data.dir, sep = " "), file = f, append = TRUE)
write(paste("Before:", dim(data), c("genes", "cells"), sep = " "), file = f, append = TRUE)
write(paste("Filter Umi Count:", filter.umi.count, "| filter features:", filter.feature.expr.cells, sep = " "),
      file = f, append = TRUE)

data <- data[keep.features, keep.cells]
dim(data)
write(paste("Filter cells and genes:", dim(data), c("genes","cells"), "left.", sep = " "), file = f, append = TRUE)


# Create SingleCellExperiment object with sparse matrix
sce <- SingleCellExperiment(assays = list(counts = data))

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
dim(sce)
write(paste("Filter after libsize, feature count, mt expr:", dim(sce), c("genes", "cells"), "left.", sep = " "), file = f, append = TRUE)


# normalize with scater function normalize using library sizes as size factors
sce <- normalize(sce)
# pdf(paste(project, ".pdf", sep = ""))
# plotHighestExprs(sce)
# plotExpression(sce, features = c("Rho", "Opn1sw", "Opn1mw", "Arr3", "Nrl", "Gnat2"))



# Load into Seurat objet
Bsn <- as.Seurat(sce, data = "counts", assay = "RNA", counts = "counts", project = project)

# normalize the data
Bsn <- NormalizeData(Bsn, normalization.method = "LogNormalize", scale.factor = 1000)
Bsn <- FindVariableFeatures(Bsn, selection.method = "vst", nfeatures = 500)

sce <- as.SingleCellExperiment(Bsn)
sce.logcounts <- logcounts(sce)

### markers used for hierarchical clustering
# rod.markers <- c("Rho", "Nt5e", "Nr2e3", "Cngb1")
# cone.markers <- c("Opn1sw", "Opn1mw", "Arr3", "Gnat2", "Pde6h", "Gnb3", "Gngt2")

# Alex approved Markers
rod.markers <- c("Rho", "Nt5e", "Nr2e3", "Gnat1", "Cngb1")
cone.markers <- c("Opn1sw", "Opn1mw", "Arr3", "Gnat2", "Pde6h", "Gngt2", "Gnb3")


### CellCycleScoring Test ###

GetCellType <- function(
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

Bsn <- GetCellType(object = Bsn, type.1 = "Cone", type.2 = "Rod", type.1.features = cone.markers,
                   type.2.features = rod.markers, name = "Cell Type")

FindAllMarkers(Bsn, only.pos = T)
table(Idents(Bsn))
Bsn$metadata

## Heatmap
markers <- c(rod.markers, cone.markers)

subset <- as.data.frame(t(as.matrix(GetAssayData(Bsn, slot = "data")[markers, ])))

ident <- Idents(Bsn)
subset$ident <- ident
cells <- rownames(subset)
subset$cells <- cells
subset <- subset %>% dplyr::arrange(ident)

rownames(subset) <- subset$cells
ident <- as.data.frame(subset$ident)
rownames(ident) <- subset$cells

subset <- select(subset, -c(ident, cells))

marker.genes <- data.frame(Marker = factor(rep(c("Rod", "Cone"), times = c(5,7))))
rownames(marker.genes) <- colnames(subset)

library(pheatmap)
pheatmap(mat = subset, cluster_cols = F,
         cluster_rows = F,
         annotation_row = ident, 
         annotation_colors = colors[1],
         show_rownames = F,
         annotation_col = marker.genes)


#check if all markers occur in the data set
markers <- intersect(rownames(sce.logcounts), markers)

subset <- sce.logcounts[markers, ]

dim(subset)

hc <- hclust(dist(t(subset)), method = "complete")

clusters <- cutree(hc, k = 5)
print(table(clusters))
Bsn <- SetIdent(Bsn, value = clusters)

all.markers <- FindAllMarkers(Bsn, logfc.threshold = 0.15)

print(all.markers)
table(Idents(Bsn))
file = "/media/data2/Daniel/R/all.markers.txt"
write(project, file = file, append = TRUE)
write.table(table(grp), file = file, append = TRUE, quote = F, col.names = T)
write.table(all.markers, file = file, append = TRUE, quote = F, col.names = T)
# }

FindMarkers(Bsn, ident.1 = 1, ident.2 = 3)


library(dendextend)





dist <- dist(as.matrix(t(test)))
t <- hclust(dist)
dend <- as.dendrogram(t)
dend <- color_branches(dend, k = 4)

grp <- cutree(dend, k = 5)

Bsn <- SetIdent(object = Bsn, value = grp)

FindAllMarkers(Bsn)

FindMarkers(Bsn, ident.1 = 5, ident.2 = 1)


n.cells <- ncol(Bsn)
n <- as.integer(n.cells/10000)
if (n.cells/n < 2000){
  n <- n+1
}

for (i in 1:n){
  if (i < n){
    
  }
  else
  j <- n*10000
}




sce <- as.SingleCellExperiment(Bsn)

sce.markers <- sce[markers, ]
rowData(sce.markers) <- rownames(sce.markers)
colnames(rowData(sce.markers)) <- "feature_symbol"

counts(sce.markers) <- as.matrix(counts(sce.markers))
logcounts(sce.markers) <- as.matrix(logcounts(sce.markers))

sce.markers <- sc3_prepare(sce.markers, kmeans_nstart = 50)
sce.markers <- sc3_calc_dists(sce.markers)


sc3_kmeans(sce, k = 2:5)


distance <-SC3::consensus_matrix(sce)

### Test hiearchical clustering ###
library(RCA)
library(WGCNA)
data <- Bsn@assays$RNA[]

test <- data[markers, ]

library(dendextend)
dist <- dist(as.matrix(t(test)))
t <- hclust(dist)
dend <- as.dendrogram(t)
dend <- color_branches(dend, k = 4)

grp <- cutree(dend, k = 5)

Bsn <- SetIdent(object = Bsn, value = grp)

FindAllMarkers(Bsn)

FindMarkers(Bsn, ident.1 = 5, ident.2 = 1)

# works

plot(dend)

construct <- dataConstruct(test)

wat <- featureConstruct(construct)


cluster <- cellClust(construct, method = "hclust", deepSplit_wgcna = 3)

sum(Bsn@assays$RNA["Bsn", ])/ncol(Bsn@assays$RNA)



# Find variable Features, top 2000
Bsn <- FindVariableFeatures(Bsn, selection.method = "vst", nfeatures = 500)


# Scale Data for PCA
Bsn <- ScaleData(Bsn)

## load into monocle objet
library(monocle)
library(garnett)
library(org.Mm.eg.db)
cds <- as.CellDataSet(Bsn)
cds <- estimateSizeFactors(cds)

marker.file.path <- system.file("extdata", 
                                "/home/daniel/master_thesis/bassoon_data/Cellranger output/markers.txt",
                                package = "garnett")
marker.file.path <- "/home/daniel/master_thesis/bassoon_data/Cellranger output/markers.txt"


marker.check <- check_markers(cds, marker.file.path, db = org.Mm.eg.db)
plot_markers(marker.check)

# train classifier

cds.classifier <- train_cell_classifier(cds = cds, 
                                        marker_file = marker.file.path, 
                                        db = org.Mm.eg.db)
# cell cycle scoring with rod and cone markers
rod.markers <- c("Rho", "Gnat1", "Nr2e3")
cone.markers <- c("Opn1sw", "Opn1mw", "Arr3", "Gnat2", "Pde6h", "Gnb3")

score <- CellCycleScoring(Bsn, s.features = rod.markers, g2m.features = cone.markers,
                          set.ident = T)


RidgePlot(Bsn, features = c("Rho", "Opn1sw", "Arr3", "Opn1mw", "Gnat1", "Gnat2"))
head(score@meta.data)

Idents(score)


# Run PCA
features = c("Rho", "Opn1sw", "Opn1mw", "Arr3", "Nt5e", "Gnat2", "Gngt2", "Pde6h", "Nr2e3", "Gnat1", "Cngb1")
features = c("Rho", "Opn1sw", "Opn1mw", "Arr3", "Gnat2", "Gnat1")

Bsn <- RunPCA(Bsn)


PCAPlot(Bsn)

VizDimLoadings(Bsn)
ElbowPlot(Bsn)

# Identify PCs to take

# pdf(paste(project, ".2.pdf", sep = ""))
# Cluster cells to 2,3,4,5 different clusters
# 3 cell clusters optimal
# check ElbowPlot for high variance pc

seq <- c(seq(0.01,0.1,0.01), seq(0.12, 0.34, 0.02))

dims <- 5:10
for (dim in dims){
  Bsn <- FindNeighbors(Bsn, dims = 1:dim, verbose = F)
  for (res in seq){
    Bsn <- FindClusters(Bsn, resolution = res, verbose = F)
    markers <- FindAllMarkers(Bsn, only.pos = T, logfc.threshold = 0.1, verbose = F)
    n.cluster <- length(levels(Idents(Bsn)))
    if (ncol(markers) > 0){
      markers <- markers %>% select(gene, p_val, avg_logFC, p_val_adj, cluster)
    }
    n.cells.cluster <- table(Idents(Bsn))
    head <- paste("Dims:", dim, "|", "resolution:", res, "|", "cluster:", n.cluster, sep = " ")
    print(head)
    if (ncol(markers) > 0){
      print(markers %>% group_by(cluster) %>% top_n(-4, p_val_adj))
    }
    print(n.cells.cluster)
    write(head, file = f, append = T)
    write.table(markers, file = f, append = TRUE, quote = F, col.names = T, row.names = F)
    write("Number of cells in each cluster:", file = f, append = TRUE)
    write.table(n.cells.cluster, file = f, append = TRUE, quote = F, col.names = F)
    write("", file = f, append = T)
  }
}


f2 <- paste0(project, "umi_count_", filter.umi.count, ".only_variable_genes.markers.txt")
write(paste("Project:", project, sep = " "), file = f2)
write(paste("Data:", data.dir, sep = " "), file = f2, append = TRUE)

Bsn <- FindNeighbors(Bsn, features = VariableFeatures(object = Bsn), verbose = F)

for (res in seq){
  Bsn <- FindClusters(Bsn, resolution = res,verbose = F)
  markers <- FindAllMarkers(Bsn, only.pos = T, logfc.threshold = 0.1, verbose = F)
  n.cluster <- length(levels(Idents(Bsn)))
  if (ncol(markers) > 0){
    markers <- markers %>% select(gene, p_val, avg_logFC, p_val_adj, cluster)
  }
  n.cells.cluster <- table(Idents(Bsn))
  head <- paste("Only Variable Features. Resolution:", res, "|", "cluster:", n.cluster, sep = " ")
  print(head)
  if (ncol(markers) > 0){
    print(markers %>% group_by(cluster) %>% top_n(-4, p_val_adj))
  }
  print(n.cells.cluster)
  write(head, file = f2, append = T)
  write.table(markers, file = f2, append = TRUE, quote = F, col.names = T, row.names = F)
  write("Number of cells in each cluster:", file = f2, append = TRUE)
  write.table(n.cells.cluster, file = f2, append = TRUE, quote = F, col.names = F)
  write("", file = f2, append = T)
}



# Run UMAP and tSNE
Bsn <- RunUMAP(Bsn, dims = dims)
Bsn <- RunTSNE(Bsn, dims = dims)

DimPlot(Bsn, reduction = "umap")
DimPlot(Bsn, reduction = "tsne")

rod.vs.cones <- FindMarkers(Bsn, ident.1 = 0, ident.2 = 2, min.pct = 0.25)
write.table(rod.vs.cones, file = paste(project, "rod_vs_cones_marker.txt", sep = "_"))

cluster.ids <- c("Rod", "undefined", "Cone")
names(cluster.ids) <- levels(Bsn)
Bsn <- RenameIdents(Bsn, cluster.ids)


DimPlot(Bsn, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
DimPlot(Bsn, reduction = "tsne", label = TRUE, pt.size = 0.1) + NoLegend()


markers <- c("Rho", "Opn1sw", "Opn1mw")
l.markers <- length(markers)

marker.table <- GetAssayData(Bsn, slot = "counts")[markers, ]
cluster <- Bsn$RNA_snn_res.0.04
marker.table <- rbind(as.matrix(marker.table), cluster)
dim(marker.table)
marker.table <- t(marker.table)
marker.df <- as.data.frame(marker.table)
marker.df$cell <- rownames(marker.df)
marker.df <- marker.df %>% gather("marker", "count", 1:l.markers);

k <- 1;
for (k in 1:5){
  z <- marker.df %>% dplyr::group_by(cluster, marker) %>% 
    dplyr::summarise(pos=length(which(count>=k)), total=n()) %>% dplyr::mutate(percent=pos/total*100);
  z$cluster <- as.factor(z$cluster);
  print(z)
  
  p <- ggplot(z, aes(x=cluster, y=percent)) + geom_bar(stat="identity", position="dodge", width=0.9, aes(fill=marker));
  p <- p + theme_bw();
  p <- p + theme(text=element_text(colour="black"), axis.title=element_text(size=12), 
                 axis.text.x=element_text(colour="black", size=12), 
                 axis.text.y=element_text(colour="black", size=12));
  p <- p + labs(y="Percent", x="") + ggtitle(paste("Cells with", k, "or more reads"));
  p <- p + scale_fill_npg();
  p <- p + theme(legend.position="bottom");
  p <- p + theme(plot.title=element_text(hjust=0.5));
  print(p);
  z <- marker.df %>% dplyr::mutate(threshold=ifelse(count>=k, 1, 0)) %>% 
    dplyr::select(-count) %>% tidyr::spread(marker, threshold);
  z$cluster <- NULL;
  z$cell <- NULL;
  u <- upset(z)
  print(u)
}
dev.off()