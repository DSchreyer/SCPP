# This is an example analysis of the sample 221934
# filter cells with less than 100 genes expressed --> do able with this sample


library(Matrix)
library(Seurat)
library(scater)
library(scran)
library(dplyr)
library(ggpubr)
library(UpSetR)
library(ggsci)
library(tidyr)
library(DropletUtils)

# Project name: Don't use special symbols for project name: /?^%$ --> Is used for file names
project = "221934_BL6J_Ex45"

# Read in sparse matrix with cellranger output
# data.dir = path to cellranger output matrix.mtx, features.tsv, barcodes.tsv
data.dir = "/media/data2/Daniel/cellranger_output/221934/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data.dir, gene.column = 2, unique.features = TRUE)

keep.features <- Matrix::rowSums(data > 0) >= 20
keep.cells <- Matrix::colSums(data > 0) >= 100

data <- data[keep.features, keep.cells]
dim(data)

# Create SingleCellExperiment object with sparse matrix
sce <- SingleCellExperiment(assays = list(counts = data))

# Quality Control with scater
is.mito <- grepl("^mt-", rownames(sce))
# quality control with caluclateQCMetrics
sce <- calculateQCMetrics(sce, exprs_values = "counts",
                          feature_controls = list(MT = is.mito))

# Remove cells with fewer than 50 expressed genes
# remove genes with fewer than 20 expressed cells
# keep.cells <- sce$total_features_by_counts >= 50
# keep.features <- nexprs(sce, byrow = TRUE) >= 20
# sce <- sce[keep.features, keep.cells]

# Drop below and above 3 median absolute deviations: libsize, mitochondrial gene count, expressed genes
libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "higher", log = T)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads = 3, type = "lower", log = T)
mito.drop <- isOutlier(sce$pct_counts_MT, nmads = 3, type = "higher")

sce <- sce[, !(libsize.drop | feature.drop | mito.drop)]
dim(sce)

# normalize with scater function normalize using library sizes as size factors
sce <- normalize(sce)
pdf(paste(project, ".pdf", sep = ""))
plotHighestExprs(sce)
plotExpression(sce, features = c("Rho", "Opn1sw", "Opn1mw", "Arr3", "Nrl", "Gnat2"))

# Load into Seurat objet
Bsn <- as.Seurat(sce, data = "counts", assay = "RNA", counts = "counts", project = project)

# normalize the data
Bsn <- NormalizeData(Bsn, normalization.method = "LogNormalize", scale.factor = 1000)

# Find variable Features, top 2000
Bsn <- FindVariableFeatures(Bsn, selection.method = "vst", nfeatures = 2000)

# Scale Data for PCA
Bsn <- ScaleData(Bsn)

# Run PCA
Bsn <- RunPCA(Bsn)

# Identify PCs to take
ElbowPlot(Bsn)
dev.off()

pdf(paste(project, ".2.pdf", sep = ""))
# Cluster cells to 2,3,4,5 different clusters
# 3 cell clusters optimal
# check ElbowPlot for high variance pc
res = 0.4
dims = 1:8
Bsn <- FindNeighbors(Bsn, dims = dims)
Bsn <- FindClusters(Bsn, resolution = res)
markers <- FindAllMarkers(Bsn, only.pos = T, logfc.threshold = 0.25)

write.table(markers, file = paste0(project, "all.markers.txt"))

# table(Bsn$RNA_snn_res.0.4)

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
  dplyr::summarise(pos=length(which(count>=k)), total=n()) %>%
  dplyr::mutate(percent=pos/total*100);
  z$cluster <- as.factor(z$cluster);
  print(z)

  p <- ggplot(z, aes(x=cluster, y=percent)) +
  geom_bar(stat="identity", position="dodge", width=0.9, aes(fill=marker));
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
