library(Seurat)
library(scater)
library(scran)
library(dplyr)
library(ggpubr)
library(UpSetR)
library(ggsci)
library(tidyr)

# Project name: Don't use special symbols for project name: /?^%$ --> Is used for file names
sample <- "221940"
path.to.cellranger.outs <- paste("/media/data2/Daniel/cellranger_output/", sample, "/outs/", sep = "/")
project <- sample
seq <- c(seq(0.01,0.1,0.01), seq(0.12, 0.4, 0.02))

f <- paste0(project, ".umi_count_", 100, ".markers.supervised.txt")
write(paste("Project:", project, sep = " "), file = f)

# Read in sparse matrix with cellranger output
# data.dir = path to cellranger output matrix.mtx, features.tsv, barcodes.tsv
data.dir = paste(path.to.cellranger.outs, "/raw_feature_bc_matrix/", sep = "/")
data <- Read10X(data.dir = data.dir, gene.column = 2, unique.features = TRUE)

data
write(paste("Data:", data.dir, sep = " "), file = f, append = TRUE)
write(paste("Before:", dim(data), c("genes", "cells"), sep = " "), file = f, append = TRUE)
write(paste("Filter Umi Count:", filter.umi.count, "| filter features:", filter.feature.expr.cells, sep = " ", ), file = f, append = TRUE)


filter.umi.count <- 100
filter.feature.expr.cells <- ncol(data)*0.005
keep.features <- Matrix::rowSums(data > 0) >= filter.feature.expr.cells
keep.cells <- Matrix::colSums(data) >= filter.umi.count

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
libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "higher", log = T)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads = 3, type = "lower", log = T)
mito.drop <- isOutlier(sce$pct_counts_MT, nmads = 3, type = "higher")

sce <- sce[, !(libsize.drop | feature.drop | mito.drop)]
dim(sce)
write(paste("Filter after libsize, feature count, mt expr:", dim(sce), c("genes", "cells"), "left.", sep = " "), file = f, append = TRUE)


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
