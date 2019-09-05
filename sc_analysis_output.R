library(scater)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(Matrix)
library(tibble)
library(dplyr)
# SingleCellExperiment | Scater
data.dir <- "/home/daniel/master_thesis/bassoon_data/Cellranger output/Aggr/raw_feature_bc_matrix/"
output.dir <- "/home/daniel/master_thesis/bassoon_data/Output/Tables_Graphs"

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
  sample.names = names(count.tables),
  generate.info = TRUE
){
  if (is.null(sample.names)){
    sample.names <- paste0("sample", as.character(seq(1, length(count.tables))))
  }
  i <- 1
  sce.list <- list()
  for (count.table in count.tables){
    sample <- sample.names[i]
    # Create SingleCellExperiment object with sparse matrix
    sce <- SingleCellExperiment(assays = list(counts = count.table))
    
    # Quality Control with scater
    is.mito <- grepl("^mt-", rownames(sce))
    # quality control with caluclateQCMetrics
    sce <- calculateQCMetrics(sce, exprs_values = "counts",
                              feature_controls = list(MT = is.mito))
    n.cells.pre <- ncol(sce)
    # Drop cells above [n] median absolute deviations: libsize, mitochondrial gene count, expressed genes
    libsize.drop <- isOutlier(sce$total_counts, nmads = MAD, type = "higher", log = T)
    feature.drop <- isOutlier(sce$total_features_by_counts, nmads = MAD, type = "higher", log = T)
    mito.drop <- isOutlier(sce$pct_counts_MT, nmads = MAD, type = "higher")
    sce <- sce[, !(libsize.drop | feature.drop | mito.drop)]
    n.cells.post <- ncol(sce)
    if (generate.info){
      if (i == 1){
        libsizes <- c()
        featuredrops <- c()
        mitodrops <- c()
        before.qc <- c()
        after.qc <- c()
        samples <- c()
        median.umi <- c()
        median.gene <- c()
        mean.umi <- c()
        mean.gene <- c()
      }
      libsizes <- c(libsizes, sum(libsize.drop))
      featuredrops <- c(featuredrops, sum(feature.drop))
      mitodrops <- c(mitodrops, sum(mito.drop))
      samples <- c(samples, sample)
      before.qc <- c(before.qc, n.cells.pre)
      after.qc <- c(after.qc, n.cells.post)
      median.umi <- c(median.umi, median(sce$total_counts))
      median.gene <- c(median.gene, median(sce$total_features_by_counts))
      mean.umi <- c(mean.umi, mean(sce$total_counts))
      mean.gene <- c(mean.gene, mean(sce$total_features_by_counts))
    }
    
    sce.list[[sample]] <- sce
    print(paste(sample, "Done"))
    i <- i + 1
  }
  qc.info <<- cbind(before.qc, after.qc, libsizes, featuredrops, mitodrops, median.umi, median.gene, mean.umi, mean.gene)
  qc.info <<- as.data.frame(qc.info)
  qc.info$samples <<- samples
  
  return(sce.list)
}
sce.list <- qcControl(count.tables = count.tables, MAD = 5, sample.names = names(count.tables))

seurat.objects$genotype <- NA
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221931"] <- "C57BL/6J_WT"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221932"] <- "C57BL/6J_WT"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221933"] <- "C57BL/6J_Bsn_mut"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221934"] <- "C57BL/6J_Bsn_mut"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221935"] <- "C57BL/6J_Bsn_mut"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221936"] <- "C57BL/6NJ_WT"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221937"] <- "C57BL/6NJ_WT"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221938"] <- "C57BL/6NJ_Bsn_KO"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221939"] <- "C57BL/6NJ_Bsn_KO"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221940"] <- "C57BL/6NJ_Bsn_KO"
## Generell sample overview with qc.info
qc.info$samples <- as.character(seq(221931, 221940))
qc.info$genotype <- c("C57BL/6J_WT", "C57BL/6J_WT", "C57BL/6J_Bsn_mut", "C57BL/6J_Bsn_mut", "C57BL/6J_Bsn_mut",
                                 "C57BL/6NJ_WT", "C57BL/6NJ_WT", "C57BL/6NJ_Bsn_KO",  "C57BL/6NJ_Bsn_KO",  "C57BL/6NJ_Bsn_KO")
qc.info$genotype <- factor(qc.info$genotype, levels = c("C57BL/6J_WT",  "C57BL/6J_Bsn_mut",  "C57BL/6NJ_WT",  "C57BL/6NJ_Bsn_KO"))


# Plot: Samples on x-axis, total cells on y - axis
# Plot: Samples on x-axis, total cells before + after qc

pdf("sample_graph_info.pdf")
# colors <- c(rep("#A4BFEB", 2), rep("#8CABBE", 3), rep("#BBA0B2", 2), rep("#9D858D", 3))
colors <- c("#A4BFEB", "#8CABBE", "#BBA0B2", "#9D858D")

ggplot(data = qc.info, aes(x = samples, y = before.qc, fill = genotype)) + 
  geom_bar(stat = "identity", colour = "black") + scale_fill_manual(values = colors) + 
  theme_minimal() + 
  scale_y_continuous(limits = c(0,30000)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.1, vjust = 1.5))
  
ggplot(data = qc.info, aes(x = samples, y = after.qc, fill = genotype)) + 
  geom_bar(stat = "identity", colour = "black") + scale_fill_manual(values = colors) + 
  theme_minimal() +
  scale_y_continuous(limits = c(0,30000)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.1, vjust = 1.5))

# Plot: Samples on x-axis, median.umi and median.gene on y-axis
ggplot(data = qc.info, aes(x = samples, y = median.umi, fill = genotype)) + 
  geom_bar(stat = "identity", colour = "black") + scale_fill_manual(values = colors) + 
  theme_minimal() + 
  scale_y_continuous(limits = c(0,200)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.1, vjust = 1.5))

ggplot(data = qc.info, aes(x = samples, y = median.gene, fill = genotype)) + 
  geom_bar(stat = "identity", , colour = "black") + scale_fill_manual(values = colors) + 
  theme_minimal() + 
  scale_y_continuous(limits = c(0,150)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.1, vjust = 1.5))




barplot(seq(1,10), col=colors)



## Seurat
input.dir <- "/home/daniel/master_thesis/bassoon_data/Output/post_ct_ident/"
output.dir <- "/home/daniel/master_thesis/bassoon_data/Output/downstream_analyses/"
setwd(output.dir)

samples <- paste0("Bsn.", seq(221931, 221940))


LoadSeuratFiles <- function(
  dir,
  log.counts = TRUE,
  samples,
  counts = "logcounts.mtx",
  feature = "features.txt",
  barcode = "barcodes.txt",
  celltype = "celltype.txt"
){
  seurat.objects <- list()
  ls <- list.files(dir)
  setwd(dir)
  first <- TRUE
  for (sample in samples){
    ind <- grep(pattern = sample, ls)
    count.file <- grep(counts, ls[ind], value = TRUE)
    feature.file <- grep(feature, ls[ind], value = TRUE)
    barcode.file <- grep(barcode, ls[ind], value = TRUE)
    celltype.file <- grep(celltype, ls[ind], value = TRUE)
    count.table <- readMM(count.file)
    features <- read.csv(feature.file, stringsAsFactors = F, header = F)
    barcodes <- read.csv(barcode.file, stringsAsFactors = F, header = F)
    celltypes <- read.csv(celltype.file, stringsAsFactors = F, header = F)
    rownames(count.table) <- features$V1
    colnames(count.table) <- barcodes$V1
    seurat.object <- CreateSeuratObject(project = sample, counts = count.table)
    seurat.object$CellType <- celltypes$V1
    Idents(seurat.object) <- seurat.object$CellType
    seurat.object <- RenameCells(seurat.object, add.cell.id = seurat.object$orig.ident)
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
seurat.objects <- LoadSeuratFiles(dir = input.dir, samples = samples)

# 8,9,10: C57BL/6NJ Bsn_KO
seurat.objects$genotype <- NA
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221931"] <- "C57BL/6J_WT"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221932"] <- "C57BL/6J_WT"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221933"] <- "C57BL/6J_Bsn_mut"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221934"] <- "C57BL/6J_Bsn_mut"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221935"] <- "C57BL/6J_Bsn_mut"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221936"] <- "C57BL/6NJ_WT"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221937"] <- "C57BL/6NJ_WT"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221938"] <- "C57BL/6NJ_Bsn_KO"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221939"] <- "C57BL/6NJ_Bsn_KO"
seurat.objects$genotype[seurat.objects$orig.ident == "Bsn.221940"] <- "C57BL/6NJ_Bsn_KO"

## Table with Number of cells, genes, median.umi, median.gene, rods, cones, unknown, proportions
seurat.objects$genotype <- factor(seurat.objects$genotype, levels = c("C57BL/6J_WT",  "C57BL/6J_Bsn_mut",  "C57BL/6NJ_WT",  "C57BL/6NJ_Bsn_KO"))
clean <- subset(seurat.objects, subset = CellType != "Unknown")


n.cells.genotype <- as.data.frame(table(seurat.objects$genotype))
colnames(n.cells.genotype) <- c("Genotype", "n.cells")

ggplot(data = n.cells.genotype, aes(x = Genotype, y = n.cells, fill = Genotype)) + 
  geom_bar(stat = "identity", colour = "black") + scale_fill_manual(values = colors) + 
  theme_minimal() + theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0, 60000)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1))

n.cells.genotype.wo.unk <- as.data.frame(table(clean$genotype))
colnames(n.cells.genotype.wo.unk) <- c("Genotype", "n.cells")


ggplot(data = n.cells.genotype.wo.unk, aes(x = Genotype, y = n.cells, fill = Genotype)) + 
  geom_bar(stat = "identity", colour = "black") + scale_fill_manual(values = colors) + 
  theme_minimal() + theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0, 60000)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1))


# Celltype plots
colors3 <- c("#FFEDDF", "#C5D86D", "#666469")


ct.table <- as.data.frame(t(table(seurat.objects$CellType, seurat.objects$genotype)))
colnames(ct.table) <- c("genotype", "celltype", "n.cells")
ct.table$celltype <- factor(ct.table$celltype, levels = c("Unknown", "Cone", "Rod"))

ggplot(ct.table, aes(x = genotype, y = n.cells, fill = celltype)) + 
  geom_bar(stat = "identity", colour = "black") + 
  scale_fill_manual(values = colors3) + theme_minimal()

ct.table.wo.unk <- dplyr::filter(ct.table, celltype != "Unknown")
colors2 <- colors3[c(2,3)]
ggplot(ct.table.wo.unk, aes(x = genotype, y = n.cells, fill = celltype)) + 
  geom_bar(stat = "identity", colour = "black") + 
  scale_fill_manual(values = colors2)  + theme_minimal()

ct.sample <- as.data.frame(t(table(seurat.objects$CellType, seurat.objects$orig.ident)))
colnames(ct.sample) <- c("sample", "celltype", "n.cells")
ct.sample$sample <- as.character(seq(221931, 221940))
ct.sample$genotype <-  c(rep("C57BL/6J_WT", 2), rep("C57BL/6J_Bsn_mut",3), 
                         rep("C57BL/6NJ_WT", 2), rep("C57BL/6NJ_Bsn_KO", 3))

ct.sample$celltype <- factor(ct.sample$celltype, levels = c("Unknown", "Cone", "Rod"))
ggplot(ct.sample,aes( x = sample, y = n.cells, fill = celltype)) + 
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = colors3) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1))
  

ct.sample.wo.unk <- dplyr::filter(ct.sample, celltype != "Unknown")
ggplot(ct.sample.wo.unk, aes(x = sample, y = n.cells, fill = celltype)) + 
  geom_bar(stat = "identity", colour = "black") + 
  scale_fill_manual(values = colors2) + theme_minimal()
dev.off()



# Create tables with information
file = "Bsn_info_table.csv"
write.table(table(seurat.objects$CellType, seurat.objects$orig.ident), file = file, sep = " ", col.names = T)
write.table(prop.table(table(seurat.objects$CellType, seurat.objects$orig.ident), margin = 2), file = file, sep = " ", append = T)
write.table(table(seurat.objects$CellType), file = file, sep = " ", append = T)
write.table(prop.table(table(seurat.objects$CellType)), file = file, sep = " ", append = T)
write.table(table(seurat.objects$CellType, seurat.objects$genotype), file = file, sep = " ", append = T)
write.table(table(seurat.objects$genotype), file = file, sep = " ", append = T)
write.table(prop.table(table(seurat.objects$CellType, seurat.objects$genotype), margin = 2), file = file, append = T ,sep = " ")
write.table(table(seurat.objects$orig.ident), file = file, append = T, sep = " ", col.names = T)

# Without Unknown cells
clean <- subset(seurat.objects, subset = CellType != "Unknown")
write.table(table(clean$CellType, clean$orig.ident), file = file, append = T, sep = " ")
write.table(prop.table(table(clean$CellType, clean$orig.ident), margin = 2), file = file, append = T, sep = " ")
write.table(table(clean$CellType), file = file, append = T, sep = " ")
write.table(prop.table(table(clean$CellType)), file = file, append = T, sep = " ")
write.table(table(clean$CellType, clean$genotype), file = file, sep = " ", append = T)
write.table(prop.table(table(clean$CellType, clean$genotype), margin = 2), file = file, append = T, sep = " ")
write.table(table(clean$orig.ident), file = file, append = T, sep = " ", col.names = T)

test <- clean@meta.data %>% group_by(orig.ident) %>% summarise(median = median(nFeature_RNA))



# sample 2: TSNE with Identified Populations
sample <- subset(seurat.objects, subset = orig.ident == "Bsn.221932")
sample <- FindVariableFeatures(sample)
sample <- ScaleData(sample)
sample <- RunPCA(sample)
sample <- RunTSNE(sample, dims = 1:6)
sample <- RunUMAP(sample, dims = 1:10)
DimPlot(sample, reduction = "tsne")
DimPlot(sample, reduction = "umap")

# generate different plots
# raw counts read in
counts <- readMM("/home/daniel/master_thesis/bassoon_data/Output/post_qc/Bsn.221932.raw.counts.mtx")
features <- read.table("/home/daniel/master_thesis/bassoon_data/Output/post_qc/Bsn.221932.features.txt")
barcodes <- read.table("/home/daniel/master_thesis/bassoon_data/Output/post_qc/Bsn.221932.barcodes.txt")
rownames(counts) <- features$V1
colnames(counts) <- barcodes$V1

sample.sce <- SingleCellExperiment(list(counts = counts))
is.mito <- grepl("^mt-", rownames(sample.sce))
sample.sce <- calculateQCMetrics(sample.sce, feature_controls = list(Mt=is.mito))

par(mfrow=c(1,1))
hist(sample.sce$total_counts, xlab="Library sizes", main="", 
     breaks=50, col="grey80", ylab="Number of cells")
hist(sample.sce$total_features_by_counts, xlab="Number of expressed genes", main="", 
     breaks=50, col="grey80", ylab="Number of cells")
hist(sample.sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=50, main="", col="grey80",xlim = c(0,100))

sample.sce.norm <- normalize(sample.sce)
genes <- c("Rho", "Opn1sw", "Opn1mw", "Gnat1", "Arr3", "Gnat2")
plotExpression(sample.sce.norm, features = genes)

plotHighestExprs(sample.sce, exprs_values = "counts")
