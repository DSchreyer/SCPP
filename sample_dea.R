# Identify Marker gene expression
library(Matrix)
library(Seurat)
library(tibble)
library(plyr)
library(dplyr)
library(xtable)
library(enrichR)
library(pheatmap)
library(RColorBrewer)
library(pheatmap)

input.dir <- "/home/daniel/master_thesis/bassoon_data/Output/new_ct_ident/"
output.dir <- "/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/"
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
setwd(output.dir)

samples <- paste0("Bsn.", seq(221931, 221940))
samples <- samples[-5]


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


# ### TRANSFER CELL LABELS TO OTHER OBJECTS
ref <- subset(seurat.objects, subset = orig.ident == "Bsn.221932")
ref <- FindVariableFeatures(ref, nfeatures = 1000)
ref <- ScaleData(ref)
ref <- RunPCA(ref)
ref <- FindNeighbors(ref, dims = 1:7)
# 0.05 works
ref <- FindClusters(ref, resolution = 0.1)
m <- FindAllMarkers(ref, only.pos = T)
# 2 is cones 0 is rods
ref <- RunTSNE(ref, dims = 1:7)

ref <- RunUMAP(ref, dims = 1:7)
pdf("/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/221932_ref_umap.pdf")
DimPlot(ref, reduction = "umap")
dev.off()

pdf("/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/221932_ref_tsne.pdf")
DimPlot(ref, reduction = "tsne")
dev.off()

seurat.object.list <- SplitObject(seurat.objects, split.by = "orig.ident")


## Transfer Anchors from ref data set
new.so.list <- list()
for (object in seurat.object.list){
  name <- unique(object$orig.ident)
  if (ncol(object) < 500){
    next
  }
  anchors <- FindTransferAnchors(ref, object)
  predictions <- TransferData(anchorset = anchors, refdata = ref$seurat_clusters, dims = 1:7)
  object <- AddMetaData(object, metadata = predictions)
  Idents(object) <- object$predicted.id
  new.so.list[[name]] <- object
}

save(new.so.list, file = "/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/post_transfer_so.Robj")
load("/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/post_transfer_so.Robj")


markers <- FindAllMarkers(new.so.list$Bsn.221932, only.pos = T)

so.list <- new.so.list
so.list$Bsn.221931 <- NULL

so <- merge(x = new.so.list$Bsn.221931, so.list)





# Combine Genotypes together
# 1,2: C57BL/6J_WT
# 3,4,5: C57BL/6J Bsn_mut
# 6,7: C57BL/6NJ WT
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

seurat.objects$ct.gt <- paste(seurat.objects$genotype, seurat.objects$CellType, sep = ".")
Idents(seurat.objects) <- seurat.objects$ct.gt


# Compare Rods and Cones to each other --> Different vulnerability
deg <- list()
deg[["wt1.rod.cone"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6J_WT.rod",
                                   ident.2 = "C57BL/6J_WT.cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)
deg[["wt2.rod.cone"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6NJ_WT.rod",
                                   ident.2 = "C57BL/6NJ_WT.cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)

deg[["mut.rod.cone"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6J_Bsn_mut.rod", 
                               ident.2 = "C57BL/6J_Bsn_mut.cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)
deg[["ko.rod.cone"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6NJ_Bsn_KO.rod", 
                                     ident.2 = "C57BL/6NJ_Bsn_KO.cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)
deg[["wt1.wt2.rod"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6J_WT.rod", 
                                ident.2 = "C57BL/6NJ_WT.rod", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)
deg[["wt1.wt2.cone"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6J_WT.cone", 
                                     ident.2 = "C57BL/6NJ_WT.cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)


# Compare rod vs. rod, cone vs. cone
# mutant strain
deg[["wt.mut.rod"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6J_WT.rod",
                                    ident.2 = "C57BL/6J_Bsn_mut.rod", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)
deg[["wt.mut.cone"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6J_WT.cone",
                                    ident.2 = "C57BL/6J_Bsn_mut.cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)

# KO strain
deg[["wt.ko.rod"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6NJ_WT.rod",
                                   ident.2 = "C57BL/6NJ_Bsn_KO.rod", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)
deg[["wt.ko.cone"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6NJ_WT.cone",
                                   ident.2 = "C57BL/6NJ_Bsn_KO.cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)


# Compare Mutant to KO
deg[["mut.ko.cone"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6J_Bsn_mut.cone", 
                               ident.2 = "C57BL/6NJ_Bsn_KO.cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)
deg[["mut.ko.rod"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6J_Bsn_mut.rod", 
                               ident.2 = "C57BL/6NJ_Bsn_KO.rod", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)

save(deg, file = "/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/deg_list.Robj")
load("/home/daniel/master_thesis/bassoon_data/Output/New_analysis_output/deg_list.Robj")


# Enrichment analysis
dbs <- listEnrichrDbs()
websiteLive <- ifelse(is.null(dbs), FALSE, TRUE)
if (websiteLive) head(dbs)

databases <- c("GO_Biological_Process_2018", "KEGG_2019_Mouse", "MGI_Mammalian_Phenotype_Level_4_2019", "Mouse_Gene_Atlas", "GO_Molecular_Function_2018")
databases <- c("KEGG_2019_Mouse")

enriched <- list()
i <- 1

deg.sep <- list()
for (table in deg){
  name <- names(deg)[i]
  table$expr <- ifelse(table$avg_logFC > 0, "up", "down")
  up <- dplyr::filter(table, expr == "up")
  down <- dplyr::filter(table, expr == "down")
  deg.sep[[paste(name, "up", sep = ".")]] <- up
  deg.sep[[paste(name, "down", sep = ".")]] <- down
  enriched.all <- enrichr(table$gene, databases)
  enriched[[paste(name, "all", sep = ".")]] <- enriched.all
  enriched.up <- enrichr(up$gene, databases)
  enriched[[paste(name, "up", sep = ".")]] <- enriched.up
  enriched.down <- enrichr(down$gene, databases)
  enriched[[paste(name, "down", sep = ".")]] <- enriched.down
  i <- i + 1
}

# Generate Heatmap with DE genes of all. LogFC as color -- only cone genes
# pdf("deg.heatmap.pdf")
gene.list <- lapply(deg, "[", , c(1,3))
cone.deg.list <- lapply(gene.list, "[", , 1)
cone.deg.list <- unique(unlist(cone.deg.list[grepl(names(cone.deg.list), pattern = "cone")], use.names = F))

deg.table <- data.frame(matrix(nrow = length(cone.deg.list)))
deg.table$gene <- cone.deg.list
new.list <- list()
new.list$gene <- deg.table
new.list <- append(new.list, gene.list)
full.table <- join_all(new.list, by = "gene", type = "left")
rownames(full.table) <- full.table$gene
colnames(full.table) <- paste0("V", seq(1, length(full.table)))
full.table <- select(full.table, 3:length(full.table))
colnames(full.table) <- names(gene.list)

colors <- colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))
full.table[is.na(full.table)] <- 0
cone.genes.deg <- pheatmap(full.table, cluster_rows = F, cluster_cols = FALSE,
                           breaks = seq(-1.5,1.5, by = 0.05), 
                           color = colors(60), fontsize_row = 5)
cone.genes.deg <- pheatmap(full.table, cluster_rows = T, cluster_cols = FALSE,
         breaks = seq(-1.5,1.5, by = 0.05), 
         color = colors(60), fontsize_row = 5)
write.table(x = full.table, file = "deg.cones.heatmap.logfc.csv", sep = ",")


# all genes
gene.list <- lapply(deg, "[", , c(1,3))
all.deg.list <- lapply(gene.list, "[", , 1)
all.deg.list <- unique(unlist(all.deg.list))

deg.table <- data.frame(matrix(nrow = length(all.deg.list)))
deg.table$gene <- all.deg.list
new.list <- list()
new.list$gene <- deg.table
new.list <- append(new.list, gene.list)
full.table <- join_all(new.list, by = "gene", type = "left")
rownames(full.table) <- full.table$gene
colnames(full.table) <- paste0("V", seq(1, length(full.table)))
full.table <- select(full.table, 3:length(full.table))
colnames(full.table) <- names(gene.list)

colors <- colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))
full.table[is.na(full.table)] <- 0
pheatmap(full.table, cluster_cols = F, cluster_rows = F,
         breaks = seq(-1.5,1.5, by = 0.05), 
         color = colors(60))
pheatmap(full.table, cluster_cols = F, cluster_rows = T,
         breaks = seq(-1.5,1.5, by = 0.05), 
         color = colors(60))
write.table(x = full.table, file = "deg.all.heatmap.logfc.csv", sep = ",")


# Generate Heatmap with DE genes -- separated up and downregulated
gene.list <- lapply(deg.sep, "[", , c(1,3))
cone.deg.list <- lapply(gene.list, "[", , 1)
cone.deg.list <- unique(unlist(cone.deg.list[grepl(names(cone.deg.list), pattern = "cone")], use.names = F))

deg.table <- data.frame(matrix(nrow = length(cone.deg.list)))
deg.table$gene <- cone.deg.list
new.list <- list()
new.list$gene <- deg.table
new.list <- append(new.list, gene.list)
full.table <- join_all(new.list, by = "gene", type = "left")
rownames(full.table) <- full.table$gene
colnames(full.table) <- paste0("V", seq(1, length(full.table)))
full.table <- select(full.table, 3:length(full.table))
colnames(full.table) <- names(gene.list)
heat.table <- ifelse(is.na(full.table), 0, 1)
pheatmap(heat.table, cluster_rows = T, cluster_cols = F,
         fontsize_row = 6)
write.table(x = heat.table, file = "deg.up.down.heatmap.csv", sep = ",")
dev.off()

pdf("pathways_heatmap.pdf")
## Pathway analysis Heatmap --> Only upregulated pathways
up.down <- grepl(names(enriched), pattern = "up|down")
new <- enriched[up.down]

kegg.list <- lapply(new, function(x){
  x <- as.data.frame(x$KEGG_2019_Mouse)
  x <- dplyr::filter(x, Adjusted.P.value < 0.05)
  return(x)
})

path.pval <- lapply(kegg.list, "[", , c(1,4))
pathways <- unique(unlist(lapply(kegg.list, "[", , 1)))
path.pval.test <- Filter(nrow, path.pval)

path.table <- data.frame(matrix(nrow = length(pathways)))
path.table$Term <- pathways
new.list <- list()
new.list$Term <- path.table
new.list <- append(new.list, path.pval.test)
full.table <- join_all(new.list, by = "Term", type = "left")
rownames(full.table) <- full.table$Term
colnames(full.table) <- paste0("V", seq(1, length(full.table)))
full.table <- select(full.table, 3:length(full.table))
colnames(full.table) <- names(path.pval.test)
heat.table <- ifelse(is.na(full.table), 0, 1)
pheatmap(heat.table, cluster_rows = T, cluster_cols = F)
write.table(x = heat.table, file = "pathways.up.down.heatmap.csv", sep = ",")

#all genes in one 
all <- grepl(names(enriched), pattern = "all")
new <- enriched[all]

kegg.list <- lapply(new, function(x){
  x <- as.data.frame(x$KEGG_2019_Mouse)
  x <- dplyr::filter(x, Adjusted.P.value < 0.05)
  return(x)
})


path.pval <- lapply(kegg.list, "[", , c(1,4))
pathways <- unique(unlist(lapply(kegg.list, "[", , 1)))
path.pval.test <- Filter(nrow, path.pval)

path.table <- data.frame(matrix(nrow = length(pathways)))
path.table$Term <- pathways
new.list <- list()
new.list$Term <- path.table
new.list <- append(new.list, path.pval.test)
full.table <- join_all(new.list, by = "Term", type = "left")
rownames(full.table) <- full.table$Term
colnames(full.table) <- paste0("V", seq(1, length(full.table)))
full.table <- select(full.table, 3:length(full.table))
colnames(full.table) <- names(path.pval.test)
heat.table <- ifelse(is.na(full.table), 0, 1)
pheatmap(heat.table, cluster_rows = FALSE, cluster_cols = FALSE)
write.table(x = heat.table, file = "pathways.all.heatmap.csv", sep = ",")
dev.off()






# Principal Component Analysis of Cones
good <- subset(subset, subset  = orig.ident == "Bsn.221932")
good <- subset(good, subset = CellType == "Cone")
good <- FindVariableFeatures(good, nfeatures = 500)
good <- ScaleData(good)
good <- RunPCA(good)
good <- RunTSNE(good)
DimPlot(good, reduction = "pca")
DimPlot(good, reduction = "tsne")
good <- FindNeighbors(good)
good <- FindClusters(good)
FindAllMarkers(good)


bad <- subset(subset, subset  = orig.ident == "Bsn.221939")
bad <- subset(bad, subset = CellType == "Cone")
bad <- FindVariableFeatures(bad, nfeatures = 500)
bad <- ScaleData(bad)
bad <- RunPCA(bad)
bad <- RunTSNE(bad)
DimPlot(bad, reduction = "pca")
DimPlot(bad, reduction = "tsne")
bad <- FindNeighbors(bad)
bad <- FindClusters(bad, resolution = 0.6)
FindAllMarkers(bad)
table(Idents(bad))



## DEG TABLES WITHOUT COMMON GENES
deg.filtered <- list()
# Upregulated genes in mut cones to wt1
table <- deg.sep$wt.mut.cone.down
genes <- table$gene
# Filter out genes, which are downregulated in wt cones
genes <- setdiff(genes, deg.sep$cone.wt.down$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["mut.cone.up"]] <- subset

# Upregulated genes in ko cones to wt2
table <- deg.sep$wt.ko.cone.down
genes <- table$gene
# Filter out genes, which are downregulated in wt cones 2
genes <- setdiff(genes, deg.sep$cone.wt2.down$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["ko.cone.up"]] <- subset

# Upregulated genes in mut cones to ko cones
table <- deg.sep$mut.ko.cone.up
genes <- table$gene
# Filter out genes, which are upregulated in wt cones
genes <- setdiff(genes, deg.sep$cone.wt.up$gene)
# Filter out genes, which are downregulated in ko vs wt
genes <- setdiff(genes, deg.sep$wt.ko.cone.up$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["mut.ko.cone.up"]] <- subset

# Upregulated genes in ko cones to mut cones
table <- deg.sep$mut.ko.cone.down
genes <- table$gene
# Filter out genes, which are upregulated in wt cones
genes <- setdiff(genes, deg.sep$cone.wt.down$gene)
# Filter out genes, which are downregulated in ko vs wt
genes <- setdiff(genes, deg.sep$wt.ko.cone.up$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["mut.ko.cone.down"]] <- subset



dbs <- "KEGG_2019_Mouse"
enr <- enrichr(deg$rod.cone.wt.1$gene, dbs)

enrichr(setdiff(deg.sep$cone.wt.down$gene, deg.sep$wt.mut.cone.down$gene), dbs)

enriched$rod.cone.wt.1.down$KEGG_2019_Mouse
enriched$rod.cone.mut.down$KEGG_2019_Mouse
enriched$wt.mut.cone.up$KEGG_2019_Mouse
enriched$wt.mut.cone.down$KEGG_2019_Mouse

write("", "enriched.csv")
i <- 1
length(names(enriched))
for (comp in enriched){
  name <- names(enriched)[i]
  j <- 1
  print(name)
  for (table in comp){
    comp.name <- names(comp)[j]
    table <- table %>% filter(Adjusted.P.value < 0.05)
    write(name, "enriched.csv", append = T)
    write(comp.name, "enriched.csv", append = T)
    write.table(table, "enriched.csv", append = T, sep = "\t", row.names = F)
    write("", "enriched.csv", append = T)
    j <- j + 1
  }  
  i <- i + 1
}



### Generate Differentially Expressed Gene Tables
i <- 1
for(table in deg){
  if (i == 1){
    names <- names(deg)
    file <- "deg.enriched.tex"
    write(x = "", file = file)
    write(x = "", file = "deg.table.csv")
  }
  name <- names[i]
  table <- table %>% select(gene, p_val_adj, avg_logFC, pct.1, pct.2) %>% filter(p_val_adj < 0.05) %>% 
    arrange(p_val_adj)
  gene <- select(table, gene)
  p.adj <- table$p_val_adj
  p.adj <- formatC(p.adj, format = "e", digits = 3)
  table <- table %>% select(avg_logFC, pct.1, pct.2)
  table$expr <- ifelse(table$avg_logFC > 0, "up", "down")
  table <- cbind(gene, p.adj, table) %>% arrange(desc(expr, p_val_adj))
  write(x = name, file = "deg.table.csv", append = T)
  write.table(x = table, file = "deg.table.csv", append = T, sep = "\t", quote = F, row.names = F)
  write(x = "", file = "deg.table.csv", append = T)
  write(print(xtable(table, digits = 3), tabular.environment= "longtable", floating = FALSE), file = file, append = T)
  for (enrich.df in enriched[[i]]){
    enrich.df <- enrich.df %>% select(Term, Overlap, P.value, Adjusted.P.value, Genes) %>%
      filter(Adjusted.P.value < 0.1)
    write(print(xtable(enrich.df, digits = 3), tabular.environment = "longtable", floating = FALSE), file = file, append = T)
  }
  i <- i + 1
}




test <- subset(seurat.objects, subset = orig.ident == "Bsn.221936")
1
test <- FindVariableFeatures(test)
test <- ScaleData(test)
test <- RunPCA(test)
test <- RunTSNE(test)
DimPlot(test, reduction = "tsne")
test <- FindNeighbors(test)
test <- FindClusters(test, resolution = 0.5)
FindAllMarkers(test)

test <- enriched$GO_Biological_Process_2018

org <- org.Mm.eg.db
keytypes(org)
go.ids <- keys(org, keytype = "GO")
genes <- keys(org, keytype = "ALIAS")


AnnotationDbi::select(org, keys = markers$gene, columns = "GO", keytype = "ALIAS")



Bsn.all$sample <- Bsn.all@meta.data$orig.ident
Bsn.all <- NormalizeData(Bsn.all)
table(Idents(Bsn.all))
table(Bsn.all$sample)
prop.table(table(Bsn.all$CellType))
table(Idents(Bsn.all), Bsn.all$sample)
prop.table(table(Idents(Bsn.all), Bsn.all$sample), margin = 2)
subset(Bsn.all, idents = c("Cone", "Rod"))
average_expressions <- AverageExpression(Bsn, return.seurat = TRUE)
CellScatter(average_expressions, cell1 = "Cone", cell2 = "Rod")


library(KEGG.db)
KEGG <- as.list(KEGGPATHID2NAME)
Idents(Bsn.all) <- Bsn.all$id.celltype
markers <- FindMarkers(Bsn.all, ident.1 = "Bsn.221933.Cone", ident.2 = "Bsn.221940.Cone") %>% rownames_to_column("gene") %>%  filter(p_val_adj < 0.05)
avg.logfc <- ifelse(markers$avg_logFC > 0, "pos", "neg")
markers$regulation <- avg.logfc
kegg.ids <- AnnotationDbi::select(org, keys = markers$gene, columns = "PATH", keytype = "ALIAS")
kegg.ids$gene <- kegg.ids$ALIAS

marker.kegg <- inner_join(kegg.ids, markers, by = "gene")


path <- c()
for (id in marker.kegg$PATH){
  if (!is.na(id)){
    path <- c(path, KEGG[[as.character(id)]])
  }
  else
    path <- c(path, NA)
}
marker.kegg$PATHNAME <- path
table(marker.kegg$PATHNAME, marker.kegg$regulation)

all.sample.ct <- levels(Idents(Bsn.all))


all.markers <- c()
markers <- c()
write("Cone transcriptome comparison between all samples", file = "/home/daniel/master_thesis/bassoon_data/Cellranger output/cone_marker_table.txt")
for (ct.1 in all.sample.ct){
  for (ct.2 in all.sample.ct){
    if (ct.1 != ct.2){
      if (str_detect(ct.1, "Cone") && str_detect(ct.2, "Cone")){
        markers <- try(FindMarkers(Bsn.all, ident.1 = ct.1, ident.2 = ct.2) %>% rownames_to_column("gene") %>% filter(p_val_adj < 0.05) %>% select(gene, avg_logFC, p_val_adj))
        write(x = ct.1, file = "/home/daniel/master_thesis/bassoon_data/Cellranger output/cone_marker_table.txt", append = T)
        write(x = ct.2, file = "/home/daniel/master_thesis/bassoon_data/Cellranger output/cone_marker_table.txt", append = T)
        try(write.table(markers, file = "/home/daniel/master_thesis/bassoon_data/Cellranger output/cone_marker_table.txt", append = T, quote = F, row.names = F))
        all.markers <- c(all.markers, ct.1, ct.2, markers)
      }
    }
  }
}

m <- rownames(markers)

Seurat.objects$Bsn.221935

for (n in m){
  pathway <- grep(n, mm_GO)
  print(n)
  for (i in pathway){
    print(unlist(mm_GO[i])[1])
  }
}

library(gskb)

data(mm_GO)
data("mm_pathway")
mm_GO
mm_pathway
AcBsn.test <- FindVariableFeatures(Bsn.all)
Bsn.test <- ScaleData(Bsn.test)
Bsn.test <- RunPCA(Bsn.test)

features <- rownames(Bsn)
barcodes <- colnames(Bsn)