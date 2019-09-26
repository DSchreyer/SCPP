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
library(ggplot2)

input.dir <- "/home/daniel/master_thesis/bassoon_data/Output/final_ct_class///"
output.dir <- "/home/daniel/master_thesis/bassoon_data/Output/final_downstream_analysis/"
analysis.out <- "/home/daniel/master_thesis/bassoon_data/Output/final_downstream_analysis/"
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
seurat.objects$CellType <- Idents(seurat.objects)

seurat.objects$ct.gt <- paste(seurat.objects$genotype, seurat.objects$CellType, sep = ".")
Idents(seurat.objects) <- seurat.objects$ct.gt

# generate useful info
table_gen_ct <- table(seurat.objects$genotype, seurat.objects$CellType)
write.table(file = paste0(analysis.out, "/table_gen_ct.csv"), x = table_gen_ct, sep = ",")


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
deg[["mut.ko.rod"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6J_Bsn_mut.rod", 
                               ident.2 = "C57BL/6NJ_Bsn_KO.rod", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)
deg[["mut.ko.cone"]] <- FindMarkers(seurat.objects, ident.1 = "C57BL/6J_Bsn_mut.cone", 
                                    ident.2 = "C57BL/6NJ_Bsn_KO.cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% dplyr::filter(p_val_adj <= 0.05)


lapply(deg, function(x){write.table(as.data.frame(x), 
                                    file = paste0(analysis.out, "/deg.table.csv"),
                                    append = T, 
                                    sep = ",")})
save(deg, file = "/home/daniel/master_thesis/bassoon_data/Output/final_downstream_analysis/deg_list.Robj")




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
head(enriched)
i <- 1
names <- names(enriched)
file <- "/home/daniel/master_thesis/bassoon_data/Output/final_downstream_analysis/pathway_analysis.csv"
write("", file = file)
for (x in enriched){
  name <- names[i]
  x <- as.data.frame(x)
  x <- dplyr::select(x, KEGG_2019_Mouse.Term, KEGG_2019_Mouse.Overlap,
                     KEGG_2019_Mouse.Adjusted.P.value, KEGG_2019_Mouse.Genes)
  x <- dplyr::filter(x, KEGG_2019_Mouse.Adjusted.P.value < 0.2)
  i <- i + 1
  write(name, file = file, append = T)
  write.table(x, file = file, append = T)
}

# Generate Heatmap with DE genes of all. LogFC as color -- only cone genes
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
cone.genes.deg.clusterd <- pheatmap(full.table, cluster_rows = T, cluster_cols = FALSE,
         breaks = seq(-1.5,1.5, by = 0.05), 
         color = colors(60), fontsize_row = 4)
ggsave("deg_cones_involved_heatmap.png", plot = cone.genes.deg.clusterd, height = 30, width = 10, units = "cm")
write.table(x = full.table, file = "/home/daniel/master_thesis/bassoon_data/Output/Tables_Graphs/deg.cones.heatmap.logfc.csv", sep = ",")


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
all.genes.deg <- pheatmap(full.table, cluster_cols = F, cluster_rows = F,
         breaks = seq(-1.5,1.5, by = 0.05), 
         color = colors(60))
all.genes.deg.clustered <- pheatmap(full.table, cluster_cols = F, 
                                    cluster_rows = T,
         breaks = seq(-1.5,1.5, by = 0.05), 
         color = colors(60), fontsize_row = 4)
ggsave("deg_all_heatmap.png", plot = all.genes.deg.clustered, height = 60, width = 25, units = "cm")





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
path.pval <- lapply(path.pval, function(x) if (nrow(x) == 0){
  x <- as.data.frame(matrix(nrow = length(pathways)))
  x$Term <- pathways
  x$Adjusted.P.value <- c(0)
  x$V1 <- NULL
  return(x)}
  else {return(x)}
)

path.table <- data.frame(matrix(nrow = length(pathways)))
path.table$Term <- pathways
new.list <- list()
new.list$Term <- path.table
new.list <- append(new.list, path.pval)
full.table <- join_all(new.list, by = "Term", type = "left")
rownames(full.table) <- full.table$Term
colnames(full.table) <- paste0("V", seq(1, length(full.table)))
full.table <- select(full.table, 3:length(full.table))
colnames(full.table) <- names(path.pval)
heat.table <- ifelse(full.table == 0 | is.na(full.table), 0, 1)
pathway.sep.heatmap <- pheatmap(heat.table, cluster_rows = T, cluster_cols = F)
ggsave(plot = pathway.sep.heatmap, file = "kegg_pathway.png")
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
path.pval <- lapply(path.pval, function(x) if (nrow(x) == 0){
  x <- as.data.frame(matrix(nrow = length(pathways)))
  x$Term <- pathways
  x$Adjusted.P.value <- c(0)
  x$V1 <- NULL
  return(x)}
  else {return(x)}
)

path.table <- data.frame(matrix(nrow = length(pathways)))
path.table$Term <- pathways
new.list <- list()
new.list$Term <- path.table
new.list <- append(new.list, path.pval.test)
full.table <- join_all(new.list, by = "Term", type = "left")
rownames(full.table) <- full.table$Term
colnames(full.table) <- paste0("V", seq(1, length(full.table)))
full.table <- select(full.table, 3:length(full.table))
colnames(full.table) <- names(path.pval)
heat.table <- ifelse(is.na(full.table), 0, 1)
pathway.all.heatmap <- pheatmap(heat.table, cluster_rows = FALSE, cluster_cols = FALSE)
write.table(x = heat.table, file = "pathways.all.heatmap.csv", sep = ",")
dev.off()

pdf("pathway_sep_heatmap.pdf")
pathway.sep.heatmap
dev.off()
pdf("pathway_all_heatmap.pdf")
pathway.all.heatmap
dev.off()
pdf("cone_deg_heatmap.pdf", height = 18)
cone.genes.deg.clusterd
dev.off()
pdf("all_deg_heatmap.pdf", height = 35)
all.genes.deg.clustered
dev.off()


## DEG TABLES WITHOUT COMMON GENES
deg.filtered <- list()

# DE genes in WT1 rod vs cones
table <- deg$wt1.rod.cone
genes <- table$gene
# Filter out genes
# genes <- setdiff(genes, deg$mut.rod.cone$gene)
# genes <- setdiff(genes, deg$wt.mut.rod$gene)
# genes <- setdiff(genes, deg$wt.mut.cone$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["wt1.rod.cone"]] <- subset

###############

# DE genes in WT2 rod vs cones
table <- deg$wt2.rod.cone
genes <- table$gene
# Filter out genes, which are upregulated in wt cones
# genes <- setdiff(genes, deg$ko.rod.cone$gene)
# genes <- setdiff(genes, deg$wt.ko.rod$gene)
# genes <- setdiff(genes, deg$wt.ko.cone$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["wt2.rod.cone"]] <- subset

###############

# DE genes in mut rod vs cones
table <- deg$mut.rod.cone
genes <- table$gene
# Filter out genes, which are upregulated in wt cones
genes <- setdiff(genes, deg$wt1.rod.cone$gene)
# genes <- setdiff(genes, deg$wt.mut.rod$gene)
# genes <- setdiff(genes, deg$wt.mut.cone$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["mut.rod.cone"]] <- subset

###############

# DE genes in ko rod vs cones
table <- deg$ko.rod.cone
genes <- table$gene
# Filter out genes, which are upregulated in wt cones
# genes <- setdiff(genes, deg$wt1.rod.cone$gene)
genes <- setdiff(genes, deg$wt2.rod.cone$gene)
# genes <- setdiff(genes, deg$wt.ko.rod$gene)
# genes <- setdiff(genes, deg$wt.ko.cone$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["ko.rod.cone"]] <- subset

###############

deg.filtered[["wt1.wt2.cone"]] <- deg$wt1.wt2.cone
deg.filtered[["wt1.wt2.rod"]] <- deg$wt1.wt2.rod


# DE genes in wt1 cone vs mut cone
table <- deg$wt.mut.cone
genes <- table$gene
# Filter out genes
# genes <- setdiff(genes, deg$wt.mut.rod$gene)
# genes <- setdiff(genes, deg$mut.rod.cone$gene)
# genes <- setdiff(genes, deg$wt1.rod.cone$gene)
# genes <- setdiff(genes, deg$wt1.wt.2.cone$gene)
# genes <- setdiff(genes, deg$wt2.rod.cone$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["wt1.mut.cone"]] <- subset

################

# DE genes in wt1 rod vs mut rod
table <- deg$wt.mut.rod
genes <- table$gene
# Filter out genes
# genes <- setdiff(genes, deg$wt.mut.cone$gene)
# genes <- setdiff(genes, deg$mut.rod.cone$gene)
# genes <- setdiff(genes, deg$wt1.rod.cone$gene)
# genes <- setdiff(genes, deg$wt1.wt.2.rod$gene)
# genes <- setdiff(genes, deg$wt2.rod.cone$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["wt1.mut.rod"]] <- subset

###############

# DE genes in wt2 vs ko cone
table <- deg$wt.ko.cone
genes <- table$gene
# Filter out genes
# genes <- setdiff(genes, deg$wt1.wt2.cone$gene)
# genes <- setdiff(genes, deg$wt.ko.rod$gene)
# genes <- setdiff(genes, deg$ko.rod.cone$gene)
# genes <- setdiff(genes, deg$wt1.rod.cone$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["wt2.ko.cone"]] <- subset


###############

# DE genes in wt2 vs ko cone
table <- deg$wt.ko.rod
genes <- table$gene
# Filter out genes
# genes <- setdiff(genes, deg$wt1.wt2.rod$gene)
# genes <- setdiff(genes, deg$wt.ko.cone$gene)
# genes <- setdiff(genes, deg$wt2.rod.cone$gene)
# genes <- setdiff(genes, deg$ko.rod.cone$gene)
# genes <- setdiff(genes, deg$wt1.rod.cone$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["wt2.ko.rod"]] <- subset


###############

# DE genes in mut cones to ko cones
table <- deg$mut.ko.cone
genes <- table$gene
# Filter out genes
genes <- setdiff(genes, deg$wt1.wt2.cone$gene)
# genes <- setdiff(genes, deg$wt.mut.rod$gene)
# genes <- setdiff(genes, deg$wt.ko.rod$gene)
# genes <- setdiff(genes, deg$wt.mut.cone$gene)
# genes <- setdiff(genes, deg$wt.ko.cone$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["mut.ko.cone"]] <- subset

###############

# DE genes in mut rods to ko rods
table <- deg$mut.ko.rod
genes <- table$gene
# Filter out genes
genes <- setdiff(genes, deg$wt1.wt2.rod$gene)
# genes <- setdiff(genes, deg$wt.mut.cone$gene)
# genes <- setdiff(genes, deg$wt.ko.cone$gene)
# genes <- setdiff(genes, deg$wt.mut.rod$gene)
# genes <- setdiff(genes, deg$wt.ko.rod$gene)
subset <- subset(table, gene %in% genes)
deg.filtered[["mut.ko.rod"]] <- subset

### Generate filtered heatmap
gene.list <- lapply(deg.filtered, "[", , c(1,3))
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
full.table[is.na(full.table)] <- 0
colors <- colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))
all.deg.filtered <- pheatmap(full.table, cluster_cols = F, 
                                    cluster_rows = T,
                                    breaks = seq(-1.5,1.5, by = 0.05), 
                                    color = colors(60), fontsize_row = 5)

mut.and.ko <- full.table %>% select(wt1.mut.cone, wt2.ko.cone, mut.ko.cone)
mut.and.ko <- mut.and.ko[rowSums(mut.and.ko) != 0, ]
mut.and.ko.heat <- pheatmap(mut.and.ko, cluster_cols = F, 
                             cluster_rows = T,
                             breaks = seq(-1.5,1.5, by = 0.05), 
                             color = colors(60), fontsize_row = 5,
                            fontsize_col = 8, angle_col = 0,
                            labels_col = c("WT cones vs. mutant cones", "WT cones vs. knockout cones", "mutant cones vs. knockout cones"))
ggsave(plot = mut.and.ko.heat, filename = "/home/daniel/master_thesis/bassoon_data/Output/Tables_Graphs/mut_ko_cone_heatmap.png")

# rod.cone <- full.table %>% select(wt1.rod.cone, wt2.rod.cone, mut.rod.cone, ko.rod.cone)
rod.cone <- full.table %>% select(mut.rod.cone, ko.rod.cone)
rod.cone <- rod.cone[rowSums(rod.cone) != 0, ]
rod.cone.heat <- pheatmap(rod.cone, cluster_cols = F, 
                            cluster_rows = T,
                            breaks = seq(-1.5,1.5, by = 0.05), 
                            color = colors(60), fontsize_row = 5, border_color = NA, 
                          labels_col = c("Bsn mutant rods vs. cones", "Bsn knockout rods vs cones"),
                          angle_col = 0)
ggsave(plot = rod.cone.heat, filename = "/home/daniel/master_thesis/bassoon_data/Output/Tables_Graphs/rod_cone_heatmap.png")

mut.ko.rod <- full.table %>% select(wt1.mut.rod, wt2.ko.rod, mut.ko.rod)
mut.ko.rod <- mut.ko.rod[rowSums(mut.ko.rod) != 0, ]
mut.ko.rod.heat <- pheatmap(mut.ko.rod, cluster_cols = F, 
                          cluster_rows = T,
                          breaks = seq(-1.5,1.5, by = 0.05), 
                          color = colors(60), fontsize_row = 6)
ggsave(plot = mut.ko.rod.heat, filename = "/home/daniel/master_thesis/bassoon_data/Output/Tables_Graphs/mut_ko_rod_heatmap.png",
       height = 40, width = 20, units = "cm")
######################

# databases <- c("GO_Biological_Process_2018")
databases <- c("KEGG_2019_Mouse")
enriched.filt <- list()
i <- 1
names <- names(enriched)

for (table in deg.filtered){
  name <- names(deg.filtered)[i]
  table$expr <- ifelse(table$avg_logFC > 0, "up", "down")
  up <- dplyr::filter(table, expr == "up")
  down <- dplyr::filter(table, expr == "down")
  deg.sep[[paste(name, "up", sep = ".")]] <- up
  deg.sep[[paste(name, "down", sep = ".")]] <- down
  enriched.filt.all <- enrichr(table$gene, databases)
  enriched.filt[[paste(name, "all", sep = ".")]] <- enriched.filt.all
  enriched.filt.up <- enrichr(up$gene, databases)
  enriched.filt[[paste(name, "up", sep = ".")]] <- enriched.filt.up
  enriched.filt.down <- enrichr(down$gene, databases)
  enriched.filt[[paste(name, "down", sep = ".")]] <- enriched.filt.down
  i <- i + 1
}
up.down <- grepl(names(enriched.filt), pattern = "up|down")
new <- enriched.filt[up.down]

kegg.list <- lapply(new, function(x){
  # x <- as.data.frame(x$GO_Biological_Process_2018)
  x <- as.data.frame(x$KEGG_2019_Mouse)
  x <- dplyr::filter(x, Adjusted.P.value < 0.2)
  return(x)
})

path.pval <- lapply(kegg.list, "[", , c(1,4))
pathways <- unique(unlist(lapply(kegg.list, "[", , 1)))
path.pval <- lapply(path.pval, function(x) if (nrow(x) == 0){
  x <- as.data.frame(matrix(nrow = length(pathways)))
  x$Term <- pathways
  x$Adjusted.P.value <- c(0)
  x$V1 <- NULL
  return(x)}
  else {return(x)}
)

path.table <- data.frame(matrix(nrow = length(pathways)))
path.table$Term <- pathways
new.list <- list()
new.list$Term <- path.table
new.list <- append(new.list, path.pval)
full.table <- join_all(new.list, by = "Term", type = "left")
rownames(full.table) <- full.table$Term
colnames(full.table) <- paste0("V", seq(1, length(full.table)))
full.table <- select(full.table, 3:length(full.table))
colnames(full.table) <- names(path.pval)
heat.table <- ifelse(full.table == 0 | is.na(full.table), 0, 1)
pathway.sep.heatmap <- pheatmap(heat.table, cluster_rows = T, cluster_cols = F, fontsize_row = 6)
ggsave(plot = pathway.sep.heatmap, file = "go_term_heatmap.png", height = 20, units = "cm")

#############################
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
m <- rownames(markers)
