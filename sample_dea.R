# Identify Marker gene expression
library(Matrix)
library(Seurat)
library(tibble)
library(dplyr)
library(xtable)
library(enrichR)

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

# sample 2: TSNE with Identified Populations
sample <- subset(seurat.objects, subset = orig.ident == "Bsn.221932")
sample <- FindVariableFeatures(sample)
sample <- ScaleData(sample)
sample <- RunPCA(sample)
sample <- RunTSNE(sample, dims = 1:6)
sample <- RunUMAP(sample, dims = 1:10)
DimPlot(sample, reduction = "tsne")
DimPlot(sample, reduction = "umap")



Idents(seurat.objects) <- seurat.objects$ct.gt

subset <- subset(seurat.objects, subset  = CellType != "Unknown")


# Compare Rods and Cones to each other --> Different vulnerability
deg <- list()
deg[["rod.cone.wt.1"]] <- FindMarkers(subset, ident.1 = "C57BL/6J_WT.Rod",
                                   ident.2 = "C57BL/6J_WT.Cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)
deg[["rod.cone.wt.2"]] <- FindMarkers(subset, ident.1 = "C57BL/6NJ_WT.Rod",
                                   ident.2 = "C57BL/6NJ_WT.Cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)

deg[["rod.cone.mut"]] <- FindMarkers(subset, ident.1 = "C57BL/6J_Bsn_mut.Rod", 
                               ident.2 = "C57BL/6J_Bsn_mut.Cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)
deg[["rod.cone.ko"]] <- FindMarkers(subset, ident.1 = "C57BL/6NJ_Bsn_KO.Rod", 
                                     ident.2 = "C57BL/6NJ_Bsn_KO.Cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)
deg[["cone.wt"]] <- FindMarkers(subset, ident.1 = "C57BL/6J_WT.Cone", 
                                  ident.2 = "C57BL/6NJ_WT.Cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)
deg[["rod.wt"]] <- FindMarkers(subset, ident.1 = "C57BL/6J_WT.Rod", 
                                ident.2 = "C57BL/6NJ_WT.Rod", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)


# Compare rod vs. rod, cone vs. cone
# mutant strain
deg[["wt.mut.rod"]] <- FindMarkers(subset, ident.1 = "C57BL/6J_WT.Rod",
                                    ident.2 = "C57BL/6J_Bsn_mut.Rod", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)
deg[["wt.mut.cone"]] <- FindMarkers(subset, ident.1 = "C57BL/6J_WT.Cone",
                                    ident.2 = "C57BL/6J_Bsn_mut.Cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)

# KO strain
deg[["wt.ko.rod"]] <- FindMarkers(subset, ident.1 = "C57BL/6NJ_WT.Rod",
                                   ident.2 = "C57BL/6NJ_Bsn_KO.Rod", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)
deg[["wt.ko.cone"]] <- FindMarkers(subset, ident.1 = "C57BL/6NJ_WT.Cone",
                                   ident.2 = "C57BL/6NJ_Bsn_KO.Cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)


# Compare Mutant to KO
deg[["mut.ko.cone"]] <- FindMarkers(subset, ident.1 = "C57BL/6J_Bsn_mut.Cone", 
                               ident.2 = "C57BL/6NJ_Bsn_KO.Cone", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)
deg[["mut.ko.rod"]] <- FindMarkers(subset, ident.1 = "C57BL/6J_Bsn_mut.Rod", 
                               ident.2 = "C57BL/6NJ_Bsn_KO.Rod", logfc.threshold = 0.1, min.pct = 0.05) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)


# Enrichment analysis
dbs <- listEnrichrDbs()
websiteLive <- ifelse(is.null(dbs), FALSE, TRUE)
if (websiteLive) head(dbs)

databases <- c("GO_Biological_Process_2018", "KEGG_2019_Mouse", "MGI_Mammalian_Phenotype_Level_4_2019", "Mouse_Gene_Atlas", "GO_Molecular_Function_2018")

enriched <- list()
i <- 1

deg.sep <- list()
for (table in deg){
  name <- names(deg)[i]
  table$expr <- ifelse(table$avg_logFC > 0, "up", "down")
  up <- filter(table, expr == "up")
  down <- filter(table, expr == "down")
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