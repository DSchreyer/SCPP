# Identify Marker gene expression
# library(Matrix)
library(Seurat)
library(tibble)
# library(stringr)
library(org.Mm.eg.db)
library(dplyr)

input.dir <- "/home/daniel/master_thesis/bassoon_data/Output/post_ct_ident/"
output.dir <- "/home/daniel/master_thesis/bassoon_data/Output/downstream_analyses/"

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
      seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
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

Idents(seurat.objects) <- seurat.objects$ct.gt

subset <- subset(seurat.objects, subset  = CellType != "Unknown")



# Compare Rods and Cones to each other --> Different vulnerability
deg.rod.cone.1 <- FindMarkers(subset, ident.1 = "C57BL/6J_WT.Rod", ident.2 = "C57BL/6J_WT.Cone", logfc.threshold = 0.1) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)
deg.rod.cone.2 <- FindMarkers(subset, ident.1 = "C57BL/6NJ_WT.Rod", ident.2 = "C57BL/6NJ_WT.Cone", logfc.threshold = 0.1) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)

# Compare Mutant and KO Cones to WTs
deg.wt.mut.cone <- FindMarkers(subset, ident.1 = "C57BL/6J_WT.Cone", ident.2 = "C57BL/6J_Bsn_mut.Cone", logfc.threshold = 0.1) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)
deg.wt.ko.cone <- FindMarkers(subset, ident.1 = "C57BL/6NJ_WT.Cone", ident.2 = "C57BL/6NJ_Bsn_KO.Cone", logfc.threshold = 0.1) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)

# Compare Mutant to KO Cones
deg.mut.ko <- FindMarkers(subset, ident.1 = "C57BL/6J_Bsn_mut.Cone", ident.2 = "C57BL/6NJ_Bsn_KO.Cone", logfc.threshold = 0.1) %>% 
  rownames_to_column("gene") %>% filter(p_val_adj <= 0.05)


# LIBRARY ENRICHR!!!

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