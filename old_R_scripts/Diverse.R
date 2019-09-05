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
