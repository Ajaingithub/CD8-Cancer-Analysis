library(Seurat)
### extracting out the CD8 T cells
## from this /diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/CD8_Genen/code/Integrated_obj.R
clusters <- c(0, 2, 3, 5, 8, 9, 13, 19, 20, 21, 22, 25)
CD8_Tcells <- subset(Tcell_obj, idents = clusters)

DefaultAssay(CD8_Tcells) <- "RNA"
CD8_Tcells <- NormalizeData(CD8_Tcells)

savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/CD8/paper_integrated/"
dir.create(paste0(savedir, "UMAP"), showWarning = FALSE, recursive = TRUE)

pdf(paste0(savedir, "UMAP/CD8_dimplot.pdf"))
DimPlot(CD8_Tcells)
dev.off()

pdf(paste0(savedir, "UMAP/CD8_celltype.pdf"))
DimPlot(CD8_Tcells, group.by = "cell_type")
dev.off()

### Making the UMAP based on the new CD8
DefaultAssay(CD8_Tcells) <- "integrated"
CD8_Tcells <- RunUMAP(CD8_Tcells, dims = 1:30, reduction = "pca")

savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/CD8/paper_integrated/"
pdf(paste(savedir, "UMAP/CD8_Tcell_subcelltype.pdf", sep = ""), width = 6, height = 6)
DimPlot(CD8_Tcells, reduction = "umap", group.by = "subtype")
dev.off()

pdf(paste(savedir, "UMAP/CD8_Tcell_celltype.pdf", sep = ""), width = 9, height = 6)
DimPlot(CD8_Tcells, reduction = "umap", group.by = "cell_type", label = TRUE)
dev.off()

Tcell_type <- grep("T cells|T helper|Th17", unique(CD8_Tcells@meta.data$cell_type), value = TRUE)

rm(plot_list)
plot_list <- list()
for (i in 1:length(Tcell_type)) {
    cellhighlight <- rownames(CD8_Tcells@meta.data[grep(Tcell_type[i], CD8_Tcells@meta.data$cell_type), ])
    p <- DimPlot(CD8_Tcells, group.by = "cell_type", cells.highlight = cellhighlight, cols.highlight = "red", label = FALSE) + ggtitle(Tcell_type[i])
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "UMAP/splitted_Tcells.pdf"))
plot_list
dev.off()

CD8_Tcells <- FindNeighbors(CD8_Tcells, dims = 1:30)

#### CD8 T cell Phenotyping
dir.create(paste0(savedir, "featureplot"), showWarnings = FALSE)
genes <- c("CD8A", "CD8B", "CD4")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- FeaturePlot(CD8_Tcells, features = genes[i], reduction = "umap", min.cutoff = 0.5) +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra)
    # scale_color_gradientn(colors = c("royalblue2", "royalblue1", "lightpink2", "lightpink2", "lightpink2", "indianred1", "indianred3", "indianred4", "firebrick4"))
}

dir.create(paste0(savedir, "featureplot"), showWarnings = FALSE)
pdf(paste0(savedir, "featureplot/CD8_CD4_2.pdf"), width = 5.5, height = 15)
plot_list[[1]] + plot_list[[2]] + plot_list[[3]]
dev.off()

SCM <- c("PTPRC", "CCR7", "SELL", "CD27", "CD28", "IL7R", "FAS", "ITGAL", "IL2RB", "CD58", "CXCR3", "TCF7", "LEF1")
Naive <- c("TCF7", "LEF1", "FOXP1", "NOSIP", "CCR7", "SELL", "PTPRC")
CM <- c("CD44", "CCR7", "PTPRC")
EM <- c("CD58", "KLRG1", "GZMB", "GZMK", "NKG7", "ZEB2", "GZMA", "CX3CR1", "IFNG")
cytotoxic <- c("PRF1", "TBX21", "LAMP1", "FGFBP2", "KLRG1", "GLNY", "ZNF683", "GZMH", "FASL", "IFNG")
Exhausted <- c("PDCD1", "LAG3", "HAVCR2", "TIM3", "TOX")
TRM <- c("CD69", "CD103", "CXCR6")
pro_exhausted <- c("CD39", "EOMES", "CXCR5", "TCF7", "LEF1", "BCL6", "ID3", "SELL", "SLAMF6")
trans_exhausted <- c("TBX21", "ZEB2", "KLF2", "RUNX1", "RUNX2", "GZMB", "IFNG", "CX3CR1")
terminal_exhausted <- c("TOX", "NFAT", "EOMES", "BLIMP1", "CTLA4", "TIGIT", "CD101", "PDCD1", "TNFSF8")
TEMRA <- c("PTPRC", "FAS", "ITGAL", "IL2RB", "CD58", "CD57")
Tregs <- c("CD25", "FOXP3", "CTLA4", "PDCD1", "CD39", "STAT5")

# Use co-expression of CD45RA/CD45RO markers with other differentiation markers to infer their state:
# If CCR7 and SELL (CD62L) are co-expressed → likely CD45RA (Naive/TEMRA cells).
# If CD44 and CD27 are co-expressed → likely CD45RO (Memory cells).

genes <- unique(c(SCM, Naive, CM, EM, cytotoxic, Exhausted, TRM, pro_exhausted, trans_exhausted, terminal_exhausted, TEMRA, Tregs))

res <- c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
res <- 1.2
DefaultAssay(CD8_Tcells) <- "integrated"
# res <- 1.4
for (i in 1:length(res)) {
    CD8_Tcells <- FindClusters(
        CD8_Tcells,
        resolution = res[i]
    )

    dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
    pdf(paste0(savedir, "UMAP/dimplot_", res[i], ".pdf"), width = 12, height = 5.5)
    print(DimPlot(CD8_Tcells, reduction = "umap", label = TRUE))
    dev.off()


    p <- DotPlot(CD8_Tcells, genes, assay = "RNA") +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
        coord_flip() +
        scale_size(range = c(1, 10)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "/dotplot/CD8_celltypes_RNA_", res[i], ".pdf"), height = 15, width = 15)
    print(p)
    dev.off()

    #### Highlighting the cells
    rm(plot_list)
    clus <- levels(CD8_Tcells@meta.data$seurat_clusters)
    plot_list <- list()
    for (j in 1:length(clus)) {
        cellhighlight <- rownames(CD8_Tcells@meta.data[grep(paste0("^", clus[j], "$"), CD8_Tcells@meta.data$seurat_clusters), ])
        p <- DimPlot(CD8_Tcells,
            group.by = "seurat_clusters", reduction = "umap",
            cells.highlight = cellhighlight, cols.highlight = "red", label = FALSE
        ) + ggtitle(clus[j])
        plot_list[[j]] <- p
    }

    pdf(paste0(savedir, "UMAP/splitted_clusters_", res[i], ".pdf"), width = 6, height = 5)
    print(plot_list)
    dev.off()
}

dir.create(paste0(savedir, "saveRDS_obj/"))

pdf(paste(savedir, "UMAP/Tcell_res0.6.pdf", sep = ""), width = 6, height = 5)
print(DimPlot(immmune_int, reduction = "umap", label = TRUE, raster = TRUE))
dev.off()


rm(plot_list)
plot_list <- list()
for (i in 1:length(Tcell_type)) {
    cellhighlight <- rownames(CD8_Tcells@meta.data[grep(Tcell_type[i], CD8_Tcells@meta.data$cell_type), ])
    p <- DimPlot(CD8_Tcells,
        group.by = "cell_type", reduction = "umap",
        cells.highlight = cellhighlight, cols.highlight = "red", label = FALSE
    ) + ggtitle(Tcell_type[i])
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "UMAP/splitted_celltypes_Tcells.pdf"))
plot_list
dev.off()

dir.create("/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/CD8/paper_integrated/Table/", showWarning = FALSE)
CD8_Tcellnames <- rownames(CD8_Tcells@meta.data)
write.table(CD8_Tcellnames, paste0(savedir, "Table/CD8_Tcellnames.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
