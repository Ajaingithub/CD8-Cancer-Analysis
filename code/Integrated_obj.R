library(Seurat)
immmune_int <- readRDS("/diazlab/data3/.abhinav/.immune/cancer_combined/project/raw_data/TICAtlas.rds")
DefaultAssay(immmune_int) <- "integrated"
# immmune_int <- RunUMAP(immmune_int, dims = 1:30, reduction = "pca")

savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/paper_integrated/"
dir.create(paste0(savedir, "UMAP"), showWarnings = F, recursive = T)
pdf(paste0(savedir, "UMAP/integrated_obj.pdf"), width = 12.4, height = 5.5)
DimPlot(immmune_int, reduction = "umap", group.by = "cell_type", label = TRUE)
dev.off()

immmune_int <- FindNeighbors(immmune_int, dims = 1:30)
immmune_int <- FindClusters(immmune_int, resolution = 0.6, algorithm = 1, verbose = FALSE)


pdf(paste(savedir, "UMAP/integrated_res0.6.pdf", sep = ""), width = 6, height = 6)
print(DimPlot(immmune_int, reduction = "umap", label = TRUE, raster = TRUE))
dev.off()

dir.create(paste(saveDir, "QC_Vln/", sep = ""), showWarnings = FALSE)
pdf(paste(saveDir, "QC_Vln/", objname, "_", Assay, "_", process, ".pdf", sep = ""), width = 9, height = 6)
print(VlnPlot(ADT_integrated_impute, QC_features))
dev.off()

#### Lets try only with the T cells since clustering is mainly happening in the myeloid or B cells
Tcells_names <- grep("T cells", levels(immmune_int@meta.data$cell_type), value = TRUE)
Tcell_type <- c(Tcells_names, "T helper cells", "Th17 cells")
Tcellnames <- rownames(immmune_int@meta.data[grep(paste0("^", Tcell_type, "$", collapse = "|"), immmune_int@meta.data$cell_type), ])
Tcell_obj <- subset(immmune_int, cells = Tcellnames)

Tcell_obj <- RunUMAP(Tcell_obj, dims = 1:30, reduction = "pca")

pdf(paste(savedir, "UMAP/Tcell_subcelltype.pdf", sep = ""), width = 6, height = 6)
DimPlot(Tcell_obj, reduction = "umap", group.by = "subtype")
dev.off()

pdf(paste(savedir, "UMAP/Tcell_celltype.pdf", sep = ""), width = 9, height = 6)
DimPlot(Tcell_obj, reduction = "umap", group.by = "cell_type", label = TRUE)
dev.off()


Tcell_type <- grep("T cells|T helper|Th17", unique(immune@meta.data$cell_type), value = TRUE)

rm(plot_list)
plot_list <- list()
for (i in 1:length(Tcell_type)) {
    cellhighlight <- rownames(immune@meta.data[grep(Tcell_type[i], immune@meta.data$cell_type), ])
    p <- DimPlot(immune, group.by = "cell_type", cells.highlight = cellhighlight, cols.highlight = "red", label = FALSE) + ggtitle(Tcell_type[i])
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "UMAP/splitted_Tcells.pdf"))
plot_list
dev.off()

Tcell_obj <- FindNeighbors(Tcell_obj, dims = 1:30)

genes <- c("CD3E", "CD3D", "CD3G", "TRAC", "TRBC", "PTPRC", "CD4", "IL7R", "FOXP3", "CXCR5", "CD8A", "CD8B", "GZMB", "PRF1", "GZMK")
# res <- c(1, 1.2, 1.4, 1.6, 1.8, 2)
DefaultAssay(Tcell_obj) <- "integrated"
res <- 1.4
for (i in 1:length(res)) {
    Tcell_obj <- FindClusters(
        Tcell_obj,
        resolution = res[i]
    )

    dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
    pdf(paste0(savedir, "UMAP/dimplot_", res[i], ".pdf"), width = 12, height = 5.5)
    print(DimPlot(Tcell_obj, reduction = "umap", label = TRUE))
    dev.off()

    p <- DotPlot(Tcell_obj, genes, assay = "RNA") +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
        coord_flip() +
        scale_size(range = c(1, 10)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "/dotplot/CD4_CD8_RNA_", res[i], ".pdf"), height = 10, width = 15)
    print(p)
    dev.off()
}

pdf(paste(savedir, "UMAP/Tcell_res0.6.pdf", sep = ""), width = 6, height = 6)
print(DimPlot(immmune_int, reduction = "umap", label = TRUE, raster = TRUE))
dev.off()


rm(plot_list)
plot_list <- list()
for (i in 1:length(Tcell_type)) {
    cellhighlight <- rownames(Tcell_obj@meta.data[grep(Tcell_type[i], Tcell_obj@meta.data$cell_type), ])
    p <- DimPlot(Tcell_obj,
        group.by = "cell_type", reduction = "umap",
        cells.highlight = cellhighlight, cols.highlight = "red", label = FALSE
    ) + ggtitle(Tcell_type[i])
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "UMAP/splitted_celltypes_Tcells.pdf"))
plot_list
dev.off()

CD4_subset <- subset(x = Tcell_obj, subset = CD4 > 1)
CD8_subset <- subset(x = Tcell_obj, subset = (CD8A > 1 | CD8B > 1))

pdf(paste0(savedir, "UMAP/CD4_Tcells.pdf"))
DimPlot(CD4_subset, group.by = "cell_type")
dev.off()

pdf(paste0(savedir, "UMAP/CD8_Tcells.pdf"))
DimPlot(CD8_subset, group.by = "cell_type")
dev.off()

other_clus <- grep(paste(c("^", 0, 2, 3, 8, 9, 13, 19, 20, 21, 22, 25, 16, 28, "$"), collapse = "$|^", sep = ""), Tcell_obj@meta.data$seurat_clusters, value = TRUE, invert = TRUE) %>% unique()

Tcells_obj@meta.data$seurat_clusters <- factor(Tcells_obj@meta.data$seurat_clusters, levels = c(0, 2, 3, 8, 9, 13, 19, 20, 21, 22, 25, 16, 28, other_clus))

genes <- c("CD4", "CD8A", "CD8B")
p <- DotPlot(Tcells_obj, genes, assay = "RNA", group.by = "seurat_clusters") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD4_CD8_RNA_celltype.pdf"), height = 3, width = 12)
print(p)
dev.off()

DefaultAssay(Tcells_obj) <- "RNA"

limiting <- c(5, 5, 3.5)
dir.create(paste0(savedir, "featureplot"), showWarnings = FALSE)
genes <- c("CD8A", "CD8B", "CD4")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- FeaturePlot(Tcells_obj, features = genes[i], reduction = "umap", min.cutoff = 0.5, max.cutoff = 5) +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra, limits = c(0.5, limiting[i]))
    # scale_color_gradientn(colors = c("royalblue2", "royalblue1", "lightpink2", "lightpink2", "lightpink2", "indianred1", "indianred3", "indianred4", "firebrick4"))
}

dir.create(paste0(savedir, "featureplot"), showWarnings = FALSE)
pdf(paste0(savedir, "featureplot/CD8_CD4_2.pdf"), width = 5.5, height = 15)
plot_list[[1]] + plot_list[[2]] + plot_list[[3]]
dev.off()

dir.create(paste0(savedir, "UMAP/"), showWarnings = FALSE)
pdf(paste0(savedir, "UMAP/Tcell_dimplot.pdf"))
DimPlot(Tcells_obj, label = TRUE, label.size = 5)
dev.off()

dir.create(paste0(savedir, "UMAP/"), showWarnings = FALSE)
pdf(paste0(savedir, "UMAP/Tcell_dimplot.pdf"))
DimPlot(Tcells_obj, label = TRUE)
dev.off()

saveRDS(Tcell_obj, "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/paper_integrated/Tcell_obj.RDS")

### extracting out the CD8 T cells
clusters <- c(0, 2, 3, 8, 9, 13, 19, 20, 21, 22, 25)
CD8_Tcells <- subset(Tcell_obj, idents = clusters)

savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/CD8/"
dir.create(paste0(savedir, "UMAP"), showWarning = FALSE, recursive = TRUE)

pdf(paste0(savedir, "UMAP/CD8_dimplot.pdf"))
DimPlot(CD8_Tcells)
dev.off()

pdf(paste0(savedir, "UMAP/CD8_celltype.pdf"))
DimPlot(CD8_Tcells, group.by = "cell_type")
dev.off()
