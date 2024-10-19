library(Seurat)
immmune_int <- readRDS("/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/raw_data/TICAtlas.rds")
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
Tcell_obj <- FindClusters(Tcell_obj, resolution = 0.8, algorithm = 1, verbose = FALSE)

pdf(paste(savedir, "UMAP/Tcell_res0.6.pdf", sep = ""), width = 6, height = 6)
print(DimPlot(immmune_int, reduction = "umap", label = TRUE, raster = TRUE))
dev.off()
