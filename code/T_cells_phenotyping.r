#### Extracting out the T cells
library(Seurat)
library(dplyr)
immune <- readRDS("/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/raw_data/TICAtlas_RNA.rds")
Tcells_names <- grep("T cells", levels(immune@meta.data$cell_type), value = TRUE)
Tcell_type <- c(Tcells_names, "T helper cells", "Th17 cells")
Tcellnames <- rownames(immune@meta.data[grep(paste0("^", Tcell_type, "$", collapse = "|"), immune@meta.data$cell_type), ])
Tcell_obj <- subset(immune, cells = Tcellnames)

genes <- c("CD4", "CD8A", "CD8B", "IFIT3", "FOXP3", "TIGIT", "IFNG", "CCL4", "IL7R", "GZMA", "PDCD1", "HAVCR2")
p <- VlnPlot(Tcell_obj, genes, group.by = "cell_type", pt.size = 0)

dir.create(paste0(savedir, "Vlnplot"), showWarnings = FALSE)
pdf(paste0(savedir, "Vlnplot/Tcell_genes.pdf"), height = 15, width = 15)
p
dev.off()

rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- FeaturePlot(immune, features = genes[i], reduction = "umap", cols = c("grey", "red"), min.cutoff = 1.5)
    # scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste0(savedir, "featureplot"), showWarnings = FALSE)
pdf(paste0(savedir, "featureplot/immune_CD8_CD4.pdf"), width = 5.5, height = 15)
plot_list[[1]] + plot_list[[2]] + plot_list[[3]]
dev.off()

genes <- c("CD8A", "CD8B", "CD4")
pdf(paste0(savedir, "featureplot/immune_CD8_CD4_min1.5.pdf"))
FeaturePlot(immune, features = genes, reduction = "umap", cols = c("grey", "red"), min.cutoff = 1.5, max.cutoff = 5)
dev.off()

##### Performing the analysis
savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/Tcells/"
dir.create(savedir, showWarnings = FALSE)

### Normalization without Integration
DefaultAssay(Tcell_obj) <- "RNA"

Tcell_obj <- NormalizeData(Tcell_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# features <- rownames(Tcell_obj)
# chunk_size <- 10000

# for (i in seq(1, length(features), by = chunk_size)) {
#     chunk <- features[i:min(i + chunk_size - 1, length(features))]
#     Tcell_obj <- ScaleData(Tcell_obj, features = chunk, do.center = TRUE, do.scale = TRUE)
# }
Tcell_obj <- ScaleData(Tcell_obj)

Tcell_obj <- FindVariableFeatures(Tcell_obj, nfeatures = 4000)
Tcell_obj <- Tcell_obj %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/Tcells/"
dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE, recursive = TRUE)
pdf(paste(savedir, "UMAP/Tcell.pdf", sep = ""))
DimPlot(Tcell_obj, reduction = "umap", group.by = "cell_type")
dev.off()

dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/Tcell_vague_patient.pdf", sep = ""))
DimPlot(Tcell_obj, reduction = "umap", group.by = "patient") + NoLegend()
dev.off()

dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/Tcell_vague_subtype.pdf", sep = ""))
DimPlot(Tcell_obj, reduction = "umap", group.by = "subtype")
dev.off()

#### For now I am using harmony because it is fast but for later would go to sctransform with CCA
### It require integration
library(harmony)
Tcell_obj <- RunHarmony(
    Tcell_obj,
    group.by.vars = "patient",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "harmony"
)

Tcell_obj <- RunUMAP(
    Tcell_obj,
    assay = "RNA",
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

p1 <- DimPlot(Tcell_obj, reduction = "harmonyumap", label = TRUE, group.by = "cell_type")
p2 <- DimPlot(Tcell_obj, reduction = "harmonyumap", label = FALSE, group.by = "patient") + NoLegend()
p3 <- DimPlot(Tcell_obj, reduction = "harmonyumap", label = TRUE, group.by = "subtype")

pdf(paste0(savedir, "UMAP/harmony_celltype_patient_subtype.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

genes <- c("CD8A", "CD8B", "CD4")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- FeaturePlot(Tcell_obj, features = genes[i], reduction = "harmonyumap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste0(savedir, "featureplot"), showWarnings = FALSE)
pdf(paste0(savedir, "featureplot/CD8_CD4.pdf"), width = 5.5, height = 15)
plot_list[[1]] + plot_list[[2]] + plot_list[[3]]
dev.off()

p <- DotPlot(Tcell_obj, genes, assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD4_CD8_SCT.pdf"), height = 5, width = 10)
print(p)
dev.off()

### Performing the clustering
# Find neighbors based on the UMAP reduction
Tcell_obj <- FindNeighbors(
    Tcell_obj,
    reduction = "harmony",
    dims = 1:30
)

### Performing the imputation for CD8A, CD8B, and CD4
library(Seurat)
library(reticulate)
library(Rmagic)
use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check

Tcell_obj <- magic(Tcell_obj, npca = 30)

DefaultAssay(Tcell_obj) <- "MAGIC_RNA"

genes <- c("CD8A", "CD8B", "CD4")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- FeaturePlot(Tcell_obj, features = genes[i], reduction = "harmonyumap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste0(savedir, "featureplot"), showWarnings = FALSE)
pdf(paste0(savedir, "featureplot/CD8_CD4_imputation.pdf"), width = 5.5, height = 15)
plot_list[[1]] + plot_list[[2]] + plot_list[[3]]
dev.off()

# Perform clustering
DefaultAssay(Tcell_obj) <- "RNA"
genes <- c("CD3E", "CD3D", "CD3G", "TRAC", "TRBC", "PTPRC", "CD4", "IL7R", "FOXP3", "CXCR5", "CD8A", "CD8B", "GZMB", "PRF1", "GZMK")

res <- c(0.2, 0.4, 0.6)
for (i in 1:length(res)) {
    Tcell_obj <- FindClusters(
        Tcell_obj,
        resolution = res[i] # You can adjust the resolution parameter based on your needs
    )
    # savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/"
    dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
    # Optional: Visualize the clusters using UMAP
    pdf(paste0(savedir, "UMAP/harmony_dinmplot_", res[i], ".pdf"), width = 5.5, height = 5.5)
    print(DimPlot(Tcell_obj, reduction = "harmonyumap", label = TRUE))
    dev.off()

    p <- DotPlot(Tcell_obj, genes, assay = "RNA") +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
        coord_flip() +
        scale_size(range = c(1, 10)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "/dotplot/CD4_CD8_RNA_", res[i], ".pdf"), height = 5, width = 15)
    print(p)
    dev.off()
}

p <- DotPlot(Tcell_obj, genes, assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD4_CD8_RNA_", res[i], ".pdf"), height = 5, width = 15)
print(p)
dev.off()

genes <- c("CD4", "CD8A", "CD8B", "GZMB", "PRF1", "GZMK")
savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/Tcells/"
p <- DotPlot(Tcell_obj, genes, assay = "RNA", group.by = "cell_type") +
    scale_color_gradientn(colors = c(
        "royalblue2", "royalblue1", "royalblue1", "lightblue2", "lightblue2",
        "lightpink2", "indianred1", "indianred3", "indianred4", "firebrick4"
    )) +
    # scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    # coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD4_CD8_RNA_only_Tcells_4.pdf"), height = 7, width = 9)
print(p)
dev.off()

rm(plot_list)
plot_list <- list()
for (i in 1:length(Tcell_type)) {
    cellhighlight <- rownames(Tcell_obj@meta.data[grep(Tcell_type[i], immune@meta.data$cell_type), ])
    p <- DimPlot(Tcell_obj,
        group.by = "cell_type", reduction = "harmonyumap",
        cells.highlight = cellhighlight, cols.highlight = "red", label = FALSE
    ) + ggtitle(Tcell_type[i])
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "UMAP/splitted_celltypes.pdf"))
plot_list
dev.off()


#### Performing integration using RPCA
savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/Tcells/"
source("/diazlab/data3/.abhinav/resources/all_scripts/R/scCITESeq_ADT_merging.R")
Tcell_obj <- ADT_merging(Tcell_obj, savedir,
    Assay = "RNA", dims = 30, numfeatures = 3000,
    process = "RPCA", objname = "Tcell_RPCA",
    sample_tree = NULL, split_by = "subtype",
    reference = NULL
)

DefaultAssay(Tcell_obj) <- "RNA"
Tcell_obj <- NormalizeData(Tcell_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell_obj <- FindVariableFeatures(Tcell_obj, nfeatures = 4000)
Tcell_obj <- ScaleData(Tcell_obj)


Tcell_obj <- Tcell_obj %>%
    RunPCA(verbose = FALSE)

library(harmony)
Tcell_obj <- RunHarmony(
    Tcell_obj,
    group.by.vars = "subtype",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "harmony"
)

Tcell_obj <- RunUMAP(
    Tcell_obj,
    assay = "RNA",
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

p1 <- DimPlot(t_obj, reduction = "harmonyumap", label = TRUE, group.by = "cell_type")
p2 <- DimPlot(t_obj, reduction = "harmonyumap", label = FALSE, group.by = "patient") + NoLegend()
p3 <- DimPlot(t_obj, reduction = "harmonyumap", label = TRUE, group.by = "subtype")

pdf(paste0(savedir, "UMAP/harmony_celltype_patient_subtype.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

Tcells_names <- grep("T cells", levels(immune@meta.data$cell_type), value = TRUE)
Tcell_type <- c(Tcells_names, "T helper cells", "Th17 cells")

rm(plot_list)
plot_list <- list()
for (i in 1:length(Tcell_type)) {
    cellhighlight <- rownames(Tcell_obj@meta.data[grep(Tcell_type[i], Tcell_obj@meta.data$cell_type), ])
    p <- DimPlot(Tcell_obj,
        group.by = "cell_type", reduction = "harmonyumap",
        cells.highlight = cellhighlight, cols.highlight = "red", label = FALSE
    ) + ggtitle(Tcell_type[i])
    plot_list[[i]] <- p
}
savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/Tcells/"
pdf(paste0(savedir, "UMAP/splitted_celltypes_harmony_Tcells.pdf"))
plot_list
dev.off()
