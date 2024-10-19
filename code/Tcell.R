library(Seurat)
library(dplyr)
library(ArchR)

immune <- readRDS("/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/raw_data/TICAtlas_RNA.rds")
metadata <- read.csv("/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/raw_data/TICAtlas_metadata.csv")
stopifnot(all(immune@meta.data$cell_type == metadata$cell_type)) ## Check if this true

celltypes <- table(metadata$cell_type) %>% as.data.frame()

savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/downstream/"
dir.create(paste0(savedir, "/Table"), showWarnings = TRUE, recursive = TRUE)
write.table(celltypes, paste0(savedir, "Table/celltypes.txt"), row.names = F, col.names = T, sep = ",", quote = F)

### Subsetting the T cells which are not subsetted as CD4 or CD8
Tcell_vague <- c("Naive T cells", "Proliferative T cells", "Regulatory T cells")

cellnames <- rownames(immune@meta.data[grep(paste0("^", Tcell_vague, "$", collapse = "|"), immune@meta.data$cell_type), ])

t_obj <- subset(immune, cells = cellnames)

### Normalization without Integration
DefaultAssay(t_obj) <- "RNA"

t_obj <- NormalizeData(t_obj, normalization.method = "LogNormalize", scale.factor = 10000)
t_obj <- ScaleData(t_obj, features = rownames(immune))

t_obj <- FindVariableFeatures(t_obj)
t_obj <- t_obj %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

# pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/EndoDimPlot_patient.pdf")
dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/Tcell_vague.pdf", sep = ""))
DimPlot(t_obj, reduction = "umap", group.by = "cell_type")
dev.off()

dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/Tcell_vague_patient.pdf", sep = ""))
DimPlot(t_obj, reduction = "umap", group.by = "patient") + NoLegend()
dev.off()

dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/Tcell_vague_subtype.pdf", sep = ""))
DimPlot(t_obj, reduction = "umap", group.by = "subtype")
dev.off()

#### For now I am using harmony because it is fast but for later would go to sctransform with CCA
### It require integration
library(harmony)
t_obj <- RunHarmony(
    t_obj,
    group.by.vars = "patient",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "harmony"
)

t_obj <- RunUMAP(
    t_obj,
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

genes <- c("CD8A", "CD8B", "CD4")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- FeaturePlot(t_obj, features = genes[i], reduction = "harmonyumap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste0(savedir, "featureplot"), showWarnings = FALSE)
pdf(paste0(savedir, "featureplot/CD8_CD4.pdf"), width = 5.5, height = 15)
plot_list[[1]] + plot_list[[2]] + plot_list[[3]]
dev.off()

p <- DotPlot(t_obj, genes, assay = "RNA") +
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
t_obj <- FindNeighbors(
    t_obj,
    reduction = "harmony",
    dims = 1:30
)

### Performing the imputation for CD8A, CD8B, and CD4
library(Seurat)
library(reticulate)
library(Rmagic)
use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check

t_obj <- magic(t_obj, npca = 30)

DefaultAssay(t_obj) <- "MAGIC_RNA"

genes <- c("CD8A", "CD8B", "CD4")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- FeaturePlot(t_obj, features = genes[i], reduction = "harmonyumap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste0(savedir, "featureplot"), showWarnings = FALSE)
pdf(paste0(savedir, "featureplot/CD8_CD4_imputation.pdf"), width = 5.5, height = 15)
plot_list[[1]] + plot_list[[2]] + plot_list[[3]]
dev.off()

# Perform clustering
DefaultAssay(t_obj) <- "RNA"

res <- c(1.2, 1.4, 1.6, 1.8, 2)
for (i in 1:length(res)) {
    t_obj <- FindClusters(
        t_obj,
        resolution = res[i] # You can adjust the resolution parameter based on your needs
    )
    # savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/"
    dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
    # Optional: Visualize the clusters using UMAP
    pdf(paste0(savedir, "UMAP/harmony_dinmplot_", res[i], ".pdf"), width = 5.5, height = 5.5)
    print(DimPlot(t_obj, reduction = "harmonyumap", label = TRUE))
    dev.off()

    p <- DotPlot(t_obj, genes, assay = "RNA") +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
        coord_flip() +
        scale_size(range = c(1, 10)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "/dotplot/CD4_CD8_RNA_", res[i], ".pdf"), height = 5, width = 10)
    print(p)
    dev.off()

    # p <- DotPlot(t_obj, genes, assay = "MAGIC_RNA") +
    #     scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    #     coord_flip() +
    #     scale_size(range = c(1, 10)) +
    #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    # dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
    # pdf(paste0(savedir, "/dotplot/CD4_CD8_MAGIC_RNA_", res[i], ".pdf"), height = 5, width = 10)
    # print(p)
    # dev.off()
}

### Extracting CD8 T cells from cluster 5 and 6
t_obj_UMAP <- t_obj@reductions$harmonyumap@cell.embeddings
stopifnot(all(rownames(t_obj@meta.data) == rownames(t_obj_UMAP)))
t_obj_UMAP_df <- as.data.frame(as.matrix(t_obj_UMAP))
stopifnot(all(rownames(t_obj@meta.data) == rownames(t_obj_UMAP_df)))
t_obj_UMAP_df$seurat_clusters <- t_obj@meta.data$seurat_clusters
CD4_CD8_transpose <- t(t_obj@assays$RNA@data[c("CD4", "CD8A", "CD8B"), ])
t_obj_UMAP_CD4_CD8 <- cbind(t_obj_UMAP_df, CD4_CD8_transpose)

### extracting cluster 5 and 6
t_obj_UMAP_CD4_CD8_clus5_6 <- t_obj_UMAP_CD4_CD8[grep(c("^5$|^6$"), t_obj_UMAP_CD4_CD8$seurat_clusters), ]
t_obj_UMAP_CD4_CD8_clus5_6[t_obj_UMAP_CD4_CD8_clus5_6$CD4 > 1, "harmonyUMAP_2"]
max(t_obj_UMAP_CD4_CD8_clus5_6[t_obj_UMAP_CD4_CD8_clus5_6$CD4 > 1, "harmonyUMAP_2"])

summary(t_obj_UMAP_CD4_CD8_clus5_6[
    t_obj_UMAP_CD4_CD8_clus5_6$CD4 > 1 &
        t_obj_UMAP_CD4_CD8_clus5_6$CD8A < 1 &
        t_obj_UMAP_CD4_CD8_clus5_6$CD8B < 1,
    "harmonyUMAP_2"
])


## Since mean and median is 2.2
### I will keep the cutoff of 2.2 the cell below -2.2 are CD4 and above are CD8
t_obj_UMAP_CD4_CD8_clus5_6$new_seurat_cluster <- t_obj_UMAP_CD4_CD8_clus5_6$seurat_clusters
t_obj_UMAP_CD4_CD8_clus5_6$new_seurat_cluster <- as.character(t_obj_UMAP_CD4_CD8_clus5_6$new_seurat_cluster) ## to remove the factors
t_obj_UMAP_CD4_CD8_clus5_6[t_obj_UMAP_CD4_CD8_clus5_6$harmonyUMAP_2 > -2.2, "new_seurat_cluster"] <- "CD8"
t_obj_UMAP_CD4_CD8_clus5_6[t_obj_UMAP_CD4_CD8_clus5_6$harmonyUMAP_2 <= -2.2, "new_seurat_cluster"] <- "CD4"
t_obj_UMAP_CD4_CD8_clus5_6$CD4_CD8_clus <- paste(t_obj_UMAP_CD4_CD8_clus5_6$new_seurat_cluster, t_obj_UMAP_CD4_CD8_clus5_6$seurat_cluster, sep = "_")

### Adding it to the seurat cluster
# t_obj@meta.data$seurat_clusters2 <- t_obj@meta.data$seurat_clusters
# t_obj@meta.data$seurat_clusters2 <- as.character(t_obj@meta.data$seurat_clusters2) ## to remove the factors
t_obj@meta.data$seurat_clusters2 <- t_obj_UMAP_CD4_CD8_clus5_6[match(rownames(t_obj@meta.data), rownames(t_obj_UMAP_CD4_CD8_clus5_6)), "CD4_CD8_clus"]
t_obj@meta.data[is.na(t_obj@meta.data$seurat_clusters2), "seurat_clusters2"] <- "clus_"
t_obj@meta.data$seurat_clusters3 <- paste(t_obj@meta.data$seurat_clusters2, t_obj@meta.data$seurat_clusters, sep = "_")

pdf(paste0(savedir, "UMAP/harmony_dimplot_1_CD4_CD8.pdf", sep = ""))
DimPlot(t_obj, reduction = "harmonyumap", group.by = "seurat_clusters3", label = TRUE)
dev.off()

p <- DotPlot(t_obj, genes, assay = "RNA", group.by = "seurat_clusters3") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    # coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD4_CD8_RNA_seurat_clus_3.pdf"), height = 7, width = 5)
print(p)
dev.off()

### Trying to identify all the CD8 T cells
### Checking CD8 T cells in all the object of the immune celltypes
savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/"
pdf(paste0(savedir, "UMAP/immune_celltypes.pdf"), width = 12, height = 5)
DimPlot(immune, reduction = "umap", group.by = "cell_type", label = TRUE)
dev.off()

genes <- c("CD3E", "CD3D", "CD3G", "TRAC", "TRBC", "PTPRC", "CD4", "IL7R", "FOXP3", "CXCR5", "CD8A", "CD8B", "GZMB", "PRF1", "GZMK")
p <- DotPlot(immune, genes, assay = "RNA", group.by = "cell_type", scale = FALSE) +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    # coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/immune_celltype_CD4_CD8_noscaling.pdf"), height = 10, width = 10)
print(p)
dev.off()

Tcell_type <- grep("T cells|T helper|Th17", unique(immune@meta.data$cell_type), value = TRUE)

rm(plot_list)
plot_list <- list()
for (i in 1:length(Tcell_type)) {
    cellhighlight <- rownames(immune@meta.data[grep(Tcell_type[i], immune@meta.data$cell_type), ])
    p <- DimPlot(immune, group.by = "cell_type", cells.highlight = cellhighlight, cols.highlight = "red", label = FALSE) + ggtitle(Tcell_type[i])
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "UMAP/splitted_celltypes.pdf"))
plot_list
dev.off()

immune@meta.data$cell_type_broad <- immune@meta.data$cell_type
immune@meta.data$cell_type_broad <- as.character(immune@meta.data$cell_type_broad)
immune@meta.data[grep("CD4 T cells|Th17 cells", immune@meta.data$cell_type_broad), "cell_type_broad"] <- "CD4"
immune@meta.data[grep("CD8 T cells", immune@meta.data$cell_type_broad), "cell_type_broad"] <- "CD8"
immune@meta.data[grep("T helper cells|Regulatory T cells|Naive T cells", immune@meta.data$cell_type_broad), "cell_type_broad"] <- "T cells"

Tcell_type <- c("CD4", "CD8", "T cells")
rm(plot_list)
plot_list <- list()
for (i in 1:length(Tcell_type)) {
    cellhighlight <- rownames(immune@meta.data[grep(Tcell_type[i], immune@meta.data$cell_type_broad), ])
    p <- DimPlot(immune, group.by = "cell_type_broad", cells.highlight = cellhighlight, cols.highlight = "red", label = FALSE) + ggtitle(Tcell_type[i])
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "UMAP/splitted_celltypes_braod.pdf"))
plot_list
dev.off()

CD4_subset <- subset(x = immune, subset = CD4 > 1.5)
CD8_subset <- subset(x = immune, subset = (CD8A > 1.5 | CD8B > 1.5))

pdf(paste0(savedir, "UMAP/CD4_subset.pdf"))
DimPlot(CD4_subset)
dev.off()

pdf(paste0(savedir, "UMAP/CD8_subset.pdf"))
DimPlot(CD8_subset)
dev.off()


source("/diazlab/data3/.abhinav/resources/all_scripts/R/scCITESeq_sctransform_V2.R")
objname <- "Tcell_obj"
Assay <- "RNA"
process <- "sctransform"
Tcell_obj_sctransformed <- sctransform_V2_integration(
    obj = Tcell_obj,
    saveDir = savedir,
    ngenes = 4000,
    regress = c("nCount_RNA"),
    dims = 30,
    Assay = Assay,
    process = process,
    objname = objname,
    split_by = "subtype",
    reference = NULL,
    sample_tree = NULL
)

Tcell_subset <- sample(rownames(Tcell_obj@meta.data), size = 70000)
Tcell_obj_subset <- subset(Tcell_obj, cells = Tcell_subset)

savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/Tcells_subset/"
dir.create(savedir, showWarnings = FALSE)

source("/diazlab/data3/.abhinav/resources/all_scripts/R/scCITESeq_sctransform_V2.R")
objname <- "Tcell_obj_subset"
Assay <- "RNA"
process <- "sctransform"
Tcell_obj_subset_sctransformed <- sctransform_V2_integration(
    obj = Tcell_obj_subset,
    saveDir = savedir,
    ngenes = 4000,
    regress = c("nCount_RNA"),
    dims = 30,
    Assay = Assay,
    process = process,
    objname = objname,
    split_by = "subtype",
    reference = NULL,
    sample_tree = NULL
)
