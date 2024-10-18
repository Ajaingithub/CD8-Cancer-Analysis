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

# Perform clustering
t_obj <- FindClusters(
    t_obj,
    resolution = 0.6 # You can adjust the resolution parameter based on your needs
)

savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/"
dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
# Optional: Visualize the clusters using UMAP
pdf(paste0(savedir, "UMAP/harmony_dinmplot.pdf"), width = 5.5, height = 5.5)
DimPlot(t_obj, reduction = "harmonyumap", label = TRUE)
dev.off()


### Performing the imputation for CD8A, CD8B, and CD4
library(Seurat)
library(reticulate)
library(Rmagic)
use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check

t_obj <- magic(t_obj, npca = 30)

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
