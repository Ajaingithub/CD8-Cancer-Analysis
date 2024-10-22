library(Seurat)
library(ggplot2)

### CD8 T cellanmes from here
immune <- readRDS("/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/raw_data/TICAtlas_RNA.rds")
savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/CD8/CCA_integrated/"
dir.create(paste0(savedir, "UMAP"), showWarning = FALSE, recursive = TRUE)

CD8_Tcellnames <- read.table("/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/CD8/paper_integrated/Table/CD8_Tcellnames.txt", header = FALSE)[, 1]
CD8_obj <- subset(immune, cells = CD8_Tcellnames)

dir.create(paste0(savedir, "QC"), showWarnings = FALSE)
CD8_obj[["percent.mt"]] <- PercentageFeatureSet(CD8_obj, pattern = "^MT\\.")
pdf(paste0(savedir, "QC/feature_mt_percent.pdf"))
VlnPlot(CD8_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf(paste0(savedir, "QC/nCount_RNA.pdf"))
VlnPlot(CD8_obj, features = c("nCount_RNA"), ncol = 1) + coord_cartesian(ylim = c(0, 50000))
dev.off()

CD8_obj <- NormalizeData(CD8_obj, normalization.method = "LogNormalize", scale.factor = 10000)
CD8_subset <- subset(CD8_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 10)
CD8_subset <- subset(CD8_subset, features = rownames(CD8_subset)[Matrix::rowSums(CD8_subset@assays$RNA@counts > 0) >= 3])

cellnames <- rownames(CD8_subset@meta.data[grep("OC", CD8_subset@meta.data$subtype, invert = TRUE), ])
CD8_subset_2 <- subset(CD8_subset, cells = cellnames)
CD8_subset_2 <- NormalizeData(CD8_subset_2)

pdf(paste0(savedir, "QC/feature_mt_percent_after_QC.pdf"))
VlnPlot(CD8_subset_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#### Performing the CCA integration
source("/diazlab/data3/.abhinav/resources/all_scripts/R/scCITESeq_sctransform_V2.R")
objname <- "CD8_subset_2"
Assay <- "RNA"
process <- "sctransform"
Tcell_obj_subset_sctransformed <- sctransform_V2_integration(
    obj = CD8_subset_2,
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

Tcell_obj_subset_sctransformed <- RunUMAP(Tcell_obj_subset_sctransformed, dims = 1:30, reduction = "pca")

pdf(paste0(savedir, "UMAP/Tcell_integrated.pdf"))
DimPlot(Tcell_obj_subset_sctransformed, reduction = "umap")
dev.off()

pdf(paste0(savedir, "UMAP/Tcell_integrated_celltype.pdf"))
DimPlot(Tcell_obj_subset_sctransformed, reduction = "umap", group.by = "cell_type")
dev.off()

pdf(paste0(savedir, "UMAP/Tcell_integrated_celltype_splitted.pdf"))
DimPlot(Tcell_obj_subset_sctransformed, reduction = "umap", group.by = "cell_type", split.by = "cell_type", ncol = 3) + NoLegend()
dev.off()

pdf(paste0(savedir, "UMAP/Tcell_integrated_cancertype.pdf"))
DimPlot(Tcell_obj_subset_sctransformed, reduction = "umap", group.by = "subtype") + NoLegend()
dev.off()

pdf(paste0(savedir, "UMAP/Tcell_integrated_cancertype_splitted.pdf"))
DimPlot(Tcell_obj_subset_sctransformed, reduction = "umap", group.by = "subtype", split.by = "subtype", ncol = 3)
dev.off()

pdf(paste0(savedir, "UMAP/Tcell_integrated_patient_splitted.pdf"))
DimPlot(Tcell_obj_subset_sctransformed, reduction = "umap", group.by = "patient", split.by = "patient", ncol = 6) + NoLegend()
dev.off()

pdf(paste0(savedir, "UMAP/Tcell_integrated_patient.pdf"))
DimPlot(Tcell_obj_subset_sctransformed, reduction = "umap", group.by = "patient") + NoLegend()
dev.off()

Tcell_obj_subset_sctransformed <- FindNeighbors(Tcell_obj_subset_sctransformed, dims = 1:30)
CD8_Tcells <- Tcell_obj_subset_sctransformed

MAGIC_assay <- CD8_imputed[["MAGIC_RNA"]]
CD8_Tcells[["MAGIC_RNA"]] <- MAGIC_assay

SCM <- c("PTPRC", "CCR7", "SELL", "CD27", "CD28", "IL7R", "FAS", "ITGAL", "IL2RB", "CD58", "CXCR3", "TCF7", "LEF1")
Naive <- c("PTPRC", "TCF7", "LEF1", "FOXP1", "NOSIP", "CCR7", "SELL", "PTPRC", "IL7R", "TGFBR1", "TGFBR2", "FOXO1", "FOXO3", "GATA3", "SP1", "GATA3")
CM <- c("PTPRC", "CCR7", "ITGAL", "SELL", "CD44", "IL7R", "CCR7", "PTPRC", "TCF7", "BACH2", "ID3", "BCL6", "FOXO1", "MYB", "STAT3", "EOMES", "ITGA4", "CD27", "CD28", "LFA3", "CD44", "ITGA4")
EM <- c("IL2RB", "CD58", "KLRG1", "GZMM", "GZMB", "IFNG", "GZMK", "NKG7", "ZEB2", "GZMA", "CX3CR1", "IFNG", "HLA-DRA", "HLA-DRB1", "CD44")
cytotoxic <- c("GZMB", "PRF1", "TBX21", "LAMP1", "FGFBP2", "KLRG1", "GNLY", "ZNF683", "GZMH", "FASL", "IFNG", "NKG7")
TEMRA <- c("PTPRC", "FAS", "KLRG1", "NKG2D", "LFA3", "B3GAT1", "LILRB1", "TYROBP", "KIR2DL1", "KIR3DL1", "KLRK1")
Exhausted <- c("PDCD1", "LAG3", "HAVCR2", "TOX")
TRM <- c("CD69", "ITGAE", "CXCR6", "ZNF683", "PDCD1", "ENTPD1", "CCR2", "HAVCR2", "LATN", "CTLA4", "PRDM1", "RUNX3", "KLF2", "TOX")
Ag_experienc <- c("ENTPD1", "ITGAE", "PDCD1", "HAVCR2", "LAG3", "TNFRSF9", "CD69")
pro_exhausted <- c(
    "TCF1", "LEF1", "BCL6", "ID3", "SELL", "MYB", "CXCR5", "GZMK", "ENTPD1", "PDCD1", "CD69", "IL7R", "CCR7", "CD127",
    "EOMES", "CXCR5", "TCF7", "LEF1", "BCL6", "ID3", "SELL", "SLAMF6", "CCR7", "LAG3"
)
trans_exhausted <- c("TBX21", "ZEB2", "KLF2", "RUNX1", "RUNX2", "ID2", "GZMB", "IFNG", "CX3CR1")
terminal_exhausted <- c(
    "TOX", "NFAT", "EOMES", "PRDM1", "CTLA4", "TIGIT", "CD101", "PDCD1", "TNFSF8", "EGR2",
    "IRF4", "BATF", "ID2", "SOSTDC1", "IL27", "IL35", "PDCD1", "TNFSF8", "CTLA4", "TIGIT",
    "HAVCR2", "TNFRSF4", "TNFRSF9", "TNFRSF18", "TNFRSF8"
)
TEMRA <- c("PTPRC", "FAS", "ITGAL", "IL2RB", "CD58", "B3GAT1")
Tregs <- c("IL2RA", "FOXP3", "CTLA4", "PDCD1", "CD39", "STAT5A", "STAT5B", "IKZF2")

# Use co-expression of CD45RA/CD45RO markers with other differentiation markers to infer their state:
# If CCR7 and SELL (CD62L) are co-expressed → likely CD45RA (Naive/TEMRA cells).
# If CD44 and CD27 are co-expressed → likely CD45RO (Memory cells).

genes <- unique(c(SCM, Naive, CM, EM, cytotoxic, Exhausted, TRM, pro_exhausted, trans_exhausted, terminal_exhausted, TEMRA, Tregs))

res <- c(0.8, 1, 1.2, 1.4)
# res <- 1.2
DefaultAssay(CD8_Tcells) <- "integrated"
# res <- 1.4
for (i in 1:length(res)) {
    CD8_Tcells <- FindClusters(
        CD8_Tcells,
        resolution = res[i]
    )

    dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
    pdf(paste0(savedir, "UMAP/dimplot_", res[i], ".pdf"), width = 6, height = 5.5)
    print(DimPlot(CD8_Tcells, reduction = "umap", label = TRUE))
    dev.off()


    p <- DotPlot(CD8_Tcells, genes, assay = "RNA") +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
        coord_flip() +
        scale_size(range = c(1, 10)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "/dotplot/CD8_celltypes_integrated_", res[i], ".pdf"), height = 19, width = 19)
    print(p)
    dev.off()

    p <- DotPlot(CD8_Tcells, genes, assay = "MAGIC_RNA") +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
        coord_flip() +
        scale_size(range = c(1, 10)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "/dotplot/CD8_celltypes_integrated_", res[i], "_MAGIC_RNA.pdf"), height = 25, width = 19)
    print(p)
    dev.off()

    #### Highlighting the cells
    rm(plot_list)
    clus <- levels(CD8_Tcells@meta.data$seurat_clusters)
    plot_list <- list()
    for (j in 1:length(clus)) {
        cellhighlight <- rownames(CD8_Tcells@meta.data[grep(paste0("^", clus[j], "$"), CD8_Tcells@meta.data$seurat_clusters), ])
        p <- DimPlot(CD8_Tcells,
            group.by = "seurat_clusters", reduction = "umap", sizes.highlight = 0.2,
            cells.highlight = cellhighlight, cols.highlight = "red", label = FALSE
        ) + ggtitle(clus[j])
        plot_list[[j]] <- p
    }

    pdf(paste0(savedir, "UMAP/splitted_clusters_", res[i], ".pdf"), width = 6, height = 5)
    print(plot_list)
    dev.off()
}

#### Imputing the scRNAseq
DefaultAssay(CD8_subset_2) <- "RNA"
CD8_subset_2 <- NormalizeData(CD8_subset_2)
library(reticulate)
library(Rmagic)
use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check
CD8_subset_2 <- magic(CD8_subset_2, npca = 30) ## imputing the RNA data as for RNA PCs are 20
DefaultAssay(CD8_subset_2) <- "MAGIC_RNA"

saveRDS(CD8_subset_2, "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/raw_data/CD8_subset_2_impute.RDS")

CD8_Tcells[["MAGIC_RNA"]] <- CD8_subset_2[["MAGIC_RNA"]]

### CD8 T cells subset celltypes
library(ArchR)
library(UCell)
library(ggplot2)
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(CD8_Tcells) <- "RNA"
rm(markers)
markers <- list()
files <- c("SCM", "Naive", "CM", "EM", "cytotoxic", "Exhausted", "TRM", "pro_exhausted", "trans_exhausted", "terminal_exhausted", "TEMRA", "Tregs")

for (i in 1:length(files)) {
    Tcellsubset <- unique(get(files[i]))
    Tcellsubset <- rownames(CD8_Tcells)[match(Tcellsubset, rownames(CD8_Tcells), nomatch = 0)]
    markers[[files[i]]] <- Tcellsubset
    # assign(files[i], Tcellsubset)
}


# pdf(paste0(savedir, "featureplot/markers_genes.pdf"))
# plot_list
# dev.off()

DefaultAssay(CD8_Tcells) <- "MAGIC_RNA"
CD8_Tcells <- AddModuleScore(CD8_Tcells, features = markers, slot = "data")
colnames(CD8_Tcells@meta.data)[18:29] <- files

rm(plot_list)
plot_list <- list()
for (i in 1:length(files)) {
    plot_list[[i]] <- FeaturePlot(CD8_Tcells, features = files[i], reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

p3 <- VlnPlot(CD8_Tcells, features = paste0(files), group.by = "integrated_snn_res.1.2", pt.size = 0)

rm(plot_list2)
plot_list2 <- list()
for (i in 1:length(files)) {
    plot_list2[[i]] <- VlnPlot(CD8_Tcells, features = paste0(files[i]), group.by = "integrated_snn_res.1.2", pt.size = 0) + geom_boxplot()
}


rm(plot_list)
plot_list <- list()
for (i in 1:length(files)) {
    plot_list[[i]] <- FeaturePlot(CD8_Tcells, features = files, reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/CD8/CCA_integrated/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/CD8_Tcell_subset_score_1.2.pdf", sep = ""), width = 20, height = 15)
print(p3)
dev.off()

dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/CD8_Tcell_subset_score_1.2_boxplot.pdf", sep = ""), width = 10, height = 5.5)
print(plot_list2)
dev.off()

dir.create(paste(savedir, "featureplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "featureplot/CD8_Tcell_subset_score_boxplot.pdf", sep = ""))
plot_list2
dev.off()

savedir <- "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/cell_phenotyping/downstream2/CD8/CCA_integrated/"
markers_1.2 <- FindAllMarkers(CD8_Tcells, assay = "RNA", group.by = "seurat_clusters")

dir.create(paste0(savedir, "Table"), showWarnings = FALSE)
write.table(markers_1.2, paste0(savedir, "Table/CD8_markers_res_1.2.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

dir.create(paste0(savedir, "saveRDS_obj"), showWarnings = FALSE)
saveRDS(CD8_Tcells, paste0(savedir, "saveRDS_obj/CD8_Tcells_imputed_integrated.RDS"))

#### Performing FAST Topics
CD8_subset_2 <- readRDS("/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/raw_data/CD8_subset_2.RDS")
counts <- t(GetAssayData(CD8_subset_2, slot = "counts"))
counts_req <- counts[, colSums(counts) > 0] ### remove gene with 0 UMIs
cellnames <- sample(rownames(counts_req), size = 30000)
counts_req_2 <- counts_req[cellnames, ]

fit10 <- fit_topic_model(counts_req_2, k = 10)
fit20 <- fit_topic_model(counts_req_2, k = 20)
fit30 <- fit_topic_model(counts_req_2, k = 30)

saveRDS(fit10, paste0(savedir, "fit10_topic.RDS"))
saveRDS(fit20, paste0(savedir, "fit20_topic.RDS"))
saveRDS(fit30, paste0(savedir, "fit30_topic.RDS"))

# Convert the log-likelihood list to a numeric vector
loglik_values <- unlist(fit_list)

# Plot the elbow curve using base R or ggplot2
plot(k_values, loglik_values,
    type = "b", xlab = "Number of Topics (k)", ylab = "Log-Likelihood",
    main = "Elbow Plot for Topic Modeling"
)

# Alternatively, using ggplot2 for a more polished plot
library(ggplot2)
df <- data.frame(k = k_values, loglik = loglik_values)
ggplot(df, aes(x = k, y = loglik)) +
    geom_line() +
    geom_point() +
    labs(title = "Elbow Plot for Topic Modeling", x = "Number of Topics (k)", y = "Log-Likelihood") +
    theme_minimal()

### Running a trajectory
suppressPackageStartupMessages({
    library(slingshot)
    library(SingleCellExperiment)
    library(RColorBrewer)
    library(scales)
    library(viridis)
    library(UpSetR)
    library(pheatmap)
    library(msigdbr)
    library(fgsea)
    library(knitr)
    library(ggplot2)
    library(gridExtra)
    library(tradeSeq)
    library(Seurat)
    # library(bioc2020trajectories)
})

# Save the objects as separate matrices for input in slingshot
## convert back to singleCellExperiment
DefaultAssay(CD8_Tcells) <- "RNA"
sce <- as.SingleCellExperiment(CD8_Tcells, assay = "RNA")

dir.create(paste0(savedir, "trajectory/UMAP"))
pdf(paste0(savedir, "trajectory/UMAP/CD8_celltype_annotated.pdf"))
DimPlot(CD8_Tcells, reduction = "umap", label = TRUE)
dev.off()

CD4_sce <- slingshot(sce,
    clusterLabels = "seurat_clusters",
    reducedDim = "UMAP", approx_points = 100,
    omega = TRUE,
    omega_scale = 1.5
)

CD4_sce1 <- slingshot(sce,
    clusterLabels = "seurat_clusters",
    reducedDim = "UMAP", approx_points = 100,
    start.clus = "6",
    omega = TRUE,
    omega_scale = 1.5
)

library(ArchR)
library(scater)
embedded_orig <- embedCurves(CD4_sce, "UMAP")

rm(plot_list)
plot_list <- list()
for (i in 1:21) {
    embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
    embedded <- data.frame(embedded$s[embedded$ord, ])
    g <- plotUMAP(CD4_sce, colour_by = paste("slingPseudotime_", i, sep = ""))
    stopifnot(all(rownames(CD8_Tcells@reductions$umap@cell.embeddings) == rownames(g$data)))
    data <- merge(CD8_Tcells@reductions$umap@cell.embeddings, g$data, by = "row.names")
    colnames(data) <- c("cellname", "UMAP_1", "UMAP_2", "X", "Y", paste("Lineage", i, sep = ""))
    p <- ggplot(data, aes_string("UMAP_1", "UMAP_2", color = paste("Lineage", i, sep = ""))) +
        geom_point(size = 0.01) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
        geom_path(data = embedded, aes(x = umap_1, y = umap_2), color = "black", size = 1.2) +
        theme_bw()
    plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir, "trajectory/UMAP/CD8_Tcells_sce_UMAP_splitted_unbias.pdf", sep = ""), width = 25, height = 20)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
    plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
    plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], plot_list[[13]], plot_list[[14]],
    plot_list[[15]], plot_list[[16]], plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]], plot_list[[21]],
    nrow = 5, ncol = 5
)
dev.off()

library(ArchR)
library(scater)
embedded_orig <- embedCurves(CD4_sce1, "UMAP")
rm(plot_list)
plot_list <- list()

for (i in 1:21) {
    embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
    embedded <- data.frame(embedded$s[embedded$ord, ])
    g <- plotUMAP(CD4_sce1, colour_by = paste("slingPseudotime_", i, sep = ""))
    stopifnot(all(rownames(CD8_Tcells@reductions$umap@cell.embeddings) == rownames(g$data)))
    data <- merge(CD8_Tcells@reductions$umap@cell.embeddings, g$data, by = "row.names")
    colnames(data) <- c("cellname", "UMAP_1", "UMAP_2", "X", "Y", paste("Lineage", i, sep = ""))
    p <- ggplot(data, aes_string("UMAP_1", "UMAP_2", color = paste("Lineage", i, sep = ""))) +
        geom_point(size = 0.01) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
        geom_path(data = embedded, aes(x = umap_1, y = umap_2), color = "black", size = 1.2) +
        theme_bw()
    plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir, "trajectory/UMAP/CD8_Tcells_sce_UMAP_splitted_start_naive.pdf", sep = ""), width = 25, height = 20)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
    plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
    plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], plot_list[[13]], plot_list[[14]],
    plot_list[[15]], plot_list[[16]], plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]], plot_list[[21]],
    nrow = 5, ncol = 5
)
dev.off()

CD4_sce7 <- slingshot(sce,
    clusterLabels = "seurat_clusters",
    reducedDim = "UMAP", approx_points = 100,
    start.clus = "7",
    omega = TRUE,
    omega_scale = 1.5
)

library(ArchR)
library(scater)
embedded_orig <- embedCurves(CD4_sce1, "UMAP")
rm(plot_list)
plot_list <- list()

for (i in 1:21) {
    embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
    embedded <- data.frame(embedded$s[embedded$ord, ])
    g <- plotUMAP(CD4_sce1, colour_by = paste("slingPseudotime_", i, sep = ""))
    stopifnot(all(rownames(CD8_Tcells@reductions$umap@cell.embeddings) == rownames(g$data)))
    data <- merge(CD8_Tcells@reductions$umap@cell.embeddings, g$data, by = "row.names")
    colnames(data) <- c("cellname", "UMAP_1", "UMAP_2", "X", "Y", paste("Lineage", i, sep = ""))
    p <- ggplot(data, aes_string("UMAP_1", "UMAP_2", color = paste("Lineage", i, sep = ""))) +
        geom_point(size = 0.01) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
        geom_path(data = embedded, aes(x = umap_1, y = umap_2), color = "black", size = 1.2) +
        theme_bw()
    plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir, "trajectory/UMAP/CD8_Tcells_sce_UMAP_splitted_start_naive.pdf", sep = ""), width = 25, height = 20)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
    plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
    plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], plot_list[[13]], plot_list[[14]],
    plot_list[[15]], plot_list[[16]], plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]], plot_list[[21]],
    nrow = 5, ncol = 5
)
dev.off()
