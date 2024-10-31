library(Seurat)
library(ggplot2)

### CD8 T cellanmes from here
immune <- readRDS("/diazlab/data3/.abhinav/.immune/cancer_combined/project/raw_data/TICAtlas_RNA.rds")
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

SCM <- c("CCR7", "SELL", "CD27", "CD28", "IL7R", "FAS", "ITGAL", "IL2RB", "CD58", "CXCR3", "TCF7", "LEF1", "BCL2")
Naive <- c("TCF7", "LEF1", "FOXP1", "NOSIP", "CCR7", "SELL", "IL7R", "TGFBR1", "TGFBR2", "FOXO1", "FOXO3", "GATA3", "SP1", "GATA3")
CM <- c("CCR7", "ITGAL", "SELL", "CD44", "IL7R", "CCR7", "TCF7", "BACH2", "ID3", "BCL6", "FOXO1", "MYB", "STAT3", "EOMES", "ITGA4", "CD27", "CD28", "LFA3", "CD44", "ITGA4")
EM <- c("IL2RB", "CD58", "KLRG1", "GZMM", "GZMB", "IFNG", "GZMK", "NKG7", "ZEB2", "GZMA", "CX3CR1", "IFNG", "HLA-DRA", "HLA-DRB1", "CD44")
cytotoxic <- c("GZMB", "PRF1", "TBX21", "LAMP1", "FGFBP2", "KLRG1", "GNLY", "ZNF683", "GZMH", "FASL", "IFNG", "NKG7")
TEMRA <- c("FAS", "KLRG1", "NKG2D", "LFA3", "B3GAT1", "LILRB1", "TYROBP", "KIR2DL1", "KIR3DL1", "KLRK1")
Exhausted <- c("PDCD1", "LAG3", "HAVCR2", "TOX")
TRM <- c("CD69", "ITGAE", "CXCR6", "ZNF683", "PDCD1", "ENTPD1", "CCR2", "HAVCR2", "LATN", "CTLA4", "PRDM1", "RUNX3", "KLF2", "TOX")
Ag_experienced <- c("ENTPD1", "ITGAE", "PDCD1", "HAVCR2", "LAG3", "TNFRSF9", "CD69")
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
TEMRA <- c("FAS", "ITGAL", "IL2RB", "CD58", "B3GAT1")
Tregs <- c("IL2RA", "FOXP3", "CTLA4", "PDCD1", "CD39", "STAT5A", "STAT5B", "IKZF2")
Activated <- c("FOS", "FOSB", "JUN", "CD31", "OX40", "CD69", "TNFRSF9", "IL2RA")
Interferon <- c("ISG15", "MX1", "MX2", "OAS1", "OAS2", "TNFSF19", "IFI44", "IFI44L")
TEMRA <- c("TYROBP", "FCER1G", "KLRB1", "KLRD1", "KLRC1", "IRF8")

# Use co-expression of CD45RA/CD45RO markers with other differentiation markers to infer their state:
# If CCR7 and SELL (CD62L) are co-expressed → likely CD45RA (Naive/TEMRA cells).
# If CD44 and CD27 are co-expressed → likely CD45RO (Memory cells).

genes <- unique(c(SCM, Naive, CM, EM, cytotoxic, Exhausted, TRM, pro_exhausted, trans_exhausted, terminal_exhausted, TEMRA, Tregs, Activated, Interferon, TEMRA))

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

savedir <- c("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/")

Tumor_reactive <- c(
    "ENTPD1", "PDCD1", "TIM3", "LAG3", "TNFRSF9",
    "LAYN", "TOX", "HAVCR2", "CXCL13", "CTLA4", "TNFRSF18", "CXCR6"
)
newly_entering <- c("CXCR4", "CCR7", "SELL", "S1PR1", "CXCR5", "CCR5", "CXCR3", "S1PR1", "GZMK", "CXCR4", "KLRG1", "EOMES")
Trm <- c("CD69", "ITGAE", "ZNF683", "GZMB", "PRF1", "IFNG", "CX3CR1", "TNF", "BHLHE40", "NR4A1")

genes <- unique(c(Tumor_reactive, Trm, newly_entering))
p <- DotPlot(CD8_Tcells_noCD4, genes, assay = "RNA", group.by = "Celltypes2") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/Ag_experience_newly_entered.pdf"), width = 8, height = 9)
print(p)
dev.off()


### Final Genes
clus_CT <- read.table("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/Table/cluster_CT.txt", header = TRUE, sep = "\t")
patterns <- clus_CT$Clusters
replacements <- clus_CT$Celltypes

CD8_Tcells_noCD4@meta.data$celltypes <- CD8_Tcells_noCD4@meta.data$seurat_clusters
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD8_Tcells_noCD4@meta.data$celltypes <- str_replace_all(CD8_Tcells_noCD4@meta.data$celltypes, pattern, replacements[i])
}

naive <- c("CCR7", "SELL", "CD27", "CD28", "IL7R", "TCF7", "LEF1", "TGFBR2", "CXCR3")
SCM <- c("CD44", "EOMES", "KLF2", "FAS", "ITGAL", "EOMES", "FAS", "ITGAL", "CXCR5", "GZMK")
CM_activated <- c("CD69", "FOS", "JUN")
CM_EM <- c("FOXP1", "NOSIP", "BCL6", "STAT3", "PRDM1", "RUNX3", "IRF4", "EGR2", "TNFRSF9", "STAT5A", "STAT5B")
new_genes <- c("ITGAL", "LY6A", "CD58", "KLRG1", "GZMB", "GZMK", "ZEB2", "HLA-DRA", "HLA-DRB1")
EM_ISG <- c("ISG15", "IFI44", "OAS1", "OAS2", "MX1", "MX2", "IFI44")
EM2 <- c("IFNG", "NKG7", "GZMA", "ZEB2", "PRF1", "TBX21")
EM_NK <- c("TYROBP", "FCER1G", "KLRB1", "KLRD1", "KLRC1", "IRF8")
TEMRA <- c("GNLY", "ZNF683", "GZMH", "FGFBP2", "LAMP1")
Proliferating <- c("MKI67", "TOP2A", "MYB")
Pre_TEX <- c("BCL2", "ID3", "STAT3", "ITGA4", "PDCD1", "LAG3", "HAVCR2", "TOX", "ITGAE", "CXCR6", "ENTPD1", "CTLA4", "ID2", "TIGIT")
Ter_TEX <- c("TNFRSF4", "TNFRSF9", "TNFRSF18", "TNFRSF8", "TNFSF8")
TRM <- c("ITGAE", "ZNF683", "ITGA1", "TRDC", "TRGC1", "TRGC2")
MAIT <- c("KLRB1", "ZBTB16")
NK <- c("NCAM1")

genes <- unique(c(naive, SCM, CM_activated, new_genes, CM_EM, EM_ISG, EM2, EM_NK, TEMRA, Pre_TEX, Ter_TEX, Proliferating, TRM, MAIT))

p <- DotPlot(CD8_Tcells_noCD4, genes, assay = "RNA", group.by = "subtype") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
pdf(paste0(savedir, "dotplot/celltypes_genes_subtypes.pdf"), width = 8, height = 20)
p
dev.off()

# CD8_Tcells_noCD4@meta.data$celltypes <- factor(CD8_Tcells_noCD4@meta.data$celltypes,
#     levels = c(
#         "Naive", "CM", "CM Activated", "CM-EM1", "CM-EM2", "EM ISG", "EM2", "EM2 NK", "TEMRA",
#         "TEMRA CD62L", "Pre TEX1", "Pre TEX2", "Trans TEX", "Term TEX", "Proliferating"
#     )
# )

# cellnames_noCD4 <- row.names(CD8_Tcells@meta.data[!is.na(CD8_Tcells@meta.data$celltypes), ])
# CD8_Tcells_noCD4 <- subset(CD8_Tcells, cells = cellnames_noCD4)
# CD8_Tcells_noCD4@meta.data$Celltypes2 <- as.character(CD8_Tcells_noCD4@meta.data$celltypes)
# CD8_Tcells_noCD4$Celltypes2 <- gsub("Naive|^CM$", "SCM-like", CD8_Tcells_noCD4$Celltypes2) %>%
#     gsub("CM Activated", "Effector Activated", .) %>%
#     gsub("Pre TEX1", "TEX1-ISG", .) %>%
#     gsub("^CM-", "", .)

patterns <- c("EM2", "SCM-like", "EM1", "Effector Activated", "TEMRA", "GD CD62L", "EM ISG", "TRM", "GD", "TEX ISG")
replacements <- c("Effector", "Stem-like", "Early Eff", "Activated", "Cytotoxic", "Cytotoxic", "ISG Eff", "TRM", "Gamma Delta", "ISG TEX")

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/Final_presentation/"
CD8_Tcells_noCD4@meta.data$celltypes2 <- CD8_Tcells_noCD4@meta.data$Celltypes2
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD8_Tcells_noCD4@meta.data$celltypes2 <- str_replace_all(CD8_Tcells_noCD4@meta.data$celltypes2, pattern, replacements[i])
}

pdf(paste0(savedir, "UMAP/CD8_Tcells_3.pdf"), width = 5, height = 5.5)
DimPlot(CD8_Tcells_noCD4, group.by = "celltypes1", reduction = "umap", label = TRUE) + NoLegend()
dev.off()

# Define your color palette (adjust colors as needed)
celltype_colors <- c(
    "Effector" = "#c237b9",
    "Stem-like" = "mediumseagreen",
    "Early Eff" = "#777ee0",
    "ISG Eff" = "#444fe3",
    "Activated" = "#f28deb",
    "Cytotoxic" = "#b87f16",
    "Gamma Delta" = "#edb141",
    "TRM" = ""
    "Pre TEX" = "#EF9A9A",
    "ISG TEX" = "#e08989",
    "Trans TEX" = "#F44336",
    "Term TEX" = "#db2a2a",
    "Proliferating" = "#BBDEFB"
)

pdf(paste0(savedir, "UMAP/CD8_Tcells_2_noLabels.pdf"), width = 5, height = 5.5)
DimPlot(CD8_Tcells_noCD4, group.by = "celltypes1", reduction = "umap", label = FALSE, cols = celltype_colors) + NoLegend()
dev.off()

pdf(paste0(savedir, "UMAP/CD8_Tcells_2_noLabels_nocol.pdf"), width = 5, height = 5.5)
DimPlot(CD8_Tcells_noCD4, group.by = "celltypes1", reduction = "umap", label = FALSE) + NoLegend()
dev.off()


Cairo::CairoPDF(paste0(savedir, "UMAP/CD8_Tcells.pdf"), width = 7, height = 5.5)
DimPlot(CD8_Tcells_noCD4, group.by = "celltypes1", reduction = "umap", pt.size = 0.000001, label = TRUE)
dev.off()

CD8_Tcells_noCD4@meta.data$Celltypes2 <- factor(CD8_Tcells_noCD4@meta.data$celltypes,
    levels = c(
        "SCM-like", "Effector Activated", "Early EM", "TRM", "EM ISG", "EM2", "GD", "GD CD62L", "TEMRA",
        "TEX ISG", "Pre TEX", "Trans TEX", "Term TEX", "Proliferating"
    )
)

p <- DotPlot(CD8_Tcells_noCD4, genes, assay = "RNA", group.by = "Celltypes2") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_celltypes_RNA_2.pdf"), height = 20, width = 10)
print(p)
dev.off()

p <- DotPlot(CD8_Tcells_noCD4, proteins, assay = "predicted_ADT", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_celltypes_protein.pdf"), height = 14, width = 8)
print(p)
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
pdf(paste0(savedir, "UMAP/celltypes_noCD4.pdf"), width = 7, height = 5.5)
DimPlot(CD8_Tcells_noCD4, reduction = "umap", group.by = "celltypes", label = TRUE)
dev.off()

### DotPlot for the markers
score <- c("SCM", "Naive", "CM", "EM", "cytotoxic", "Exhausted", "TRM", "pro_exhausted", "trans_exhausted", "terminal_exhausted", "TEMRA", "Tregs")

p <- DotPlot(CD8_Tcells_noCD4, score, group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_celltypes_scores.pdf"), height = 8, width = 9)
print(p)
dev.off()

#### Identifying the Tissue Resident Memory (TRM), Tumor entering and Tumore reactive T cells
Trm <- c("CD69", "ITGAE", "ZNF683", "CD101", "CCR2", "HAVCR2", "LAYN", "CTLA4", "CX3CR1", "CXCR6", "PRDM1", "RUNX3", "KLF2", "TOX")
reactive <- c("CD39", "CD103", "PDCD1", "TIM3", "LAG3", "TNFRSF9", "LAYN", "TOX", "GZMB", "PRF1", "TNFRSF18")
newly <- c("CCR7", "CD27", "CD28", "SELL", "CXCR5", "CCR5", "CXCR3", "CCR7", "S1PR1", "MKI67")

genes <- unique(c(Trm, reactive, newly))

p <- DotPlot(CD8_Tcells_noCD4, genes, group.by = "celltypes", assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_celltypes_Trm_reactive_newly.pdf"), height = 12, width = 9)
print(p)
dev.off()

Trm <- c("CD69", "ITGAE", "ZNF683", "ITGA1", "CD101", "CCR2", "CX3CR1", "CXCR6", "PRDM1", "RUNX3", "BHLHE40", "NR4A1", "CD44", "IL2RB")
reactive <- c("ENTPD1", "ITGAE", "PDCD1", "HAVCR2", "LAG3", "TNFRSF9", "LAYN", "TOX", "GZMB", "PRF1", "TNFRSF18", "GZMA", "TIGIT", "CXCL13")
newly <- c("CCR7", "CD27", "CD28", "SELL", "CXCR5", "CCR5", "CXCR3", "CCR7", "S1PR1", "MKI67", "FCGR3A")

genes <- unique(c(Trm, reactive, newly))

p <- DotPlot(CD8_Tcells_noCD4, genes, group.by = "seurat_clusters", assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_clus_Trm_reactive_newly.pdf"), height = 12, width = 15)
print(p)
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
genes <- c("GPR183", "LMNA", "ANXA1", "CXCR4", "TNF", "IFNG")

p <- DotPlot(CD8_Tcells_noCD4, genes, group.by = "Celltypes2", assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_IFNG_TNF_imputed.pdf"), height = 7, width = 15)
print(p)
dev.off()

genes <- c("CALM3", "CALM2", "CALM1")

p <- DotPlot(Azizi_obj, genes, group.by = "Celltypes2", assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_GO_topic3_Azizi.pdf"), height = 5, width = 8)
print(p)
dev.off()

BC <- c("SELL", "TCF7", "LEF1", "CD44", "KLF2", "PRDM1", "ETS1", "JUN", "TBX21", "PRF1")
BCC <- c("HSP90AB1", "HSPD1", "HSPA4", "HSP90AA1", "DNAJB6")
CM <- c("EOMES", "SLAMF6", "SLAMF7", "GZMK", "TOX", "LAG3", "HAVCR2")
UM <- c("PDCD1", "CD27", "ENTPD1")
CRC <- c("KLRB1", "GNLY", "GZMB", "IRGAE", "TNFRSF18")
EA <- c("IFITM1", "IFITM2", "IFITM3", "IFNG", "EGR1")
PDAC <- c("MKI67", "TOP2A", "TUBB")
RCC <- c("NKG7", "GZMH", "KLRG1", "ITGB1")
HCC <- c("IRF1", "KLRF1", "CCL3", "CCL4")

genes <- unique(c(BC, BCC, CM, UM, CRC, EA, PDAC, RCC, HCC))
cancerlevel <- c("BC", "BCC", "CM", "UM", "CRC", "EA", "PDAC", "RCC", "HCC", "ICC", "NSCLC", "SCC")
CD8_Tcells_noCD4@meta.data$subtype <- factor(CD8_Tcells_noCD4@meta.data$subtype, levels = cancerlevel)
p <- DotPlot(CD8_Tcells_noCD4, genes, group.by = "subtype", assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_differential_exp_subtype.pdf"), height = 12, width = 8)
print(p)
dev.off()

### Gene SCore for these
library(ArchR)
library(UCell)
library(ggplot2)
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(CD8_Tcells) <- "RNA"
rm(markers)
markers <- list()
files <- c("Trm", "reactive", "newly")

for (i in 1:length(files)) {
    Tcellsubset <- unique(get(files[i]))
    Tcellsubset <- rownames(CD8_Tcells)[match(Tcellsubset, rownames(CD8_Tcells), nomatch = 0)]
    markers[[files[i]]] <- Tcellsubset
    # assign(files[i], Tcellsubset)
}

DefaultAssay(CD8_Tcells) <- "MAGIC_RNA"
CD8_Tcells <- AddModuleScore(CD8_Tcells, features = markers, slot = "data", name = files)

files <- c("Trm1", "reactive2", "newly3")

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

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/CD8_Tcell_TRM_reactive_newly_1.2_2.pdf", sep = ""), width = 20, height = 15)
print(p3)
dev.off()

dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/CD8_Tcell_TRM_reactive_newly_boxplot_2.pdf", sep = ""), width = 10, height = 5.5)
print(plot_list2)
dev.off()

dir.create(paste(savedir, "featureplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "featureplot/CD8_Tcell_Tcell_TRM_reactive_2.pdf", sep = ""))
plot_list
dev.off()

p <- DotPlot(CD8_Tcells, files, group.by = "seurat_clusters") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_clus_Trm_reactive_newly_gene_score.pdf"), height = 5, width = 12)
print(p)
dev.off()

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
Stem_like <- c(
    "CCR7", "SELL", "CD27", "CD28", "IL7R", "FAS", "ITGAL", "IL2RB", "CD58", "CXCR3", "TCF7", "LEF1", "BCL2",
    "TCF7", "LEF1", "FOXP1", "NOSIP", "CCR7", "SELL", "IL7R", "TGFBR1", "TGFBR2", "FOXO1", "FOXO3", "GATA3", "SP1", "GATA3",
    "CCR7", "ITGAL", "SELL", "CD44", "IL7R", "CCR7", "TCF7", "BACH2", "ID3", "BCL6", "FOXO1", "MYB", "STAT3", "EOMES", "ITGA4", "CD27", "CD28", "LFA3", "CD44", "ITGA4"
)
EM <- c("IL2RB", "CD58", "KLRG1", "GZMM", "GZMB", "IFNG", "GZMK", "NKG7", "ZEB2", "GZMA", "CX3CR1", "IFNG", "HLA-DRA", "HLA-DRB1")
cytotoxic <- c(
    "GZMB", "PRF1", "TBX21", "LAMP1", "FGFBP2", "KLRG1", "GNLY", "ZNF683", "GZMH", "FASL", "IFNG", "NKG7",
    "KLRG1", "NKG2D", "LFA3", "B3GAT1", "LILRB1", "TYROBP", "KIR2DL1", "KIR3DL1", "KLRK1", "IFNG", "ITGA4", "CXCR6"
)
Exhausted <- c("PDCD1", "LAG3", "HAVCR2", "TOX", "ENTPD1", "ID2", "CTLA4", "TNFRSF4", "TNFRSF18", "TNFSF8")
# TRM <- c("CD69", "ITGAE", "CXCR6", "ZNF683", "PDCD1", "ENTPD1", "CCR2", "HAVCR2", "LATN", "CTLA4", "PRDM1", "RUNX3", "KLF2", "TOX")
# Ag_experienced <- c("ENTPD1", "ITGAE", "PDCD1", "HAVCR2", "LAG3", "TNFRSF9", "CD69")
# pro_exhausted <- c(
#     "TCF1", "LEF1", "BCL6", "ID3", "SELL", "MYB", "CXCR5", "GZMK", "ENTPD1", "PDCD1", "CD69", "IL7R", "CCR7", "CD127",
#     "EOMES", "CXCR5", "TCF7", "LEF1", "BCL6", "ID3", "SELL", "SLAMF6", "CCR7", "LAG3"
# )
# trans_exhausted <- c("TBX21", "ZEB2", "KLF2", "RUNX1", "RUNX2", "ID2", "GZMB", "IFNG", "CX3CR1")
# terminal_exhausted <- c(
#     "TOX", "NFAT", "EOMES", "PRDM1", "CTLA4", "TIGIT", "CD101", "PDCD1", "TNFSF8", "EGR2",
#     "IRF4", "BATF", "ID2", "SOSTDC1", "IL27", "IL35", "PDCD1", "TNFSF8", "CTLA4", "TIGIT",
#     "HAVCR2", "TNFRSF4", "TNFRSF9", "TNFRSF18", "TNFRSF8"
# )
# TEMRA <- c("FAS", "ITGAL", "IL2RB", "CD58", "B3GAT1")
# Tregs <- c("IL2RA", "FOXP3", "CTLA4", "PDCD1", "CD39", "STAT5A", "STAT5B", "IKZF2")
# Activated <- c("FOS", "FOSB", "JUN", "CD31", "OX40", "CD69", "TNFRSF9", "IL2RA")
# Interferon <- c("ISG15", "MX1", "MX2", "OAS1", "OAS2", "TNFSF19", "IFI44", "IFI44L")
# TEMRA <- c("TYROBP", "FCER1G", "KLRB1", "KLRD1", "KLRC1", "IRF8")


library(ArchR)
library(UCell)
library(ggplot2)
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(CD8_Tcells) <- "RNA"
rm(markers)
markers <- list()
files <- c("Stem_like", "EM", "cytotoxic", "Exhausted")

for (i in 1:length(files)) {
    Tcellsubset <- unique(get(files[i]))
    Tcellsubset <- rownames(CD8_Tcells)[match(Tcellsubset, rownames(CD8_Tcells), nomatch = 0)]
    markers[[files[i]]] <- Tcellsubset
    # assign(files[i], Tcellsubset)
}


# pdf(paste0(savedir, "featureplot/markers_genes.pdf"))
# plot_list
# dev.off()

DefaultAssay(CD8_Tcells_noCD4) <- "MAGIC_RNA"
CD8_Tcells_noCD4 <- AddModuleScore(CD8_Tcells_noCD4, features = markers, slot = "data")
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
CD8_subset <- readRDS("/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/raw_data/CD8_subset_2.RDS")
counts <- t(GetAssayData(CD8_subset, slot = "counts"))
counts_req <- counts[, colSums(counts) > 0] ### remove gene with 0 UMIs
cellnames <- sample(rownames(counts_req), size = 30000)
counts_req_2 <- counts_req[cellnames, ]

fit10 <- fit_topic_model(counts_req_2, k = 10)
fit20 <- fit_topic_model(counts_req_2, k = 20)
fit30 <- fit_topic_model(counts_req_2, k = 30)

saveRDS(fit10, paste0(savedir, "fit10_topic.RDS"))
saveRDS(fit20, paste0(savedir, "fit20_topic.RDS"))
saveRDS(fit30, paste0(savedir, "fit30_topic.RDS"))

cellnames <- rownames(fit10$L)
CD8_subset <- subset(CD8_Tcells, cells = cellnames)

fit10_L <- fit10$L[match(rownames(CD8_subset@meta.data), rownames(fit10$L)), ]
all(rownames(fit10_L) == rownames(CD8_subset@meta.data))

fit10_L_mat <- as.data.frame(as.matrix(fit10_L))

for (i in 1:10) {
    CD8_subset@meta.data[, paste0("k", i)] <- fit10_L_mat[, paste0("k", i)]
}

fit10_L <- fit10$L[match(rownames(CD8_subset@meta.data), rownames(fit10$L)), ]
all(rownames(fit10_L) == rownames(CD8_subset@meta.data))

fit10_L_mat <- as.data.frame(as.matrix(fit10_L))

for (i in 1:10) {
    CD8_subset@meta.data[, paste0("k", i)] <- fit10_L_mat[, paste0("k", i)]
}

df <- CD8_subset@meta.data[, paste0("k", 1:10)]

for (i in 1:nrow(df)) {
    df[i, "max_k1"] <- colnames(df)[(which.max(df[i, ]))]
}

all(rownames(df) == rownames(CD8_subset@meta.data))
CD8_subset@meta.data$max_k1 <- df$max_k1

k1_CT <- as.matrix(table(CD8_subset@meta.data$max_k1, CD8_subset@meta.data$subtype))

k1_CT2 <- k1_CT
total_sum <- colSums(k1_CT)
for (i in 1:ncol(k1_CT)) {
    k1_CT2[, i] <- k1_CT[, i] / total_sum[i]
}

### Sanity Check
colSums(k1_CT2)

library(ComplexHeatmap)
p <- Heatmap(as.matrix(k1_CT2), cluster_rows = FALSE, cluster_columns = FALSE)

pdf(paste0(savedir, "fit10_heatmap_subtype.pdf"))
p
dev.off()

source("/diazlab/data3/.abhinav/resources/all_scripts/R/express_cell_front1.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(proteins)) {
    p <- featureplot_front(CD8_subset, "k1",
        reduction = "umap", x = "umap_1", y = "umap_2", size = 0.1
    ) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "featureplot/topic_k1.pdf"))
p
dev.off()

metadata <- as.data.frame(CD8_subset@meta.data)
metadata_arrnge <- metadata[match(rownames(fit10$L), rownames(metadata)), ]
metadata_arrnge[is.na(metadata_arrnge$celltypes), "celltypes"] <- "Unknown"
metadata_arrnge$celltypes <- as.factor(metadata_arrnge$celltypes)

topic_colors <- c(
    "skyblue", "forestgreen", "darkmagenta", "dodgerblue",
    "gold", "darkorange", "red", "yellow", "magenta", "black"
)
p <- structure_plot(fit10,
    colors = topic_colors, topics = 1:10, gap = 25,
    grouping = metadata_arrnge$celltypes
)

pdf(paste0(savedir, "Topic_modelling/stucture.pdf"), width = 14, height = 5.5)
p
dev.off()

p <- structure_plot(fit20,
    topics = 1:20, gap = 25,
    grouping = metadata_arrnge$celltypes
)

pdf(paste0(savedir, "Topic_modelling/stucture_20.pdf"), width = 14, height = 5.5)
p
dev.off()

p <- structure_plot(fit20, topics = 1:20)

pdf(paste0(savedir, "Topic_modelling/stucture_20_nolabels.pdf"), width = 14, height = 5.5)
p
dev.off()

# Identifying clusters from the topic proportions
# One simple strategy that often works well is to first identify clusters in the topic proportions matrix, 𝐿. Here we use 𝑘-means to identify clusters, but other methods could be used.

set.seed(1)
pca <- prcomp(fit20$L)$x
clusters <- kmeans(pca, centers = 12, iter.max = 100)$cluster
summary(factor(clusters))

structure_plot(fit20, topics = 1:12, gap = 25, grouping = clusters) +
    theme(axis.text.x = element_blank())

pdf(paste0(savedir, "Topic_modelling/stucture_20_with_clusters.pdf"), width = 14, height = 5.5)
p
dev.off()

### fit30
p <- structure_plot(fit30, topics = 1:30, grouping = metadata_arrnge$celltypes)

pdf(paste0(savedir, "Topic_modelling/stucture_30.pdf"), width = 14, height = 5.5)
p
dev.off()

df <- as.data.frame(rbind(fit10$F["PDCD1", ], fit10$F["LAG3", ], fit10$F["HAVCR2", ], fit10$F["TOX", ]))

max_info <- function(row) {
    max_value <- max(row)
    max_col <- names(row)[which.max(row)]
    return(c(max_value, max_col))
}

result <- t(apply(df, 1, max_info))

# Specify the columns to check
columns_to_check <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10")

df <- CD8_subset@meta.data[, columns_to_check]

for (i in 1:nrow(df)) {
    df[i, "max_k1"] <- colnames(df)[(which.max(df[i, ]))]
}

all(rownames(df) == row.names(CD8_subset@meta.data))
CD8_subset@meta.data[, "max_k1"] <- df[, "max_k1"]

pdf(paste0(savedir, "UMAP/topic_fit10_clsuters.pdf"))
DimPlot(CD8_subset, group.by = "max_k1")
dev.off()

#### Topic modelling new
savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/Topic_modelling/"
fit10 <- readRDS(paste0(savedir, "fit10_topic.RDS"))
fit20 <- readRDS(paste0(savedir, "fit20_topic.RDS"))
fit30 <- readRDS(paste0(savedir, "fit30_topic.RDS"))

df <- as.data.frame(fit10$F)
df_2 <- df[grep("^RPL|^RPS|^MT\\.", rownames(df), invert = TRUE), ]

library(dplyr)
# Function to get top 50 IDs for a given column
get_top_ids <- function(data, column) {
    data %>%
        arrange(desc(!!sym(column))) %>%
        slice_head(n = 20) %>%
        select(ID = 1, Value = !!sym(column)) %>%
        mutate(Column_Name = column)
}

# Get top IDs for each value column
top_ids_list <- lapply(names(df_2)[-1], function(col) {
    get_top_ids(df_2, col)
})

# Combine results into a single dataframe
top_ids_combined <- do.call(rbind, top_ids_list)

genes <- rownames(top_ids_combined)
h <- DoHeatmap(CD8_Tcells, genes, group.by = "Celltypes2")
pdf(paste0("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/trajectory/heatmap_30000_celltypes.pdf"))
h
dev.off()

h <- DoHeatmap(CD8_Tcells_noCD4, genes, group.by = "subtype")
pdf(paste0("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/trajectory/heatmap_subtype_2.pdf"))
h
dev.off()

CD8_Tcells_noCD4 <- ScaleData(CD8_Tcells_noCD4, features = genes)

h <- DoHeatmap(CD8_Tcells_noCD4, genes, group.by = "Study")
pdf(paste0("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/Topic_modelling/heatmap_study_2.pdf"))
h
dev.off()

### Since Topic modelling showing the batch effect we are trying to perform the topic modelling with each Study and find the transition state
## Taking Yost SCC and BCC separately and performing Topic modeling, Li, and Azizi

### Azizi
Azizi_cells <- rownames(CD8_Tcells_noCD4@meta.data[grep("Azizi", CD8_Tcells_noCD4@meta.data$Study), ])
Azizi_obj <- subset(CD8_Tcells_noCD4, cells = Azizi_cells)
DefaultAssay(Azizi_obj) <- "RNA"
counts <- t(GetAssayData(Azizi_obj, slot = "counts"))
counts_req <- counts[, colSums(counts) > 0] ### remove gene with 0 UMIs
saveRDS(counts_req, "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/Topic_modelling/Azizi_TM/counts_Azizi.RDS")

### Yost BCC
CD8_Tcells_noCD4@meta.data$Study_subtype <- paste(CD8_Tcells_noCD4@meta.data$Study, CD8_Tcells_noCD4@meta.data$subtype, sep = "_")
Yost_BCC_cells <- rownames(CD8_Tcells_noCD4@meta.data[grep("Yost_BCC", CD8_Tcells_noCD4@meta.data$Study_subtype), ])
Yost_BCC_obj <- subset(CD8_Tcells_noCD4, cells = Yost_BCC_cells)
DefaultAssay(Yost_BCC_obj) <- "RNA"
counts <- t(GetAssayData(Yost_BCC_obj, slot = "counts"))
counts_req <- counts[, colSums(counts) > 0] ### remove gene with 0 UMIs
saveRDS(counts_req, "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/Topic_modelling/Yost_BCC_TM/counts_Yost_BCC.RDS")

### Yost SCC
# CD8_Tcells_noCD4@meta.data$Study_subtype <- paste(CD8_Tcells_noCD4@meta.data$Study, CD8_Tcells_noCD4@meta.data$subtype, sep ="_")
Yost_SCC_cells <- rownames(CD8_Tcells_noCD4@meta.data[grep("Yost_SCC", CD8_Tcells_noCD4@meta.data$Study_subtype), ])
Yost_SCC_obj <- subset(CD8_Tcells_noCD4, cells = Yost_SCC_cells)
DefaultAssay(Yost_SCC_obj) <- "RNA"
counts <- t(GetAssayData(Yost_SCC_obj, slot = "counts"))
counts_req <- counts[, colSums(counts) > 0] ### remove gene with 0 UMIs
saveRDS(counts_req, "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/Topic_modelling/Yost_SCC_TM/counts_Yost_SCC.RDS")

### Li CM
# CD8_Tcells_noCD4@meta.data$Study_subtype <- paste(CD8_Tcells_noCD4@meta.data$Study, CD8_Tcells_noCD4@meta.data$subtype, sep ="_")
Li_CM_cells <- rownames(CD8_Tcells_noCD4@meta.data[grep("Li_CM", CD8_Tcells_noCD4@meta.data$Study_subtype), ])
Li_CM_obj <- subset(CD8_Tcells_noCD4, cells = Li_CM_cells)
DefaultAssay(Li_CM_obj) <- "RNA"
counts <- t(GetAssayData(Li_CM_obj, slot = "counts"))
counts_req <- counts[, colSums(counts) > 0] ### remove gene with 0 UMIs
saveRDS(counts_req, "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/Topic_modelling/Li_CM_TM/counts_Li_CM.RDS")

#### Slimiing down the seurat object
CD8_Tcells_noCD4 <- readRDS("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/saveRDS_obj/CD8_Tcells_noCD4.RDS")
CD8_Tcells_noCD4@meta.data$patient_subtype <- paste(CD8_Tcells_noCD4@meta.data$patient, CD8_Tcells_noCD4@meta.data$subtype, sep = ":")

Total_cellnumber <- table(CD8_Tcells_noCD4@meta.data$patient_subtype) %>% as.data.frame()
Total_cellnumber_patient <- Total_cellnumber[Total_cellnumber$Freq > 50, ] ### minimum number of cells cutoff
Total_cellnumber_patient$patients <- gsub(":.*.", "", Total_cellnumber_patient$Var1)
Total_cellnumber_patient$cancer_type <- gsub(".*.:", "", Total_cellnumber_patient$Var1)

#### Extracting out the cells for the samples > 50 cells
cellnames <- rownames(CD8_Tcells_noCD4@meta.data[grep(paste("^", Total_cellnumber_patient$patients, "$", sep = "", collapse = "|"), CD8_Tcells_noCD4@meta.data$patient), ])
CD8_Tcells_noCD4_QC <- subset(CD8_Tcells_noCD4, cells = cellnames)

set.seed(7)
# downsampled_cells <- sample(rownames(CD8_Tcells_noCD4_QC@meta.data), size = 30000)
# write.table(downsampled_cells, paste0(savedir, "Table/downsample_cells.txt"), quote = F, row.names = F, col.names = F, sep = "\t")

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
downsampled_cells <- read.table(paste0(savedir, "Table/downsample_cells.txt"), header = FALSE)[, 1] ### For reproducibility
CD8_Tcells_noCD4_QC_down <- subset(CD8_Tcells_noCD4_QC, cells = downsampled_cells)

DefaultAssay(CD8_Tcells_noCD4_QC_down) <- "RNA"
# Identify genes that are greater then 3 UMIs
genes_to_keep <- rowSums(CD8_Tcells_noCD4_QC_down@assays$RNA@counts) >= 3
CD8_Tcells_noCD4_QC_down <- subset(CD8_Tcells_noCD4_QC_down, features = rownames(CD8_Tcells_noCD4_QC_down)[genes_to_keep])
CD8_Tcells_noCD4_QC_down_RNA <- DietSeurat(CD8_Tcells_noCD4_QC_down, assays = "RNA", dimreducs = c("pca", "umap"), graphs = "integrated_nn")


saveRDS(CD8_Tcells_noCD4_QC_down_RNA, "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/saveRDS_obj/CD8_Tcells_noCD4_QC_down_RNA.RDS")

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
DefaultAssay(CD8_Tcells_noCD4) <- "RNA"
sce <- as.SingleCellExperiment(CD8_Tcells_noCD4, assay = "RNA")

dir.create(paste0(savedir, "trajectory/UMAP"), showWarnings = FALSE)
pdf(paste0(savedir, "trajectory/UMAP/CD8_celltype_annotated.pdf"))
DimPlot(CD8_Tcells_noCD4, reduction = "umap", label = TRUE)
dev.off()

### celltypes
CD4_sce <- slingshot(sce,
    clusterLabels = "Celltypes2",
    reducedDim = "UMAP", approx_points = 100,
    omega = TRUE,
    omega_scale = 1.5
)

CD4_sce1 <- slingshot(sce,
    clusterLabels = "Celltypes2",
    reducedDim = "UMAP", approx_points = 100,
    start.clus = "SCM-like",
    omega = TRUE,
    omega_scale = 1.5
)

### clusters
CD4_sce <- slingshot(sce,
    clusterLabels = "Celltypes2",
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
for (i in 1:6) {
    embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
    embedded <- data.frame(embedded$s[embedded$ord, ])
    g <- plotUMAP(CD4_sce, colour_by = paste("slingPseudotime_", i, sep = ""))
    stopifnot(all(rownames(CD8_Tcells_noCD4@reductions$umap@cell.embeddings) == rownames(g$data)))
    data <- merge(CD8_Tcells_noCD4@reductions$umap@cell.embeddings, g$data, by = "row.names")
    colnames(data) <- c("cellname", "UMAP_1", "UMAP_2", "X", "Y", paste("Lineage", i, sep = ""))
    p <- ggplot(data, aes_string("UMAP_1", "UMAP_2", color = paste("Lineage", i, sep = ""))) +
        geom_point(size = 0.01) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
        geom_path(data = embedded, aes(x = umap_1, y = umap_2), color = "black", size = 1.2) +
        theme_bw()
    plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir, "trajectory/UMAP/CD8_Tcells_sce_UMAP_splitted_unbias_celltypes.pdf", sep = ""), width = 12.5, height = 9)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
    plot_list[[5]], plot_list[[6]],
    nrow = 2, ncol = 3
)
dev.off()

library(ArchR)
library(scater)
embedded_orig <- embedCurves(CD4_sce1, "UMAP")
rm(plot_list)
plot_list <- list()

for (i in 1:20) {
    embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
    embedded <- data.frame(embedded$s[embedded$ord, ])
    g <- plotUMAP(CD4_sce1, colour_by = paste("slingPseudotime_", i, sep = ""))
    stopifnot(all(rownames(CD8_Tcells_noCD4@reductions$umap@cell.embeddings) == rownames(g$data)))
    data <- merge(CD8_Tcells_noCD4@reductions$umap@cell.embeddings, g$data, by = "row.names")
    colnames(data) <- c("cellname", "UMAP_1", "UMAP_2", "X", "Y", paste("Lineage", i, sep = ""))
    p <- ggplot(data, aes_string("UMAP_1", "UMAP_2", color = paste("Lineage", i, sep = ""))) +
        geom_point(size = 0.01) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
        geom_path(data = embedded, aes(x = umap_1, y = umap_2), color = "black", size = 1.2) +
        theme_bw()
    plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir, "trajectory/UMAP/CD8_Tcells_sce_UMAP_splitted_start_naive_celltypes.pdf", sep = ""), width = 12.5, height = 9)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
    plot_list[[5]], plot_list[[6]],
    nrow = 2, ncol = 3
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
embedded_orig <- embedCurves(CD4_sce7, "UMAP")
rm(plot_list)
plot_list <- list()

for (i in 1:21) {
    embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
    embedded <- data.frame(embedded$s[embedded$ord, ])
    g <- plotUMAP(CD4_sce7, colour_by = paste("slingPseudotime_", i, sep = ""))
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

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
require(gridExtra)
pdf(paste(savedir, "trajectory/UMAP/CD8_Tcells_sce_UMAP_splitted_start_cluster6.pdf", sep = ""), width = 25, height = 20)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
    plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
    plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], plot_list[[13]], plot_list[[14]],
    plot_list[[15]], plot_list[[16]], plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]],
    nrow = 5, ncol = 5
)
dev.off()

#### TradeSeq workflow
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)

set.seed(7)
downsampled_cells <- sample(rownames(CD8_Tcells_noCD4_QC@meta.data), size = 30000)
write.table(downsampled_cells, paste0(savedir, "Table/downsample_cells.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
downsampled_cells <- read.table(paste0(savedir, "Table/downsample_cells.txt"), header = FALSE)[, 1] ### For reproducibility
CD8_Tcells_noCD4_QC_down <- subset(CD8_Tcells_noCD4_QC, cells = downsampled_cells)
pdf(savedir)

DefaultAssay(CD8_Tcells_noCD4_QC_down) <- "RNA"
sce <- as.SingleCellExperiment(CD8_Tcells_noCD4_QC_down, assay = "RNA")

CD4_sce1 <- slingshot(sce,
    clusterLabels = "seurat_clusters",
    reducedDim = "UMAP", approx_points = 100,
    start.clus = "6",
    omega = TRUE,
    omega_scale = 1.5
)

CD4_sce_naive <- slingshot(sce,
    clusterLabels = "Celltypes2",
    reducedDim = "UMAP", approx_points = 100,
    start.clus = "SCM-like",
    omega = TRUE,
    omega_scale = 1.5
)

counts <- as.matrix(CD8_Tcells_noCD4_QC_down@assays$RNA@counts)
CD4_sds <- SlingshotDataSet(CD4_sce1)
CD4_sds_naive <- SlingshotDataSet(CD4_sce_naive)

sce <- fitGAM(counts = counts, sds = CD4_sds)
sce_naive <- fitGAM(counts = counts, sds = CD4_sds_naive)

plotGeneCount(curves, filt_counts, clusters = clustering, models = sce)

set.seed(5)
icMat <- evaluateK(
    counts = counts, sds = CD4_sds, k = 3:10,
    nGenes = 500, verbose = T
)

#### Running TradeSeq on individual Studies i.e. Yost_BCC and Yost_SCC
CD8_Tcells_noCD4@meta.data[grep("Yost_BCC", CD8_Tcells_noCD4@meta.data$Study_subtype), ] %>% rownames() -> Yost_BCC_cellnames
Yost_BCC_obj <- subset(CD8_Tcells_noCD4, cells = Yost_BCC_cellnames)
DefaultAssay(Yost_BCC_obj) <- "RNA"
Yost_BCC_obj_Diet <- DietSeurat(Yost_BCC_obj, assays = "RNA", dimreducs = c("pca", "umap"), graphs = "integrated_nn")
Yost_BCC_obj_Diet_feat <- subset(Yost_BCC_obj_Diet, features = rownames(Yost_BCC_obj_Diet)[Matrix::rowSums(Yost_BCC_obj_Diet@assays$RNA@counts > 0) >= 3])
Yost_BCC_obj_Diet_feat <- NormalizeData(Yost_BCC_obj_Diet_feat)

saveRDS(
    Yost_BCC_obj_Diet_feat,
    "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/saveRDS_obj/Yost_BCC_obj_Diet_feat.RDS"
)

CD8_Tcells_noCD4@meta.data[grep("Yost_SCC", CD8_Tcells_noCD4@meta.data$Study_subtype), ] %>% rownames() -> Yost_SCC_cellnames
Yost_SCC_obj <- subset(CD8_Tcells_noCD4, cells = Yost_SCC_cellnames)
DefaultAssay(Yost_SCC_obj) <- "RNA"
Yost_SCC_obj_Diet <- DietSeurat(Yost_SCC_obj, assays = "RNA", dimreducs = c("pca", "umap"), graphs = "integrated_nn")
Yost_SCC_obj_Diet_feat <- subset(Yost_SCC_obj_Diet, features = rownames(Yost_SCC_obj_Diet)[Matrix::rowSums(Yost_SCC_obj_Diet@assays$RNA@counts > 0) >= 3])
Yost_SCC_obj_Diet_feat <- NormalizeData(Yost_SCC_obj_Diet_feat)

saveRDS(
    Yost_SCC_obj_Diet_feat,
    "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/saveRDS_obj/Yost_SCC_obj_Diet_feat.RDS"
)

clus_CT <- read.table("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/Table/cluster_CT.txt", header = TRUE, sep = "\t")
patterns <- clus_CT$Clusters
replacements <- clus_CT$Celltypes

CD8_Tcells@meta.data$celltypes <- CD8_Tcells@meta.data$seurat_clusters
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD8_Tcells@meta.data$celltypes <- str_replace_all(CD8_Tcells@meta.data$celltypes, pattern, replacements[i])
}

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
pdf(paste0(savedir, "UMAP/celltypes.pdf"), width = 7, height = 5.5)
DimPlot(CD8_Tcells, reduction = "umap", group.by = "celltypes", label = TRUE)
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
pdf(paste0(savedir, "UMAP/Celltypes2.pdf"), width = 7, height = 5.5)
DimPlot(CD8_Tcells_noCD4, reduction = "umap", group.by = "Celltypes2", label = TRUE)
dev.off()

#### Performing reference mapping
# update object and then convert Assay to SCTAssay
CD8_Ag <- readRDS("/diazlab/data3/.abhinav/.immune/cancer_combined/project/resources/GSE275633_CD8_Antigen_BEAM_T.RDS")
CD8_Ag <- UpdateSeuratObject(CD8_Ag)
CD8_Ag[["integrated"]] <- as(object = CD8_Ag[["integrated"]], Class = "SCTAssay")

CD8_Ag <- RunUMAP(CD8_Ag,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_",
    return.model = TRUE
)

DefaultAssay(CD8_Ag) <- "integrated"
sample.anchors <- FindTransferAnchors(
    reference = CD8_Ag,
    query = CD8_Tcells,
    dims = 1:30,
    reference.reduction = "pca",
    normalization.method = "SCT",
    recompute.residuals = FALSE
)

CD8_Tcells <- MapQuery(
    anchorset = sample.anchors,
    query = CD8_Tcells,
    reference = CD8_Ag,
    refdata = list(
        celltype.l1 = "celltypes",
        predicted_ADT = "ADT"
    ),
    reference.reduction = "pca",
    reduction.model = "wnn.umap"
)

### Got the predicted ADT data also
proteins <- rownames(CD8_Tcells@assays$predicted_ADT$data)
p <- DotPlot(CD8_Tcells, proteins, assay = "predicted_ADT") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_celltypes_integrated_protein.pdf"), height = 11, width = 19)
print(p)
dev.off()

# CD45RA_RO <- c("CD45RA-protein", "CD45RO-protein")
source("/diazlab/data3/.abhinav/resources/all_scripts/R/express_cell_front1.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(proteins)) {
    p <- featureplot_front(CD8_Tcells, proteins[i],
        reduction = "umap", x = "umap_1", y = "umap_2", size = 0.1
    ) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "featureplot/All_proteins_Ag.pdf"))
plot_list
dev.off()

source("/diazlab/data3/.abhinav/resources/all_scripts/R/express_cell_front1.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(proteins)) {
    p <- featureplot_front(CD8_Tcells, proteins[i],
        reduction = "umap", x = "umap_1", y = "umap_2", size = 0.1
    ) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "featureplot/All_proteins_Bulk.pdf"))
plot_list
dev.off()

#### Performing the imputation
DefaultAssay(CD8_Tcells) <- "predicted_ADT"
# CD8_Tcells <- NormalizeData(CD8_Tcells)
library(reticulate)
library(Rmagic)
use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check
CD8_Tcells <- magic(CD8_Tcells, npca = 10) ## imputing the RNA data as for RNA PCs are 20
DefaultAssay(CD8_Tcells) <- "MAGIC_predicted_ADT"

p <- DotPlot(CD8_Tcells, proteins, assay = "MAGIC_predicted_ADT") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_celltypes_integrated_protein_imputed_Ag.pdf"), height = 11, width = 19)
print(p)
dev.off()

source("/diazlab/data3/.abhinav/resources/all_scripts/R/express_cell_front1.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(proteins)) {
    p <- featureplot_front(CD8_Tcells, proteins[i],
        reduction = "umap", x = "umap_1", y = "umap_2", size = 0.1
    ) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "featureplot/All_proteins_Ag_imputed.pdf"))
plot_list
dev.off()

saveRDS(CD8_subset_2, "/diazlab/data3/.abhinav/.industry/Genentech/panel_discussion/raw_data/CD8_subset_2_impute.RDS")


pdf(paste0(savedir, "UMAP/ADT_predicted_Ag.pdf"))
DimPlot(CD8_Tcells, reduction = "umap", group.by = "predicted.celltype.l1")
dev.off()

pdf(paste0(savedir, "UMAP/ADT_predicted_Ag_split.pdf"), width = 12, height = 12)
DimPlot(CD8_Tcells, reduction = "umap", group.by = "predicted.celltype.l1", split.by = "predicted.celltype.l1", ncol = 3)
dev.off()

### Performing the same thing with the CD8 Bulk
CD8_Bulk <- readRDS("/diazlab/data3/.abhinav/.immune/cancer_combined/project/resources/GSE275633_CD8_bulk_BEAM_T.RDS")
CD8_Bulk <- UpdateSeuratObject(CD8_Bulk)
CD8_Bulk[["integrated"]] <- as(object = CD8_Bulk[["integrated"]], Class = "SCTAssay")

CD8_Bulk <- RunUMAP(CD8_Bulk, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = TRUE)

DefaultAssay(CD8_Tcells) <- "integrated"
DefaultAssay(CD8_Bulk) <- "integrated"

sample.anchors <- FindTransferAnchors(
    reference = CD8_Bulk,
    query = CD8_Tcells,
    dims = 1:30,
    reference.reduction = "pca",
    normalization.method = "SCT",
    recompute.residuals = TRUE
)

CD8_Tcells <- MapQuery(
    anchorset = sample.anchors,
    query = CD8_Tcells,
    reference = CD8_Bulk,
    refdata = list(
        celltype.l1 = "celltypes",
        predicted_ADT = "ADT"
    ),
    reference.reduction = "pca",
    reduction.model = "wnn.umap"
)

proteins <- rownames(CD8_Tcells@assays$predicted_ADT$data)
p <- DotPlot(CD8_Tcells, proteins, assay = "predicted_ADT") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_celltypes_integrated_protein_bulk.pdf"), height = 11, width = 19)
print(p)
dev.off()


#### Removing CD4 T cells FOXP3 and Th17 and making a UMAP again
cellnames_noCD4 <- rownames(CD8_Tcells@meta.data[!is.na((CD8_Tcells@meta.data$celltypes)), ])
CD8_Tcells_noCD4 <- subset(CD8_Tcells, cells = cellnames_noCD4)

#### Downloaded the metadata with the publication version
new_metadata <- read.csv("/diazlab/data3/.abhinav/.immune/cancer_combined/project/raw_data/TICAtlas_metadata_new.csv")
new_metadata$new_cellnames <- gsub(paste("^", unique(new_metadata$patient), "_", collapse = "|", sep = ""), "", rownames(new_metadata))

match(rownames(CD8_Tcells_noCD4@meta.data), new_metadata$new_cellnames)
CD8_Tcells_noCD4@meta.data$paper_celltypes <- new_metadata[match(rownames(CD8_Tcells_noCD4@meta.data), new_metadata$new_cellnames), "lv2_annot"]

cells <- rownames(CD8_Tcells_noCD4@meta.data[!is.na(CD8_Tcells_noCD4@meta.data$paper_celltypes), ])
CD8_Tcells_noCD4_subset <- subset(CD8_Tcells_noCD4, cells = cells)

pdf(paste0(savedir, "UMAP/paper_celltypes.pdf"))
DimPlot(CD8_Tcells_noCD4_subset, group.by = "paper_celltypes", reduction = "umap")
dev.off()


### Got the predicted ADT data also
proteins <- rownames(CD8_Tcells@assays$predicted_ADT$data)
p <- DotPlot(CD8_Tcells, proteins, assay = "predicted_ADT") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_celltypes_integrated_protein.pdf"), height = 11, width = 19)
print(p)
dev.off()

#### Downstream Analysis
#####
### Checking how many samples would left if I removed samples with less than 50 cells

### Adding the information of the metadata to the object
sample_info <- read.table("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/Table/sample_info.txt",
    header = TRUE, sep = "\t"
)

patterns <- sample_info$SampleID
replacements <- sample_info$Study
CD8_Tcells_noCD4@meta.data$Study <- CD8_Tcells_noCD4@meta.data$patient
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD8_Tcells_noCD4@meta.data$Study <- str_replace_all(CD8_Tcells_noCD4@meta.data$Study, pattern, replacements[i])
}

CD8_Tcells_noCD4@meta.data$patient_subtype <- paste(CD8_Tcells_noCD4@meta.data$patient, CD8_Tcells_noCD4@meta.data$subtype, sep = ":")

Total_cellnumber <- table(CD8_Tcells_noCD4@meta.data$patient_subtype) %>% as.data.frame()
Total_cellnumber_patient <- Total_cellnumber[Total_cellnumber$Freq > 50, ] ### minimum number of cells cutoff
Total_cellnumber_patient$patients <- gsub(":.*.", "", Total_cellnumber_patient$Var1)
Total_cellnumber_patient$cancer_type <- gsub(".*.:", "", Total_cellnumber_patient$Var1)

#### Extracting out the cells for the samples > 50 cells
cellnames <- rownames(CD8_Tcells_noCD4@meta.data[grep(paste("^", Total_cellnumber_patient$patients, "$", sep = "", collapse = "|"), CD8_Tcells_noCD4@meta.data$patient), ])
CD8_Tcells_noCD4_QC <- subset(CD8_Tcells_noCD4, cells = cellnames)

## Quantitative Differences
savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
celltypes_subtype <- table(CD8_Tcells_noCD4_QC@meta.data$subtype, CD8_Tcells_noCD4_QC@meta.data$celltypes2)
write.table(celltypes_subtype, paste0(savedir, "Table/celltypes_subtype_QCed_2.txt"),
    row.names = T, col.names = T, sep = "\t", quote = F
)

celltypes_subtype <- read.table(paste0(savedir, "Table/celltypes_subtype.txt"), header = TRUE, sep = "\t")
n_cells <- celltypes_subtype
n_cells_sum <- as.vector(rowSums(celltypes_subtype))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
    df[j, ] <- (n_cells[j, ] / n_cells_sum[j]) * 100
}

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("Cancertype", "celltype", "percentage")
library(ggplot2)
# create a dataset
# Stacked + percent

p <- ggplot(df_melted, aes(fill = celltype, y = percentage, x = Cancertype)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white"
    )) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(paste0(savedir, "Table//marker_celltype_Qced_2.pdf"), width = 9, height = 6.5)
p
dev.off()


### Performing based on each individuals
library(ggpubr)
library(rstatix)
library(stringr)

patient_celltypes <- table(CD8_Tcells_noCD4_QC@meta.data$patient_subtype, CD8_Tcells_noCD4_QC@meta.data$celltypes2)
write.table(patient_celltypes, paste0(savedir, "Table/patient_celltypes_2.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
patient_celltypes <- read.table(paste0(savedir, "Table/patient_celltypes_2.txt"), header = TRUE, sep = "\t")
n_cells <- patient_celltypes
n_cells_sum <- as.vector(rowSums(patient_celltypes))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
    df[j, ] <- (n_cells[j, ] / n_cells_sum[j]) * 100
}

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
df_melted$patient_id <- gsub(":.*.", "", df_melted$sample)
df_melted$cancer_type <- gsub(".*.:", "", df_melted$sample)

colnames(df_melted)[1:3] <- c("sampleid:cacertype", "celltype", "percentage")
write.table(df_melted, paste0(savedir, "Table/df_melted_2.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

### Checking whether to use parametric or non-parametric test
# Load necessary libraries
normality_results <- aggregate(percentage ~ celltype, data = df_melted, FUN = function(x) shapiro.test(x)$p.value)

# Print Shapiro-Wilk test results
print("Shapiro-Wilk Test P-values by Cell Type:")
print(normality_results)

write.table(normality_results, paste0(savedir, "Table/normality_test_shapiro_wilk_2.txt"))

# Visualize the distribution of percentages by cell type
ggplot(df_melted, aes(x = percentage, fill = celltype)) +
    geom_histogram(binwidth = 5, position = "dodge") +
    labs(title = "Distribution of Percentages by Cell Type", x = "Percentage", y = "Frequency")

# Q-Q plot for one of the cell types (e.g., SCM.like)
ggplot(df_melted[df_melted$celltype == "SCM.like", ], aes(sample = percentage)) +
    geom_qq() +
    geom_qq_line() +
    labs(title = "Q-Q Plot for SCM.like Cell Type")


# Perform Kruskal-Wallis test to assess if there are any differences
kruskal_results <- df_melted %>%
    group_by(celltype) %>%
    kruskal_test(percentage ~ cancer_type) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")


# Perform Dunn’s test for pairwise comparisons within each cell type
stat.test <- df_melted %>%
    group_by(celltype) %>%
    dunn_test(percentage ~ cancer_type, p.adjust.method = "bonferroni")

write.table(stat.test, paste(savedir, "Table/kruskal_dunn_test_2.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
write.table(kruskal_results, paste(savedir, "Table/kruskal_results_2.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)

# Filter to retain only significant comparisons (raw p-value < 0.05)
stat.test2 <- stat.test %>% filter(p.adj < 0.1)

# Add x and y positions for p-value annotations
stat.test2 <- stat.test2 %>%
    add_xy_position(x = "celltype", dodge = 0.8)

# Create the faceted boxplot, split by cell type

bxp <- ggboxplot(
    df_melted,
    x = "cancer_type", y = "percentage",
    color = "cancer_type",
    palette = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white"
    )
) +
    facet_wrap(~celltype, scales = "free_y", labeller = "label_both", nrow = 2) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        strip.placement = "outside", # Moves facet labels outside
        panel.spacing.y = unit(1, "lines"), # Adds space between rows
        axis.ticks.length.y = unit(0.25, "cm"), # Ensure all facets show ticks
        strip.text = element_text(size = 16, face = "bold"),
        # Legend customization
        legend.text = element_text(size = 14),   # Legend text size
        legend.title = element_text(size = 16, face = "bold"),  # Legend title size
        legend.key.size = unit(1.5, "lines")     # Increases legend key size
    )

bxp2 <- ggboxplot(
    df_melted,
    x = "celltype", y = "percentage",
    color = "celltype",
    palette = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white"
    )
) +
    facet_wrap(~cancer_type, scales = "free_y", labeller = "label_both", nrow = 2) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside", # Moves facet labels outside
        panel.spacing.y = unit(1, "lines"), # Adds space between rows
        axis.ticks.length.y = unit(0.25, "cm"), # Ensure all facets show ticks
    )


# # Create boxplots
# bxp <- ggboxplot(
#     df_melted,
#     x = "celltype", y = "percentage",
#     color = "cancer_type", palette = c(
#         "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
#         "#8c564b", "#e377c2", "#b5bd61", "#17becf",
#         "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
#         "yellow", "black", "blanchedalmond", "blue", "white"
#     )
# ) +
#     facet_wrap(~celltype, scales = "free_y") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add significant p-values to each facet
bxp_with_pvals <- bxp + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0.01
) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "Significant Comparisons of Cancer Types within Cell Types")

pdf(paste(savedir, "Table/krushkal_wallis_test_between_cancer_type.pdf", sep = ""), width = 15, height = 20)
bxp_with_pvals
dev.off()

pdf(paste(savedir, "Table/krushkal_wallis_test_between_cancer_type_nopvalues_2.pdf", sep = ""), width = 25, height = 10)
bxp
dev.off()

pdf(paste(savedir, "Table/krushkal_wallis_test_between_celltype_nopvalues.pdf", sep = ""), width = 25, height = 10)
bxp2
dev.off()

# Add p-values to the plot
stat.test <- stat.test %>%
    add_xy_position(x = "cancer_type")

bxp_with_pvals <- bxp + stat_pvalue_manual(
    stat.test,
    label = "p.adj", tip.length = 0.01
) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "Comparison of Cancer Types within Cell Types")

# Save the plot as PDF
output_file <- "/mnt/data/kruskal_wallis_boxplot.pdf"
ggsave(output_file, plot = bxp_with_pvals, width = 10, height = 6)

# Display the plot
print(bxp_with_pvals)


bxp2 <- bxp + geom_dotplot(
    aes(fill = timepts, color = timepts),
    trim = FALSE,
    binaxis = "y", stackdir = "center", dotsize = 0.15,
    position = position_dodge(0.8)
) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))

# Add p-values onto the box plots
stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) + theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.grid.major = element_line(colour = "white")
)

# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir, "Table/t_test_paired_primary_recurrent_samples_bxplot_two_sided_celltype2.pdf", sep = ""), width = 8, height = 6)
bxp4
dev.off()


# Performing Qualittative differences using single cell RNAseq in each celltype
CD8_Tcells_noCD4_QC <- NormalizeData(CD8_Tcells_noCD4_QC)

subtypes <- unique(CD8_Tcells_noCD4_QC@meta.data$subtype)
celltypes <- levels(CD8_Tcells_noCD4_QC@meta.data$Celltypes2)
dir.create(paste(savedir, "scdifferential/", sep = ""), recursive = TRUE, showWarnings = FALSE)

for (k in 9:length(celltypes)) {
    print(k)
    celltype_cells <- rownames(CD8_Tcells_noCD4_QC@meta.data[grep(paste0("^", celltypes[k], "$"), CD8_Tcells_noCD4_QC@meta.data$Celltypes2), ])
    obj_subset <- subset(CD8_Tcells_noCD4_QC, cell = celltype_cells)
    stopifnot(unique(obj_subset@meta.data$Celltypes2) == celltypes[k])
    obj_subset <- NormalizeData(obj_subset)
    dir.create(paste(savedir, "scdifferential/", celltypes[k], sep = ""), recursive = TRUE, showWarnings = FALSE)
    for (i in 1:length(subtypes)) {
        dir.create(paste(savedir, "scdifferential/", celltypes[k], "/", subtypes[i], sep = ""), recursive = TRUE, showWarnings = FALSE)
        # for (j in i:(length(subtypes) - 1)) {
        for (j in 1:length(subtypes)) {
            try({
                if (i != j) {
                    scRNA_differential <- FindMarkers(obj_subset, ident.1 = subtypes[i], ident.2 = subtypes[j], group.by = "subtype")
                    write.table(scRNA_differential,
                        paste0(savedir, "scdifferential/", celltypes[k], "/", subtypes[i], "/", celltypes[k], "_", subtypes[i], "_vs_", subtypes[j], ".txt"),
                        quote = F, sep = "\t"
                    )
                }
            })
        }
    }
}

# Performing for all the cells
for (i in 1:length(subtypes)) {
    dir.create(paste(savedir, "scdifferential/all/", subtypes[i], sep = ""), recursive = TRUE, showWarnings = FALSE)

    for (j in 1:length(subtypes)) {
        try({
            if (i != j) {
                scRNA_differential <- FindMarkers(CD8_Tcells_noCD4_QC, ident.1 = subtypes[i], ident.2 = subtypes[j], group.by = "subtype")
                write.table(scRNA_differential,
                    paste0(savedir, "scdifferential/all/", subtypes[i], "/all_", subtypes[i], "_vs_", subtypes[j], ".txt"),
                    quote = F, sep = "\t"
                )
            }
        })
    }
}

#  awk '{if ($4 > 0.1 && $5 > 0.1) print $0}' SCM-like_NSCLC_vs_BC.txt | grep -v "RP"

#### Finding the common genes
### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
filepath <- list.files("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/scdifferential/all/BC/", pattern = ".txt", full.names = TRUE)
filename <- gsub(".txt", "", basename(filepath))
# filename2 <- gsub(".*._","",filename)

for (i in 1:length(filepath)) {
    PTmarkers <- read.table(filepath[i], header = T, sep = "\t")
    PTmarkers <- PTmarkers[PTmarkers$pct.1 > 0.05 & PTmarkers$pct.2 > 0.05, ]
    highmarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC > 0.5, ]) ## Putting the cutoff < 0.05 and avg_log2FC > 1
    assign(paste(filename[i], "_sig_high", sep = ""), highmarkers)
    lowmarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC < -0.5, ])
    assign(paste(filename[i], "_sig_low", sep = ""), lowmarkers)
}

filename2 <- ls(pattern = "_sig_high")

H <- list(
    "BC_BCC_high" = get(filename2[1]),
    "BC_CM_high" = get(filename2[2]),
    "BC_CRC_high" = get(filename2[3]),
    "BC_EA_high" = get(filename2[4]),
    "BC_HCC_high" = get(filename2[5]),
    "BC_ICC_high" = get(filename2[6]),
    "BC_NSCLC_high" = get(filename2[7]),
    "BC_PDAC_high" = get(filename2[8]),
    "BC_RCC_high" = get(filename2[9]),
    "BC_SCC_high" = get(filename2[10]),
    "BC_UM_high" = get(filename2[11])
)

# create customised venn diagram
# library(ggvenn)
# a <- ggvenn(H, show_elements = FALSE, stroke_color = "Black", stroke_linetype = "solid")
# savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/"
# dir.create(paste(savedir, "Table/", sep = ""), showWarnings = FALSE)
# pdf(paste(savedir, "Table/MAST_recurrent_high_venny.pdf", sep = ""))
# a
# dev.off()

### Making an upset plot
n <- max(
    length(get(filename2[1])), length(get(filename2[2])),
    length(get(filename2[3])), length(get(filename2[4])),
    length(get(filename2[5])), length(get(filename2[6])),
    length(get(filename2[7])), length(get(filename2[8])),
    length(get(filename2[9])), length(get(filename2[10])),
    length(get(filename2[11]))
)

length(all_BC_vs_BCC_sig_high) <- n
length(all_BC_vs_CM_sig_high) <- n
length(all_BC_vs_CRC_sig_high) <- n
length(all_BC_vs_EA_sig_high) <- n
length(all_BC_vs_HCC_sig_high) <- n
length(all_BC_vs_ICC_sig_high) <- n
length(all_BC_vs_NSCLC_sig_high) <- n
length(all_BC_vs_PDAC_sig_high) <- n
length(all_BC_vs_RCC_sig_high) <- n
length(all_BC_vs_SCC_sig_high) <- n
length(all_BC_vs_UM_sig_high) <- n

df <- as.data.frame(cbind(
    (all_BC_vs_BCC_sig_high), (all_BC_vs_CM_sig_high),
    (all_BC_vs_CRC_sig_high), (all_BC_vs_EA_sig_high),
    (all_BC_vs_HCC_sig_high), (all_BC_vs_ICC_sig_high),
    (all_BC_vs_NSCLC_sig_high), (all_BC_vs_PDAC_sig_high),
    (all_BC_vs_RCC_sig_high), (all_BC_vs_SCC_sig_high),
    (all_BC_vs_UM_sig_high)
))

library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- c(
    ("all_BC_vs_BCC_sig_high"), ("all_BC_vs_CM_sig_high"),
    ("all_BC_vs_CRC_sig_high"), ("all_BC_vs_EA_sig_high"),
    ("all_BC_vs_HCC_sig_high"), ("all_BC_vs_ICC_sig_high"),
    ("all_BC_vs_NSCLC_sig_high"), ("all_BC_vs_PDAC_sig_high"),
    ("all_BC_vs_RCC_sig_high"), ("all_BC_vs_SCC_sig_high"),
    ("all_BC_vs_UM_sig_high")
)

pdf(paste(savedir, "scdifferential/all/BC/combined_all_upset_plot.pdf", sep = ""), width = 20, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
    "BC_BCC_high" = get(filename2[1]),
    "BC_CM_high" = get(filename2[2]),
    "BC_CRC_high" = get(filename2[3]),
    "BC_EA_high" = get(filename2[4]),
    "BC_HCC_high" = get(filename2[5]),
    "BC_ICC_high" = get(filename2[6]),
    "BC_NSCLC_high" = get(filename2[7]),
    "BC_PDAC_high" = get(filename2[8]),
    "BC_RCC_high" = get(filename2[9]),
    "BC_SCC_high" = get(filename2[10]),
    "BC_UM_high" = get(filename2[11])
)

df2 <- data.frame(gene = unique(unlist(list_filter)))

library(dplyr)
df1 <- lapply(list_filter, function(x) {
    data.frame(gene = x)
}) %>%
    bind_rows(.id = "path")

df_int <- lapply(df2$gene, function(x) {
    # pull the name of the intersections
    intersection <- df1 %>%
        dplyr::filter(gene == x) %>%
        arrange(path) %>%
        pull("path") %>%
        paste0(collapse = "|")

    # build the dataframe
    data.frame(gene = x, int = intersection)
}) %>%
    bind_rows()

library(dplyr)
final <- df_int %>%
    group_by(int) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n))

df_int$int2 <- gsub("\\|", "__", df_int$int)

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "_high")
all_sample <- df_int[(df_int$sample_num > 6), 1]
morethan_3 <- df_int[(df_int$sample_num >= 6), 1]

write.table(df_int, paste(savedir, "Table/MAST_recurrent_gene_intersect.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "Table/MAST_recurrent_intersected_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

write.table(all_sample, paste(savedir, "Table/morethan_6.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(morethan_3, paste(savedir, "Table/morethan_orequal_3.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")


file_full_path <- list.files("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/scdifferential/all/", full.names = T)
cancername <- basename(file_full_path)
source("/diazlab/data3/.abhinav/.immune/cancer_combined/project/CD8_Genen/code/upset_plot.R")
for (i in 1:length(file_full_path)) {
    df_int <- upset_R(file_full_path[i])
    assign(cancername[i], df_int)
}

### Identified interesting genes
genes <- unique(read.table("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/scdifferential/all/combined/differential_genes")[, 1])

p <- DotPlot(CD8_Tcells_noCD4, genes, assay = "RNA", group.by = "subtype") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/"
pdf(paste0(savedir, "dotplot/differential_genes_subtypes.pdf"), width = 8, height = 20)
p
dev.off()

#### Performing with the linent cutoff for logFC == 0
library(UpSetR)
library(dplyr)
filename <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table",
    pattern = "_recurrent_vs_primary_markers_MAST_nocutoff.txt", full.names = TRUE
) %>%
    grep("_Endo_|_HSP_", ., value = TRUE, invert = TRUE)
markersname <- gsub("_recurrent_vs_primary_markers_MAST_nocutoff.txt", "", basename(filename))

# extracting the primary and recurrent markers
for (i in 1:length(markersname)) {
    PTmarkers <- read.table(filename[i], sep = "\t", header = TRUE)
    recurrentmarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC > 0, ]) ## Putting the cutoff < 0.05 and avg_log2FC > 1
    assign(paste(markersname[i], "_markers_MAST_recurrent", sep = ""), recurrentmarkers)
    primarymarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC < 0, ])
    assign(paste(markersname[i], "_markers_MAST_primary", sep = ""), primarymarkers)
}

recurrentname <- ls(pattern = "_MAST_recurrent")

H <- list(
    "PT_7WYPEC3Q_MAST_recurrent" = PT_7WYPEC3Q_markers_MAST_recurrent,
    "PT_9S6WMQ92_MAST_recurrent" = PT_9S6WMQ92_markers_MAST_recurrent,
    "PT_H3WWDMW9_MAST_recurrent" = PT_H3WWDMW9_markers_MAST_recurrent,
    "PT_XA98HG1C_MAST_recurrent" = PT_XA98HG1C_markers_MAST_recurrent
)

# create customised venn diagram
library(ggvenn)
a <- ggvenn(H, show_elements = FALSE, stroke_color = "Black", stroke_linetype = "solid")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/"
dir.create(paste(savedir, "Table/", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "Table/MAST_recurrent_high_venny_lfc_0.pdf", sep = ""))
a
dev.off()

### Making an upset plot
n <- max(
    length(PT_7WYPEC3Q_markers_MAST_recurrent), length(PT_9S6WMQ92_markers_MAST_recurrent),
    length(PT_H3WWDMW9_markers_MAST_recurrent), length(PT_XA98HG1C_markers_MAST_recurrent)
)

length(PT_7WYPEC3Q_markers_MAST_recurrent) <- n
length(PT_9S6WMQ92_markers_MAST_recurrent) <- n
length(PT_H3WWDMW9_markers_MAST_recurrent) <- n
length(PT_XA98HG1C_markers_MAST_recurrent) <- n

df <- as.data.frame(cbind(
    (PT_7WYPEC3Q_markers_MAST_recurrent), (PT_9S6WMQ92_markers_MAST_recurrent),
    (PT_H3WWDMW9_markers_MAST_recurrent), (PT_XA98HG1C_markers_MAST_recurrent)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- c("PT_7WYPEC3Q", "PT_9S6WMQ92", "PT_H3WWDMW9", "PT_XA98HG1C")

pdf(paste(savedir, "Table/MAST_recurrent_markers_upset_plot_lfc_0.pdf", sep = ""), width = 18, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
    "PT_7WYPEC3Q" = PT_7WYPEC3Q_markers_MAST_recurrent,
    "PT_9S6WMQ92" = PT_9S6WMQ92_markers_MAST_recurrent,
    "PT_H3WWDMW9" = PT_H3WWDMW9_markers_MAST_recurrent,
    "PT_XA98HG1C" = PT_XA98HG1C_markers_MAST_recurrent
)

df2 <- data.frame(gene = unique(unlist(list_filter)))

df1 <- lapply(list_filter, function(x) {
    data.frame(gene = x)
}) %>%
    bind_rows(.id = "path")

df_int <- lapply(df2$gene, function(x) {
    # pull the name of the intersections
    intersection <- df1 %>%
        dplyr::filter(gene == x) %>%
        arrange(path) %>%
        pull("path") %>%
        paste0(collapse = "|")

    # build the dataframe
    data.frame(gene = x, int = intersection)
}) %>%
    bind_rows()

library(dplyr)
library(stringr)
final <- df_int %>%
    group_by(int) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n))

df_int$int2 <- gsub("\\|", "__", df_int$int)

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
all_sample <- df_int[(df_int$sample_num > 3), 1]
morethan_3 <- df_int[(df_int$sample_num >= 3), 1]

write.table(df_int, paste(savedir, "Table/MAST_recurrent_gene_intersect_lfc_0.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "Table/MAST_recurrent_intersected_gene_number_lfc_0.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

write.table(all_sample, paste(savedir, "Table/MAST_recurrent_all_samples_lfc_0.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(morethan_3, paste(savedir, "Table/MAST_recurrent_morethan_3_lfc_0.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")


#### TCR data for trajectory
BC09_expr <- Read10X("/diazlab/data3/.abhinav/.immune/cancer_combined/project/raw_data/TCR_study/Azizi_study/seurat_obj/GSM3148580_BC09_TUMOR1/")
BC09_object <- CreateSeuratObject(counts = BC09_expr)
BC09_tum1_TCR <- read.csv("/diazlab/data3/.abhinav/.immune/cancer_combined/project/raw_data/TCR_study/Azizi_study/seurat_obj/GSM3148580_BC09_TUMOR1/GSM3148580_BC09_TUMOR1_filtered_contig_annotations.csv")

df2 <- df2[match(rownames(df1), rownames(df2), nomatch = 0), ]
df1 <- df1[match(rownames(df2), rownames(df1), nomatch = 0), ]
stopifnot(rownames(df2) == rownames(df1))
for (i in 1:ncol(df1)) {
    for (j in 1:ncol(df2)) {
        if (all(df1[, i] == df2[, j]) == TRUE) {
            print(paste(i, "_", j))
        }
    }
}

matches <- outer(
    1:ncol(df1), 1:ncol(df2),
    Vectorize(function(i, j) all(df1[, i] == df2[, j]))
)
which(matches, arr.ind = TRUE) # Gives the matching pairs


Azizi_obj <- NormalizeData(Azizi_obj)
Azizi_obj <- FindVariableFeatures(Azizi_obj, nfeatures = 3000)

BC09_object <- NormalizeData(BC09_object)
BC09_object <- FindVariableFeatures(BC09_object, nfeatures = 3000)

anchors <- FindTransferAnchors(
    reference = Azizi_obj,
    query = BC09_object,
    #   normalization.method = "SCT",
    reference.reduction = "pca",
    dims = 1:30
)

BC09_object <- MapQuery(
    anchorset = anchors,
    query = BC09_object,
    reference = Azizi_obj,
    reference.reduction = "pca",
    reduction.model = "umap"
)

#### Aziz Data does not work. Using Yost TCR which BCC and SCC
BCC_TCR <- read.table("/diazlab/data3/.abhinav/.immune/cancer_combined/project/raw_data/TCR_study/Yost_BCC_SCC/GSE123813_bcc_tcr.txt", header = TRUE)
SCC_TCR <- read.table("/diazlab/data3/.abhinav/.immune/cancer_combined/project/raw_data/TCR_study/Yost_BCC_SCC/GSE123813_scc_tcr.txt", header = TRUE)
BCC_SCC_TCR <- rbind(BCC_TCR, SCC_TCR)

### Since it has TRA1, TRB1
### Making it for TRA1, TRA2, TRB1, TRB2

# Function to extract TRA1, TRA2, TRB1, TRB2 columns from sequences
extract_tcr_columns <- function(tcr_vector) {
    # Initialize empty lists to store the columns
    TRA1nt <- TRA2nt <- TRB1nt <- TRB2nt <- vector("character", length(tcr_vector))

    # Loop through each element in the input vector
    for (i in seq_along(tcr_vector)) {
        # Split each entry by ';' to get individual chains
        chains <- unlist(strsplit(tcr_vector[i], ";"))

        # Track the number of TRA and TRB chains found
        tra_count <- trb_count <- 0

        # Assign chains to the appropriate variables
        for (chain in chains) {
            if (grepl("^TRA:", chain)) {
                tra_count <- tra_count + 1
                if (tra_count == 1) {
                    TRA1nt[i] <- sub("^TRA:", "", chain)
                } else if (tra_count == 2) {
                    TRA2nt[i] <- sub("^TRA:", "", chain)
                }
            } else if (grepl("^TRB:", chain)) {
                trb_count <- trb_count + 1
                if (trb_count == 1) {
                    TRB1nt[i] <- sub("^TRB:", "", chain)
                } else if (trb_count == 2) {
                    TRB2nt[i] <- sub("^TRB:", "", chain)
                }
            }
        }
    }

    # Combine into a data frame
    data.frame(TRA1nt, TRA2nt, TRB1nt, TRB2nt, stringsAsFactors = FALSE)
}

tcr_nt_df <- extract_tcr_columns(BCC_SCC_TCR$cdr3s_nt)

### Amino Acid
# Function to extract TRA1, TRA2, TRB1, TRB2 columns from sequences
extract_tcraa_columns <- function(tcr_vector) {
    # Initialize empty lists to store the columns
    TRA1aa <- TRA2aa <- TRB1aa <- TRB2aa <- vector("character", length(tcr_vector))

    # Loop through each elemeaa in the input vector
    for (i in seq_along(tcr_vector)) {
        # Split each eaary by ';' to get individual chains
        chains <- unlist(strsplit(tcr_vector[i], ";"))

        # Track the number of TRA and TRB chains found
        tra_count <- trb_count <- 0

        # Assign chains to the appropriate variables
        for (chain in chains) {
            if (grepl("^TRA:", chain)) {
                tra_count <- tra_count + 1
                if (tra_count == 1) {
                    TRA1aa[i] <- sub("^TRA:", "", chain)
                } else if (tra_count == 2) {
                    TRA2aa[i] <- sub("^TRA:", "", chain)
                }
            } else if (grepl("^TRB:", chain)) {
                trb_count <- trb_count + 1
                if (trb_count == 1) {
                    TRB1aa[i] <- sub("^TRB:", "", chain)
                } else if (trb_count == 2) {
                    TRB2aa[i] <- sub("^TRB:", "", chain)
                }
            }
        }
    }

    # Combine into a data frame
    data.frame(TRA1aa, TRA2aa, TRB1aa, TRB2aa, stringsAsFactors = FALSE)
}

tcr_aa_df <- extract_tcraa_columns(BCC_SCC_TCR$cdr3s_aa)
BCC_SCC_TCR_combined <- cbind(BCC_SCC_TCR, tcr_nt_df, tcr_aa_df)
BCC_SCC_TCR_combined$cellbarcode <- rownames(BCC_SCC_TCR_combined)

CD8_Tcells_noCD4@meta.data$cellbarcode <- rownames(CD8_Tcells_noCD4@meta.data)
CD8_Tcells_noCD4@meta.data$cellbarcode <- gsub("_16$", "", CD8_Tcells_noCD4@meta.data$cellbarcode)

# combined_df <- merge(CD8_Tcells_noCD4@meta.data, BCC_SCC_TCR_combined, by = "cellbarcode", all.x = TRUE)

#### Add each column
CD8_Tcells_noCD4@meta.data$cellbarcode2 <- BCC_SCC_TCR_combined[match(CD8_Tcells_noCD4@meta.data$cellbarcode, BCC_SCC_TCR_combined$cellbarcode), "cellbarcode"] ## for sanity check
testing <- CD8_Tcells_noCD4@meta.data[!is.na(CD8_Tcells_noCD4@meta.data$cellbarcode2), ]
stopifnot(all(testing$cellbarcode == testing$cellbarcode2))

columntokeep <- colnames(BCC_SCC_TCR_combined)[1:10]

for (i in 1:length(columntokeep)) {
    CD8_Tcells_noCD4@meta.data[, columntokeep[i]] <- BCC_SCC_TCR_combined[match(CD8_Tcells_noCD4@meta.data$cellbarcode, BCC_SCC_TCR_combined$cellbarcode), columntokeep[i]]
}

### extracting the required columns for the TCR analysis
required_df <- CD8_Tcells_noCD4@meta.data[, c(columntokeep, "seurat_clusters", "celltypes2")]
require_df_noNA <- na.omit(required_df)

library(ggplot2)
library(Seurat)
library(dplyr)
library(reshape2)
library(ArchR)
library(stringr)

total_celltype <- length(unique(require_df_noNA$Celltypes2))
celltypes <- unique(require_df_noNA$Celltypes2)
df <- data.frame(matrix(ncol = length(celltypes), nrow = length(celltypes)))
colnames(df) <- celltypes
rownames(df) <- celltypes

object_TCR_cluster <- table(require_df_noNA[, "TRB1aa"], require_df_noNA[, "Celltypes2"])
# object_TCR_cluster <- object_TCR_cluster[grep("^_$", rownames(object_TCR_cluster), invert = TRUE), ]
### Remove the first since it is empty
object_TCR_cluster_2 <- object_TCR_cluster[-1, ]
object_TCR_cluster_2[object_TCR_cluster_2 > 0] <- 1

for (i in 1:length(df)) {
    for (j in 1:length(df)) {
        df[i, j] <- as.data.frame(table(object_TCR_cluster[, i] + object_TCR_cluster[, j] == 2)[2])[, 1]
    }
}

# df[is.na(df)] <- 0
df$cluster <- rownames(df)
hm <- melt(df)
colnames(hm) <- c("Celltypes1", "Celltypes2", "Clones")

hm$Celltypes1 <- factor(hm$Celltypes1, levels = c(
    "SCM-like", "EM1", "Effector Activated", "EM ISG", "EM2", "GD",
    "GD CD62L", "TRM", "TEX ISG", "TEMRA", "Pre TEX", "Trans TEX", "Term TEX", "Proliferating"
))
hm$Celltypes2 <- factor(hm$Celltypes2, levels = c(
    "SCM-like", "EM1", "Effector Activated", "EM ISG", "EM2", "GD",
    "GD CD62L", "TRM", "TEX ISG", "TEMRA", "Pre TEX", "Trans TEX", "Term TEX", "Proliferating"
))

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/TCR/"
library(ggplot2)
p <- ggplot(hm, aes(x = Celltypes1, y = Celltypes2, fill = Clones)) +
    geom_tile(color = "black") +
    theme_bw() +
    coord_equal() +
    scale_fill_distiller(palette = "RdYlBu", direction = 1) +
    scale_fill_gradientn(limits = c(0, 1371), colours = c(ArchRPalettes$solarExtra)) +
    # guides(fill=F) + # removing legend for `fill`
    labs(title = "Celltypes") + # using a title instead
    geom_text(aes(label = Clones), color = "black") + # printing values +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

dir.create(paste0(savedir, "clonal_sharing"), showWarnings = FALSE)
pdf(paste(savedir, "/", "clonal_sharing/BCC_SCC_TCR_number.pdf", sep = ""), width = 10, height = 8)
print(p)
dev.off()

write.table(df, paste(savedir, "/", "clonal_sharing/BCC_SCC_TCR_number.txt", sep = ""),
    quote = F, row.names = T, col.names = T
)

### Seurat Clusters
### Using Bhattacharya cofficient ####
# from this paper https://www.sciencedirect.com/science/article/pii/S0092867421013337?via%3Dihub#bib28
require_df_noNA$TRA1aa_TRA2aa_TRB1aa_TRB2aa <- paste(require_df_noNA$TRA1aa, require_df_noNA$TRA2aa, require_df_noNA$TRB1aa, require_df_noNA$TRB2aa, sep = "_")

object_TCR_cluster <- table(require_df_noNA[, "TRB1aa"], require_df_noNA[, "Celltypes2"])
object_TCR_cluster <- object_TCR_cluster[grep("^_$", rownames(object_TCR_cluster), invert = TRUE), ]
### Remove the first since it is empty
object_TCR_cluster <- object_TCR_cluster[-1, ]

rm(df)
df <- data.frame(matrix(ncol = ncol(object_TCR_cluster), nrow = ncol(object_TCR_cluster)))
colnames(df) <- colnames(object_TCR_cluster)
rownames(df) <- colnames(object_TCR_cluster)


# test <- head(object_TCR_cluster,20) %>% tail(10) ## this is more interesting
column_sum <- colSums(object_TCR_cluster) %>% as.vector()
vec <- c()

for (i in 1:length(df)) {
    for (j in 1:length(df)) {
        req <- object_TCR_cluster[object_TCR_cluster[, i] + object_TCR_cluster[, j] > 1, ]
        if (length(req) > 0) {
            if (length(req) > length(column_sum)) {
                req2 <- req[, c(i, j)]
                vec <- c()
                for (k in 1:nrow(req2)) {
                    vec[k] <- sqrt((req2[k, 1] / column_sum[i]) * (req2[k, 2] / column_sum[j]))
                }
            } else if (length(req) == length(column_sum)) {
                req_df <- as.data.frame(req)
                req2 <- req_df[c(i, j), ]
                vec <- c()
                vec <- sqrt((req2[1] / column_sum[i]) * (req2[2] / column_sum[j]))
            }
            df[i, j] <- sum(vec)
        }
    }
}

# df[is.na(df)] <- 0

df$cluster <- rownames(df)
hm <- melt(df)
colnames(hm) <- c("Celltypes1", "Celltypes2", "TRSS")
hm$TRSS <- round(hm$TRSS, digits = 2)

hm$Celltypes1 <- factor(hm$Celltypes1, levels = c(
    "SCM-like", "EM1", "Effector Activated", "EM ISG", "EM2", "GD",
    "GD CD62L", "TRM", "TEX ISG", "TEMRA", "Pre TEX", "Trans TEX", "Term TEX", "Proliferating"
))
hm$Celltypes2 <- factor(hm$Celltypes2, levels = c(
    "SCM-like", "EM1", "Effector Activated", "EM ISG", "EM2", "GD",
    "GD CD62L", "TRM", "TEX ISG", "TEMRA", "Pre TEX", "Trans TEX", "Term TEX", "Proliferating"
))

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/TCR/"
library(ggplot2)
p <- ggplot(hm, aes(x = Celltypes1, y = Celltypes2, fill = TRSS)) +
    geom_tile(color = "black") +
    theme_bw() +
    coord_equal() +
    scale_fill_distiller(palette = "RdYlBu", direction = 1) +
    scale_fill_gradientn(limits = c(0, 0.5), colours = c(ArchRPalettes$solarExtra)) +
    # guides(fill=F) + # removing legend for `fill`
    labs(title = "Bhattacharya Coefficient Celltypes") + # using a title instead
    geom_text(aes(label = TRSS), color = "black") + # printing values +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

dir.create(paste0(savedir, "clonal_sharing"), showWarnings = FALSE)
pdf(paste(savedir, "/", "clonal_sharing/BCC_SCC_TRA1_TRA2_TRB1_TRB2_Bhattacharya.pdf", sep = ""), width = 10, height = 8)
print(p)
dev.off()

write.table(df, paste(savedir, "/", "clonal_sharing/BCC_SCC_TCRA1_TCRA2_TRAB1_TRAB2_number_Bhattacharya.txt", sep = ""),
    quote = F, row.names = T, col.names = T
)

### Diversity
size <- nrow(require_df_noNA)
celltypes <- unique(require_df_noNA$Celltypes2)
sha_ent <- rep(NA, length(celltypes) + 1)
TCR_col <- "TRA1aa_TRA2aa_TRB1aa_TRB2aa"


df <- as.data.frame(table(require_df_noNA[, TCR_col]))
# Replace empty strings with NA
df[df == ""] <- NA
df <- na.omit(df)
clone_size <- df[, 2]
total_sum <- sum(clone_size)
clone_size_prop <- (clone_size / total_sum)
sha_ent[1] <- -(sum(clone_size_prop * (log2(clone_size_prop))))

cluster_num <- levels(celltypes)

for (i in 1:length(cluster_num)) {
    tryCatch(
        {
            require_df_noNA_subset_clus <- require_df_noNA[grep(paste("^", cluster_num[i], "$", sep = ""), require_df_noNA[, "Celltypes2"]), ]
            df <- as.data.frame(table(require_df_noNA_subset_clus[, TCR_col]))
            df[df == ""] <- NA
            df <- na.omit(df)
            clone_size <- df[, 2]
            total_sum <- sum(clone_size)
            clone_size_prop <- (clone_size / total_sum)
            sha_ent[i + 1] <- -(sum(clone_size_prop * (log2(clone_size_prop))))
            # print(-(sum(clone_size_prop * (log2(clone_size_prop)))))
        },
        error = function(e) {
            # Error handling: If an error occurs, assign NA to the result
            sha_ent[i + 1] <- NA
        }
    )
}


diversity <- data.frame("Celltypes" = c("all", as.character(cluster_num)), "Shannon_TRB1" = sha_ent, )
diversity$Shannon_TRA1_TRA2_TRB1_TRB2 <- sha_ent

### Simpson Diversity
size <- nrow(require_df_noNA)
celltypes <- levels(require_df_noNA$Celltypes2)
sha_ent <- rep(NA, length(celltypes) + 1)
TCR_col <- "TRA1aa_TRB1aa"

df <- as.data.frame(table(require_df_noNA[, TCR_col]))
# Replace empty strings with NA
df[df == ""] <- NA
df <- na.omit(df)
clone_size <- df[, 2]

total_sum <- sum(clone_size)
clone_size_prop <- (clone_size / total_sum)^2
sha_ent[1] <- sum(clone_size_prop)

cluster_num <- levels(require_df_noNA$Celltypes2)

for (i in 1:length(cluster_num)) {
    tryCatch(
        {
            require_df_noNA_subset_clus <- require_df_noNA[grep(paste("^", cluster_num[i], "$", sep = ""), require_df_noNA[, "Celltypes2"]), ]
            df <- as.data.frame(table(require_df_noNA_subset_clus[, TCR_col]))
            df[df == ""] <- NA
            df <- na.omit(df)
            clone_size <- df[, 2]
            total_sum <- sum(clone_size)
            clone_size_prop <- (clone_size / total_sum)^2
            sha_ent[i + 1] <- sum(clone_size_prop)
        },
        error = function(e) {
            # Error handling: If an error occurs, assign NA to the result
            sha_ent[i + 1] <- NA
        }
    )
}

diversity$Inv_Simp_TRA1_TRB1aa <- 1 / sha_ent

write.table(diversity,
    "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/TCR/diversity.txt",
    quote = F, row.names = F, col.names = T, sep = "\t"
)

### Clonality
rm(sha_ent)
TCR_col <- "TRA1aa_TRA2aa_TRB1aa_TRB2aa"
clus_col <- "Celltypes2"
sha_ent <- rep(NA, length(levels(require_df_noNA$Celltypes2)) + 1)

df <- as.data.frame(table(require_df_noNA[, TCR_col]))
df[df == ""] <- NA
df <- na.omit(df)
clone_size <- df[, 2]

total_sum <- sum(clone_size)
clone_size_prop <- (clone_size / total_sum)

sorted_proportions <- sort(clone_size_prop)
n <- length(sorted_proportions)
denominator <- n * (sum(sorted_proportions))

rm(combine)
combine <- vector()

for (j in 1:length(sorted_proportions)) {
    combine[j] <- (2 * j - n - 1) * sorted_proportions[j]
}

sha_ent[1] <- sum(combine) / denominator

cluster_num <- levels(require_df_noNA$Celltypes2)

for (i in 1:length(cluster_num)) {
    tryCatch(
        {
            require_df_noNA_subset_clus <- require_df_noNA[grep(paste("^", cluster_num[i], "$", sep = ""), require_df_noNA[, clus_col]), ]
            df <- as.data.frame(table(require_df_noNA_subset_clus[, TCR_col]))
            df[df == ""] <- NA
            df <- na.omit(df)
            clone_size <- df[, 2]
            total_sum <- sum(clone_size)
            clone_size_prop <- (clone_size / total_sum)
            sorted_proportions <- sort(clone_size_prop)
            n <- length(sorted_proportions)
            denominator <- n * (sum(sorted_proportions))

            rm(combine)
            combine <- vector()

            for (j in 1:length(sorted_proportions)) {
                combine[j] <- (2 * j - n - 1) * sorted_proportions[j]
            }
            sha_ent[i + 1] <- sum(combine) / denominator
        },
        error = function(e) {
            # Error handling: If an error occurs, assign NA to the result
            sha_ent[i + 1] <- NA
        }
    )
}

# clonality <- data.frame("Celltypes" = c("all", as.character(cluster_num)), "Gini_Coef_TRB1" = sha_ent)
clonality$Gini_Coef_TRA1aa_TRA2aa_TRB1aa_TRB2aa <- sha_ent

write.table(clonality, "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/TCR/clonaity.txt",
    quote = F,
    row.names = F, col.names = T, sep = "\t"
)

#### Making the clones to show sharing
dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)

object_TCR_cluster_2[rowSums(object_TCR_cluster_2) > 10, ]
clones <- c("CASGLAGYNEQFF", "CATSDLGQGVGNEQFF")
clones <- c("CASRRILAGGPEGTQYF", "CASSLGLVQPQHF", "CASSLGTVNTEAFF")
clones <- c("CASSLGQGGSLWANVLTF")
clones <- c("CSARDRGYGNTIYF", "CASSYSGTNMNTEAFF")


rm(plot_list)
plot_list <- list()
for (i in 1:length(clones)) {
    cellnames <- rownames(CD8_Tcells_noCD4@meta.data[grep(clones[i], CD8_Tcells_noCD4@meta.data$TRB1aa), ])
    plot_list[[i]] <- DimPlot(CD8_Tcells_noCD4, cells.highlight = cellnames, reduction = "umap", sizes.highlight = 1) + ggtitle(clones[i])
}

pdf(paste0(savedir, "UMAP/Sharing.pdf"))
plot_list
dev.off()

### Constructing a Phylogenetic Tree
# Load libraries
library(ape) # For phylogenetic tree construction and visualization
library(seqinr) # For handling sequence data
library(phangorn) # For building and analyzing phylogenetic trees
tcr_sequences <- require_df_noNA[, "TRB1aa"]
names(tcr_sequences) <- paste0("TCR_", seq_along(tcr_sequences))

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/TCR/phylogenetic/"
# Write sequences to a FASTA file
write.fasta(
    sequences = as.list(tcr_sequences),
    names = names(tcr_sequences),
    file.out = paste0(savedir, "tcr_sequences.fasta")
)

#### Alignment
# conda activate daisy
# mafft --auto tcr_sequences.fasta > aligned_tcr.fasta

aligned <- read.alignment("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/TCR/phylogenetic/aligned_tcr.fasta", format = "fasta")

dist_matrix <- dist.alignment(aligned, "identity")
nj_tree <- nj(dist_matrix)
rooted_tree <- root(nj_tree, outgroup = "TCR_1", resolve.root = TRUE)
plot(rooted_tree, main = "TCR Phylogenetic Tree", cex = 0.8)

##### Performing Analysis on Seurat clusters
### Using Bhattacharya cofficient ####
# from this paper https://www.sciencedirect.com/science/article/pii/S0092867421013337?via%3Dihub#bib28

require_df_noNA$TRA1aa_TRA2aa_TRB1aa_TRB2aa <- paste(require_df_noNA$TRA1aa, require_df_noNA$TRA2aa, require_df_noNA$TRB1aa, require_df_noNA$TRB2aa, sep = "_")
clus_col <- "seurat_clusters"
object_TCR_cluster <- table(require_df_noNA[, "TRB1aa"], require_df_noNA[, clus_col])
object_TCR_cluster <- object_TCR_cluster[grep("^_$", rownames(object_TCR_cluster), invert = TRUE), ]
### Remove the first since it is empty
object_TCR_cluster <- object_TCR_cluster[-1, ]

rm(df)
df <- data.frame(matrix(ncol = ncol(object_TCR_cluster), nrow = ncol(object_TCR_cluster)))
colnames(df) <- colnames(object_TCR_cluster)
rownames(df) <- colnames(object_TCR_cluster)


# test <- head(object_TCR_cluster,20) %>% tail(10) ## this is more interesting
column_sum <- colSums(object_TCR_cluster) %>% as.vector()
vec <- c()

for (i in 1:length(df)) {
    for (j in 1:length(df)) {
        req <- object_TCR_cluster[object_TCR_cluster[, i] + object_TCR_cluster[, j] > 1, ]
        if (length(req) > 0) {
            if (length(req) > length(column_sum)) {
                req2 <- req[, c(i, j)]
                vec <- c()
                for (k in 1:nrow(req2)) {
                    vec[k] <- sqrt((req2[k, 1] / column_sum[i]) * (req2[k, 2] / column_sum[j]))
                }
            } else if (length(req) == length(column_sum)) {
                req_df <- as.data.frame(req)
                req2 <- req_df[c(i, j), ]
                vec <- c()
                vec <- sqrt((req2[1] / column_sum[i]) * (req2[2] / column_sum[j]))
            }
            df[i, j] <- sum(vec)
        }
    }
}

df <- df[-c(10, 25), ] #### Since they are Nan

df$cluster <- rownames(df)
hm <- melt(df)
colnames(hm) <- c("Celltypes1", "Celltypes2", "TRSS")
hm$TRSS <- round(hm$TRSS, digits = 2)

hm$Celltypes1 <- factor(hm$Celltypes1, levels = c(0:37))
hm$Celltypes2 <- factor(hm$Celltypes2, levels = c(0:37))

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/TCR/"
library(ggplot2)
p <- ggplot(hm, aes(x = Celltypes1, y = Celltypes2, fill = TRSS)) +
    geom_tile(color = "black") +
    theme_bw() +
    coord_equal() +
    scale_fill_distiller(palette = "RdYlBu", direction = 1) +
    scale_fill_gradientn(limits = c(0, 0.6), colours = c(ArchRPalettes$solarExtra)) +
    # guides(fill=F) + # removing legend for `fill`
    labs(title = "Bhattacharya Coefficient Celltypes") + # using a title instead
    geom_text(aes(label = TRSS), color = "black") + # printing values +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

dir.create(paste0(savedir, "clonal_sharing"), showWarnings = FALSE)
pdf(paste(savedir, "/", "clonal_sharing/BCC_SCC_TRB1_Bhattacharya_clusters.pdf", sep = ""), width = 16, height = 15)
print(p)
dev.off()

write.table(df, paste(savedir, "/", "clonal_sharing/BCC_SCC_TCRA1_TCRA2_TRAB1_TRAB2_number_Bhattacharya.txt", sep = ""),
    quote = F, row.names = T, col.names = T
)

### Diversity
size <- nrow(require_df_noNA)
celltypes <- unique(require_df_noNA$Celltypes2)
sha_ent <- rep(NA, length(celltypes) + 1)
TCR_col <- "TRA1aa_TRA2aa_TRB1aa_TRB2aa"


df <- as.data.frame(table(require_df_noNA[, TCR_col]))
# Replace empty strings with NA
df[df == ""] <- NA
df <- na.omit(df)
clone_size <- df[, 2]
total_sum <- sum(clone_size)
clone_size_prop <- (clone_size / total_sum)
sha_ent[1] <- -(sum(clone_size_prop * (log2(clone_size_prop))))

cluster_num <- levels(celltypes)

for (i in 1:length(cluster_num)) {
    tryCatch(
        {
            require_df_noNA_subset_clus <- require_df_noNA[grep(paste("^", cluster_num[i], "$", sep = ""), require_df_noNA[, "Celltypes2"]), ]
            df <- as.data.frame(table(require_df_noNA_subset_clus[, TCR_col]))
            df[df == ""] <- NA
            df <- na.omit(df)
            clone_size <- df[, 2]
            total_sum <- sum(clone_size)
            clone_size_prop <- (clone_size / total_sum)
            sha_ent[i + 1] <- -(sum(clone_size_prop * (log2(clone_size_prop))))
            # print(-(sum(clone_size_prop * (log2(clone_size_prop)))))
        },
        error = function(e) {
            # Error handling: If an error occurs, assign NA to the result
            sha_ent[i + 1] <- NA
        }
    )
}


diversity <- data.frame("Celltypes" = c("all", as.character(cluster_num)), "Shannon_TRB1" = sha_ent, )
diversity$Shannon_TRA1_TRA2_TRB1_TRB2 <- sha_ent

### Simpson Diversity
size <- nrow(require_df_noNA)
celltypes <- levels(require_df_noNA$Celltypes2)
sha_ent <- rep(NA, length(celltypes) + 1)
TCR_col <- "TRA1aa_TRB1aa"

df <- as.data.frame(table(require_df_noNA[, TCR_col]))
# Replace empty strings with NA
df[df == ""] <- NA
df <- na.omit(df)
clone_size <- df[, 2]

total_sum <- sum(clone_size)
clone_size_prop <- (clone_size / total_sum)^2
sha_ent[1] <- sum(clone_size_prop)

cluster_num <- levels(require_df_noNA$Celltypes2)

for (i in 1:length(cluster_num)) {
    tryCatch(
        {
            require_df_noNA_subset_clus <- require_df_noNA[grep(paste("^", cluster_num[i], "$", sep = ""), require_df_noNA[, "Celltypes2"]), ]
            df <- as.data.frame(table(require_df_noNA_subset_clus[, TCR_col]))
            df[df == ""] <- NA
            df <- na.omit(df)
            clone_size <- df[, 2]
            total_sum <- sum(clone_size)
            clone_size_prop <- (clone_size / total_sum)^2
            sha_ent[i + 1] <- sum(clone_size_prop)
        },
        error = function(e) {
            # Error handling: If an error occurs, assign NA to the result
            sha_ent[i + 1] <- NA
        }
    )
}

diversity$Inv_Simp_TRA1_TRB1aa <- 1 / sha_ent

write.table(diversity,
    "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/TCR/diversity.txt",
    quote = F, row.names = F, col.names = T, sep = "\t"
)

### Clonality
rm(sha_ent)
TCR_col <- "TRA1aa_TRA2aa_TRB1aa_TRB2aa"
clus_col <- "Celltypes2"
sha_ent <- rep(NA, length(levels(require_df_noNA$Celltypes2)) + 1)

df <- as.data.frame(table(require_df_noNA[, TCR_col]))
df[df == ""] <- NA
df <- na.omit(df)
clone_size <- df[, 2]

total_sum <- sum(clone_size)
clone_size_prop <- (clone_size / total_sum)

sorted_proportions <- sort(clone_size_prop)
n <- length(sorted_proportions)
denominator <- n * (sum(sorted_proportions))

rm(combine)
combine <- vector()

for (j in 1:length(sorted_proportions)) {
    combine[j] <- (2 * j - n - 1) * sorted_proportions[j]
}

sha_ent[1] <- sum(combine) / denominator

cluster_num <- levels(require_df_noNA$Celltypes2)

for (i in 1:length(cluster_num)) {
    tryCatch(
        {
            require_df_noNA_subset_clus <- require_df_noNA[grep(paste("^", cluster_num[i], "$", sep = ""), require_df_noNA[, clus_col]), ]
            df <- as.data.frame(table(require_df_noNA_subset_clus[, TCR_col]))
            df[df == ""] <- NA
            df <- na.omit(df)
            clone_size <- df[, 2]
            total_sum <- sum(clone_size)
            clone_size_prop <- (clone_size / total_sum)
            sorted_proportions <- sort(clone_size_prop)
            n <- length(sorted_proportions)
            denominator <- n * (sum(sorted_proportions))

            rm(combine)
            combine <- vector()

            for (j in 1:length(sorted_proportions)) {
                combine[j] <- (2 * j - n - 1) * sorted_proportions[j]
            }
            sha_ent[i + 1] <- sum(combine) / denominator
        },
        error = function(e) {
            # Error handling: If an error occurs, assign NA to the result
            sha_ent[i + 1] <- NA
        }
    )
}

# clonality <- data.frame("Celltypes" = c("all", as.character(cluster_num)), "Gini_Coef_TRB1" = sha_ent)
clonality$Gini_Coef_TRA1aa_TRA2aa_TRB1aa_TRB2aa <- sha_ent

write.table(clonality, "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/TCR/clonaity.txt",
    quote = F,
    row.names = F, col.names = T, sep = "\t"
)
