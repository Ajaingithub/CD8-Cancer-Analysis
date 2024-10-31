### Final Presentation Figures
library(Seurat)
library(ArchR)
library(ggplot2)
library(dplyr)

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/Final_presentation/"
immune <- readRDS("/diazlab/data3/.abhinav/.immune/cancer_combined/project/raw_data/TICAtlas_RNA.rds")
Tcells <- readRDS("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/paper_integrated/Tcell_obj.RDS")

DefaultAssay(immune) <- "RNA"
genes <- c("CD3E", "CD3D", "CD3G", "CD4", "CD8A", "CD8B")

p <- DotPlot(immune, genes, assay = "RNA", group.by = "cell_type") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    # coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD4_CD8_RNA.pdf"), height = 10, width = 8)
print(p)
dev.off()

### CD8 T cells subset celltypes
Stem_like <- c(
    "CCR7", "CD27", "CD28", "IL7R", "FAS", "TCF7", "LEF1", "BCL2", "CXCR3", "TGFBR1", "TGFBR2", "ITGAL"
)
Effector <- c(
    "IL2RB", "CD58", "KLRG1", "GZMM", "GZMB", "IFNG", "GZMK",
    "NKG7", "ZEB2", "GZMA", "CX3CR1", "IFNG", "NOSIP", "FOXP1",
    "ZEB2", "BCL6", "EOMES", "PRDM1", "RUNX3", "ITGA4", "CXCR6"
)

Cytotoxic <- c(
    "GZMB", "PRF1", "TBX21", "LAMP1", "FGFBP2", "KLRG1", "GNLY", "ZNF683", "GZMH", "FASL", "IFNG", "NKG7",
    "KLRG1", "NKG2D", "LFA3", "B3GAT1", "LILRB1", "TYROBP", "KIR2DL1", "KIR3DL1", "KLRK1", "IFNG", "ITGA4", "CXCR6"
)
Exhausted <- c("PDCD1", "LAG3", "HAVCR2", "TOX", "ENTPD1", "ID2", "CTLA4", "TNFRSF4", "TNFRSF18", "TNFSF8")

library(ArchR)
library(UCell)
library(ggplot2)
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(CD8_Tcells) <- "RNA"
rm(markers)
markers <- list()
files <- c("Stem_like", "Effector")

for (i in 1:length(files)) {
    Tcellsubset <- unique(get(files[i]))
    Tcellsubset <- rownames(CD8_Tcells)[match(Tcellsubset, rownames(CD8_Tcells), nomatch = 0)]
    markers[[files[i]]] <- Tcellsubset
    # assign(files[i], Tcellsubset)
}

DefaultAssay(CD8_Tcells_noCD4) <- "MAGIC_RNA"
CD8_Tcells_noCD4 <- AddModuleScore(CD8_Tcells_noCD4, features = markers, slot = "data")
colnames(CD8_Tcells@meta.data)[50:53] <- files

CD8_Tcells_noCD4@meta.data$celltypes1 <- factor(CD8_Tcells_noCD4@meta.data$celltypes1,
    levels = c("Stem-like", "Activated", "Early Eff", "ISG Eff", "Effector", "Gamma Delta", "Cytotoxic", "ISG TEX", "Pre TEX", "Trans TEX", "Term TEX", "Proliferating")
)

files <- c("Cytotoxic", "Exhaustion", "Stem-like", "Effector")
colnames(CD8_Tcells_noCD4@meta.data)[52:55] <- files
p <- DotPlot(CD8_Tcells_noCD4, c("Stem-like", "Effector", "Cytotoxic", "Exhaustion"), assay = "RNA", group.by = "celltypes1") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD4_CD8_RNA_2.pdf"), height = 4, width = 7.5)
print(p)
dev.off()

genes <- c(
    "CCR7", "TCF7", "LEF1", "IL7R",
    "CD69", "FOS", "JUN",
    "RUNX3", "PRDM1", "ZEB2",
    "ISG15", "OAS1", "MX1",
    "IFNG", "NKG7", "EGR2",
    "KLRC1", "TRDC", "TRGC1",
    "FGFBP2", "GNLY", "TBX21", "PRF1",
    "PDCD1", "LAG3", "TOX", "ENTPD1",
    "MKI67", "TOP2A"
)

p <- DotPlot(CD8_Tcells_noCD4, genes, assay = "RNA", group.by = "celltypes1") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_Tcell_phenotyping_genes.pdf"), height = 9, width = 7)
print(p)
dev.off()

### TRM
patterns <- c(
    "Effector", "Stem-like", "Early Eff", "Activated", "Cytotoxic",
    "Term TEX", "Proliferating", "ISG Eff", "Trans TEX", "Pre TEX", "Gamma Delta", "ISG TEX"
)

replacements <- c(
    "Acute Ag", "Newly Entering", "Newly Entering", "Newly Entering", "Cytotoxic",
    "Chronic Ag", "Proliferating", "ISG Eff", "Chronic Ag", "Chronic Ag", "Gamma Delta", "Chronic Ag"
)

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/Final_presentation/"
CD8_Tcells_noCD4@meta.data$celltypes3 <- CD8_Tcells_noCD4@meta.data$celltypes2
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD8_Tcells_noCD4@meta.data$celltypes3 <- str_replace_all(CD8_Tcells_noCD4@meta.data$celltypes3, pattern, replacements[i])
}

pdf(paste0(savedir, "UMAP/CD8_Tcells_3.pdf"), width = 5, height = 5.5)
DimPlot(CD8_Tcells_noCD4, group.by = "celltypes3", reduction = "umap", label = TRUE) + NoLegend()
dev.off()

# Define your color palette (adjust colors as needed)
celltype_colors <- c(
    "Newly Entering" = "mediumseagreen",
    "ISG Eff" = "gray87",
    "Acute Ag" = "#b87f16",
    "Gamma Delta" = "gray87",
    "Cytotoxic" = "gray87",
    "TRM" = "#16b8b0",
    "Chronic Ag" = "#db2a2a",
    "Proliferating" = "gray87"
)

pdf(paste0(savedir, "UMAP/CD8_Tcells_3_noLabels.pdf"), width = 5, height = 5.5)
DimPlot(CD8_Tcells_noCD4, group.by = "celltypes3", reduction = "umap", label = FALSE, cols = celltype_colors) + NoLegend()
dev.off()

#### Identifying the Tissue Resident Memory (TRM), Tumor entering and Tumore reactive T cells
Trm <- c("CD69", "ITGAE", "ZNF683", "TNF")
Ag_expr <- c("CXCR6", "ENTPD1", "ITGAE", "PDCD1", "TIM3", "LAG3", "TNFRSF9", "LAYN", "TOX", "GZMB", "PRF1", "TNFRSF18")
Newly <- c("CXCR4", "CCR7", "CD27", "CD28", "CXCR5", "CCR5", "CCR7", "S1PR1")

library(ArchR)
library(UCell)
library(ggplot2)
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
rm(markers)
markers <- list()
files <- c("Newly", "Trm", "Ag_expr")

for (i in 1:length(files)) {
    Tcellsubset <- unique(get(files[i]))
    Tcellsubset <- rownames(CD8_Tcells)[match(Tcellsubset, rownames(CD8_Tcells), nomatch = 0)]
    markers[[files[i]]] <- Tcellsubset
    # assign(files[i], Tcellsubset)
}

DefaultAssay(CD8_Tcells_noCD4) <- "MAGIC_RNA"
CD8_Tcells_noCD4 <- AddModuleScore(CD8_Tcells_noCD4, features = markers, slot = "data")
colnames(CD8_Tcells_noCD4@meta.data)[63:65] <- c("Newly_Tumor", "Tissue_Resident", "Antigen Experienced")

CD8_Tcells_noCD4@meta.data$celltypes2 <- factor(CD8_Tcells_noCD4@meta.data$celltypes2,
    levels = c("Stem-like", "Activated", "Early Eff", "ISG Eff", "Effector", "Gamma Delta", "TRM", "Cytotoxic", "ISG TEX", "Pre TEX", "Trans TEX", "Term TEX", "Proliferating")
)

# files <- c("Cytotoxic", "Exhaustion", "Stem-like", "Effector")
# colnames(CD8_Tcells_noCD4@meta.data)[52:55] <- files
p <- DotPlot(CD8_Tcells_noCD4, features = c("Newly_Tumor", "Tissue_Resident", "Antigen Experienced"), assay = "RNA", group.by = "celltypes2") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_TRM_Ag_expr_Newly_2.pdf"), height = 3.5, width = 8)
print(p)
dev.off()

genes <- c(
    "CXCR4", "CCR7", "S1PR1", "GZMK",
    "CD69", "ITGAE", "ZNF683", "PRF1", "IFNG", "TNF",
    "CXCR6", "TNFRSF18", "CXCL13", "HACR2", "TOX"
)

p <- DotPlot(CD8_Tcells_noCD4, features = genes, assay = "RNA", group.by = "celltypes2") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/CD8_TRM_genes.pdf"), height = 6, width = 7.5)
print(p)
dev.off()

# Gene Number
cancer_genes <- read.table("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/scdifferential/all/cancer_gene_number.txt", header = TRUE)

# Load necessary library
library(ggplot2)

p <- ggplot(cancer_genes, aes(x = Cancer, y = DEGs, fill = Cancer)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white"
    )) +
    labs(
        title = "Number of DEGs per Cancer Type",
        x = "Cancer Type", y = "Number of DEGs"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

pdf("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/scdifferential/all/DEGs.pdf", width = 7, height = 5)
p
dev.off()


### Seurat Clusters
### Using Bhattacharya cofficient ####
# from this paper https://www.sciencedirect.com/science/article/pii/S0092867421013337?via%3Dihub#bib28
patterns <- c("EM2", "SCM-like", "EM1", "Effector Activated", "TEMRA", "GD CD62L", "EM ISG", "TRM", "GD", "TEX ISG")
replacements <- c("Effector", "Stem-like", "Early Eff", "Activated", "Cytotoxic", "Cytotoxic", "ISG Eff", "TRM", "Gamma Delta", "ISG TEX")

savedir <- "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/Final_presentation/"
CD8_Tcells_noCD4@meta.data$celltypes2 <- CD8_Tcells_noCD4@meta.data$Celltypes2
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD8_Tcells_noCD4@meta.data$celltypes2 <- str_replace_all(CD8_Tcells_noCD4@meta.data$celltypes2, pattern, replacements[i])
}

require_df_noNA_2 <- require_df_noNA[grep("Gamma Delta", require_df_noNA$celltypes2, invert = TRUE), ]
require_df_noNA <- require_df_noNA_2
require_df_noNA$TRA1aa_TRA2aa_TRB1aa_TRB2aa <- paste(require_df_noNA$TRA1aa, require_df_noNA$TRA2aa, require_df_noNA$TRB1aa, require_df_noNA$TRB2aa, sep = "_")

object_TCR_cluster <- table(require_df_noNA[, "TRB1aa"], require_df_noNA[, "celltypes2"])
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
    c(
        "Stem-like", "Activated", "Early Eff", "Effector", "ISG Eff",
        "TRM", "Cytotoxic", "Pre TEX", "Trans TEX", "Term TEX",
        "ISG TEX", "Proliferating"
    )
))
hm$Celltypes2 <- factor(hm$Celltypes2, levels = c(
    "Stem-like", "Activated", "Early Eff", "Effector", "ISG Eff",
    "TRM", "Cytotoxic", "Pre TEX", "Trans TEX", "Term TEX",
    "ISG TEX", "Proliferating"
))

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
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20)
    )

dir.create(paste0(savedir, "clonal_sharing"), showWarnings = FALSE)
pdf(paste(savedir, "/", "clonal_sharing/BCC_SCC_TRB1aa_Bhattacharya_Celltypes2.pdf", sep = ""), width = 10, height = 8)
print(p)
dev.off()

write.table(df, paste(savedir, "/", "clonal_sharing/BCC_SCC_TRB1aa_Bhattacharya_Celltypes2.txt", sep = ""),
    quote = F, row.names = T, col.names = T
)

merged_data <- readRDS("CZI_TS_Fibroblast.rds")
DefaultAssay(merged_data) <- "RNA"
merged_data <- NormalizeData(merged_data, normalization.method = "LogNormalize", scale.factor = 10000)
merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged_data)
merged_data <- ScaleData(merged_data)
merged_data <- RunPCA(merged_data)
merged_data <- FindNeighbors(merged_data, dims = 1:10)
merged_data <- FindClusters(merged_data, resolution = c(5, 3, 1, 0.9, 0.7, 0.5, 0.3, 0.1, 0.09, 0.07, 0.05, 0.03, 0.02, 0.01, 0.008, 0.006, 0.004, 0.002, 0.0001, 0.0005))
Idents(merged_data) <- "tissue_cell_ontology"
FibroMarkers <- FindAllMarkers(merged_data)
write.table(FibroMarkers, file = "CZI_TS_FibroMarkers.txt", sep = "\t")
saveRDS(merged_data, file = "CZI_TS_Fibroblast_done.rds")

library(Seurat)
merged_data <- readRDS("CZI_TS_Tcell_Muscle.rds")
DefaultAssay(merged_data) <- "RNA"
merged_data <- NormalizeData(merged_data, normalization.method = "LogNormalize", scale.factor = 10000)
merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 2000)
merged_data <- ScaleData(merged_data)
merged_data <- RunPCA(merged_data)
merged_data <- FindNeighbors(merged_data, dims = 1:10)
merged_data <- FindClusters(merged_data, resolution = c(5, 3, 1, 0.9, 0.7, 0.5, 0.3, 0.1, 0.09, 0.07, 0.05, 0.03, 0.02, 0.01, 0.008, 0.006, 0.004, 0.002, 0.0001, 0.0005))
Idents(merged_data) <- "tissue_cell_ontology"
FibroMarkers <- FindAllMarkers(merged_data)
write.table(FibroMarkers, file = "CZI_TS_TcellMarkers.txt", sep = "\t")
saveRDS(merged_data, file = "CZI_TS_Tcell_done.rds")

library(Seurat)
merged_data <- readRDS("CZI_TS_Endo2_Endothelial.rds")
DefaultAssay(merged_data) <- "RNA"
merged_data <- NormalizeData(merged_data, normalization.method = "LogNormalize", scale.factor = 10000)
merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 2000)
merged_data <- ScaleData(merged_data)
merged_data <- RunPCA(merged_data)
merged_data <- FindNeighbors(merged_data, dims = 1:10)
merged_data <- FindClusters(merged_data, resolution = c(5, 3, 1, 0.9, 0.7, 0.5, 0.3, 0.1, 0.09, 0.07, 0.05, 0.03, 0.02, 0.01, 0.008, 0.006, 0.004, 0.002, 0.0001, 0.0005))
Idents(merged_data) <- "tissue_cell_ontology"
FibroMarkers <- FindAllMarkers(merged_data)
write.table(FibroMarkers, file = "CZI_TS_EndoEndo2Markers.txt", sep = "\t")
saveRDS(merged_data, file = "CZI_TS_Endothelial_done.rds")

library(Seurat)
merged_data <- readRDS("CZI_TS_SMC_Muscle.rds")
DefaultAssay(merged_data) <- "RNA"
merged_data <- NormalizeData(merged_data, normalization.method = "LogNormalize", scale.factor = 10000)
merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged_data)
merged_data <- ScaleData(merged_data)
merged_data <- RunPCA(merged_data)
merged_data <- FindNeighbors(merged_data, dims = 1:10)
merged_data <- FindClusters(merged_data, resolution = c(5, 3, 1, 0.9, 0.7, 0.5, 0.3, 0.1, 0.09, 0.07, 0.05, 0.03, 0.02, 0.01, 0.008, 0.006, 0.004, 0.002, 0.0001, 0.0005))
Idents(merged_data) <- "tissue_cell_ontology"
FibroMarkers <- FindAllMarkers(merged_data)
write.table(FibroMarkers, file = "CZI_TS_SMCMarkers.txt", sep = "\t")
saveRDS(merged_data, file = "CZI_TS_SMC_done.rds")
