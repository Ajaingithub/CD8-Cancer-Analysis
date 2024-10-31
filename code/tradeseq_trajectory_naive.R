# Load required libraries
library(Seurat)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(dplyr)

# Read the Seurat object and convert it to SingleCellExperiment
CD8_Tcells_noCD4_QC_down_RNA <- readRDS(
    "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/saveRDS_obj/CD8_Tcells_noCD4_QC_down_RNA.RDS"
)
sce <- as.SingleCellExperiment(CD8_Tcells_noCD4_QC_down_RNA, assay = "RNA")

# Apply Slingshot for trajectory inference
CD4_sce1 <- slingshot(
    sce,
    clusterLabels = "Celltypes2",
    reducedDim = "UMAP",
    approx_points = 100, # Controls smoothing of the trajectory
    start.clus = "SCM-like", # Start trajectory from cluster 6
    omega = TRUE, # Include cluster connectivity
    omega_scale = 1.5 # Adjust the weight of edges
)

saveRDS(
    CD4_sce1,
    "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/trajectory/SCM-like_start_slingshot.RDS"
)

# Prepare count matrix and SlingshotDataSet
counts <- as.matrix(CD8_Tcells_noCD4_QC_down_RNA@assays$RNA@counts)
CD4_sds <- SlingshotDataSet(CD4_sce1)

# Fit the generalized additive model (GAM) using tradeSeq
sce_gam <- fitGAM(counts = counts, sds = CD4_sds)

# Save the fitted GAM model to an RDS file
saveRDS(
    sce_gam,
    "/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/trajectory/sce_SCM_fit_GAM.RDS"
)

# Extract necessary data for plotting
curves <- CD4_sds@curves # Ensure the curves object is defined properly
clustering <- colData(sce)$Celltypes2 # Use cluster labels from metadata
filt_counts <- counts[rownames(sce), ] # Filtered counts based on features in sce

# Plot the gene counts along the trajectory
pdf(paste0("/diazlab/data3/.abhinav/.immune/cancer_combined/project/cell_phenotyping/downstream2/CD8/CCA_integrated/trajectory/CD8_SCM_genes_count.pdf"))
plotGeneCount(curves, filt_counts, clusters = clustering, models = sce_gam)
dev.off()
