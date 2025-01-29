# ==========================================================
# MDS Patient Single-cell Multimodal Analysis Pipeline
# ==========================================================

# ==========================================================
# 1. Overview Section
# ==========================================================

# ----------------------------------------------------------
# Title and Description
# ----------------------------------------------------------
# Analysis of pre/post-transplant MDS patient samples integrating:
# - RNA sequencing
# - TCR repertoire
# - Peptide dextramer (dex) libraries
# - CITE-seq hashing (to distinguish dex+ and dex- cells)

# ----------------------------------------------------------
# Patient Information
# ----------------------------------------------------------
# SRSF2_9:  MDS-EB1 (15-month post-transplant)
#           7.3K dex(+) + 8.7K dex(-) cells
# SRSF2_10: MDS-EB1 (pre-transplant)
#           5K dex(+) + 11K dex(-) cells
# Dextramer: 46 Test, 1 Positive, 3 Negative controls

# ----------------------------------------------------------
# Project Goals
# ----------------------------------------------------------
# - Identify TCR clonotypes binding peptide neoantigens (dex+ cells).
# - Characterize immune changes pre- and post-transplant.

# ==========================================================
# 2. Setup 
# ==========================================================

# ----------------------------------------------------------
# Library Loading
# ----------------------------------------------------------
library(tidyverse)      # Data manipulation and visualization
library(Seurat)         # Single-cell analysis
library(patchwork)      # Plot composition
library(ggplot2)        # Plotting
library(RColorBrewer)   # Color palettes
library(scRepertoire)   # TCR analysis
library(gridExtra)      # Plot arrangement
library(ggrepel)        # Text labels
library(cowplot)        # Plot composition

# ----------------------------------------------------------
# Global Variables
# ----------------------------------------------------------

threshold_umis <- 3
threshold_mito <- 0.10
threshold_nCount_RNA <- 500
threshold_nFeature_RNA <- 250

# ----------------------------------------------------------
# Directory Paths
# ----------------------------------------------------------

# Load the helper functions
source("helper_functions/data_loading_helpers.R") # Different modality data loadings and integrations
source("helper_functions/cite_data_processing.R")  # Define the Dex+ population using cite-seq manual thresholding
source("helper_functions/plotting_functions.R") # Plotting functions
source("helper_functions/dextramer_peptide_scoring.R") # Z-score analysis for peptide clonotype pairs

# Set directories for each patient
set_directories <- function(patient_id, base_dir) {
  list(
    dir_gex = file.path(base_dir, paste0(patient_id, "/CellRangerGex_results")),
    dir_dex = file.path(base_dir, paste0(patient_id, "_dextramer_count/umi_count")),
    dir_TCR = file.path(base_dir, paste0(patient_id, "_TCR_VDJ/CellRangerVdj_results")),
    dir_CITE = file.path(base_dir, paste0(patient_id, "_hash_count/umi_count"))
  )
}

# Define base directory and patient IDs
base_dir <- "/Users/ecrosse/Desktop/"
dirs_SRSF2_9 <- set_directories("data_for_edie_third_batch_january/WJK-2859_SRSF2_9", base_dir)
dirs_SRSF2_10 <- set_directories("dextramer_data_for_edie_january_part_2/WJK-2864_SRSF2_10", base_dir)

# ==========================================================
# 1. Load Expression, CITE, Dextramer, and TCR Data into Seurat Objects
# ==========================================================

# Load and integrate data for SRSF2_9
adata <- load_expression_data(dirs_SRSF2_9$dir_gex)
adata <- load_dextramer_data(dirs_SRSF2_9$dir_dex, adata)
adata <- load_cite_data(dirs_SRSF2_9$dir_CITE, adata)
adata <- load_tcr_data(dirs_SRSF2_9$dir_TCR, adata)

# Load and integrate data for SRSF2_10
adata2 <- load_expression_data(dirs_SRSF2_10$dir_gex)
adata2 <- load_dextramer_data(dirs_SRSF2_10$dir_dex, adata2)
adata2 <- load_cite_data(dirs_SRSF2_10$dir_CITE, adata2)
adata2 <- load_tcr_data(dirs_SRSF2_10$dir_TCR, adata2)

# ==========================================================
# 2. Seurat Object Pre-Processing 
# ==========================================================

# ----------------------------------------------------------
# Seurat Merged Object Creation
# ----------------------------------------------------------
fss <- merge(adata, y = adata2, add.cell.ids = c("SRSF2_9", "SRSF2_10"), project = "MDS_pre_post_transplant")

# ----------------------------------------------------------
# Data Filtering, Normalization and Dataset Integration
# ----------------------------------------------------------

# Filter out low-quality cells
fss$mitoRatio <- PercentageFeatureSet(object = fss, pattern = "^MT-") / 100
fss <- subset(x = fss, 
              subset = (nCount_RNA >= threshold_nCount_RNA) & 
                       (nFeature_RNA >= threshold_nFeature_RNA) & 
                       (mitoRatio < threshold_mito))

# Normalize the data
fss <- NormalizeData(fss)

# Find variable features
fss <- FindVariableFeatures(fss, selection.method = "vst", nfeatures = 4000)

VFs <- VariableFeatures(fss)

# Filter out genes starting with 'TR'
VFs_filtered <- VFs[!grepl("^TR", VFs)]

# Update the Seurat object with the filtered variable features list
VariableFeatures(fss) <- VFs_filtered

all.genes <- rownames(fss)

# Scale the data
fss <- ScaleData(fss, features = all.genes)

# Run PCA
fss <- RunPCA(fss)

# Run CCA integration
fss <- IntegrateLayers(object = fss, method = CCAIntegration, 
                       orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)

# re-join layers after integration
fss[["RNA"]] <- JoinLayers(fss[["RNA"]])

# ==========================================================
# 3. CITE-seq Analysis - To Determine Cut-off for Dex+ Cells
# ==========================================================

# Normalize CITE-seq data
fss <- NormalizeData(fss, assay = "CITE", normalization.method = "CLR")

# Add sample information to the metadata
fss$sample <- ifelse(grepl("SRSF2_9", colnames(fss)), "SRSF2_9", "SRSF2_10")

# Subset the fss object into two samples to analyze separately
fss_SRSF2_9 <- subset(fss, subset = sample == "SRSF2_9")
fss_SRSF2_10 <- subset(fss, subset = sample == "SRSF2_10")

# Generate the density plot for SRSF2_9
cite_data_SRSF2_9 <- generate_cite_density_plot(
  seurat_obj = fss_SRSF2_9,
  assay_name = "CITE",
  title_suffix = "SRSF2_9"
)

# View the density plot
print(cite_data_SRSF2_9$plot)

# Apply the threshold for SRSF2_9
fss_SRSF2_9 <- apply_cite_threshold(
  seurat_obj = fss_SRSF2_9,
  assay_name = "CITE",
  hash_threshold = 1.5 # Determined by the above density plot
)

# Repeat for SRSF2_10
cite_data_SRSF2_10 <- generate_cite_density_plot(
  seurat_obj = fss_SRSF2_10,
  assay_name = "CITE",
  title_suffix = "SRSF2_10"
)

# View the density plot for SRSF2_10
print(cite_data_SRSF2_10$plot)

# Apply the threshold for SRSF2_10
fss_SRSF2_10 <- apply_cite_threshold(
  seurat_obj = fss_SRSF2_10,
  assay_name = "CITE",
  hash_threshold = 1 # Determined by the above density plot
)

# Subset on the dex+ (i.e. Cite antibody negative) cells

# For SRSF2_9
Idents(fss_SRSF2_9) <- "manual_hash_dmux"
fss_SRSF2_9_dex <- subset(fss_SRSF2_9, idents = "Negative")

# For SRSF2_10
Idents(fss_SRSF2_10) <- "manual_hash_dmux"
fss_SRSF2_10_dex <- subset(fss_SRSF2_10, idents = "Negative")

# ==========================================================
# 4. Comparison of clonotypes in Dex+ and Dex- cells
# ==========================================================

# Apply the function to both samples
plot_SRSF2_9 <- plot_clonotype_proportions(fss_SRSF2_9)
plot_SRSF2_10 <- plot_clonotype_proportions(fss_SRSF2_10)

# Check out the plots
print(plot_SRSF2_9)
print(plot_SRSF2_10)

### Here we can see that for SRSF2_9 clonotype 16 and 17 are enriched in 
### Dex+ cells indicating that they recognize specific neoantigens.

# ----------------------------------------------------------
# 5. Dextramer neoantigen peptide analysis
# ----------------------------------------------------------

# Perform z-score analysis on the dex+ cells to identify 
# predicted peptide clonotype pairs

result_SRSF2_9 <- perform_z_score_analysis(fss_SRSF2_9_dex)
result_SRSF2_10 <- perform_z_score_analysis(fss_SRSF2_10_dex)

## Here we can see that for SRSF2_9 clonotype 16 is predicted to 
## recognize the peptide "SRSF2-31" (RHOT2-5) with a high z-score.

# ==========================================================
# 6. Dimension reduction and UMAP plotting
# ==========================================================

fss <- FindNeighbors(fss, reduction = "integrated.cca", dims = 1:30)
fss <- FindClusters(fss, resolution = 0.5)

# UMAP
fss <- RunUMAP(fss, dims = 1:30, reduction = "integrated.cca")

# Define palettes
pal1 <- c("#A2D47C", "#C1D375", "#317FE5", "#13741F", "#396185", "#7FB285",
          "#FFD166", "#E09540", "#6C91C2", "#805D93", "#BFC2C6", "#FFA9AD",
          "#FFD7D5", "#699684")

pal2 <- c("#699684", "#FFD166", "#7FB285", "#317FE5", "#6C91C2", "#805D93",
          "#E09540", "#396185", "#FFA9AD", "#396185")

# UMAP plotting
DimPlot(fss, group.by="ident", cols = pal1, reduction = "umap", pt.size = 0.1, 
        label.size = 8, label = TRUE)

# Remove contaminating populations
FeaturePlot(fss, features = c("CD4", "ITGAM", "PF4", "CD8A"))

fss <- subset(fss, subset = seurat_clusters %in% c(4, 9, 12), invert = TRUE)

fss <- FindNeighbors(fss, reduction = "integrated.cca", dims = 1:30)
fss <- FindClusters(fss, resolution = 0.5)

fss <- RunUMAP(fss, dims = 1:30, reduction = "integrated.cca")

DimPlot(fss, group.by="ident", cols = pal2, reduction = "umap", pt.size = 0.1, 
        label.size = 8, label = TRUE)


DefaultAssay(fss) <- 'CITE'

# ==========================================================
# 6. Clonotype 16 plotting
# ==========================================================

# Define TCR directories
tcr_dir_9 <- dirs_SRSF2_9$dir_TCR
tcr_dir_10 <- dirs_SRSF2_10$dir_TCR

# Integrate TCR data into fss
fss <- integrate_tcr_into_fss(fss, tcr_dir_9, tcr_dir_10, threshold_umis = 3)

# Add sample ID prefix to clonotype IDs
fss@meta.data$raw_clonotype_id <- ifelse(
  grepl("^SRSF2_9", rownames(fss@meta.data)),
  paste0("SRSF2_9_", fss@meta.data$raw_clonotype_id),
  paste0("SRSF2_10_", fss@meta.data$raw_clonotype_id)
)

# Add binary indicator for clonotype 16
fss[["clonotype16"]] <- ifelse(
  fss@meta.data$raw_clonotype_id == "SRSF2_9_clonotype16",
  1, 0
)

# Plot UMAP highlighting clonotype 16
umap_plot <- generate_clonotype_umap_plot(fss, "clonotype16")
print(umap_plot)