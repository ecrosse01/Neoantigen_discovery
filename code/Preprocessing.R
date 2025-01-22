# ==========================================================
# 1. Overview and Patient Information
# ==========================================================

# Overview
# This pipeline performs multimodal single-cell analysis of patient data including:
# - Transcriptomics
# - TCR clonotype profiling
# - Dextramer dCODEs
# - Hashing CITE-seq to determine Dex+ and Dex- fractions

# Patient Information
# SRSF2_9 Fresh (MDS-EB1; 15 months post-transplant): 7.3K Dex(+) + 8.7K Dex(-) cells
# SRSF2_10 Fresh (MDS-EB1; pre-transplant, s/p 2 cycles DEC): 5K Dex(+) + 11K Dex(-) cells

# Dextramer counts:
# - 46 Test
# - 1 Positive
# - 3 Negative

# ==========================================================
# 2. Setup: Required Libraries and Directory Paths
# ==========================================================

# Required Libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(scRepertoire)
library(gridExtra)
library(ggrepel)
library(cowplot)

# Directory Paths
## SRSF2_9 Paths
dir_gex <- "/path/to/SRSF2_9/CellRangerGex_results"
dir_dex <- "/path/to/SRSF2_9_dextramer_count/umi_count"
dir_TCR <- "/path/to/SRSF2_9_TCR_VDJ/CellRangerVdj_results"
dir_CITE <- "/path/to/SRSF2_9_hash_count/umi_count"

## SRSF2_10 Paths
dir_gex2 <- "/path/to/SRSF2_10/CellRangerGex_results"
dir_dex2 <- "/path/to/SRSF2_10_dextramer_count/umi_count"
dir_TCR2 <- "/path/to/SRSF2_10_TCR_VDJ/CellRangerVdj_results"
dir_CITE2 <- "/path/to/SRSF2_10_hash_count/umi_count"

# Define output directories
local_home <- "/path/to/output/data"
plots_dir <- "/path/to/output/plots"

# ==========================================================
# 3. Data Loading and Initialization
# ==========================================================

# Load Expression Data
df <- Read10X(data.dir = file.path(dir_gex, "filtered_feature_bc_matrix/"))
df2 <- Read10X(data.dir = file.path(dir_gex2, "filtered_feature_bc_matrix/"))

# Initialize Seurat Objects
adata <- CreateSeuratObject(counts = df, project = "dextramer_pilot")
adata2 <- CreateSeuratObject(counts = df2, project = "dextramer_pilot")

# ==========================================================
# 4. CITE-seq Data Integration
# ==========================================================

# Load and process SRSF2_9 CITE-seq data
df <- Read10X(data.dir = dir_CITE)
CITE_assay <- CreateAssayObject(counts = df)
new_colnames_CITE <- paste0(colnames(CITE_assay), "-1")
colnames(CITE_assay@counts) <- new_colnames_CITE
colnames(CITE_assay@data) <- new_colnames_CITE

common_barcodes <- intersect(rownames(adata@meta.data), colnames(CITE_assay@counts))
adata <- subset(adata, cells = common_barcodes)
CITE_assay_common <- subset(CITE_assay, cells = common_barcodes)
adata[["CITE"]] <- CITE_assay_common

# Repeat for SRSF2_10 CITE-seq data
df <- Read10X(data.dir = dir_CITE2)
CITE_assay <- CreateAssayObject(counts = df)
new_colnames_CITE <- paste0(colnames(CITE_assay), "-1")
colnames(CITE_assay@counts) <- new_colnames_CITE
colnames(CITE_assay@data) <- new_colnames_CITE

common_barcodes <- intersect(rownames(adata2@meta.data), colnames(CITE_assay@counts))
adata2 <- subset(adata2, cells = common_barcodes)
CITE_assay_common <- subset(CITE_assay, cells = common_barcodes)
adata2[["CITE"]] <- CITE_assay_common

# Merge Seurat Objects
fss <- merge(adata, y = adata2, add.cell.ids = c("SRSF2_9", "SRSF2_10"), project = "MDS_pre_post_transplant")

# ==========================================================
# 5. Quality Control and Normalization
# ==========================================================

# Calculate mitochondrial ratio
fss$mitoRatio <- PercentageFeatureSet(object = fss, pattern = "^MT-")
fss$mitoRatio <- fss@meta.data$mitoRatio / 100

# Filter low-quality cells
fss <- subset(x = fss, 
              subset = (nCount_RNA >= 500) & 
                (nFeature_RNA >= 250) & 
                (mitoRatio < 0.10))

# Normalize and identify variable features
fss <- NormalizeData(fss)
fss <- FindVariableFeatures(fss, selection.method = "vst", nfeatures = 4000)

# Exclude TCR genes from variable features
VFs <- VariableFeatures(fss)
VFs_filtered <- VFs[!grepl("^TR", VFs)]
VariableFeatures(fss) <- VFs_filtered

# ==========================================================
# 6. Dimensionality Reduction and Clustering
# ==========================================================

# Scaling, PCA, and Integration
all.genes <- rownames(fss)
fss <- ScaleData(fss, features = all.genes)
fss <- RunPCA(fss)

fss <- IntegrateLayers(object = fss, method = CCAIntegration, 
                       orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)

# Clustering and UMAP
fss <- FindNeighbors(fss, reduction = "integrated.cca", dims = 1:30)
fss <- FindClusters(fss, resolution = 0.5)
fss <- RunUMAP(fss, dims = 1:30, reduction = "integrated.cca")

# ==========================================================
# 7. Visualization 
# ==========================================================

# UMAP Visualization
pal1 <- c("#A2D47C", "#C1D375", "#317FE5", "#13741F", "#396185", "#7FB285",
          "#FFD166", "#E09540", "#6C91C2", "#805D93", "#BFC2C6", "#FFA9AD",
          "#FFD7D5", "#699684")
DimPlot(fss, group.by = "ident", cols = pal1, reduction = "umap", pt.size = 0.1, label.size = 8, label = TRUE)

# Save Processed Object
saveRDS(fss, file = file.path(local_home, "SRSF2_9_10_merge_seurat_processed.rds"))

# ==========================================================
# 8. Removing Contaminating Cells
# ==========================================================

FeaturePlot(fss, features = c("CD4", "ITGAM", "PF4", "CD8A"))
fss <- subset(fss, subset = seurat_clusters %in% c(4, 9, 12), invert = TRUE)

# Recompute neighbors, clusters, and UMAP
fss <- FindNeighbors(fss, reduction = "integrated.cca", dims = 1:30)
fss <- FindClusters(fss, resolution = 0.5)
fss <- RunUMAP(fss, dims = 1:30, reduction = "integrated.cca")

DimPlot(fss, group.by = "ident", cols = pal1, reduction = "umap", pt.size = 0.1, label.size = 8, label = TRUE)

# Save Updated Object
saveRDS(fss, file = file.path(local_home, "SRSF2_9_10_merge_seurat_processed_SUBSET.rds"))

# ==========================================================
# 9. CITE-seq Analysis
# ==========================================================

# CLR Normalization
DefaultAssay(fss) <- 'CITE'
fss <- NormalizeData(fss, assay = "CITE", normalization.method = "CLR")

# Density Plot
data <- GetAssayData(fss, assay = "CITE", slot = "data")
total_reads <- Matrix::colSums(data)
ggplot(mapping = aes(x = total_reads)) + 
  geom_density(fill = "blue", alpha = 0.5) +
  labs(x = "Total Reads", y = "Density", title = "Density Plot of Reads for CITE Assay") +
  theme_minimal()

# Hash Counts
hash_counts <- fss@assays$CITE@data["C0251-GTCAACTCTTTAGCG", ]
fss <- AddMetaData(fss, metadata = hash_counts, col.name = "HashCounts")
fss@meta.data$manual_hash_dmux <- ifelse(fss@meta.data$HashCounts > 1, "CITE_hash", "Negative")
FeaturePlot(fss, features = "Hash")

# Save Final Object
saveRDS(fss, file = file.path(local_home, "SRSF2_9_10_final.rds"))
