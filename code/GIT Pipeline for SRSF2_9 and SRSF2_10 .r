# ==========================================================
# MDS Patient Single-cell Multimodal Analysis Pipeline
# ==========================================================

# ==========================================================
# 1. Overview Section
# ==========================================================

# ----------------------------------------------------------
# Title and Description
# ----------------------------------------------------------
# Analysis of pre/post-transplant MDS patient samples using:
# - RNA sequencing
# - TCR repertoire
# - Dextramer binding
# - CITE-seq hashing

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
# - Characterize immune changes pre- and post-transplant.
# - Identify clonotypes linked to dextramer binding.
# - Integrate multimodal data for comprehensive insights.

# ==========================================================
# 2. Setup Section
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
# Directory Paths
# ----------------------------------------------------------
#set directory variables SRSF2_9
dir_gex="/Users/ecrosse/Desktop/data_for_edie_third_batch_january/WJK-2859_SRSF2_9/CellRangerGex_results"
dir_dex="/Users/ecrosse/Desktop/data_for_edie_third_batch_january/WJK-2859_SRSF2_9_dextramer_count/umi_count"
dir_TCR="/Users/ecrosse/Desktop/data_for_edie_third_batch_january/WJK-2859_SRSF2_9_TCR_VDJ/CellRangerVdj_results"
dir_CITE="/Users/ecrosse/Desktop/data_for_edie_third_batch_january/WJK-2859_SRSF2_9_hash_count/umi_count/"
#set directory variables SRSF2_10
dir_gex2="/Users/ecrosse/Desktop/dextramer_data_for_edie_january_part_2/WJK-2864_SRSF2_10/CellRangerGex_results"
dir_dex2="/Users/ecrosse/Desktop/dextramer_data_for_edie_january_part_2/WJK-2864_SRSF2_10_dextramer_count/umi_count"
dir_TCR2="/Users/ecrosse/Desktop/dextramer_data_for_edie_january_part_2/WJK-2864_SRSF2_10_TCR_VDJ/CellRangerVdj_results/"
dir_CITE2="/Users/ecrosse/Desktop/dextramer_data_for_edie_january_part_2/WJK-2864_SRSF2_10_hash_count/umi_count/"

local_home="/Users/ecrosse/projects/gene_expression/data/human/2023/abdel_wahab-bradley.dextramer/SRSF2"
plots_dir="/Users/ecrosse/OneDrive - Fred Hutchinson Cancer Research Center/1.Projects/data/2023/abdel-wahab-bradley.dextramers/SRSF2_9_10/Update_plots/"


# ----------------------------------------------------------
# Global Variables
# ----------------------------------------------------------
# Define global thresholds and parameters used throughout analysis.
threshold_umis <- 3
threshold_mito <- 0.10
threshold_nCount_RNA <- 500
threshold_nFeature_RNA <- 250

# ==========================================================
# 3. Data Loading Section
# ==========================================================

# ----------------------------------------------------------
# Expression Data
# ----------------------------------------------------------
df <- Read10X(data.dir = file.path(dir_gex, "filtered_feature_bc_matrix/"))
df2 <- Read10X(data.dir = file.path(dir_gex2, "filtered_feature_bc_matrix/"))

adata <- CreateSeuratObject(counts = df, project = "dextramer_pilot")
adata2 <- CreateSeuratObject(counts = df2, project = "dextramer_pilot")

# ----------------------------------------------------------
# Dextramer Data
# ----------------------------------------------------------
df_dex <- Read10X(data.dir = dir_dex)
dex_assay <- CreateAssayObject(counts = df_dex)
new_colnames_dex <- paste0(colnames(dex_assay), "-1")
colnames(dex_assay@counts) <- new_colnames_dex
colnames(dex_assay@data) <- new_colnames_dex
common_barcodes_dex <- intersect(rownames(adata@meta.data), colnames(dex_assay@counts))
adata_common_dex <- subset(adata, cells = common_barcodes_dex)
dex_assay_common <- subset(dex_assay, cells = common_barcodes_dex)
adata <- adata_common_dex
adata[["DCODE"]] <- dex_assay_common

df_dex2 <- Read10X(data.dir = dir_dex2)
dex_assay2 <- CreateAssayObject(counts = df_dex2)
new_colnames_dex2 <- paste0(colnames(dex_assay2), "-1")
colnames(dex_assay2@counts) <- new_colnames_dex2
colnames(dex_assay2@data) <- new_colnames_dex2
common_barcodes_dex2 <- intersect(rownames(adata2@meta.data), colnames(dex_assay2@counts))
adata2_common_dex <- subset(adata2, cells = common_barcodes_dex2)
dex_assay2_common <- subset(dex_assay2, cells = common_barcodes_dex2)
adata2 <- adata2_common_dex
adata2[["DCODE"]] <- dex_assay2_common

# ----------------------------------------------------------
# CITE-seq Data
# ----------------------------------------------------------
df_CITE <- Read10X(data.dir = dir_CITE)
CITE_assay <- CreateAssayObject(counts = df_CITE)
new_colnames_CITE <- paste0(colnames(CITE_assay), "-1")
colnames(CITE_assay@counts) <- new_colnames_CITE
colnames(CITE_assay@data) <- new_colnames_CITE
common_barcodes_CITE <- intersect(rownames(adata@meta.data), colnames(CITE_assay@counts))
adata_common_CITE <- subset(adata, cells = common_barcodes_CITE)
CITE_assay_common <- subset(CITE_assay, cells = common_barcodes_CITE)
adata <- adata_common_CITE
adata[["CITE"]] <- CITE_assay_common

df_CITE2 <- Read10X(data.dir = dir_CITE2)
CITE_assay2 <- CreateAssayObject(counts = df_CITE2)
new_colnames_CITE2 <- paste0(colnames(CITE_assay2), "-1")
colnames(CITE_assay2@counts) <- new_colnames_CITE2
colnames(CITE_assay2@data) <- new_colnames_CITE2
common_barcodes_CITE2 <- intersect(rownames(adata2@meta.data), colnames(CITE_assay2@counts))
adata2_common_CITE <- subset(adata2, cells = common_barcodes_CITE2)
CITE_assay2_common <- subset(CITE_assay2, cells = common_barcodes_CITE2)
adata2 <- adata2_common_CITE
adata2[["CITE"]] <- CITE_assay2_common

# ----------------------------------------------------------
# TCR Data
# ----------------------------------------------------------
VDJ <- read_csv(file.path(dir_TCR, "filtered_contig_annotations.csv"))
tcr <- VDJ %>% 
  filter(high_confidence == TRUE, 
         chain %in% c("TRA", "TRB"), 
         productive == TRUE,
         umis >= threshold_umis) %>%  
  select(barcode, chain, raw_clonotype_id, 
         v_gene, d_gene, j_gene, c_gene, cdr1_nt, 
         cdr2_nt, cdr3_nt) 

tcr <- tcr %>%
  group_by(raw_clonotype_id) %>%
  filter(any(chain == 'TRA') & any(chain == 'TRB')) %>%
  ungroup()

clono <- read_csv(file.path(dir_TCR, "clonotypes.csv"))
clono <- clono %>%
  select(clonotype_id, cdr3s_aa, cdr3s_nt) %>%
  rename(raw_clonotype_id = "clonotype_id")

tcr <- tcr %>% left_join(clono, by = "raw_clonotype_id")
clonotype_info <- tcr %>% 
  select(barcode, raw_clonotype_id) 
clonotype_info <- clonotype_info[!duplicated(clonotype_info$barcode), ]
clonotype_info <- as.data.frame(clonotype_info)
row.names(clonotype_info) <- clonotype_info$barcode
clonotype_info$barcode <- NULL
adata <- AddMetaData(adata, metadata = clonotype_info)

VDJ2 <- read_csv(file.path(dir_TCR2, "filtered_contig_annotations.csv"))
tcr2 <- VDJ2 %>% 
  filter(high_confidence == TRUE, 
         chain %in% c("TRA", "TRB"), 
         productive == TRUE,
         umis >= threshold_umis) %>%  
  select(barcode, chain, raw_clonotype_id, 
         v_gene, d_gene, j_gene, c_gene, cdr1_nt, 
         cdr2_nt, cdr3_nt) 

tcr2 <- tcr2 %>%
  group_by(raw_clonotype_id) %>%
  filter(any(chain == 'TRA') & any(chain == 'TRB')) %>%
  ungroup()

clono2 <- read_csv(file.path(dir_TCR2, "clonotypes.csv"))
clono2 <- clono2 %>%
  select(clonotype_id, cdr3s_aa, cdr3s_nt) %>%
  rename(raw_clonotype_id = "clonotype_id")

tcr2 <- tcr2 %>% left_join(clono2, by = "raw_clonotype_id")
clonotype_info2 <- tcr2 %>% 
  select(barcode, raw_clonotype_id) 
clonotype_info2 <- clonotype_info2[!duplicated(clonotype_info2$barcode), ]
clonotype_info2 <- as.data.frame(clonotype_info2)
row.names(clonotype_info2) <- clonotype_info2$barcode
clonotype_info2$barcode <- NULL
adata2 <- AddMetaData(adata2, metadata = clonotype_info2)

# ==========================================================
# 4. Data Processing Section
# ==========================================================

# ----------------------------------------------------------
# Seurat Object Creation
# ----------------------------------------------------------
fss <- merge(adata, y = adata2, add.cell.ids = c("SRSF2_9", "SRSF2_10"), project = "MDS_pre_post_transplant")

# ----------------------------------------------------------
# Data Filtering and Normalization
# ----------------------------------------------------------
fss$mitoRatio <- PercentageFeatureSet(object = fss, pattern = "^MT-") / 100
fss <- subset(x = fss, 
              subset = (nCount_RNA >= threshold_nCount_RNA) & 
                       (nFeature_RNA >= threshold_nFeature_RNA) & 
                       (mitoRatio < threshold_mito))
fss <- NormalizeData(fss)
VariableFeatures(fss) <- VariableFeatures(fss)[!grepl("^TR", VariableFeatures(fss))]

# ----------------------------------------------------------
# Integration Steps
# ----------------------------------------------------------
fss <- RunPCA(fss)
fss <- IntegrateLayers(object = fss, method = CCAIntegration, 
                       orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)

# ==========================================================
# 5. Analysis Section
# ==========================================================

# ----------------------------------------------------------
# CITE-seq Analysis
# ----------------------------------------------------------
fss <- NormalizeData(fss, assay = "CITE", normalization.method = "CLR")

# ----------------------------------------------------------
# Dextramer Analysis
# ----------------------------------------------------------
perform_z_score_analysis <- function(fss_dex) {
  DCODES <- rownames(fss_dex@assays$DCODE@counts)
  raw_counts <- fss_dex@assays$DCODE@counts %>% as.data.frame()
  raw_counts_long <- raw_counts %>%
    as_tibble(rownames = "dCODE") %>%
    gather(cell, raw_read, -dCODE)
  raw_counts_long <- raw_counts_long %>%
    mutate(log2_read = log2(raw_read + 1))
  combined_data <- raw_counts_long %>%
    inner_join(fss_dex@meta.data %>% as_tibble(rownames = "cell"), by = "cell")
  result <- purrr::map_df(c(3, 4, 5, 6, 7), function(threshold) {
    clono_counts <- combined_data %>%
      filter(log2_read >= threshold) %>%
      filter(!dCODE %in% c("SRSF2-18", "unmapped")) %>%
      group_by(raw_clonotype_id, dCODE) %>%
      summarise(Count = n_distinct(cell), .groups = 'drop')
    sum_counts <- clono_counts %>%
      group_by(raw_clonotype_id) %>%
      summarise(total_count = sum(Count), .groups = 'drop')
    merged_df <- left_join(clono_counts, sum_counts, by = "raw_clonotype_id") %>%
      filter(total_count > 10)
    merged_df <- merged_df %>%
      mutate(norm_count = Count / total_count * 1000)
    mean_std_df <- merged_df %>%
      group_by(dCODE) %>%
      summarise(mean = mean(norm_count, na.rm = TRUE), std = sd(norm_count, na.rm = TRUE), .groups = 'drop')
    merged_df <- left_join(merged_df, mean_std_df, by = "dCODE", suffix = c("", "_norm"))
    merged_df <- merged_df %>%
      mutate(z_score = (norm_count - mean) / std)
    overrepresented <- merged_df %>%
      filter(z_score > 1, Count > 10)
    overrepresented <- overrepresented %>% add_column(threshold = threshold)
    return(overrepresented)
  })
  result <- result %>% 
    filter(!str_detect(dCODE, "^[PNu]")) %>%
    arrange(desc(z_score))
  result <- result %>%
    mutate(method = "z_score stat test")
  return(result)
}
result_SRSF2_9 <- perform_z_score_analysis(fss_SRSF2_9)
result_SRSF2_10 <- perform_z_score_analysis(fss_SRSF2_10)

# ----------------------------------------------------------
# TCR Analysis
# ----------------------------------------------------------
# Perform TCR clonotype analysis and statistical tests.
result_SRSF2_9 <- perform_z_score_analysis(fss_SRSF2_9_dex)
result_SRSF2_10 <- perform_z_score_analysis(fss_SRSF2_10_dex)

## Sample Clonotype 16 is a hit for dCODE SRSF2-31 (peptide 
## RHOT2-5)

# ----------------------------------------------------------
# Comparison of clonotypes in Dex+ and Dex- cells
# ----------------------------------------------------------

