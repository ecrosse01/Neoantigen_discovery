# ==========================================================
# Data Loading Helper Functions
# ==========================================================

## This script contains helper functions for loading and integrating the expression,
## dextramer, CITE-seq, and TCR data 
## It also QCs and filters the TCR data 

# 1. Load RNA Expression Data
load_expression_data <- function(data_dir, project_name = "dextramer_pilot") {
  df <- Read10X(data.dir = file.path(data_dir, "filtered_feature_bc_matrix/"))
  CreateSeuratObject(counts = df, project = project_name)
}

# 2. Load and Integrate Dextramer Data
load_dextramer_data <- function(data_dir, seurat_obj) {
  df_dex <- Read10X(data.dir = data_dir)
  dex_assay <- CreateAssayObject(counts = df_dex)
  
  # Add "-1" suffix to barcodes
  new_colnames <- paste0(colnames(dex_assay), "-1")
  colnames(dex_assay@counts) <- new_colnames
  colnames(dex_assay@data) <- new_colnames
  
  # Find common barcodes
  common_barcodes <- intersect(rownames(seurat_obj@meta.data), colnames(dex_assay@counts))
  
  # Subset both objects for common barcodes
  seurat_obj_common <- subset(seurat_obj, cells = common_barcodes)
  dex_assay_common <- subset(dex_assay, cells = common_barcodes)
  
  # Add dextramer assay to Seurat object
  seurat_obj_common[["DCODE"]] <- dex_assay_common
  
  return(seurat_obj_common)
}

# 3. Load and Integrate CITE-seq Data
load_cite_data <- function(data_dir, seurat_obj) {
  df_cite <- Read10X(data.dir = data_dir)
  cite_assay <- CreateAssayObject(counts = df_cite)
  
  # Add "-1" suffix to barcodes
  new_colnames <- paste0(colnames(cite_assay), "-1")
  colnames(cite_assay@counts) <- new_colnames
  colnames(cite_assay@data) <- new_colnames
  
  # Find common barcodes
  common_barcodes <- intersect(rownames(seurat_obj@meta.data), colnames(cite_assay@counts))
  
  # Subset both objects for common barcodes
  seurat_obj_common <- subset(seurat_obj, cells = common_barcodes)
  cite_assay_common <- subset(cite_assay, cells = common_barcodes)
  
  # Add CITE-seq assay to Seurat object
  seurat_obj_common[["CITE"]] <- cite_assay_common
  
  return(seurat_obj_common)
}

# 4. Load and Process TCR Data
load_tcr_data <- function(tcr_dir, seurat_obj, threshold_umis = 3) {
  # Read TCR annotations
  VDJ <- read_csv(file.path(tcr_dir, "filtered_contig_annotations.csv"))
  tcr <- VDJ %>%
    filter(
      high_confidence == TRUE,
      chain %in% c("TRA", "TRB"),
      productive == TRUE,
      umis >= threshold_umis
    ) %>%
    select(barcode, chain, raw_clonotype_id, v_gene, d_gene, j_gene, c_gene, cdr1_nt, cdr2_nt, cdr3_nt) %>%
    group_by(raw_clonotype_id) %>%
    filter(any(chain == 'TRA') & any(chain == 'TRB')) %>%
    ungroup()
  
  # Read clonotypes
  clono <- read_csv(file.path(tcr_dir, "clonotypes.csv")) %>%
    select(clonotype_id, cdr3s_aa, cdr3s_nt) %>%
    rename(raw_clonotype_id = clonotype_id)
  
  # Merge TCR and clonotype data
  tcr <- tcr %>% left_join(clono, by = "raw_clonotype_id")
  
  # Create metadata for clonotypes
  clonotype_info <- tcr %>%
    select(barcode, raw_clonotype_id) %>%
    distinct() %>%
    column_to_rownames("barcode")
  
  # Add clonotype metadata to Seurat object
  seurat_obj <- AddMetaData(seurat_obj, metadata = clonotype_info)
  
  return(seurat_obj)
}
