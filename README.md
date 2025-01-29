# MDS Patient Peptide-TCR Analysis Pipeline

**Predicting Immunogenic Peptides and Their TCR Interactions in MDS Using Single-Cell Multimodal Data**

## Overview

This repository contains a single-cell multimodal analysis pipeline for studying pre- and post-transplant Myelodysplastic Syndrome (MDS) patient samples. The pipeline integrates:

- **RNA sequencing** (Gene expression analysis)
- **TCR repertoire** (Clonotype identification and tracking)
- **Dextramer binding** (Neoantigen-specific immune response detection)
- **CITE-seq hashing** (Discriminate between Dextramer + and Dextramer - (CITE hashed) cells)

It enables the identification of clonotypes associated with dextramer binding and tracks immune changes during MDS progression and treatment.

---

## Relevant Publication

**Reference Paper**:  
[**"Title of Paper"**](https://doi.org/xxxxx)  
**Authors**: 
**Abstract**:  

---

## Repository Contents

**`code/neoantigen_pipeline.R`** – Main pipeline in full
**`code/helper_functions/`** – Helper functions for the pipeline
**`data/`** - Placeholder directory for NCBI downloaded data
**`README.md`** – Documentation for setup and execution  
**`LICENSE`** – Repository license  

---

## Prerequisites

Before running the pipeline, ensure you have the following installed:

### Required Software
- **R (\u22654.0)**
- **Seurat (\u22654.0)**
- **tidyverse**
- **scRepertoire**
- **ggplot2**
- **patchwork**
- **gridExtra**
- **ggrepel**
---

## Tutorial: Running the Analysis

### 1. Clone the Repository

```bash
git clone https://github.com/YourUsername/MDS-Multimodal-Analysis.git
cd MDS-Multimodal-Analysis
```

### 2. Download the Data

The raw sequencing data is available at **NCBI GEO**:  
**[Download Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXX)**

Once downloaded, place the data in the `data/` directory:

```bash
mkdir data
mv path_to_downloaded_data/* data/
```

### 3. Start the Analysis

#### Load libraries

```r
library(tidyverse)      # Data manipulation and visualization
library(Seurat)         # Single-cell analysis
library(patchwork)      # Plot composition
library(ggplot2)        # Plotting
library(RColorBrewer)   # Color palettes
library(scRepertoire)   # TCR analysis
library(gridExtra)      # Plot arrangement
library(ggrepel)        # Text labels
library(cowplot)        # Plot composition
```

#### Set global variables

```r
threshold_umis <- 3
threshold_mito <- 0.10
threshold_nCount_RNA <- 500
threshold_nFeature_RNA <- 250
```

#### Load helper functions and set directories

```r
# Load the helper functions
source("code/helper_functions/data_loading_helpers.R") # Different modality data loadings and integrations
source("code/helper_functions/cite_data_processing.R")  # Define the Dex+ population using cite-seq manual thresholding
source("code/helper_functions/plotting_functions.R") # Plotting functions
source("code/helper_functions/dextramer_peptide_scoring.R") # Z-score analysis for peptide clonotype pairs

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
base_dir <- "/Neoantigen_discovery/data/"
dirs_SRSF2_9 <- set_directories("path/to/dir", base_dir)
dirs_SRSF2_10 <- set_directories("path/to/dir", base_dir)

```
#### Load and integrate expression and CITE-Seq data modalaties

```r
# Load data for SRSF2_9 (gene expression and CITE-seq)
adata <- load_expression_data(dirs_SRSF2_9$dir_gex)
adata <- load_cite_data(dirs_SRSF2_9$dir_CITE, adata)

# Load data for SRSF2_10 (gene expression and CITE-seq)
adata2 <- load_expression_data(dirs_SRSF2_10$dir_gex)
adata2 <- load_cite_data(dirs_SRSF2_10$dir_CITE, adata2)

```
#### Pre-process the data

These steps will merge the samples, filter out poor quality cells, normalize and scale the data and run PCA. Finally, it uses CCA integration to fully integrate the two samples.

```r

# merge the two samples datasets
fss <- merge(adata, y = adata2, add.cell.ids = c("SRSF2_9", "SRSF2_10"), project = "MDS_pre_post_transplant")

# filter out poor quality cells
fss$mitoRatio <- PercentageFeatureSet(object = fss, pattern = "^MT-") / 100
fss <- subset(x = fss, 
              subset = (nCount_RNA >= threshold_nCount_RNA) & 
                       (nFeature_RNA >= threshold_nFeature_RNA) & 
                       (mitoRatio < threshold_mito))
fss <- NormalizeData(fss)
fss <- FindVariableFeatures(fss, selection.method = "vst", nfeatures = 4000)

VFs <- VariableFeatures(fss)

# Filter out genes starting with 'TR' (T cell receptor genes)
VFs_filtered <- VFs[!grepl("^TR", VFs)]

# Update the Seurat object with the filtered variable features list
VariableFeatures(fss) <- VFs_filtered

all.genes <- rownames(fss)
fss <- ScaleData(fss, features = all.genes)

fss <- RunPCA(fss)

fss <- IntegrateLayers(object = fss, method = CCAIntegration, 
                       orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)

# re-join layers after integration
fss[["RNA"]] <- JoinLayers(fss[["RNA"]])

```

#### CITE-seq Analysis - To Define the Dex+ Populations

This section uses manual thresholding of the CITE-Seq hashing antibody data to define the negative population which contains the Dextramer + cells.

```r
# Normalize CITE-seq data
fss <- NormalizeData(fss, assay = "CITE", normalization.method = "CLR")

# Add sample information to the metadata
fss$sample <- ifelse(grepl("SRSF2_9", colnames(fss)), "SRSF2_9", "SRSF2_10")

# Subset the fss object into two samples
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

```


