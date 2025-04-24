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
[**"Mis-splicing-derived neoantigens and cognate TCRs in splicing factor mutant leukemias"**](https://doi.org/10.1016/j.cell.2025.03.047)  

**Authors**: WonJun Kim, Edie I. Crosse, Emma De Neef, Inaki Etxeberria, Erich Y. Sabio, Eric Wang, Jan Philipp Bewersdorf, Kuan-Ting Lin, Sydney X. Lu, Andrea Belleville, Nina Fox, Cynthia Castro, Pu Zhang, Takeshi Fujino, Jennifer Lewis, Jahan Rahman, Beatrice Zhang, Jacob H. Winick, Alexander M. Lewis, Robert F. Stanley, Susan DeWolf, Brigita Meskauskaite Urben, Meril Takizawa, Tobias Krause, Henrik Molina, Ronan Chaligne, Priya Koppikar, Jeffrey Molldrem, Mathieu Gigoux, Taha Merghoub, Anthony Daniyan, Smita S. Chandran, Benjamin D. Greenbaum, Christopher A. Klebanoff, Robert K. Bradley, and Omar Abdel-Wahab

**Abstract**: Mutations in RNA splicing factors are prevalent across cancers and generate recurrently mis-spliced mRNA isoforms. Here, we identified a series of bona fide neoantigens translated from highly stereotyped splicing alterations promoted by neomorphic, leukemia-associated somatic splicing machinery mutations. We utilized feature-barcoded peptide-major histocompatibility complex (MHC) dextramers to isolate neoantigen-reactive T cell receptors (TCRs) from healthy donors, patients with active myeloid malignancy, and following curative allogeneic stem cell transplant. Neoantigen-reactive CD8+ T cells were present in the blood of patients with active cancer and had a distinct phenotype from virus-reactive T cells with evidence of impaired cytotoxic function. T cells engineered with TCRs recognizing SRSF2 mutant-induced neoantigens arising from mis-splicing events in CLK3 and RHOT2 resulted in specific recognition and cytotoxicity of SRSF2-mutant leukemia. These data identify recurrent RNA mis-splicing events as sources of actionable public neoantigens in myeloid leukemias and provide proof of concept for genetically redirecting T cells to recognize these targets. 

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

The raw sequencing data is available at **NCBI GEO:GSE268157**:  
**[Download Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268157)**

#### Download the following sample files:

| Sample ID    | Description                             |
|--------------|-----------------------------------------|
| GSM8286626   | WJK-2859_SRSF2_9_RNA (SRSF2 Mutant)     |
| GSM8286627   | WJK-2859_SRSF2_9_CITE (SRSF2 Mutant)    |
| GSM8286628   | WJK-2859_SRSF2_9_TCR_VDJ (SRSF2 Mutant) |
| GSM8286629   | WJK-2864_SRSF2_10_RNA (SRSF2 Mutant)    |
| GSM8286630   | WJK-2864_SRSF2_10_CITE (SRSF2 Mutant)   |
| GSM8286631   | WJK-2864_SRSF2_10_TCR_VDJ (SRSF2 Mutant)|

SRSF2_9 = MDS patient 1, 15 months post-transplant
SRSF2_10 = MDS patient 1, pre-transplant

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

### Predict which peptides are neoantigens and the TCRs that bind them

First we will load and integrate the dextramer and TCR data for the two samples

```r

# Load dextramer and TCR data for SRSF2_9

fss_SRSF2_9 <- load_dextramer_data(
  data_dir = dirs_SRSF2_9$dir_dex,
  seurat_obj = fss_SRSF2_9,
  sample_prefix = "SRSF2_9"
)

fss_SRSF2_9 <- load_tcr_data(
  tcr_dir = dirs_SRSF2_9$dir_TCR,
  seurat_obj = fss_SRSF2_9,
  sample_prefix = "SRSF2_9"
)

# Load dextramer and TCR data for SRSF2_10
# Load dextramer data for SRSF2_10
fss_SRSF2_10 <- load_dextramer_data(
  data_dir = dirs_SRSF2_10$dir_dex,
  seurat_obj = fss_SRSF2_10,
  sample_prefix = "SRSF2_10"
)

fss_SRSF2_10 <- load_tcr_data(
  tcr_dir = dirs_SRSF2_10$dir_TCR,
  seurat_obj = fss_SRSF2_10,
  sample_prefix = "SRSF2_10"
)

```

Next create new seurat objects with just the dex+ cells using the previously defined CITE-Seq negative populations

```r

# Subset on the dex+ cells (i.e. CITE_hash negative)

#For SRSF2_9
Idents(fss_SRSF2_9) <- "manual_hash_dmux"
fss_SRSF2_9_dex <- subset(fss_SRSF2_9, idents = "Negative")

# For SRSF2_10
Idents(fss_SRSF2_10) <- "manual_hash_dmux"
fss_SRSF2_10_dex <- subset(fss_SRSF2_10, idents = "Negative")

```

Plot the TCR clonotype representations in the dextramer negative and positive populations to identify overrepresented clonotypes in the dextramer positive cells.

```r

# Apply the function to both samples
plot_SRSF2_9 <- plot_clonotype_proportions(fss_SRSF2_9)
plot_SRSF2_10 <- plot_clonotype_proportions(fss_SRSF2_10)

# Check out the plots
print(plot_SRSF2_9)
print(plot_SRSF2_10)

```

You can see that clonotypes 16 and 17 are overrepresented in the Dex+ cells in sample SRSF2_9 predictive of binding peptides.

The next functions before z-scoring analysis to predict TCR / neoantigen pairs in the Dex+ cells

```r

result_SRSF2_9 <- perform_z_score_analysis(fss_SRSF2_9_dex)
result_SRSF2_10 <- perform_z_score_analysis(fss_SRSF2_10_dex)

```

Here we can see that for SRSF2_9 clonotype 16 is predicted to recognize the peptide "SRSF2-31" (RHOT5-1).

### Visualize the integrated data with a UMAP and map on clonotype 16

These steps are for dimension reduction, UMAP visualization and filtering of contaminating populations

```r

fss <- FindNeighbors(fss, reduction = "integrated.cca", dims = 1:30)
fss <- FindClusters(fss, resolution = 0.5)

# UMAP
fss <- RunUMAP(fss, dims = 1:30, reduction = "integrated.cca")

# Define palette
pal1 <- c("#A2D47C", "#C1D375", "#317FE5", "#13741F", "#396185", "#7FB285",
          "#FFD166", "#E09540", "#6C91C2", "#805D93", "#BFC2C6", "#FFA9AD",
          "#FFD7D5", "#699684")

pal2 <- c("#699684", "#FFD166", "#7FB285", "#317FE5", "#6C91C2", "#805D93",
          "#E09540", "#396185", "#FFA9AD", "#396185")

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


```

Finally, we plot clonotype 16

```r

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

```




