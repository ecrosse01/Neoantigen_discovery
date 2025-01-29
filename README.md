# MDS Patient Peptide-TCR Analysis Pipeline

**Predicting Immunogenic Peptides and Their TCR Interactions in MDS Using Single-Cell Multimodal Data**

## Overview

This repository contains a single-cell multimodal analysis pipeline for studying pre- and post-transplant Myelodysplastic Syndrome (MDS) patient samples. The pipeline integrates:

- **RNA sequencing** (Gene expression analysis)
- **TCR repertoire** (Clonotype identification and tracking)
- **Dextramer binding** (Neoantigen-specific immune response detection)
- **CITE-seq hashing** (Surface protein profiling)

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

