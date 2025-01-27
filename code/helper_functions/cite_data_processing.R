# Manual setting of the CITE hashing antibody thresholds to identify Dex+ (AB-) cells 

# Function to generate and return the density plot
generate_cite_density_plot <- function(seurat_obj, assay_name = "CITE", title_suffix) {
  # Extract data for the CITE assay
  data <- GetAssayData(seurat_obj, assay = assay_name, slot = "data")
  
  # Calculate total reads
  total_reads <- Matrix::colSums(data)
  
  # Create the density plot
  density_plot <- ggplot(mapping = aes(x = total_reads)) + 
    geom_density(fill = "blue", alpha = 0.5) +
    labs(
      x = "Total Reads",
      y = "Density",
      title = paste("Density Plot of Reads for CITE Assay (", title_suffix, ")", sep = "")
    ) +
    theme_minimal()
  
  # Return the density plot for inspection
  return(list(plot = density_plot, total_reads = total_reads))
}


# Function to process CITE-seq data with a given hash threshold
apply_cite_threshold <- function(seurat_obj, assay_name = "CITE", hash_threshold) {
  # Extract 'Hash' counts for each cell
  hash_counts <- seurat_obj@assays[[assay_name]]@data["Hash", ]
  
  # Add these counts to the metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = hash_counts, col.name = "HashCounts")
  
  # Add a new column based on the hash threshold
  seurat_obj@meta.data$manual_hash_dmux <- ifelse(
    seurat_obj@meta.data$HashCounts > hash_threshold, "CITE_hash", "Negative"
  )
  
  # Verify the metadata addition
  print(head(seurat_obj@meta.data))
  print(table(seurat_obj@meta.data$manual_hash_dmux))
  
  # Return the updated Seurat object
  return(seurat_obj)
}
