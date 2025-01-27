## PLotting functions

# Function to plot clonotype proportions for a given sample
plot_clonotype_proportions <- function(fss_sample) {
  # Subset fss to keep only cells labeled as "CITE_hash"
  fss_CITE_hash <- subset(fss_sample, subset = manual_hash_dmux == "CITE_hash")
  
  # Subset fss to keep only cells labeled as "Negative"
  fss_Negative <- subset(fss_sample, subset = manual_hash_dmux == "Negative")
  
  # Calculate clonotype proportions for Negative
  clonotype_counts_neg <- table(fss_Negative@meta.data$raw_clonotype_id)
  clonotype_df_neg <- as.data.frame(clonotype_counts_neg)
  colnames(clonotype_df_neg) <- c("clonotype", "count")
  clonotype_df_neg$group <- "Negative"
  clonotype_df_neg$proportion <- clonotype_df_neg$count / sum(clonotype_df_neg$count)
  
  # Calculate clonotype proportions for CITE_hash
  clonotype_counts_cite <- table(fss_CITE_hash@meta.data$raw_clonotype_id)
  clonotype_df_cite <- as.data.frame(clonotype_counts_cite)
  colnames(clonotype_df_cite) <- c("clonotype", "count")
  clonotype_df_cite$group <- "CITE_hash"
  clonotype_df_cite$proportion <- clonotype_df_cite$count / sum(clonotype_df_cite$count)
  
  # Combine the two dataframes
  combined_df <- rbind(clonotype_df_neg, clonotype_df_cite)
  
  # Filter for clonotypes 1 to 50
  combined_df <- combined_df[combined_df$clonotype %in% paste0("clonotype", 1:50),]
  
  # Extract just the numeric part of the clonotype for plotting
  combined_df$clonotype_num <- as.numeric(gsub("clonotype", "", combined_df$clonotype))
  
  # Plotting the proportions
  ggplot(combined_df, aes(x = clonotype_num, y = proportion, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c("Negative" = "blue", "CITE_hash" = "red")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Clonotype ID number") +
    ylab("Proportion of total cells") +
    ggtitle("Proportion of Clonotypes 1 to 50 in Negative vs. CITE_hash") +
    scale_x_continuous(breaks = 1:50) # Display numbers 1-50 on the x-axis
}

