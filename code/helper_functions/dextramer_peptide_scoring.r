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