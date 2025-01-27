# set directory paths

set_directories <- function(patient_id, base_dir) {
  list(
    dir_gex = file.path(base_dir, paste0(patient_id, "/CellRangerGex_results")),
    dir_dex = file.path(base_dir, paste0(patient_id, "_dextramer_count/umi_count")),
    dir_TCR = file.path(base_dir, paste0(patient_id, "_TCR_VDJ/CellRangerVdj_results")),
    dir_CITE = file.path(base_dir, paste0(patient_id, "_hash_count/umi_count"))
  )
}
