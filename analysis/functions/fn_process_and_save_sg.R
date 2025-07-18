################################################################################
# Custom made function to process and save subgroup subsets and redefine cox end date in each subset
################################################################################
fn_process_and_save_sg <- function(name, subgroup, df, studyend_date) {
  subset <- df %>%
    filter(!!subgroup) %>%
    mutate(cox_date_severecovid = pmin(
      out_date_severecovid_afterlandmark,
      out_date_noncoviddeath_afterlandmark,
      cens_date_ltfu_afterlandmark,
      max_fup_date,
      studyend_date,
      na.rm = TRUE
    ))
  arrow::write_feather(
    subset,
    here::here("output", "data", "sg", paste0(name, ".arrow"))
  )
  rm(subset)
}