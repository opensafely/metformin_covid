################################################################################
# Custom made function to assess the completeness criteria
################################################################################
fn_completeness_criteria_midpoint6 <- function(data_processed){
  n_before_exclusion_processing <-
    data_processed %>%
    nrow()
  n_not_alive_at_landmark <-
    data_processed %>%
    filter(qa_bin_was_alive == FALSE) %>%
    nrow()
  n_nor_female_or_male <-
    data_processed %>%
    filter(qa_bin_is_female_or_male == FALSE) %>%
    nrow()
  n_no_imd <-
    data_processed %>%
    filter(qa_bin_known_imd == FALSE) %>%
    nrow()
  n_no_region <-
    data_processed %>%
    filter(is.na(cov_cat_region)) %>%
    nrow()
  n_not_registered <-
    data_processed %>%
    filter(qa_bin_was_registered == FALSE) %>%
    nrow()
  
  n_after_exclusion_processing <-
    data_processed %>%
    filter(qa_bin_was_alive == TRUE) %>%
    filter(qa_bin_is_female_or_male == TRUE) %>%
    filter(qa_bin_known_imd == TRUE) %>%
    filter(!is.na(cov_cat_region)) %>%
    filter(qa_bin_was_registered == TRUE) %>%
    nrow()
  
  # Output incl. redaction
  out_midpoint6 <- tibble(
      n_before_exclusion_processing_midpoint6 = fn_roundmid_any(n_before_exclusion_processing, threshold),
      n_not_alive_at_landmark_midpoint6 = fn_roundmid_any(n_not_alive_at_landmark, threshold),
      n_nor_female_or_male_midpoint6 = fn_roundmid_any(n_nor_female_or_male, threshold),
      n_no_imd_midpoint6 = fn_roundmid_any(n_no_imd, threshold),
      n_no_region_midpoint6 = fn_roundmid_any(n_no_region, threshold),
      n_not_registered_midpoint6 = fn_roundmid_any(n_not_registered, threshold),
      n_after_exclusion_processing_midpoint6 = fn_roundmid_any(n_after_exclusion_processing, threshold),
    )
  
  return(out_midpoint6)
}