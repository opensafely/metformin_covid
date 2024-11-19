################################################################################
# Custom made function to assess the completeness criteria
################################################################################
fn_completeness_criteria <- function(data_processed){
  n_before_exclusion_processing <-
    data_processed %>%
    nrow()
  n_is_alive <-
    data_processed %>%
    filter(qa_bin_was_alive == TRUE) %>%
    nrow()
  n_is_female_or_male <-
    data_processed %>%
    filter(qa_bin_is_female_or_male == TRUE) %>%
    nrow()
  n_has_imd <-
    data_processed %>%
    filter(qa_bin_known_imd == TRUE) %>%
    nrow()
  n_has_region <-
    data_processed %>%
    filter(!is.na(cov_cat_region)) %>%
    nrow()
  n_is_registered <-
    data_processed %>%
    filter(qa_bin_was_registered == TRUE) %>%
    nrow()
  
  n_after_exclusion_processing <-
    data_processed %>%
    filter(qa_bin_was_alive == TRUE) %>%
    filter(qa_bin_is_female_or_male == TRUE) %>%
    filter(qa_bin_known_imd == TRUE) %>%
    filter(!is.na(cov_cat_region)) %>%
    filter(qa_bin_was_registered == TRUE) %>%
    nrow()
  
  out <- tibble(n_before_exclusion_processing,
                n_is_alive,
                n_is_female_or_male,
                n_has_imd,
                n_has_region,
                n_is_registered,
                n_after_exclusion_processing
  )
  
}