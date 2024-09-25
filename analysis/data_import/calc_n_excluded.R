################################################################################
# A custom made function to run the eligibility criteria
# Based on and Credits to https://github.com/opensafely/pax-non-users/tree/2dbf044472efdcfeb86f8fc2c8eea222e7eefe32
################################################################################
calc_n_excluded <- function(data_processed){
  n_before_exclusion_processing <-
    data_processed %>%
    nrow()
  # completeness criteria
  n_is_alive <-
    data_processed %>%
    filter(qa_bin_was_alive == TRUE & (qa_date_of_death > baseline_date | is.na(qa_date_of_death))) %>% # additional condition since "qa_bin_was_alive == TRUE" may not cover all (e.g. pos test came out after death)
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
  # inclusion criteria
  n_has_adult_age <-
    data_processed %>%
    filter(qa_bin_was_adult == TRUE) %>%
    nrow()
  n_has_t2dm <-
    data_processed %>%
    filter(cov_bin_t2dm == TRUE) %>%
    nrow()
  n_has_covid_infection <-
    data_processed %>%
    filter(!is.na(baseline_date)) %>%
    nrow()
  # exclusion criteria (coded positively)
  n_is_inhospital <-
    data_processed %>%
    filter(cov_bin_hosp_baseline == TRUE) %>%
    nrow()
  n_already_metfin <-
    data_processed %>%
    filter(cov_bin_metfin_before_baseline == TRUE) %>%
    nrow()
  n_metfin_allergy <-
    data_processed %>%
    filter(cov_bin_metfin_allergy == TRUE) %>%
    nrow()
  n_ckd_history <-
    data_processed %>%
    filter(cov_bin_ckd_45 == TRUE) %>%
    nrow()
  n_cliver_history <-
    data_processed %>%
    filter(cov_bin_liver_cirrhosis == TRUE) %>%
    nrow()
  n_intmefin <-
    data_processed %>%
    filter(cov_bin_metfin_interaction == TRUE) %>%
    nrow()
  n_longcovid_history <-
    data_processed %>%
    filter(cov_bin_long_covid == TRUE) %>%
    nrow()

  n_after_exclusion_processing <-
    data_processed %>%
    # completeness criteria
    filter(qa_bin_was_alive == TRUE) %>%
    filter(qa_bin_is_female_or_male == TRUE) %>%
    filter(qa_bin_known_imd == TRUE) %>%
    filter(!is.na(cov_cat_region)) %>%
    filter(qa_bin_was_registered == TRUE) %>%
    # inclusion criteria
    filter(qa_bin_was_adult == TRUE) %>%
    filter(cov_bin_t2dm == TRUE) %>%
    filter(!is.na(baseline_date)) %>%
    # exclusion criteria
    filter(cov_bin_hosp_baseline == FALSE) %>% # FALSE includes missing in an ehrQL logical boolean
    filter(cov_bin_metfin_before_baseline == FALSE) %>%
    filter(cov_bin_metfin_allergy == FALSE) %>%
    filter(cov_bin_ckd_45 == FALSE) %>%
    filter(cov_bin_liver_cirrhosis == FALSE) %>%
    filter(cov_bin_metfin_interaction == FALSE) %>%
    filter(cov_bin_long_covid == FALSE) %>%
    nrow()

  out <- tibble(n_before_exclusion_processing,
                # completeness criteria
                n_is_alive,
                n_is_female_or_male,
                n_has_imd,
                n_has_region,
                n_is_registered,
                # inclusion criteria
                n_has_adult_age,
                n_has_t2dm,
                n_has_covid_infection,
                # exclusion criteria
                n_is_inhospital,
                n_already_metfin,
                n_metfin_allergy,
                n_ckd_history,
                n_cliver_history,
                n_intmefin,
                n_longcovid_history,
                n_after_exclusion_processing
                )
}
