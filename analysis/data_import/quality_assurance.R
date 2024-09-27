################################################################################
# A custom made function to run the quality assurance criteria
################################################################################
quality_assurance <- function(data_processed){
  n_before_qa_processing <-
    data_processed %>%
    nrow()
  n_yob_missing <- # Rule 1: Year of birth is missing
    data_processed %>%
    filter(is.na(qa_num_birth_year)) %>%
    nrow()
  n_yob_after_yod <- # Rule 2: Year of birth is after year of death
    data_processed %>%
    filter(qa_num_birth_year > year(qa_date_of_death)
      # (!is.na(qa_date_of_death) & !is.na(qa_num_birth_year)) &
      ) %>%
    nrow()
  n_yob_beforeNHS_aftertoday <- # Rule 3: Year of birth predates NHS established year or year of birth exceeds current date
    data_processed %>%
    filter(qa_num_birth_year < 1793 | qa_num_birth_year > year(Sys.Date())) %>%
    nrow()
  n_dob_invalid <- # Rule 4: Date of death is on or before 1/1/1900 (and not NULL) or after current date (and not NULL)
    data_processed %>%
    filter((qa_date_of_death <= as.Date("1900-01-01")) | (qa_date_of_death > Sys.Date())) %>%
    nrow()
  n_preg_men <- # Rule 5: Pregnancy/birth codes for men
    data_processed %>%
    filter(qa_bin_pregnancy == TRUE & cov_cat_sex == "Male") %>%
    nrow()
  n_hrt_men <- # Rule 6: HRT or COCP meds for men
    data_processed %>%
    filter(((cov_cat_sex == "Male" & qa_bin_hrt == TRUE) | (cov_cat_sex == "Male" & qa_bin_cocp == TRUE))) %>%
    nrow()
  n_prost_women <- # Rule 7: Prostate cancer codes for women
    data_processed %>%
    filter((qa_bin_prostate_cancer == TRUE & cov_cat_sex == "Female")) %>%
    nrow()

  n_after_qa_processing <-
    data_processed %>%
    filter(!is.na(qa_num_birth_year)) %>%
    filter(is.na(qa_date_of_death) | (qa_num_birth_year <= year(qa_date_of_death))) %>%
    filter(qa_num_birth_year >= 1793 & qa_num_birth_year <= year(Sys.Date())) %>%
    filter((qa_date_of_death > as.Date("1900-01-01")) | (qa_date_of_death < Sys.Date()) | is.na(qa_date_of_death)) %>%
    filter((cov_cat_sex == "Female" | is.na(cov_cat_sex)) | (cov_cat_sex == "Male" & (qa_bin_pregnancy == FALSE))) %>% # FALSE includes missing in a ehrQL logical
    filter((cov_cat_sex == "Female" | is.na(cov_cat_sex)) | (cov_cat_sex == "Male" & (qa_bin_hrt == FALSE)) | (cov_cat_sex == "Male" & (qa_bin_cocp == FALSE))) %>%
    filter((cov_cat_sex == "Male" | is.na(cov_cat_sex)) | (cov_cat_sex == "Female" & (qa_bin_prostate_cancer == FALSE))) %>%
    nrow()

  out <- tibble(n_before_qa_processing,
                n_yob_missing,
                n_yob_after_yod,
                n_yob_beforeNHS_aftertoday,
                n_dob_invalid,
                n_preg_men,
                n_hrt_men,
                n_prost_women,
                n_after_qa_processing
  )
}
