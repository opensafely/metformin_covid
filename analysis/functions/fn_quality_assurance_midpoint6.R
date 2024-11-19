################################################################################
# A custom made function to run the quality assurance criteria
################################################################################
fn_quality_assurance_midpoint6 <- function(data_processed){
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
  
  # Output incl. redaction
  out_midpoint6 <- tibble(
      n_before_qa_processing_midpoint6 = fn_roundmid_any(n_before_qa_processing, threshold),
      n_yob_missing_midpoint6 = fn_roundmid_any(n_yob_missing, threshold),
      n_yob_after_yod_midpoint6 = fn_roundmid_any(n_yob_after_yod, threshold),
      n_yob_beforeNHS_aftertoday_midpoint6 = fn_roundmid_any(n_yob_beforeNHS_aftertoday, threshold),
      n_dob_invalid_midpoint6 = fn_roundmid_any(n_dob_invalid, threshold),
      n_preg_men_midpoint6 = fn_roundmid_any(n_preg_men, threshold),
      n_hrt_men_midpoint6d = fn_roundmid_any(n_hrt_men, threshold),
      n_prost_women_midpoint6 = fn_roundmid_any(n_prost_women, threshold),
      n_after_qa_processing_midpoint6 = fn_roundmid_any(n_after_qa_processing, threshold),
    )

  return(out_midpoint6)
}
