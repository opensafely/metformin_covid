################################################################################
# A custom made function to assess the quality assurance criteria
################################################################################
fn_quality_assurance_midpoint6 <- function(data_processed, study_dates, threshold){
  
  # Extract study end date for Rule 3
  studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")
  
  # Apply the rules
  data_processed <- data_processed %>%
    mutate(
      # Rule 1: Year of birth is missing
      yob_missing = is.na(qa_num_birth_year), 
      # Rule 2: Year of birth is after year of death
      yob_after_yod = !is.na(qa_num_birth_year) & qa_num_birth_year > year(qa_date_of_death),
      # Rule 3: Year of birth predates NHS established year or year of birth exceeds study end date
      yob_beforeNHS_afterstudyend = qa_num_birth_year < 1838 | (qa_num_birth_year > year(studyend_date)),
      # Rule 4: Date of death is on or before 1/1/1900 (and not NULL) or after current date (and not NULL)
      dod_invalid = (qa_date_of_death <= as.Date("1900-01-01") | qa_date_of_death > Sys.Date()),
      # Rule 5: Pregnancy/birth codes for men
      preg_men = qa_bin_pregnancy == TRUE & cov_cat_sex == "Male",
      # Rule 6: HRT or COCP meds for men
      hrt_men = cov_cat_sex == "Male" & (qa_bin_hrt == TRUE | qa_bin_cocp == TRUE),
      # Rule 7: Prostate cancer codes for women
      prost_women = cov_cat_sex == "Female" & qa_bin_prostate_cancer == TRUE
    )
  
  # Count the rules
  count <- data_processed %>%
    summarise(
      n_yob_missing = sum(yob_missing, na.rm = TRUE),
      n_yob_after_yod = sum(yob_after_yod, na.rm = TRUE),
      n_yob_beforeNHS_afterstudyend = sum(yob_beforeNHS_afterstudyend, na.rm = TRUE),
      n_dod_invalid = sum(dod_invalid, na.rm = TRUE),
      n_preg_men = sum(preg_men, na.rm = TRUE),
      n_hrt_men = sum(hrt_men, na.rm = TRUE),
      n_prost_women = sum(prost_women, na.rm = TRUE)
    )
  
  # Filter
  data_filtered <- data_processed %>% # Output 1: filtered data
    filter(
      (!yob_missing | is.na(yob_missing)),
      (!yob_after_yod | is.na(yob_after_yod)),
      (!yob_beforeNHS_afterstudyend | is.na(yob_beforeNHS_afterstudyend)),
      (!dod_invalid | is.na(dod_invalid)),
      (!preg_men | is.na(preg_men)),
      (!hrt_men | is.na(hrt_men)),
      (!prost_women | is.na(prost_women))
    )
  
  n_after_qa_processing <- nrow(data_filtered)

  # Output 2: Count for flowchart, incl. redaction and including pre/post processing counts
  out_midpoint6 <- tibble(
    n_before_qa_processing_midpoint6 = fn_roundmid_any(nrow(data_processed), threshold),
    n_yob_missing_midpoint6 = fn_roundmid_any(count$n_yob_missing, threshold),
    n_yob_after_yod_midpoint6 = fn_roundmid_any(count$n_yob_after_yod, threshold),
    n_yob_beforeNHS_afterstudyend_midpoint6 = fn_roundmid_any(count$n_yob_beforeNHS_afterstudyend, threshold),
    n_dod_invalid_midpoint6 = fn_roundmid_any(count$n_dod_invalid, threshold),
    n_preg_men_midpoint6 = fn_roundmid_any(count$n_preg_men, threshold),
    n_hrt_men_midpoint6 = fn_roundmid_any(count$n_hrt_men, threshold),
    n_prost_women_midpoint6 = fn_roundmid_any(count$n_prost_women, threshold),
    n_after_qa_processing_midpoint6 = fn_roundmid_any(n_after_qa_processing, threshold)
  ) %>% 
    # pivot (for easier data review in L4)
    pivot_longer(
      cols = everything(),
      names_to = "Variable",
      values_to = "Value"
    )
  
  # Return both outputs as a list
  return(list(
    n_qa_excluded_midpoint6 = out_midpoint6,
    data_processed = data_filtered
  ))
  
}