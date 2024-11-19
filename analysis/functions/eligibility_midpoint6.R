################################################################################
# Two custom made function to a) assess (-> for flow chart) and b) apply the eligibility criteria
################################################################################

######
# a) Assess eligibility criteria for flow chart presentation
######
fn_elig_criteria_midpoint6 <- function(data_processed, study_dates, years_in_days) {
  # Extract dates from the study_dates list
  landmark_date <- study_dates$landmark_date
  mid2018_date <- study_dates$mid2018_date
  
  # Initial number of rows before exclusion criteria
  n_before_exclusion_processing <-
    data_processed %>%
    nrow()
  
  # eligibility criteria 1) T2DM diagnosis 6m prior landmark, but after the looped mid-year date
  n_no_t2dm_or_outofwindow <- # no T2DM diagnosis or out of window
    data_processed %>%
    filter(is.na(elig_date_t2dm) | elig_date_t2dm > landmark_date - days(183) | elig_date_t2dm < mid2018_date - days(years_in_days)) %>%
    nrow()
  
  # eligibility criteria 2) no prior metformin use, defined as any metformin prescription on or prior to the T2DM diagnosis, i.e., did not start metformin or only after T2DM diagnosis
  n_prior_metfin <- # prior use
    data_processed %>%
    filter(exp_date_metfin_first <= elig_date_t2dm) %>%
    nrow()

  # eligibility criteria 3) no known hypersensitivity and / or intolerance to metformin, on or prior to the T2DM diagnosis
  n_metfin_allergy <- # prior metformin allergy
    data_processed %>%
    filter(elig_date_metfin_allergy <= elig_date_t2dm) %>%
    nrow()

  # eligibility criteria 4) no clinical history of moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; chronic stage 4/5) on or prior to the T2DM diagnosis
  n_ckd45 <- # prior CKD 4/5
    data_processed %>%
    filter(elig_date_ckd_45 <= elig_date_t2dm) %>%
    nrow()
  
  # eligibility criteria 5) no clinical history of advance decompensated liver cirrhosis, on or prior to the T2DM diagnosis
  n_cirrhosis <- # prior liver cirrhosis
    data_processed %>%
    filter(elig_date_liver_cirrhosis <= elig_date_t2dm) %>%
    nrow()
  
  # eligibility criteria 6) no use of the following medications in the last 14 days, on or prior to the T2DM diagnosis: Cimetidine, hydroxychloroquine, dolutegravir, patiromer, ranolazine, Monoamine Oxide Inhibitors (Phenelzine, Tranylcypromine, Selegiline, moclobemide), Sotalol, Clonidine, Methyldopa, Prazosin, Doxazosin
  n_interaction <- # prior drug with interaction risk with metfin, in window 14 days before eligibility
    data_processed %>%
    filter(elig_date_metfin_interaction <= elig_date_t2dm & elig_date_metfin_interaction > elig_date_t2dm - days(14)) %>%
    nrow()

  n_after_exclusion_processing <-
    data_processed %>%
    # eligibilty criteria 1) T2DM status in defined period before landmark
    filter(elig_date_t2dm < landmark_date - days(183) &
             elig_date_t2dm > mid2018_date - days(years_in_days)) %>% 
    # eligibilty criteria 2) no prior metformin use, defined as any metformin prescription on or prior to the T2DM diagnosis, i.e., did not start metformin or only after T2DM diagnosis
    filter(exp_date_metfin_first > elig_date_t2dm | is.na(exp_date_metfin_first)) %>%
    # eligibilty criteria 3) no known hypersensitivity and / or intolerance to metformin, on or prior to the T2DM diagnosis
    filter(elig_date_metfin_allergy > elig_date_t2dm | is.na(elig_date_metfin_allergy)) %>%
    # eligibilty criteria 4) no clinical history of moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; chronic stage 4/5) on or prior to the T2DM diagnosis
    filter(elig_date_ckd_45 > elig_date_t2dm | is.na(elig_date_ckd_45)) %>%
    # eligibilty criteria 5) no clinical history of advance decompensated liver cirrhosis, on or prior to the T2DM diagnosis
    filter(elig_date_liver_cirrhosis > elig_date_t2dm | is.na(elig_date_liver_cirrhosis)) %>%
    # eligibilty criteria 6) no use of the following medications in the last 14 days, on or prior to the T2DM diagnosis: Cimetidine, hydroxychloroquine, dolutegravir, patiromer, ranolazine, Monoamine Oxide Inhibitors (Phenelzine, Tranylcypromine, Selegiline, moclobemide), Sotalol, Clonidine, Methyldopa, Prazosin, Doxazosin
    filter(elig_date_metfin_interaction > elig_date_t2dm | is.na(elig_date_metfin_interaction) | elig_date_metfin_interaction < elig_date_t2dm - days(14)) %>% 
    nrow()
  
  # Output
  # out <- tibble(
  #   n_before_exclusion_processing,
  #   n_no_t2dm_or_outofwindow,
  #   n_prior_metfin,
  #   n_metfin_allergy,
  #   n_ckd45,
  #   n_cirrhosis,
  #   n_interaction,
  #   n_after_exclusion_processing
  # )
  # 
  # return(out)
  
  # Output incl. redaction
  out_midpoint6 <- tibble(
      n_before_exclusion_processing_midpoint6 = fn_roundmid_any(n_before_exclusion_processing, threshold),
      n_no_t2dm_or_outofwindow_midpoint6 = fn_roundmid_any(n_no_t2dm_or_outofwindow, threshold),
      n_prior_metfin_midpoint6 = fn_roundmid_any(n_prior_metfin, threshold),
      n_metfin_allergy_midpoint6 = fn_roundmid_any(n_metfin_allergy, threshold),
      n_ckd45_midpoint6 = fn_roundmid_any(n_ckd45, threshold),
      n_cirrhosis_midpoint6 = fn_roundmid_any(n_cirrhosis, threshold),
      n_interaction_midpoint6 = fn_roundmid_any(n_interaction, threshold),
      n_after_exclusion_processing_midpoint6 = fn_roundmid_any(n_after_exclusion_processing, threshold),
    )

  return(out_midpoint6)
}


######
# b) Apply eligibility criteria
######
fn_apply_elig_criteria <- function(data_processed, study_dates, years_in_days) {
  # Extract dates from the study_dates list
  landmark_date <- study_dates$landmark_date
  mid2018_date <- study_dates$mid2018_date
  
  data_processed <- data_processed %>%
    # eligibilty criteria 1) T2DM diagnosis 6m prior landmark, but after the looped mid-year date
    filter(elig_date_t2dm < landmark_date - days(183) &
             elig_date_t2dm > mid2018_date - days(years_in_days)) %>% 
    # eligibilty criteria 2) no prior metformin use, defined as any metformin prescription on or prior to the T2DM diagnosis, i.e., did not start metformin or only after T2DM diagnosis
    filter(exp_date_metfin_first > elig_date_t2dm | is.na(exp_date_metfin_first)) %>%
    # eligibilty criteria 3) no known hypersensitivity and / or intolerance to metformin, on or prior to the T2DM diagnosis
    filter(elig_date_metfin_allergy > elig_date_t2dm | is.na(elig_date_metfin_allergy)) %>%
    # eligibilty criteria 4) no clinical history of moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; chronic stage 4/5) on or prior to the T2DM diagnosis
    filter(elig_date_ckd_45 > elig_date_t2dm | is.na(elig_date_ckd_45)) %>%
    # eligibilty criteria 5) no clinical history of advance decompensated liver cirrhosis, on or prior to the T2DM diagnosis
    filter(elig_date_liver_cirrhosis > elig_date_t2dm | is.na(elig_date_liver_cirrhosis)) %>%
    # eligibilty criteria 6) no use of the following medications in the last 14 days, on or prior to the T2DM diagnosis: Cimetidine, hydroxychloroquine, dolutegravir, patiromer, ranolazine, Monoamine Oxide Inhibitors (Phenelzine, Tranylcypromine, Selegiline, moclobemide), Sotalol, Clonidine, Methyldopa, Prazosin, Doxazosin
    filter(elig_date_metfin_interaction > elig_date_t2dm | is.na(elig_date_metfin_interaction) | elig_date_metfin_interaction < elig_date_t2dm - days(14))
}
