################################################################################
# Two custom made function to a) assess (-> for flow chart) and b) apply the eligibility criteria
################################################################################
fn_elig_criteria <- function(data_processed, study_dates, years_in_days) {
  # Initial number of rows before exclusion criteria
  n_before_exclusion_processing <- nrow(data_processed)
  
  # eligibilty criteria 1) T2DM diagnosis 6m prior landmark, but after the looped mid-year date
  data_processed_t2dm <- data_processed %>%
    filter(cov_date_t2dm < landmark_date - days(183) &
             cov_date_t2dm > mid2018_date - days(years_in_days))
  n_t2dm <- nrow(data_processed_t2dm)
  
  # eligibilty criteria 2) no prior metformin use, defined as any metformin prescription on or prior to the T2DM diagnosis, i.e., did not start metformin or only after T2DM diagnosis
  data_processed_elig2 <- data_processed_t2dm %>%
    filter(exp_date_metfin_first > cov_date_t2dm | is.na(exp_date_metfin_first))
  n_prior_metfin <- nrow(data_processed_elig2)
  
  # eligibilty criteria 3) no known hypersensitivity and / or intolerance to metformin, on or prior to the T2DM diagnosis
  data_processed_elig3 <- data_processed_t2dm %>%
    filter(elig_date_metfin_allergy > cov_date_t2dm | is.na(elig_date_metfin_allergy))
  n_metfin_allergy <- nrow(data_processed_elig3)
  
  # eligibilty criteria 4) no clinical history of moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; chronic stage 4/5) on or prior to the T2DM diagnosis
  data_processed_elig4 <- data_processed_t2dm %>%
    filter(elig_date_ckd_45 > cov_date_t2dm | is.na(elig_date_ckd_45))
  n_ckd45 <- nrow(data_processed_elig4)
  
  # eligibilty criteria 5) no clinical history of advance decompensated liver cirrhosis, on or prior to the T2DM diagnosis
  data_processed_elig5 <- data_processed_t2dm %>%
    filter(elig_date_liver_cirrhosis > cov_date_t2dm | is.na(elig_date_liver_cirrhosis))
  n_cirrhosis <- nrow(data_processed_elig5)
  
  # eligibilty criteria 6) no use of the following medications in the last 14 days, on or prior to the T2DM diagnosis: Cimetidine, hydroxychloroquine, dolutegravir, patiromer, ranolazine, Monoamine Oxide Inhibitors (Phenelzine, Tranylcypromine, Selegiline, moclobemide), Sotalol, Clonidine, Methyldopa, Prazosin, Doxazosin
  data_processed_elig6 <- data_processed_t2dm %>%
    filter(elig_date_metfin_interaction > cov_date_t2dm | 
             is.na(elig_date_metfin_interaction) | 
             elig_date_metfin_interaction < cov_date_t2dm - days(14))
  n_interaction <- nrow(data_processed_elig6)

  n_after_exclusion_processing <-
    data_processed %>%
    # eligibilty criteria 1) T2DM status in defined period before landmark
    filter(cov_date_t2dm < landmark_date - days(183) &
             cov_date_t2dm > mid2018_date - days(years_in_days)) %>% 
    # eligibilty criteria 2) no prior metformin use, defined as any metformin prescription on or prior to the T2DM diagnosis, i.e., did not start metformin or only after T2DM diagnosis
    filter(exp_date_metfin_first > cov_date_t2dm | is.na(exp_date_metfin_first)) %>%
    # eligibilty criteria 3) no known hypersensitivity and / or intolerance to metformin, on or prior to the T2DM diagnosis
    filter(elig_date_metfin_allergy > cov_date_t2dm | is.na(elig_date_metfin_allergy)) %>%
    # eligibilty criteria 4) no clinical history of moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; chronic stage 4/5) on or prior to the T2DM diagnosis
    filter(elig_date_ckd_45 > cov_date_t2dm | is.na(elig_date_ckd_45)) %>%
    # eligibilty criteria 5) no clinical history of advance decompensated liver cirrhosis, on or prior to the T2DM diagnosis
    filter(elig_date_liver_cirrhosis > cov_date_t2dm | is.na(elig_date_liver_cirrhosis)) %>%
    # eligibilty criteria 6) no use of the following medications in the last 14 days, on or prior to the T2DM diagnosis: Cimetidine, hydroxychloroquine, dolutegravir, patiromer, ranolazine, Monoamine Oxide Inhibitors (Phenelzine, Tranylcypromine, Selegiline, moclobemide), Sotalol, Clonidine, Methyldopa, Prazosin, Doxazosin
    filter(elig_date_metfin_interaction > cov_date_t2dm | is.na(elig_date_metfin_interaction) | elig_date_metfin_interaction < cov_date_t2dm - days(14)) %>% 
    nrow()
  
  out <- tibble(
    n_before_exclusion_processing,
    n_t2dm,
    n_prior_metfin,
    n_metfin_allergy,
    n_ckd45,
    n_cirrhosis,
    n_interaction,
    n_after_exclusion_processing
  )
  
  return(out)
}

fn_apply_elig_criteria <- function(data_processed, study_dates, years_in_days) {
  data_processed <- data_processed %>%
    # eligibilty criteria 1) T2DM diagnosis 6m prior landmark, but after the looped mid-year date
    filter(cov_date_t2dm < landmark_date - days(183) &
             cov_date_t2dm > mid2018_date - days(years_in_days)) %>% 
    # eligibilty criteria 2) no prior metformin use, defined as any metformin prescription on or prior to the T2DM diagnosis, i.e., did not start metformin or only after T2DM diagnosis
    filter(exp_date_metfin_first > cov_date_t2dm | is.na(exp_date_metfin_first)) %>%
    # eligibilty criteria 3) no known hypersensitivity and / or intolerance to metformin, on or prior to the T2DM diagnosis
    filter(elig_date_metfin_allergy > cov_date_t2dm | is.na(elig_date_metfin_allergy)) %>%
    # eligibilty criteria 4) no clinical history of moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; chronic stage 4/5) on or prior to the T2DM diagnosis
    filter(elig_date_ckd_45 > cov_date_t2dm | is.na(elig_date_ckd_45)) %>%
    # eligibilty criteria 5) no clinical history of advance decompensated liver cirrhosis, on or prior to the T2DM diagnosis
    filter(elig_date_liver_cirrhosis > cov_date_t2dm | is.na(elig_date_liver_cirrhosis)) %>%
    # eligibilty criteria 6) no use of the following medications in the last 14 days, on or prior to the T2DM diagnosis: Cimetidine, hydroxychloroquine, dolutegravir, patiromer, ranolazine, Monoamine Oxide Inhibitors (Phenelzine, Tranylcypromine, Selegiline, moclobemide), Sotalol, Clonidine, Methyldopa, Prazosin, Doxazosin
    filter(elig_date_metfin_interaction > cov_date_t2dm | is.na(elig_date_metfin_interaction) | elig_date_metfin_interaction < cov_date_t2dm - days(14))
}