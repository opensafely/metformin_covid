################################################################################
# Custom made function to assess the eligibility criteria
################################################################################
fn_elig_criteria_midpoint6 <- function(data_processed, study_dates, years_in_days) {
  # Extract dates from the study_dates list
  landmark_date <- as.Date(study_dates$landmark_date, format = "%Y-%m-%d")
  mid2018_date <- as.Date(study_dates$mid2018_date, format = "%Y-%m-%d")
  
  # Apply the eligibility criteria (Aged ≥ 18 years already applied in completeness criteria)
  data_processed <- data_processed %>%
    mutate(
      # Exclusion 1: no T2DM diagnosis or out of window
      no_t2dm_or_outofwindow = is.na(elig_date_t2dm) | 
        (elig_date_t2dm >= landmark_date - days(183) & elig_date_t2dm < landmark_date) | 
        (elig_date_t2dm < mid2018_date - days(years_in_days)),
      # Exclusion 2: metformin use prior to T2DM diagnosis
      prior_metfin = exp_date_metfin_first <= elig_date_t2dm,
      # Exclusion 3: metformin allergy prior to T2DM diagnosis
      prior_metfin_allergy = elig_date_metfin_allergy <= elig_date_t2dm,
      # Exclusion 4: CKD 4/5 prior to T2DM diagnosis
      prior_ckd45 = elig_date_ckd_45 <= elig_date_t2dm,
      # Exclusion 5: liver cirrhosis prior to T2DM diagnosis
      prior_cirrhosis = elig_date_liver_cirrhosis <= elig_date_t2dm,
      # Exclusion 6: prior drug with interaction risk with metfin, in 14 days window before T2DM diagnosis
      prior_interaction = (elig_date_metfin_interaction <= elig_date_t2dm) & (elig_date_metfin_interaction > elig_date_t2dm - days(14))
    )
  
  # Count the criteria
  count <- data_processed %>%
    summarise(
      n_no_t2dm_or_outofwindow = sum(no_t2dm_or_outofwindow, na.rm = TRUE),
      n_prior_metfin = sum(prior_metfin, na.rm = TRUE),
      n_prior_metfin_allergy = sum(prior_metfin_allergy, na.rm = TRUE),
      n_prior_ckd45 = sum(prior_ckd45, na.rm = TRUE),
      n_prior_cirrhosis = sum(prior_cirrhosis, na.rm = TRUE),
      n_prior_interaction = sum(prior_interaction, na.rm = TRUE)
    )
  
  # Filter
  data_filtered <- data_processed %>% # Output 1: filtered data
    filter(
      (!no_t2dm_or_outofwindow | is.na(no_t2dm_or_outofwindow)),
      (!prior_metfin | is.na(prior_metfin)),
      (!prior_metfin_allergy | is.na(prior_metfin_allergy)),
      (!prior_ckd45 | is.na(prior_ckd45)),
      (!prior_cirrhosis | is.na(prior_cirrhosis)),
      (!prior_interaction | is.na(prior_interaction))
    )
  
  n_after_exclusion_processing <- nrow(data_filtered)
  
  # Output 2: Count for flowchart, incl. redaction and including pre/post processing counts
  out_midpoint6 <- tibble(
    n_before_exclusion_processing_midpoint6 = fn_roundmid_any(nrow(data_processed), threshold),
    n_no_t2dm_or_outofwindow_midpoint6 = fn_roundmid_any(count$n_no_t2dm_or_outofwindow, threshold),
    n_prior_metfin_midpoint6 = fn_roundmid_any(count$n_prior_metfin, threshold),
    n_prior_metfin_allergy_midpoint6 = fn_roundmid_any(count$n_prior_metfin_allergy, threshold),
    n_prior_ckd45_midpoint6 = fn_roundmid_any(count$n_prior_ckd45, threshold),
    n_prior_cirrhosis_midpoint6 = fn_roundmid_any(count$n_prior_cirrhosis, threshold),
    n_prior_interaction_midpoint6 = fn_roundmid_any(count$n_prior_interaction, threshold),
    n_after_exclusion_processing_midpoint6 = fn_roundmid_any(n_after_exclusion_processing, threshold)
  )
  
  # Return both outputs as a list
  return(list(
    n_elig_excluded_midpoint6 = out_midpoint6,
    data_processed = data_filtered
  ))
  
}