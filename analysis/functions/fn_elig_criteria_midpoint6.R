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
        #(elig_date_t2dm >= landmark_date - days(183) & elig_date_t2dm < landmark_date) | # allow for T2DM diagnoses in 6m before landmark!
        (elig_date_t2dm < mid2018_date - days(years_in_days)),
      # Exclusion 2: metformin use prior to T2DM diagnosis
      prior_metfin = exp_date_metfin_first < elig_date_t2dm, # but don't count those who initiated on day of diagnosis & codelist allows for metformin including combo treatment
      # Exclusion 3: metformin allergy prior to or on T2DM diagnosis
      prior_metfin_allergy = elig_date_metfin_allergy_first <= elig_date_t2dm, # count those diagnosed with allergy on day of diagnosis
      # Exclusion 4: CKD 4/5 prior to or on T2DM diagnosis
      prior_ckd45 = elig_date_ckd_45_first <= elig_date_t2dm, # count those diagnosed with ckd on day of diagnosis
      # Exclusion 5: liver cirrhosis prior to or on T2DM diagnosis
      prior_cirrhosis = elig_date_liver_cirrhosis_first <= elig_date_t2dm, # count those diagnosed with cirrhosis on day of diagnosis
      # Exclusion 6: prior drug with interaction risk with metfin, in 14 days window before or on T2DM diagnosis day
      prior_interaction = (elig_date_metfin_interaction_first <= elig_date_t2dm) & (elig_date_metfin_interaction_first >= elig_date_t2dm - days(14))
    )
  
  # Re-Apply the time-updated eligibility criteria again at landmark (now: use _last variables, see codebook!)
  data_processed <- data_processed %>%
    mutate(
      # Exclusion 3: metformin allergy prior to or on landmark
      prior_metfin_allergy_landmark = (elig_date_metfin_allergy_first > mid2018_date - days(years_in_days)) # don't count those diagnosed with allergy on day of diagnosis again
                                      & elig_date_metfin_allergy_first <= landmark_date,
      # Exclusion 4: CKD 4/5 prior to or on landmark
      prior_ckd45_landmark = (elig_date_ckd_45_first > mid2018_date - days(years_in_days)) 
                              & elig_date_ckd_45_first <= landmark_date,
      # Exclusion 5: liver cirrhosis prior to or on landmark
      prior_cirrhosis_landmark = (elig_date_liver_cirrhosis_first > mid2018_date - days(years_in_days)) 
                                  & elig_date_liver_cirrhosis_first <= landmark_date,
      # Exclusion 6: prior drug with interaction risk with metfin, in 14 days window prior to or on landmark
      prior_interaction_landmark = (elig_date_metfin_interaction_last <= landmark_date) & (elig_date_metfin_interaction_last >= landmark_date - days(14)),
      # Exclusion 7: died prior to landmark
      prior_death_landmark = (qa_date_of_death > mid2018_date - days(years_in_days))
                              & qa_date_of_death <= landmark_date,
      # Exclusion 8: LTFU prior to landmark
      prior_ltfu_landmark = (out_date_dereg_mid2018 > mid2018_date - days(years_in_days))
                            & out_date_dereg_mid2018 <= landmark_date,
    )

  # Filter 1: main inclusion criteria: T2DM diagnosis
  data_filtered_T2DM <- data_processed %>%
    filter(
      (!no_t2dm_or_outofwindow | is.na(no_t2dm_or_outofwindow))
    )
  
  n_t2dm <- nrow(data_filtered_T2DM)
  
  # Among these, count the exclusion criteria
  count <- data_filtered_T2DM %>%
    summarise(
      n_prior_metfin = sum(prior_metfin, na.rm = TRUE),
      n_prior_metfin_allergy = sum(prior_metfin_allergy, na.rm = TRUE),
      n_prior_ckd45 = sum(prior_ckd45, na.rm = TRUE),
      n_prior_cirrhosis = sum(prior_cirrhosis, na.rm = TRUE),
      n_prior_interaction = sum(prior_interaction, na.rm = TRUE),
      n_prior_metfin_allergy_landmark = sum(prior_metfin_allergy_landmark, na.rm = TRUE),
      n_prior_ckd45_landmark = sum(prior_ckd45_landmark, na.rm = TRUE),
      n_prior_cirrhosis_landmark = sum(prior_cirrhosis_landmark, na.rm = TRUE),
      n_prior_interaction_landmark = sum(prior_interaction_landmark, na.rm = TRUE),
      n_prior_death_landmark = sum(prior_death_landmark, na.rm = TRUE),
      n_prior_ltfu_landmark = sum(prior_ltfu_landmark, na.rm = TRUE)
    )
  
  # Filter 2: apply inclusion & all exclusion criteria
  data_filtered <- data_processed %>% # Output 1: filtered data
    filter(
      (!no_t2dm_or_outofwindow | is.na(no_t2dm_or_outofwindow)),
      (!prior_metfin | is.na(prior_metfin)),
      (!prior_metfin_allergy | is.na(prior_metfin_allergy)),
      (!prior_ckd45 | is.na(prior_ckd45)),
      (!prior_cirrhosis | is.na(prior_cirrhosis)),
      (!prior_interaction | is.na(prior_interaction)),
      (!prior_metfin_allergy_landmark | is.na(prior_metfin_allergy_landmark)),
      (!prior_ckd45_landmark | is.na(prior_ckd45_landmark)),
      (!prior_cirrhosis_landmark | is.na(prior_cirrhosis_landmark)),
      (!prior_interaction_landmark | is.na(prior_interaction_landmark)),
      (!prior_death_landmark | is.na(prior_death_landmark)),
      (!prior_ltfu_landmark | is.na(prior_ltfu_landmark))
    )
  
  n_after_exclusion_processing <- nrow(data_filtered)
  
  # Output 2: Count for flowchart, without redaction, and including pre/post processing counts
  out <- tibble(
    n_before_exclusion_processing = nrow(data_processed),
    n_t2dm = n_t2dm, # counted among all data_processed
    n_prior_metfin = count$n_prior_metfin, # counted only among n_t2dm !
    n_prior_metfin_allergy = count$n_prior_metfin_allergy, # counted only among n_t2dm !
    n_prior_ckd45 = count$n_prior_ckd45, # counted only among n_t2dm !
    n_prior_cirrhosis = count$n_prior_cirrhosis, # counted only among n_t2dm !
    n_prior_interaction = count$n_prior_interaction, # counted only among n_t2dm !
    n_prior_metfin_allergy_landmark = count$n_prior_metfin_allergy_landmark, # counted only among n_t2dm !
    n_prior_ckd45_landmark = count$n_prior_ckd45_landmark, # counted only among n_t2dm !
    n_prior_cirrhosis_landmark = count$n_prior_cirrhosis_landmark, # counted only among n_t2dm !
    n_prior_interaction_landmark = count$n_prior_interaction_landmark, # counted only among n_t2dm !
    n_prior_death_landmark = count$n_prior_death_landmark, # counted only among n_t2dm !
    n_prior_ltfu_landmark = count$n_prior_ltfu_landmark, # counted only among n_t2dm !
    n_after_exclusion_processing = n_after_exclusion_processing
  )
  
  # Output 3: Count for flowchart, with redaction, and including pre/post processing counts
  out_midpoint6 <- tibble(
    n_before_exclusion_processing_midpoint6 = fn_roundmid_any(nrow(data_processed), threshold),
    n_t2dm_midpoint6 = fn_roundmid_any(n_t2dm, threshold), # counted among all data_processed
    n_prior_metfin_midpoint6 = fn_roundmid_any(count$n_prior_metfin, threshold), # counted only among n_t2dm !
    n_prior_metfin_allergy_midpoint6 = fn_roundmid_any(count$n_prior_metfin_allergy, threshold), # counted only among n_t2dm !
    n_prior_ckd45_midpoint6 = fn_roundmid_any(count$n_prior_ckd45, threshold), # counted only among n_t2dm !
    n_prior_cirrhosis_midpoint6 = fn_roundmid_any(count$n_prior_cirrhosis, threshold), # counted only among n_t2dm !
    n_prior_interaction_midpoint6 = fn_roundmid_any(count$n_prior_interaction, threshold), # counted only among n_t2dm !
    n_prior_metfin_allergy_landmark_midpoint6 = fn_roundmid_any(count$n_prior_metfin_allergy_landmark, threshold), # counted only among n_t2dm !
    n_prior_ckd45_landmark_midpoint6 = fn_roundmid_any(count$n_prior_ckd45_landmark, threshold), # counted only among n_t2dm !
    n_prior_cirrhosis_landmark_midpoint6 = fn_roundmid_any(count$n_prior_cirrhosis_landmark, threshold), # counted only among n_t2dm !
    n_prior_interaction_landmark_midpoint6 = fn_roundmid_any(count$n_prior_interaction_landmark, threshold), # counted only among n_t2dm !
    n_prior_death_landmark_midpoint6 = fn_roundmid_any(count$n_prior_death_landmark, threshold), # counted only among n_t2dm !
    n_prior_ltfu_landmark_midpoint6 = fn_roundmid_any(count$n_prior_ltfu_landmark, threshold), # counted only among n_t2dm !
    n_after_exclusion_processing_midpoint6 = fn_roundmid_any(n_after_exclusion_processing, threshold)
  )
  
  # Return outputs as a list
  return(list(
    n_elig_excluded = out,
    n_elig_excluded_midpoint6 = out_midpoint6,
    data_processed = data_filtered
  ))
  
}