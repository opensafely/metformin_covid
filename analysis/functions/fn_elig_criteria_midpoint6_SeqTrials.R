################################################################################
# Custom made function to assess the eligibility criteria
################################################################################
fn_elig_criteria_midpoint6_SeqTrials <- function(data_processed, study_dates, months_in_days) {
  # Extract dates from the study_dates list
  pandemicstart_date <- as.Date(study_dates$pandemicstart_date, format = "%Y-%m-%d")
  mid2018_date <- as.Date(study_dates$mid2018_date, format = "%Y-%m-%d")
    
    # Apply the main eligibility criteria
    data_processed <- data_processed %>%
      mutate(
       # Exclusion 1: no T2DM diagnosis a diagnosis more than 6 months prior to pandemic
        no_t2dm_or_outofwindow = is.na(elig_date_t2dm) # non-diabetics already excluded in completeness criteria...so, here it's only those with a diagnosis more than 6 months ago
        | (elig_date_t2dm < pandemicstart_date - days(months_in_days)))
    
    # Filter 1: main eligibility criteria: T2DM diagnosis
    data_filtered_T2DM <- data_processed %>%
      filter((!no_t2dm_or_outofwindow | is.na(no_t2dm_or_outofwindow)))
    
    n_t2dm <- nrow(data_filtered_T2DM)
    
    # Apply the other eligibility criteria (Aged ≥ 18 years already applied in completeness criteria)
    data_filtered_T2DM <- data_filtered_T2DM %>%
      mutate(
        # Exclusion 2: metformin use prior to pandemic start (or any other antidiabetic!)
        prior_metfin = (!is.na(exp_date_metfin_first)
                        & exp_date_metfin_first < pandemicstart_date), # but don't count those who initiated on day of diagnosis & codelist allows for metformin including combo treatment
        prior_sulfo_mono = (!is.na(exp_date_sulfo_first)
                            & exp_date_sulfo_first < pandemicstart_date),
        prior_dpp4_mono = (!is.na(exp_date_dpp4_first)
                           & exp_date_dpp4_first < pandemicstart_date),
        prior_tzd_mono = (!is.na(exp_date_tzd_first)
                          & exp_date_tzd_first < pandemicstart_date),
        prior_sglt2_mono = (!is.na(exp_date_sglt2_first)
                            & exp_date_sglt2_first < pandemicstart_date),
        prior_glp1_mono = (!is.na(exp_date_glp1_first)
                           & exp_date_glp1_first < pandemicstart_date),
        prior_megli_mono = (!is.na(exp_date_megli_first)
                            & exp_date_megli_first < pandemicstart_date),
        prior_agi_mono = (!is.na(exp_date_agi_first)
                          & exp_date_agi_first < pandemicstart_date),
        prior_insulin_mono = (!is.na(exp_date_insulin_first)
                              & exp_date_insulin_first < pandemicstart_date),
        # Exclusion 3: metformin allergy prior to or on pandemic start
        prior_metfin_allergy = (!is.na(elig_date_metfin_allergy_first) 
                                & elig_date_metfin_allergy_first <= pandemicstart_date), # count those diagnosed with allergy on day of diagnosis
        # Exclusion 4: CKD 4/5 prior to or on pandemic start
        prior_ckd45 = (!is.na(elig_date_ckd_45_first) 
                       & elig_date_ckd_45_first <= pandemicstart_date), # count those diagnosed with ckd on day of diagnosis
        # Exclusion 5: liver cirrhosis prior to or on pandemic start
        prior_cirrhosis = (!is.na(elig_date_liver_cirrhosis_first) 
                           & elig_date_liver_cirrhosis_first <= pandemicstart_date), # count those diagnosed with cirrhosis on day of diagnosis
        # Exclusion 6: prior drug with interaction risk with metfin, in 14 days window before or on pandemic start
        prior_interaction = (!is.na(elig_date_metfin_interaction_first) 
                             & (elig_date_metfin_interaction_first <= pandemicstart_date) & (elig_date_metfin_interaction_first >= pandemicstart_date - days(14))),
        # Exclusion 7: Very high HbA1c (>75), prior or on T2DM diagnosis (max 2 years back)
        prior_high_hba1c = (!is.na(cov_cat_hba1c_mmol_mol_b)
                            & cov_cat_hba1c_mmol_mol_b == "above 75")
      )
    
    # Among these, count the exclusion criteria
    count <- data_filtered_T2DM %>%
      summarise(
        n_prior_metfin = sum(prior_metfin, na.rm = TRUE),
        n_prior_sulfo_mono = sum(prior_sulfo_mono, na.rm = TRUE),
        n_prior_dpp4_mono = sum(prior_dpp4_mono, na.rm = TRUE),
        n_prior_tzd_mono = sum(prior_tzd_mono, na.rm = TRUE),
        n_prior_sglt2_mono = sum(prior_sglt2_mono, na.rm = TRUE),
        n_prior_glp1_mono = sum(prior_glp1_mono, na.rm = TRUE),
        n_prior_megli_mono = sum(prior_megli_mono, na.rm = TRUE),
        n_prior_agi_mono = sum(prior_agi_mono, na.rm = TRUE),
        n_prior_insulin_mono = sum(prior_insulin_mono, na.rm = TRUE),
        n_prior_metfin_allergy = sum(prior_metfin_allergy, na.rm = TRUE),
        n_prior_ckd45 = sum(prior_ckd45, na.rm = TRUE),
        n_prior_cirrhosis = sum(prior_cirrhosis, na.rm = TRUE),
        n_prior_interaction = sum(prior_interaction, na.rm = TRUE),
        n_prior_high_hba1c = sum(prior_high_hba1c, na.rm = TRUE)
      )
    
    # Filter 2: apply all exclusion criteria (inclusion criteria T2DM applied above)
    data_filtered <- data_filtered_T2DM %>% # Output 1: filtered data
      filter(
        (!prior_metfin | is.na(prior_metfin)),
        (!prior_sulfo_mono | is.na(prior_sulfo_mono)),
        (!prior_dpp4_mono | is.na(prior_dpp4_mono)),
        (!prior_tzd_mono | is.na(prior_tzd_mono)),
        (!prior_sglt2_mono | is.na(prior_sglt2_mono)),
        (!prior_glp1_mono | is.na(prior_glp1_mono)),
        (!prior_megli_mono | is.na(prior_megli_mono)),
        (!prior_agi_mono | is.na(prior_agi_mono)),
        (!prior_insulin_mono | is.na(prior_insulin_mono)),
        (!prior_metfin_allergy | is.na(prior_metfin_allergy)),
        (!prior_ckd45 | is.na(prior_ckd45)),
        (!prior_cirrhosis | is.na(prior_cirrhosis)),
        (!prior_interaction | is.na(prior_interaction)),
        (!prior_high_hba1c | is.na(prior_high_hba1c))
      )
    
    n_after_exclusion_processing <- nrow(data_filtered)
    
    # Output 2: Count for flowchart, without redaction, and including pre/post processing counts
    out <- tibble(
      n_before_exclusion_processing = nrow(data_processed),
      n_t2dm = n_t2dm, # counted among all data_processed
      n_prior_metfin = count$n_prior_metfin, # counted only among n_t2dm
      n_prior_sulfo_mono = count$n_prior_sulfo_mono, # counted only among n_t2dm
      n_prior_dpp4_mono = count$n_prior_dpp4_mono, # counted only among n_t2dm
      n_prior_tzd_mono = count$n_prior_tzd_mono, # counted only among n_t2dm
      n_prior_sglt2_mono = count$n_prior_sglt2_mono, # counted only among n_t2dm
      n_prior_glp1_mono = count$n_prior_glp1_mono, # counted only among n_t2dm
      n_prior_megli_mono = count$n_prior_megli_mono, # counted only among n_t2dm
      n_prior_agi_mono = count$n_prior_agi_mono, # counted only among n_t2dm
      n_prior_insulin_mono = count$n_prior_insulin_mono, # counted only among n_t2dm
      n_prior_metfin_allergy = count$n_prior_metfin_allergy, # counted only among n_t2dm
      n_prior_ckd45 = count$n_prior_ckd45, # counted only among n_t2dm
      n_prior_cirrhosis = count$n_prior_cirrhosis, # counted only among n_t2dm
      n_prior_interaction = count$n_prior_interaction, # counted only among n_t2dm
      n_prior_high_hba1c = count$n_prior_high_hba1c, # counted only among n_t2dm
      n_after_exclusion_processing = n_after_exclusion_processing
    ) %>% 
      # pivot (for easier data review in L4)
      pivot_longer(
        cols = everything(),
        names_to = "Variable",
        values_to = "Value"
      )
    
    # Output 3: Count for flowchart, with redaction, and including pre/post processing counts
    labels <- c(
      n_before_exclusion_processing_midpoint6 = "Before applying eligibility criteria",
      n_t2dm_midpoint6 = "Total with T2DM between mid2018-mid2019",
      n_prior_metfin_midpoint6 = "Prior metformin use",
      n_prior_sulfo_mono_midpoint6 = "Prior sulfonylurea",
      n_prior_dpp4_mono_midpoint6 = "Prior DPP-4 inhibitor",
      n_prior_tzd_mono_midpoint6 = "Prior Thiazolidinedione",
      n_prior_sglt2_mono_midpoint6 = "Prior SGLT2 inhibitor",
      n_prior_glp1_mono_midpoint6 = "Prior GLP-1 receptor agonist",
      n_prior_megli_mono_midpoint6 = "Prior meglitinide",
      n_prior_agi_mono_midpoint6 = "Prior alpha-glucosidase inhibitor",
      n_prior_insulin_mono_midpoint6 = "Prior insulin",
      n_prior_metfin_allergy_midpoint6 = "Prior metformin allergy",
      n_prior_ckd45_midpoint6 = "Prior CKD stage 4/5",
      n_prior_cirrhosis_midpoint6 = "Prior cirrhosis",
      n_prior_interaction_midpoint6 = "Use of drugs with interaction potential in past 14 days",
      n_prior_high_hba1c_midpoint6 = "HbA1c > 75 mmol/mol",
      n_after_exclusion_processing_midpoint6 = "After applying eligibility criteria"
    )
    out_midpoint6 <- tibble(
      n_before_exclusion_processing_midpoint6 = fn_roundmid_any(nrow(data_processed), threshold),
      n_t2dm_midpoint6 = fn_roundmid_any(n_t2dm, threshold), # counted among all data_processed
      n_prior_metfin_midpoint6 = fn_roundmid_any(count$n_prior_metfin, threshold), # counted only among n_t2dm
      n_prior_sulfo_mono_midpoint6 = fn_roundmid_any(count$n_prior_sulfo_mono, threshold), # counted only among n_t2dm
      n_prior_dpp4_mono_midpoint6 = fn_roundmid_any(count$n_prior_dpp4_mono, threshold), # counted only among n_t2dm
      n_prior_tzd_mono_midpoint6 = fn_roundmid_any(count$n_prior_tzd_mono, threshold), # counted only among n_t2dm
      n_prior_sglt2_mono_midpoint6 = fn_roundmid_any(count$n_prior_sglt2_mono, threshold), # counted only among n_t2dm
      n_prior_glp1_mono_midpoint6 = fn_roundmid_any(count$n_prior_glp1_mono, threshold), # counted only among n_t2dm
      n_prior_megli_mono_midpoint6 = fn_roundmid_any(count$n_prior_megli_mono, threshold), # counted only among n_t2dm
      n_prior_agi_mono_midpoint6 = fn_roundmid_any(count$n_prior_agi_mono, threshold), # counted only among n_t2dm
      n_prior_insulin_mono_midpoint6 = fn_roundmid_any(count$n_prior_insulin_mono, threshold), # counted only among n_t2dm
      n_prior_metfin_allergy_midpoint6 = fn_roundmid_any(count$n_prior_metfin_allergy, threshold), # counted only among n_t2dm
      n_prior_ckd45_midpoint6 = fn_roundmid_any(count$n_prior_ckd45, threshold), # counted only among n_t2dm
      n_prior_cirrhosis_midpoint6 = fn_roundmid_any(count$n_prior_cirrhosis, threshold), # counted only among n_t2dm
      n_prior_interaction_midpoint6 = fn_roundmid_any(count$n_prior_interaction, threshold), # counted only among n_t2dm
      n_prior_high_hba1c_midpoint6 = fn_roundmid_any(count$n_prior_high_hba1c, threshold), # counted only among n_t2dm
      n_after_exclusion_processing_midpoint6 = fn_roundmid_any(n_after_exclusion_processing, threshold)
    ) %>% 
      # pivot (for easier data review in L4)
      pivot_longer(
        cols = everything(),
        names_to = "Variable",
        values_to = "Value"
      ) %>%
      mutate(Variable = labels[Variable])  # Replace variable names with labels
    
    # Return outputs as a list
    return(list(
      n_elig_excluded = out,
      n_elig_excluded_midpoint6 = out_midpoint6,
      data_processed = data_filtered
    ))
  
}