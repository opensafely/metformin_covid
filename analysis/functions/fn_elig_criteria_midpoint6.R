################################################################################
# Custom made function to assess the eligibility criteria
################################################################################
fn_elig_criteria_midpoint6 <- function(data_processed, study_dates, years_in_days) {
  # Extract dates from the study_dates list
  pandemicstart_date <- as.Date(study_dates$pandemicstart_date, format = "%Y-%m-%d")
  mid2018_date <- as.Date(study_dates$mid2018_date, format = "%Y-%m-%d")
  mid2019_date <- as.Date(study_dates$mid2019_date, format = "%Y-%m-%d")
    
    # Apply the main eligibility criteria (Aged ≥ 18 years already applied in completeness criteria)
    data_processed <- data_processed %>%
      mutate(
        # Exclusion 1: no T2DM diagnosis or out of window
        no_t2dm_or_outofwindow = is.na(elig_date_t2dm) 
        | (elig_date_t2dm > mid2019_date) # in 6m window before pandemic start (or later, but should not occur, see dataset definition)
        | (elig_date_t2dm < mid2018_date - days(years_in_days))) # older than mid2018
    
    # Filter 1: main eligibility criteria: T2DM diagnosis
    data_filtered_T2DM <- data_processed %>%
      filter((!no_t2dm_or_outofwindow | is.na(no_t2dm_or_outofwindow)))
    
    n_t2dm <- nrow(data_filtered_T2DM)
    
    # Apply the other eligibility criteria (Aged ≥ 18 years & ≤ 85 years already applied in completeness criteria)
    data_filtered_T2DM <- data_filtered_T2DM %>%
      mutate(
        # Extended exclusions: 
        prior_hosp = elig_bin_hosp == TRUE,
        prior_carehome = elig_bin_carehome_status == "True",
        prior_palliative = (!is.na(elig_date_palliative)
                            & elig_date_palliative <= elig_date_t2dm),
        # Exclusion 2: metformin use prior to T2DM diagnosis (or any other antidiabetic!)
        prior_metfin = (!is.na(elig_date_metfin) # includes metformin combo, but update to elig_date_metfin later when new dataset available
                        & elig_date_metfin < elig_date_t2dm), # but don't count those who initiated on day of diagnosis & codelist allows for metformin including combo treatment
        prior_sulfo_mono = (!is.na(elig_date_sulfo)
                            & elig_date_sulfo < elig_date_t2dm),
        prior_dpp4_mono = (!is.na(elig_date_dpp4)
                           & elig_date_dpp4 < elig_date_t2dm),
        prior_tzd_mono = (!is.na(elig_date_tzd)
                          & elig_date_tzd < elig_date_t2dm),
        prior_sglt2_mono = (!is.na(elig_date_sglt2)
                            & elig_date_sglt2 < elig_date_t2dm),
        prior_glp1_mono = (!is.na(elig_date_glp1)
                           & elig_date_glp1 < elig_date_t2dm),
        prior_megli_mono = (!is.na(elig_date_megli)
                            & elig_date_megli < elig_date_t2dm),
        prior_agi_mono = (!is.na(elig_date_agi)
                          & elig_date_agi < elig_date_t2dm),
        prior_insulin_mono = (!is.na(elig_date_insulin)
                              & elig_date_insulin < elig_date_t2dm),
        # Exclusion 3: metformin allergy prior to or on T2DM diagnosis
        prior_metfin_allergy = (!is.na(elig_date_metfin_allergy_first) 
                                & elig_date_metfin_allergy_first <= elig_date_t2dm), # count those diagnosed with allergy on day of diagnosis
        # Exclusion 4: CKD 4/5 prior to or on T2DM diagnosis
        prior_ckd45 = (!is.na(elig_date_ckd_45_first) 
                       & elig_date_ckd_45_first <= elig_date_t2dm), # count those diagnosed with ckd on day of diagnosis
        # Exclusion 5: liver cirrhosis prior to or on T2DM diagnosis
        prior_cirrhosis = (!is.na(elig_date_liver_cirrhosis_first) 
                           & elig_date_liver_cirrhosis_first <= elig_date_t2dm), # count those diagnosed with cirrhosis on day of diagnosis
        # Exclusion 6: prior drug with interaction risk with metfin, in 14 days window before or on T2DM diagnosis day
        prior_interaction = (!is.na(elig_date_metfin_interaction_last) 
                             & (elig_date_metfin_interaction_last <= elig_date_t2dm) & (elig_date_metfin_interaction_last >= elig_date_t2dm - days(14))),
        # Exclusion 7: Very high HbA1c (>75), prior or on T2DM diagnosis (max 2 years back)
        prior_high_hba1c = (!is.na(cov_cat_hba1c_mmol_mol)
                            & cov_cat_hba1c_mmol_mol == "above 75")
      )
    
    # Re-apply the time-updated eligibility criteria again at landmark -> NO! Those who will be dead & LTFU by landmark, they will be out by design, that's ok to kick them out at this stage ("effect is conditional on being alive/in care at landmark"), but the other criteria should not be re-applied to respect alignment of eligibility and treatment assignment as best as possible (not re-apply elig criteria after treatment assignment which may happen anytime until landmark)
    # data_filtered_T2DM <- data_filtered_T2DM %>%
    #   mutate(
    #     # # Exclusion 8: metformin allergy prior to or on landmark, but after elig_date_t2dm (baseline) since all elig_date_t2dm in data_filtered_T2DM are in eligible time window (mid2018-2019)
    #     # prior_metfin_allergy_landmark = (!is.na(elig_date_metfin_allergy_first) 
    #     #                                  & elig_date_metfin_allergy_first > elig_date_t2dm # don't count those diagnosed with allergy on day of diagnosis again
    #     #                                  & elig_date_metfin_allergy_first <= elig_date_t2dm + days(183)),
    #     # # Exclusion 9: CKD 4/5 prior to or on landmark
    #     # prior_ckd45_landmark = (!is.na(elig_date_ckd_45_first) 
    #     #                         & elig_date_ckd_45_first > elig_date_t2dm
    #     #                         & elig_date_ckd_45_first <= elig_date_t2dm + days(183)),
    #     # # Exclusion 10: liver cirrhosis prior to or on landmark
    #     # prior_cirrhosis_landmark = (!is.na(elig_date_liver_cirrhosis_first) 
    #     #                             & elig_date_liver_cirrhosis_first > elig_date_t2dm 
    #     #                             & elig_date_liver_cirrhosis_first <= elig_date_t2dm + days(183)),
    #     # # Exclusion 11: prior drug with interaction risk with metfin, in 14 days window prior to or on landmark
    #     # prior_interaction_landmark = (!is.na(elig_date_metfin_interaction_landmark_last) 
    #     #                               & elig_date_metfin_interaction_landmark_last <= elig_date_t2dm + days(183) 
    #     #                               & elig_date_metfin_interaction_landmark_last >= elig_date_t2dm + days(169)),
    #     # Exclusion 12: died prior to landmark
    #     # prior_death_landmark = (!is.na(qa_date_of_death) 
    #     #                         & qa_date_of_death > elig_date_t2dm
    #     #                         & qa_date_of_death <= elig_date_t2dm + days(183)),
    #     # # Exclusion 13: LTFU prior to landmark
    #     # prior_ltfu_landmark = (!is.na(cens_date_dereg) 
    #     #                        & cens_date_dereg > elig_date_t2dm
    #     #                        & cens_date_dereg <= elig_date_t2dm + days(183))
    #     # # Exclusion 14: Very high HbA1c (>75), prior to landmark
    #     # prior_high_hba1c_landmark = (!is.na(elig_cat_hba1c_landmark_mmol_mol)
    #     #                     & elig_cat_hba1c_landmark_mmol_mol == "above 75")
    #   )
    
    # Among these, count the exclusion criteria
    count <- data_filtered_T2DM %>%
      summarise(
        n_prior_hosp = sum(prior_hosp, na.rm = TRUE),
        n_prior_carehome = sum(prior_carehome, na.rm = TRUE),
        n_prior_palliative = sum(prior_palliative, na.rm = TRUE),
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
        n_prior_high_hba1c = sum(prior_high_hba1c, na.rm = TRUE),
        # n_prior_metfin_allergy_landmark = sum(prior_metfin_allergy_landmark, na.rm = TRUE),
        # n_prior_ckd45_landmark = sum(prior_ckd45_landmark, na.rm = TRUE),
        # n_prior_cirrhosis_landmark = sum(prior_cirrhosis_landmark, na.rm = TRUE),
        # n_prior_interaction_landmark = sum(prior_interaction_landmark, na.rm = TRUE),
        # n_prior_death_landmark = sum(prior_death_landmark, na.rm = TRUE),
        # n_prior_ltfu_landmark = sum(prior_ltfu_landmark, na.rm = TRUE)
        # n_prior_high_hba1c_landmark = sum(prior_high_hba1c_landmark, na.rm = TRUE)
      )
    
    # Filter 2: apply all exclusion criteria (inclusion criteria T2DM applied above)
    data_filtered <- data_filtered_T2DM %>% # Output 1: filtered data
      filter(
        (!prior_hosp | is.na(prior_hosp)),
        (!prior_carehome | is.na(prior_carehome)),
        (!prior_palliative | is.na(prior_palliative)),
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
        (!prior_high_hba1c | is.na(prior_high_hba1c)),
        # (!prior_metfin_allergy_landmark | is.na(prior_metfin_allergy_landmark)),
        # (!prior_ckd45_landmark | is.na(prior_ckd45_landmark)),
        # (!prior_cirrhosis_landmark | is.na(prior_cirrhosis_landmark)),
        # (!prior_interaction_landmark | is.na(prior_interaction_landmark)),
        # (!prior_death_landmark | is.na(prior_death_landmark)),
        # (!prior_ltfu_landmark | is.na(prior_ltfu_landmark))
        # (!prior_high_hba1c_landmark | is.na(prior_high_hba1c_landmark))
      )
    
    n_after_exclusion_processing <- nrow(data_filtered)
    
    # Output 2: Count for flowchart, without redaction, and including pre/post processing counts
    out <- tibble(
      n_before_exclusion_processing = nrow(data_processed),
      n_t2dm = n_t2dm, # counted among all data_processed
      n_prior_hosp = count$n_prior_hosp, # counted only among n_t2dm
      n_prior_carehome = count$n_prior_carehome, # counted only among n_t2dm
      n_prior_palliative = count$n_prior_palliative, # counted only among n_t2dm
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
      # n_prior_metfin_allergy_landmark = count$n_prior_metfin_allergy_landmark, # counted only among n_t2dm
      # n_prior_ckd45_landmark = count$n_prior_ckd45_landmark, # counted only among n_t2dm
      # n_prior_cirrhosis_landmark = count$n_prior_cirrhosis_landmark, # counted only among n_t2dm
      # n_prior_interaction_landmark = count$n_prior_interaction_landmark, # counted only among n_t2dm
      # n_prior_death_landmark = count$n_prior_death_landmark, # counted only among n_t2dm
      # n_prior_ltfu_landmark = count$n_prior_ltfu_landmark, # counted only among n_t2dm
      # n_prior_high_hba1c_landmark = count$n_prior_high_hba1c_landmark, # counted only among n_t2dm
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
      n_prior_hosp_midpoint6 = "Hospitalized",
      n_prior_carehome_midpoint6 = "In care home",
      n_prior_palliative_midpoint6 = "In palliative care",
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
      # n_prior_metfin_allergy_landmark_midpoint6 = "New metformin allergy between baseline and landmark",
      # n_prior_ckd45_landmark_midpoint6 = "New CKD stage 4/5 between baseline and landmark",
      # n_prior_cirrhosis_landmark_midpoint6 = "New cirrhosis between baseline and landmark",
      # n_prior_interaction_landmark_midpoint6 = "New use of drugs with interaction potential in past 14 days",
      # n_prior_death_landmark_midpoint6 = "Death between baseline and landmark",
      # n_prior_ltfu_landmark_midpoint6 = "Deregistration between baseline and landmark",
      # n_prior_high_hba1c_landmark_midpoint6 = "HbA1c > 75 mmol/mol between baseline and landmark",
      n_after_exclusion_processing_midpoint6 = "After applying eligibility criteria"
    )
    out_midpoint6 <- tibble(
      n_before_exclusion_processing_midpoint6 = fn_roundmid_any(nrow(data_processed), threshold),
      n_t2dm_midpoint6 = fn_roundmid_any(n_t2dm, threshold), # counted among all data_processed
      n_prior_hosp_midpoint6 = fn_roundmid_any(count$n_prior_hosp, threshold), # counted only among n_t2dm
      n_prior_carehome_midpoint6 = fn_roundmid_any(count$n_prior_carehome, threshold), # counted only among n_t2dm
      n_prior_palliative_midpoint6 = fn_roundmid_any(count$n_prior_palliative, threshold), # counted only among n_t2dm
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
      # n_prior_metfin_allergy_landmark_midpoint6 = fn_roundmid_any(count$n_prior_metfin_allergy_landmark, threshold), # counted only among n_t2dm
      # n_prior_ckd45_landmark_midpoint6 = fn_roundmid_any(count$n_prior_ckd45_landmark, threshold), # counted only among n_t2dm
      # n_prior_cirrhosis_landmark_midpoint6 = fn_roundmid_any(count$n_prior_cirrhosis_landmark, threshold), # counted only among n_t2dm
      # n_prior_interaction_landmark_midpoint6 = fn_roundmid_any(count$n_prior_interaction_landmark, threshold), # counted only among n_t2dm
      # n_prior_death_landmark_midpoint6 = fn_roundmid_any(count$n_prior_death_landmark, threshold), # counted only among n_t2dm
      # n_prior_ltfu_landmark_midpoint6 = fn_roundmid_any(count$n_prior_ltfu_landmark, threshold), # counted only among n_t2dm
      # n_prior_high_hba1c_landmark_midpoint6 = fn_roundmid_any(count$n_prior_high_hba1c_landmark, threshold), # counted only among n_t2dm
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