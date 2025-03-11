################################################################################
# Custom made function to assess the completeness criteria
################################################################################
fn_completeness_criteria_midpoint6 <- function(data_processed, threshold){

  # Since elig_date_t2dm is now the baseline date and comes from diabetes-algo, establish this as main denominator first
  data_processed <- data_processed %>%
    mutate(
      # Exclusion 1: no T2DM diagnosis
      no_t2dm = is.na(elig_date_t2dm))

  # Filter 1: main eligibility criteria: T2DM diagnosis
  data_filtered_T2DM <- data_processed %>%
    filter((!no_t2dm | is.na(no_t2dm)))
  
  n_t2dm <- nrow(data_filtered_T2DM)
  
  # Apply the rules
  data_filtered_T2DM <- data_filtered_T2DM %>%
    mutate(
      # Rule 1: not alive at elig_date_t2dm
      not_alive_at_baseline = qa_bin_was_alive == FALSE | is.na(qa_bin_was_alive), 
      # Rule 2: not adult at elig_date_t2dm
      not_adult_at_baseline = qa_bin_was_adult == FALSE | is.na(qa_bin_was_adult), 
      # Rule 3: nor female or male
      nor_female_or_male = qa_bin_is_female_or_male == FALSE | is.na(qa_bin_is_female_or_male),
      # Rule 4: no imd at elig_date_t2dm
      no_imd = qa_bin_known_imd == FALSE | cov_cat_deprivation_5 == "unknown",
      # Rule 5: no region at elig_date_t2dm
      no_region = is.na(cov_cat_region),
      # Rule 6: not registered elig_date_t2dm
      not_registered = qa_bin_was_registered == FALSE | is.na(qa_bin_was_registered)
    )
  
  # Count the rules
  count <- data_filtered_T2DM %>%
    summarise(
      n_not_alive_at_baseline = sum(not_alive_at_baseline, na.rm = TRUE),
      n_not_adult_at_baseline = sum(not_adult_at_baseline, na.rm = TRUE),
      n_nor_female_or_male = sum(nor_female_or_male, na.rm = TRUE),
      n_no_imd = sum(no_imd, na.rm = TRUE),
      n_no_region = sum(no_region, na.rm = TRUE),
      n_not_registered = sum(not_registered, na.rm = TRUE)
    )
  
  # Output 1: filtered data
  data_filtered <- data_filtered_T2DM %>% # Output 1: filtered data
    filter(
      (!not_alive_at_baseline | is.na(not_alive_at_baseline)),
      (!not_adult_at_baseline | is.na(not_adult_at_baseline)),
      (!nor_female_or_male | is.na(nor_female_or_male)),
      (!no_imd | is.na(no_imd)),
      (!no_region | is.na(no_region)),
      (!not_registered | is.na(not_registered))
    )
  
  n_after_exclusion_processing <- nrow(data_filtered)
  
  # Output 2: Count for flowchart, without redaction, including pre/post processing counts
  labels <- c(
    n_before_exclusion_processing = "Before applying the completeness criteria",
    n_t2dm = "Total with T2DM",
    n_not_alive_at_baseline = "Not alive at baseline",
    n_not_adult_at_baseline = "Not adult at baseline",
    n_nor_female_or_male = "Neither female nor male",
    n_no_imd = "No deprivation information",
    n_no_region = "No region information",
    n_not_registered = "Not registered with TPP for past >= 12 months",
    n_after_exclusion_processing = "After applying the completeness criteria"
  )
  out <- tibble(
    n_before_exclusion_processing = nrow(data_processed),
    n_t2dm = n_t2dm, # counted among all data_processed
    n_not_alive_at_baseline = count$n_not_alive_at_baseline, # counted only among all with a T2DM diagnosis
    n_not_adult_at_baseline = count$n_not_adult_at_baseline, # counted only among all with a T2DM diagnosis
    n_nor_female_or_male = count$n_nor_female_or_male, # counted only among all with a T2DM diagnosis
    n_no_imd = count$n_no_imd, # counted only among all with a T2DM diagnosis
    n_no_region = count$n_no_region, # counted only among all with a T2DM diagnosis
    n_not_registered = count$n_not_registered, # counted only among all with a T2DM diagnosis
    n_after_exclusion_processing = n_after_exclusion_processing
  ) %>% 
    # pivot (for easier data review in L4)
    pivot_longer(
      cols = everything(),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(Variable = labels[Variable])  # Replace variable names with labels
  
  # Output 3: Count for flowchart, with redaction, including pre/post processing counts
  labels <- c(
    n_before_exclusion_processing_midpoint6 = "Before applying the completeness criteria",
    n_t2dm_midpoint6 = "Total with T2DM",
    n_not_alive_at_baseline_midpoint6 = "Not alive at baseline",
    n_not_adult_at_baseline_midpoint6 = "Not adult at baseline",
    n_nor_female_or_male_midpoint6 = "Neither female nor male",
    n_no_imd_midpoint6 = "No deprivation information",
    n_no_region_midpoint6 = "No region information",
    n_not_registered_midpoint6 = "Not registered with TPP for past >= 12 months",
    n_after_exclusion_processing_midpoint6 = "After applying the completeness criteria"
  )
  out_midpoint6 <- tibble(
    n_before_exclusion_processing_midpoint6 = fn_roundmid_any(nrow(data_processed), threshold),
    n_t2dm_midpoint6 = fn_roundmid_any(n_t2dm, threshold),
    n_not_alive_at_baseline_midpoint6 = fn_roundmid_any(count$n_not_alive_at_baseline, threshold),
    n_not_adult_at_baseline_midpoint6 = fn_roundmid_any(count$n_not_adult_at_baseline, threshold),
    n_nor_female_or_male_midpoint6 = fn_roundmid_any(count$n_nor_female_or_male, threshold),
    n_no_imd_midpoint6 = fn_roundmid_any(count$n_no_imd, threshold),
    n_no_region_midpoint6 = fn_roundmid_any(count$n_no_region, threshold),
    n_not_registered_midpoint6 = fn_roundmid_any(count$n_not_registered, threshold),
    n_after_exclusion_processing_midpoint6 = fn_roundmid_any(n_after_exclusion_processing, threshold)
  ) %>% 
    pivot_longer(
      cols = everything(),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(Variable = labels[Variable])  # Replace variable names with labels
  
  # Return all output as a list
  return(list(
    n_completeness_excluded = out,
    n_completeness_excluded_midpoint6 = out_midpoint6,
    data_processed = data_filtered
  ))
  
}