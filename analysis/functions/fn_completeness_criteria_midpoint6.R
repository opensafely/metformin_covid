################################################################################
# Custom made function to assess the completeness criteria
################################################################################
fn_completeness_criteria_midpoint6 <- function(data_processed, threshold){

  # Apply the rules
  data_processed <- data_processed %>%
    mutate(
      # Rule 1: not alive at landmark
      not_alive_at_landmark = qa_bin_was_alive == FALSE, 
      # Rule 2: not adult at landmark
      not_adult_at_landmark = qa_bin_was_adult == FALSE, 
      # Rule 3: nor female or male
      nor_female_or_male = qa_bin_is_female_or_male == FALSE,
      # Rule 4: no imd
      no_imd = qa_bin_known_imd == FALSE,
      # Rule 5: no region
      no_region = is.na(cov_cat_region),
      # Rule 6: not registered
      not_registered = qa_bin_was_registered == FALSE,
      # Additional conditions to search for deaths and LTFU between mid2018 and landmark
      not_alive_mid2018 = qa_bin_was_alive_mid2018 == FALSE, 
      not_registered_mid2018 = qa_bin_was_registered_mid2018 == FALSE
    )
  
  # Count the rules
  count <- data_processed %>%
    summarise(
      n_not_alive_at_landmark = sum(not_alive_at_landmark, na.rm = TRUE),
      n_not_adult_at_landmark = sum(not_adult_at_landmark, na.rm = TRUE),
      n_nor_female_or_male = sum(nor_female_or_male, na.rm = TRUE),
      n_no_imd = sum(no_imd, na.rm = TRUE),
      n_no_region = sum(no_region, na.rm = TRUE),
      n_not_registered = sum(not_registered, na.rm = TRUE),
      n_not_alive_mid2018 = sum(not_alive_mid2018, na.rm = TRUE),
      n_not_registered_mid2018 = sum(not_registered_mid2018, na.rm = TRUE)
    )
  
  # Filter
  data_filtered <- data_processed %>% # Output 1: filtered data: CAVE: Being alive and registration based on mid2018, not landmark!
    filter(
      #!not_alive_at_landmark,
      !not_adult_at_landmark,
      !nor_female_or_male,
      !no_imd,
      !is.na(no_region),
      #!not_registered,
      !not_alive_mid2018,
      !not_registered_mid2018
    )
  
  n_after_exclusion_processing <- nrow(data_filtered)
  
  # Output 2: Count for flowchart, without redaction, including pre/post processing counts
  out <- tibble(
    n_before_exclusion_processing = nrow(data_processed),
    #n_not_alive_at_landmark = count$n_not_alive_at_landmark,
    n_not_adult_at_landmark = count$n_not_adult_at_landmark,
    n_nor_female_or_male = count$n_nor_female_or_male,
    n_no_imd = count$n_no_imd,
    n_no_region = count$n_no_region,
    #n_not_registered = count$n_not_registered,
    n_not_alive_mid2018 = count$n_not_alive_mid2018,
    n_not_registered_mid2018 = count$n_not_registered_mid2018,
    n_after_exclusion_processing = n_after_exclusion_processing
  )
  
  # Output 3: Count for flowchart, with redaction, including pre/post processing counts
  out_midpoint6 <- tibble(
    n_before_exclusion_processing_midpoint6 = fn_roundmid_any(nrow(data_processed), threshold),
    #n_not_alive_at_landmark_midpoint6 = fn_roundmid_any(count$n_not_alive_at_landmark, threshold),
    n_not_adult_at_landmark_midpoint6 = fn_roundmid_any(count$n_not_adult_at_landmark, threshold),
    n_nor_female_or_male_midpoint6 = fn_roundmid_any(count$n_nor_female_or_male, threshold),
    n_no_imd_midpoint6 = fn_roundmid_any(count$n_no_imd, threshold),
    n_no_region_midpoint6 = fn_roundmid_any(count$n_no_region, threshold),
    #n_not_registered_midpoint6 = fn_roundmid_any(count$n_not_registered, threshold),
    n_not_alive_mid2018_midpoint6 = fn_roundmid_any(count$n_not_alive_mid2018, threshold),
    n_not_registered_mid2018_midpoint6 = fn_roundmid_any(count$n_not_registered_mid2018, threshold),
    n_after_exclusion_processing_midpoint6 = fn_roundmid_any(n_after_exclusion_processing, threshold)
  )
  
  # Return both outputs as a list
  return(list(
    n_completeness_excluded = out,
    n_completeness_excluded_midpoint6 = out_midpoint6,
    data_processed = data_filtered
  ))
  
}