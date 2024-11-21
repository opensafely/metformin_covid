################################################################################
## This script does the following:
# 1. Import/extract feather dataset from OpenSAFELY 
# 2. Basic type formatting of variables -> fn_extract_data.R()
# 3. Process some covariates and apply the diabetes algorithm -> fn_diabetes_algorithm()
# 4. Evaluate/apply the quality criteria -> fn_quality_assurance_midpoint6()
# 5. Evaluate/apply the completeness criteria: -> fn_completeness_criteria_midpoint6()
# 6. Evaluate/apply the eligibility criteria: -> fn_elig_criteria_midpoint6() | fn_apply_elig_criteria()
# (for now: just to double-check: Assign treatment and main outcome)
## Save the output: data_processed and the 1-row tables for the flow chart
################################################################################

################################################################################
# 0.0 Import libraries + functions
################################################################################
library('arrow')
library('readr')
library('here')
library('lubridate')
library('dplyr')
library('tidyr')
library('purrr')
library('forcats')
library('ggplot2')
## Import custom user functions
source(here::here("analysis", "functions", "fn_extract_data.R"))
source(here::here("analysis", "functions", "utility.R"))
source(here::here("analysis", "functions", "fn_diabetes_algorithm.R"))
source(here::here("analysis", "functions", "fn_quality_assurance_midpoint6.R"))
source(here::here("analysis", "functions", "fn_completeness_criteria_midpoint6.R"))
source(here::here("analysis", "functions", "eligibility_midpoint6.R"))

################################################################################
# 0.1 Create directories for output
################################################################################
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_properties"))

################################################################################
# 0.2 Import command-line arguments
################################################################################
args <- commandArgs(trailingOnly=TRUE)
study_dates <-
  jsonlite::read_json(path = here::here("output", "study_dates.json")) %>%
  map(as.Date)

################################################################################
# 0.3 Define redaction threshold
################################################################################
threshold <- 6

################################################################################
# 1 Import data
################################################################################
input_filename <- "dataset.arrow"

################################################################################
# 2 Reformat the imported data
################################################################################
data_extracted <- fn_extract_data(input_filename)

################################################################################
# 3 Process the data and apply diabetes algorithm
################################################################################
data_extracted <- data_extracted %>%
  mutate(
    # POPULATION/DEMOGRAPHIC ----
    cov_cat_age = cut(
      cov_num_age,
      breaks = c(18, 40, 60, 80, 110),
      labels = c("18-39", "40-59", "60-79", "80+"),
      right = FALSE),
    
    cov_cat_sex = fn_case_when(
      cov_cat_sex == "female" ~ "Female",
      cov_cat_sex == "male" ~ "Male",
      TRUE ~ NA_character_),
    
    cov_cat_region = fct_collapse(
      cov_cat_region,
      `East of England` = "East",
      `London` = "London",
      `Midlands` = c("West Midlands", "East Midlands"),
      `North East and Yorkshire` = c("Yorkshire and The Humber", "North East"),
      `North West` = "North West",
      `South East` = "South East",
      `South West` = "South West"),
    
    cov_cat_rural_urban = fn_case_when(
      cov_cat_rural_urban %in% c(1,2) ~ "Urban conurbation",
      cov_cat_rural_urban %in% c(3,4) ~ "Urban city or town",
      cov_cat_rural_urban %in% c(5,6,7,8) ~ "Rural town or village",
      TRUE ~ NA_character_),
    
    cov_cat_stp = as.factor(cov_cat_stp),
    
    cov_cat_ethnicity = fn_case_when(
      cov_cat_ethnicity == "1" ~ "White",
      cov_cat_ethnicity == "4" ~ "Black",
      cov_cat_ethnicity == "3" ~ "South Asian",
      cov_cat_ethnicity == "2" ~ "Mixed",
      cov_cat_ethnicity == "5" ~ "Other",
      cov_cat_ethnicity == "0" ~ "Unknown",
      TRUE ~ NA_character_),
    
    # Finalize smoking status
    cov_cat_smoking_status = fn_case_when(
      cov_cat_smoking_status == "S" ~ "Smoker",
      cov_cat_smoking_status == "E" ~ "Ever",
      cov_cat_smoking_status == "N" ~ "Never",
      cov_cat_smoking_status == "M" ~ "Unknown",
      TRUE ~ NA_character_),
    
    # Create one history of obesity variable
    cov_bin_obesity = fn_case_when(
      cov_bin_obesity == TRUE | cov_cat_bmi_groups == "Obese (>30)" ~ "Obese (>30)",
      TRUE ~ NA_character_),
    
    # Remove biologically implausible numerical BMI values
    cov_num_bmi = case_when(cov_num_bmi > 70 | cov_num_bmi < 12 ~ NA_real_, TRUE ~ cov_num_bmi),
    
    # TC/HDL ratio values: remove biologically implausible values: https://doi.org/10.1093/ije/dyz099
    ## remove TC < 1.75 or > 20
    ## remove HDL < 0.4 or > 5
    ## remove ratios < 1 or > 50
    tmp_cov_num_cholesterol = replace(tmp_cov_num_cholesterol, tmp_cov_num_cholesterol < 1.75 | tmp_cov_num_cholesterol > 20, NA),
    tmp_cov_num_hdl_cholesterol = replace(tmp_cov_num_hdl_cholesterol, tmp_cov_num_hdl_cholesterol < 0.4 | tmp_cov_num_hdl_cholesterol > 5, NA),
    cov_num_tc_hdl_ratio = tmp_cov_num_cholesterol / tmp_cov_num_hdl_cholesterol,
    cov_num_tc_hdl_ratio = replace(cov_num_tc_hdl_ratio, cov_num_tc_hdl_ratio > 50 | cov_num_tc_hdl_ratio < 1, NA),
  
    )

# apply diabetes algorithm and delete all helper variables (tmp & step) at the end
data_processed <- fn_diabetes_algorithm(data_extracted)

################################################################################
# 4 Apply the quality criteria
################################################################################
n_qa_excluded_midpoint6 <- fn_quality_assurance_midpoint6(data_processed) # 
data_processed <- data_processed %>%
  filter(!is.na(qa_num_birth_year)) %>%
  filter(is.na(qa_date_of_death) | (qa_num_birth_year <= year(qa_date_of_death))) %>%
  filter(qa_num_birth_year >= 1793 & qa_num_birth_year <= year(Sys.Date())) %>%
  filter((qa_date_of_death > as.Date("1900-01-01")) | (qa_date_of_death < Sys.Date()) | is.na(qa_date_of_death)) %>%
  filter((cov_cat_sex == "Female" | is.na(cov_cat_sex)) | (cov_cat_sex == "Male" & (qa_bin_pregnancy == FALSE))) %>% # FALSE includes missing
  filter((cov_cat_sex == "Female" | is.na(cov_cat_sex)) | (cov_cat_sex == "Male" & (qa_bin_hrt == FALSE)) | (cov_cat_sex == "Male" & (qa_bin_cocp == FALSE))) %>%
  filter((cov_cat_sex == "Male" | is.na(cov_cat_sex)) | (cov_cat_sex == "Female" & (qa_bin_prostate_cancer == FALSE)))

################################################################################
# 5 Apply the completeness criteria
################################################################################
n_completeness_excluded_midpoint6 <- fn_completeness_criteria_midpoint6(data_processed)
data_processed <- data_processed %>%
  filter(qa_bin_was_alive == TRUE) %>%
  filter(qa_bin_is_female_or_male == TRUE) %>%
  filter(qa_bin_known_imd == TRUE) %>%
  filter(!is.na(cov_cat_region)) %>%
  filter(qa_bin_was_registered == TRUE)

################################################################################
# 6 Apply the eligibility criteria
################################################################################
# Our primary eligibility window to define incident T2DM is mid2018-mid2019, but maybe we may want to extend the window until max. mid2013 later on => use function with loop that can be mapped to other windows
# years_in_days <- c(0, 366, 731, 1096, 1461, 1827) # define study window (mid_years until 2013)
# assign eligibility flow chart 
# n_elig_excluded_midpoint6 <- fn_elig_criteria_midpoint6(data_processed, study_dates, years_in_days = 0) # for flow chart
# assign eligibility flow chart function to all windows
n_elig_excluded_all_windows_midpoint6 <-
  map(.x = list(0, 366, 731, 1096, 1461, 1827), # define study window (mid_years 2018 until 2013)
      .f = ~ fn_elig_criteria_midpoint6(data_processed, study_dates, years_in_days = .x))
names(n_elig_excluded_all_windows_midpoint6) <- c("elig_mid2018_midpoint6", "elig_mid2017_midpoint6", "elig_mid2016_midpoint6", "elig_mid2015_midpoint6", "elig_mid2014_midpoint6", "elig_mid2013_midpoint6")

# apply eligibility criteria to define final dataset; for now, only for primary window (mid2018-mid2019 -> years_in_days = 0) 
# data_processed <- fn_apply_elig_criteria(data_processed, study_dates, years_in_days = 0)
# apply eligibility criteria function to all windows
data_processed_all_windows <-
  map(.x = list(0, 366, 731, 1096, 1461, 1827),
      .f = ~ fn_apply_elig_criteria(data_processed, study_dates, years_in_days = .x))
names(data_processed_all_windows) <- c("elig_mid2018", "elig_mid2017", "elig_mid2016", "elig_mid2015", "elig_mid2014", "elig_mid2013")

################################################################################
# 7 Double-check feasibility: Assign treatment/exposure and main outcome
################################################################################
# assign treatment/exposure
n_metfin_severecovid_midpoint6 <- map(
  .x = data_processed_all_windows,
  .f = ~ .x %>% 
    mutate(exp_bin_treat = case_when(
      exp_date_metfin_first <= study_dates$landmark_date ~ 1, # 1 if started/treated/exposed
      is.na(exp_date_metfin_first) ~ 0,                     # 0 if not started/treated/exposed until landmark
      TRUE ~ NA_real_)) %>% 
    mutate(out_bin_severecovid = case_when(
      out_date_covid19_severe > study_dates$landmark_date ~ 1, # 1 if severe covid outcome
      is.na(out_date_covid19_severe) ~ 0,                     # 0 if no severe covid outcome
      TRUE ~ NA_real_)) %>% 
    summarise(
      n_metfin_by_landmark_midpoint6 = fn_roundmid_any(sum(exp_bin_treat, na.rm = TRUE), threshold), 
      n_severeCOVID_midpoint6 = fn_roundmid_any(sum(out_bin_severecovid, na.rm = TRUE), threshold))
)
names(n_metfin_severecovid_midpoint6) <- c("treat_outcome_mid2018_midpoint6", "treat_outcome_mid2017_midpoint6", "treat_outcome_mid2016_midpoint6", "treat_outcome_mid2015_midpoint6", "treat_outcome_mid2014_midpoint6", "treat_outcome_mid2013_midpoint6")

################################################################################
# 8 Save output
################################################################################
# the data
write_rds(data_processed, here::here("output", "data", "data_processed.rds"))
# flow chart quality assurance
write.csv(n_qa_excluded_midpoint6, file = here::here("output", "data_properties", "n_qa_excluded_midpoint6.csv"))
# flow chart completeness criteria
write.csv(n_completeness_excluded_midpoint6, file = here::here("output", "data_properties", "n_completeness_excluded_midpoint6.csv"))
# flow chart eligibility criteria
purrr::walk2(
  .x = n_elig_excluded_all_windows_midpoint6, 
  .y = paste0(names(n_elig_excluded_all_windows_midpoint6), ".csv"),
  .f = ~ write.csv(.x, 
                   file = here::here("output", "data_properties", .y), 
                   row.names = FALSE)
)
# Just to double-check re feasibility: Assign treatment/exposure and main outcome to above data frames going back 6 years
purrr::walk2(
  .x = n_metfin_severecovid_midpoint6, 
  .y = paste0(names(n_metfin_severecovid_midpoint6), ".csv"),
  .f = ~ write.csv(.x, 
                   file = here::here("output", "data_properties", .y), 
                   row.names = FALSE)
)
