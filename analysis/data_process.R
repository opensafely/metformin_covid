################################################################################
# This script does the following:
# 1. Import/extract feather dataset from OpenSAFELY including basic type formatting of variables (extract_data())
# 2. Process the data
# 3. Save the processed data
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
data_extracted <- fn_extract_data(input_filename)

################################################################################
# 2 Process the data, and apply diabetes algorithm
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
# 3 Apply quality criteria
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
# 4 Apply eligibility criteria
################################################################################
# 4a. Completeness criteria
n_completeness_excluded_midpoint6 <- fn_completeness_criteria_midpoint6(data_processed)
data_processed <- data_processed %>%
  filter(qa_bin_was_alive == TRUE) %>%
  filter(qa_bin_is_female_or_male == TRUE) %>%
  filter(qa_bin_known_imd == TRUE) %>%
  filter(!is.na(cov_cat_region)) %>%
  filter(qa_bin_was_registered == TRUE)

# 4b. Eligibility criteria
# Our primary eligibility window to define incident T2DM is mid2018-mid2019, but maybe we may want to extend the window until max. mid2013 later on => use function with loop that can be mapped to other windows
# years_in_days <- c(0, 366, 731, 1096, 1461, 1827) # define study window (mid_years until 2013)

# assign eligibility flow chart 
n_elig_excluded_midpoint6 <- fn_elig_criteria_midpoint6(data_processed, study_dates, years_in_days = 0) # for flow chart
# assign eligibility flow chart function to all windows
n_elig_excluded_all_windows_midpoint6 <-
  map(.x = list(0, 366, 731, 1096, 1461, 1827), # define study window (mid_years 2018 until 2013)
      .f = ~ fn_elig_criteria_midpoint6(data_processed, study_dates, years_in_days = .x))
names(n_elig_excluded_all_windows_midpoint6) <- c("mid2018", "mid2017", "mid2016", "mid2015", "mid2014", "mid2013")

# apply eligibility criteria to define final dataset; for now, only for primary window (mid2018-mid2019 -> years_in_days = 0) 
data_processed <- fn_apply_elig_criteria(data_processed, study_dates, years_in_days = 0)
# apply eligibility criteria function to all windows
# data_processed_all_windows <-
#   map(.x = list(0, 366, 731, 1096, 1461, 1827),
#       .f = ~ fn_apply_elig_criteria(data_processed, study_dates, years_in_days = .x))
# names(data_processed_all_windows) <- c("mid2018", "mid2017", "mid2016", "mid2015", "mid2014", "mid2013")

################################################################################
# 5 Assign treatment/exposure 
################################################################################
# assign treatment/exposure; for now, only for primary window (mid2018-mid2019 -> data_processed), but can extend/map to the other windows
data_processed <- data_processed %>% 
  mutate(exp_bin_treat = case_when(exp_date_metfin_first <= study_dates$landmark_date ~ 1, # 1 if started/treated/exposed
                                   is.na(exp_date_metfin_first) ~ 0, # 0 if not started/treated/exposed until landmark
                                   TRUE ~ NA_real_)) %>% 
  mutate(tb_T2DMdiag_metfin = case_when(exp_bin_treat == 1 ~ as.numeric(difftime(exp_date_metfin_first, elig_date_t2dm, units = "days")),
                                        TRUE ~ NA_real_))



################################################################################
# Save output
################################################################################
write_rds(data_processed, here::here("output", "data", "data_processed.rds"))
