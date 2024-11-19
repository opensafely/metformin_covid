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
source(here::here("analysis", "functions", "fn_quality_assurance.R"))

################################################################################
# 0.1 Create directories for output
################################################################################
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_properties"))

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
n_qa_excluded <- fn_quality_assurance(data_processed) # 
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
n_excluded <- calc_n_excluded(data_processed$grace10)

data_processed_g10 <- data_processed_g10 %>%
  # completeness criteria
  filter(qa_bin_was_alive == TRUE & (qa_date_of_death > baseline_date | is.na(qa_date_of_death))) %>% # additional condition since "qa_bin_was_alive == TRUE" may not cover all (e.g. pos test came out after death)
  filter(qa_bin_is_female_or_male == TRUE) %>%
  filter(qa_bin_known_imd == TRUE) %>%
  filter(!is.na(cov_cat_region)) %>%
  filter(qa_bin_was_registered == TRUE) %>%
  # inclusion criteria
  filter(qa_bin_was_adult == TRUE) %>%
  filter(cov_bin_t2dm == TRUE) %>%
  filter(!is.na(baseline_date)) %>%
  # exclusion criteria
  filter(cov_bin_hosp_baseline == FALSE) %>% # FALSE includes missing in a ehrQL logical
  filter(cov_bin_metfin_before_baseline == FALSE) %>%
  filter(cov_bin_metfin_allergy == FALSE) %>%
  filter(cov_bin_ckd_45 == FALSE) %>%
  filter(cov_bin_liver_cirrhosis == FALSE) %>%
  filter(cov_bin_metfin_interaction == FALSE) %>%
  filter(cov_bin_long_covid == FALSE)

data_processed <-
  map(.x = data_processed,
      .f = ~ .x %>%
        # completeness criteria
        filter(qa_bin_was_alive == TRUE & (qa_date_of_death > baseline_date | is.na(qa_date_of_death))) %>% # additional condition since "qa_bin_was_alive == TRUE" may not cover all (e.g. pos test came out after death)
        filter(qa_bin_is_female_or_male == TRUE) %>%
        filter(qa_bin_known_imd == TRUE) %>%
        filter(!is.na(cov_cat_region)) %>%
        filter(qa_bin_was_registered == TRUE) %>%
        # inclusion criteria
        filter(qa_bin_was_adult == TRUE) %>%
        filter(cov_bin_t2dm == TRUE) %>%
        filter(!is.na(baseline_date)) %>%
        # exclusion criteria
        filter(cov_bin_hosp_baseline == FALSE) %>% # FALSE includes missing in a ehrQL logical
        filter(cov_bin_metfin_before_baseline == FALSE) %>%
        filter(cov_bin_metfin_allergy == FALSE) %>%
        filter(cov_bin_ckd_45 == FALSE) %>%
        filter(cov_bin_liver_cirrhosis == FALSE) %>%
        filter(cov_bin_metfin_interaction == FALSE) %>%
        filter(cov_bin_long_covid == FALSE)
  )



################################################################################
# 4 Apply eligibility criteria
################################################################################

################################################################################
# Save output
################################################################################
write_rds(data_processed, here::here("output", "data", "data_processed.rds"))
