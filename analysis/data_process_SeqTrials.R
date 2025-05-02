####
## This script does the following:
# 1. Import/extract feather dataset from OpenSAFELY
# 2. Basic type formatting of variables -> fn_extract_data.R()
# 3. Process covariates
# 4. Modify dummy data if run locally
# 5. Evaluate/apply the quality assurance criteria -> fn_quality_assurance
# 6. Evaluate/apply the completeness criteria: -> fn_completeness_criteria
# 7. Evaluate/apply the eligibility criteria: -> fn_elig_criteria
# 8. Restrict the dataset
# 9. Save the full dataset, the restricted dataset, the codebook and all flowchart tables
####

# Import libraries and user functions -------------------------------------
library('arrow')
library('readr')
library('here')
library('lubridate')
library('dplyr')
library('tidyr')
library('purrr')
library('forcats')
library('jsonlite')
library('skimr')
library('splines')
source(here::here("analysis", "functions", "fn_extract_data.R"))
source(here::here("analysis", "functions", "utility.R"))
source(here::here("analysis", "functions", "fn_quality_assurance_midpoint6.R"))
source(here::here("analysis", "functions", "fn_completeness_criteria_midpoint6.R"))
source(here::here("analysis", "functions", "fn_elig_criteria_midpoint6_SeqTrials.R"))


# Create directories for output -------------------------------------------
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_description_seqtrials"))


# Import dates ------------------------------------------------------------
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")
pandemicstart_date <- as.Date(study_dates$pandemicstart_date, format = "%Y-%m-%d")


# Define redaction threshold ----------------------------------------------
threshold <- 6


# Import the dataset and pre-process --------------------------------------
input_filename <- "dataset.arrow"
data_extracted <- fn_extract_data(input_filename)


# Process the data --------------------------------------------------------
data_processed <- data_extracted %>%
  mutate(
    cov_cat_age = cut(
      cov_num_age,
      breaks = c(18, 40, 60, 80, 110),
      labels = c("18-39", "40-59", "60-79", "80 or older"),
      right = FALSE),
    
    cov_cat_sex = fn_case_when(
      cov_cat_sex == "female" ~ "Female",
      cov_cat_sex == "male" ~ "Male",
      TRUE ~ NA_character_), # will have no missing, by definition (excluded)
    
    cov_cat_deprivation_5 = fn_case_when(
      cov_cat_deprivation_5 == "1 (most deprived)" ~ "1 (most deprived)",
      cov_cat_deprivation_5 == "2" ~ "2",
      cov_cat_deprivation_5 == "3" ~ "3",
      cov_cat_deprivation_5 == "4" ~ "4",
      cov_cat_deprivation_5 == "5 (least deprived)" ~ "5 (least deprived)",
      TRUE ~ NA_character_), # will have no missing, by definition (excluded)
    
    cov_cat_ethnicity = fn_case_when(
      cov_cat_ethnicity == "White" ~ "White",
      cov_cat_ethnicity == "Asian" ~ "Asian",
      cov_cat_ethnicity == "Black" ~ "Black",
      cov_cat_ethnicity == "Mixed" ~ "Mixed",
      cov_cat_ethnicity == "Other" ~ "Other",
      cov_cat_ethnicity == "Unknown" ~ "Unknown",
      TRUE ~ NA_character_), # will have no missing
    
    strat_cat_region = fn_case_when(
      strat_cat_region == "East" ~ "East",
      strat_cat_region == "London" ~ "London",
      strat_cat_region %in% c("West Midlands", "East Midlands") ~ "Midlands",
      strat_cat_region %in% c("Yorkshire and The Humber", "North East") ~ "North East and Yorkshire",
      strat_cat_region == "North West" ~ "North West",
      strat_cat_region == "South East" ~ "South East",
      strat_cat_region == "South West" ~ "South West",
      TRUE ~ NA_character_), # will have no missing, by definition (excluded)
    
    cov_cat_rural_urban = fn_case_when(
      cov_cat_rural_urban %in% c(1,2) ~ "Urban conurbation",
      cov_cat_rural_urban %in% c(3,4) ~ "Urban city or town",
      cov_cat_rural_urban %in% c(5,6,7,8) ~ "Rural town or village",
      TRUE ~ NA_character_), # will have no missing
    
    cov_cat_smoking_status = fn_case_when(
      cov_cat_smoking_status == "S" ~ "Smoker",
      cov_cat_smoking_status == "E" | cov_cat_smoking_status == "M" ~ "Ever" , # collapsed due to very few Unknown (M), to avoid variable exclusion in cox model
      cov_cat_smoking_status == "N" ~ "Never",
      TRUE ~ NA_character_), # will have no missing
    
    cov_bin_obesity_b = cov_bin_obesity_b == TRUE | cov_cat_bmi_groups_b == "Obese (>30)",
    
    # TC/HDL ratio values: remove biologically implausible values: https://doi.org/10.1093/ije/dyz099
    ## remove TC < 1.75 or > 20; remove HDL < 0.4 or > 5; remove ratios < 1 or > 50
    tmp_cov_num_cholesterol_b = replace(tmp_cov_num_cholesterol_b, tmp_cov_num_cholesterol_b < 1.75 | tmp_cov_num_cholesterol_b > 20, NA_real_),
    tmp_cov_num_hdl_cholesterol_b = replace(tmp_cov_num_hdl_cholesterol_b, tmp_cov_num_hdl_cholesterol_b < 0.4 | tmp_cov_num_hdl_cholesterol_b > 5, NA_real_),
    cov_num_tc_hdl_ratio_b = tmp_cov_num_cholesterol_b / tmp_cov_num_hdl_cholesterol_b,
    cov_num_tc_hdl_ratio_b = replace(cov_num_tc_hdl_ratio_b, cov_num_tc_hdl_ratio_b > 50 | cov_num_tc_hdl_ratio_b < 1, NA_real_),
    
    # TC/HDL ratio categories: https://www.urmc.rochester.edu/encyclopedia/content?ContentTypeID=167&ContentID=lipid_panel_hdl_ratio#:~:text=Most%20healthcare%20providers%20want%20the,1%20is%20considered%20very%20good.
    cov_cat_tc_hdl_ratio_b = cut(
      cov_num_tc_hdl_ratio_b,
      breaks = c(1, 3.5, 5.11, 50), # 50 is upper limit, see above -> NA
      labels = c("below 3.5:1" ,"3.5:1 to 5:1", "above 5:1"),
      right = FALSE),
    cov_cat_tc_hdl_ratio_b = case_when(is.na(cov_num_tc_hdl_ratio_b) ~ factor("Unknown", 
                                                                          levels = c("below 3.5:1", "3.5:1 to 5:1", "above 5:1", "Unknown")), TRUE ~ cov_cat_tc_hdl_ratio_b),
    
    # HbA1c categories: https://www.southtees.nhs.uk/resources/the-hba1c-test/
    ## remove HbA1c > 120; remove HbA1c below 0
    cov_num_hba1c_mmol_mol_b = replace(cov_num_hba1c_mmol_mol_b, cov_num_hba1c_mmol_mol_b < 0.00 | cov_num_hba1c_mmol_mol_b > 120.00, NA_real_),
    cov_cat_hba1c_mmol_mol_b = cut(
      cov_num_hba1c_mmol_mol_b,
      breaks = c(0, 42, 59, 76, 120), # 120 is upper limit, above NA
      labels = c("below 42" ,"42-58", "59-75", "above 75"),
      right = FALSE),
    cov_cat_hba1c_mmol_mol_b = case_when(is.na(cov_cat_hba1c_mmol_mol_b) ~ factor("Unknown", 
                                                                              levels = c("below 42", "42-58", "59-75", "above 75", "Unknown")), TRUE ~ cov_cat_hba1c_mmol_mol_b)
    )


# Clarify mortality outcomes ----------------------------------------------
data_processed <- data_processed %>%
  mutate(
    out_bin_death = (!is.na(qa_date_of_death) & qa_date_of_death > pandemicstart_date) & is.na(out_date_covid_death),
    out_date_death = case_when(out_bin_death == TRUE ~ qa_date_of_death, 
                               TRUE ~ as.Date(NA))
  )


# Modify dummy data -------------------------------------------------------
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  message("Running locally, adapt dummy data")
  source("analysis/modify_dummy_data_SeqTrials.R")
  message("Dummy data successfully modified")
}


# Apply the quality assurance criteria ------------------------------------
qa <- fn_quality_assurance_midpoint6(data_processed, study_dates, threshold)
n_qa_excluded_midpoint6 <- qa$n_qa_excluded_midpoint6
data_processed <- qa$data_processed


# Apply the completeness criteria -----------------------------------------
completeness <- fn_completeness_criteria_midpoint6(data_processed, threshold)
n_completeness_excluded_midpoint6 <- completeness$n_completeness_excluded_midpoint6
data_processed <- completeness$data_processed


# Apply the eligibility criteria ------------------------------------------
# Our primary eligibility window to define incident T2DM is mid2018-mid2019, but maybe we may want to extend the window until max. mid2013 later on 
# => use function with loop that can be mapped to other windows, depending on "years_in_days" input
eligibility <- fn_elig_criteria_midpoint6_SeqTrials(data_processed, study_dates, years_in_days = 0)
n_elig_excluded <- eligibility$n_elig_excluded
n_elig_excluded_midpoint6 <- eligibility$n_elig_excluded_midpoint6
data_processed <- eligibility$data_processed


# Save and inspect full processed dataset ---------------------------------
data_processed_full <- data_processed
data_processed_full_desc <- skim(data_processed_full)
write.csv(data_processed_full_desc, file = here::here("output", "data_description_seqtrials", "data_processed_full_desc_seqtrials.csv")) # for L4 reviewing only, not for release
arrow::write_feather(data_processed_full, here::here("output", "data", "data_processed_full_seqtrials.arrow"))


# Restrict the dataset for pipeline onwards -------------------------------
data_processed <- data_processed %>% 
  select(patient_id, 
         elig_date_t2dm, 
         qa_date_of_death, # To identify other deaths, e.g. between eligibility and landmark
         starts_with("elig_"), # Eligibility
         starts_with("exp_"), # Exposures
         starts_with("cov_"), # Covariates
         starts_with("strat_"), # Stratification variable
         starts_with("out_"), # Outcomes
         starts_with("cens_") # Censoring event
  )


# Save aggregate output and restricted processed dataset -------------------
# flow chart quality assurance
write.csv(n_qa_excluded_midpoint6, file = here::here("output", "data_description_seqtrials", "n_qa_excluded_seqtrials_midpoint6.csv"))
# flow chart completeness criteria
write.csv(n_completeness_excluded_midpoint6, file = here::here("output", "data_description_seqtrials", "n_completeness_excluded_seqtrials_midpoint6.csv"))
# flow chart eligibility criteria
write.csv(n_elig_excluded_midpoint6, file = here::here("output", "data_description_seqtrials", "n_elig_excluded_seqtrials_midpoint6.csv"))
write.csv(n_elig_excluded, file = here::here("output", "data_description_seqtrials", "n_elig_excluded_seqtrials.csv"))
# final (restricted) dataset
arrow::write_feather(data_processed, here::here("output", "data", "data_processed_seqtrials.arrow"))