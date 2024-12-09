################################################################################
## This script does the following:
# 1. Import/extract feather dataset from OpenSAFELY 
# 2. Basic type formatting of variables -> fn_extract_data.R()
# 3. Process some covariates and apply the diabetes algorithm -> fn_diabetes_algorithm()
# 4. Evaluate/apply the quality assurance criteria -> fn_quality_assurance_midpoint6()
# 5. Evaluate/apply the completeness criteria: -> fn_completeness_criteria_midpoint6()
# 6. Evaluate/apply the eligibility criteria: -> fn_elig_criteria_midpoint6()
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
## Import custom user functions and meta-dates
source(here::here("analysis", "functions", "fn_extract_data.R"))
source(here::here("analysis", "functions", "utility.R"))
source(here::here("analysis", "functions", "fn_diabetes_algorithm.R"))
source(here::here("analysis", "functions", "fn_quality_assurance_midpoint6.R"))
source(here::here("analysis", "functions", "fn_completeness_criteria_midpoint6.R"))
source(here::here("analysis", "functions", "fn_elig_criteria_midpoint6.R"))

################################################################################
# 0.1 Create directories for output
################################################################################
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_properties"))

################################################################################
# 0.2 Import command-line arguments and dates # to be adapted at a later stage
################################################################################
args <- commandArgs(trailingOnly=TRUE)
# study_dates <-
#    jsonlite::read_json(path = here::here("output", "study_dates.json")) %>%
#    map(as.Date)
source(here::here("analysis", "metadates.R"))
# Convert the meta-dates into Date objects
study_dates <- lapply(study_dates, function(x) as.Date(x))

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
    tmp_cov_num_cholesterol = replace(tmp_cov_num_cholesterol, tmp_cov_num_cholesterol < 1.75 | tmp_cov_num_cholesterol > 20, NA_real_),
    tmp_cov_num_hdl_cholesterol = replace(tmp_cov_num_hdl_cholesterol, tmp_cov_num_hdl_cholesterol < 0.4 | tmp_cov_num_hdl_cholesterol > 5, NA_real_),
    cov_num_tc_hdl_ratio = tmp_cov_num_cholesterol / tmp_cov_num_hdl_cholesterol,
    cov_num_tc_hdl_ratio = replace(cov_num_tc_hdl_ratio, cov_num_tc_hdl_ratio > 50 | cov_num_tc_hdl_ratio < 1, NA_real_),
    )

# apply diabetes algorithm and delete all helper variables (tmp & step) at the end
data_processed <- fn_diabetes_algorithm(data_extracted)

################################################################################
# 4 Apply the quality assurance criteria
################################################################################
qa <- fn_quality_assurance_midpoint6(data_processed, study_dates, threshold)
n_qa_excluded_midpoint6 <- qa$n_qa_excluded_midpoint6
data_processed <- qa$data_processed

################################################################################
# 5 Apply the completeness criteria
################################################################################
completeness <- fn_completeness_criteria_midpoint6(data_processed, threshold)
n_completeness_excluded_midpoint6 <- completeness$n_completeness_excluded_midpoint6
data_processed <- completeness$data_processed

################################################################################
# 6 Apply the eligibility criteria
################################################################################
# Our primary eligibility window to define incident T2DM is mid2018-mid2019, but maybe we may want to extend the window until max. mid2013 later on => use function with loop that can be mapped to other windows

# eligibility <- fn_elig_criteria_midpoint6(data_processed, study_dates, years_in_days = 0) # for flow chart, CAVE: midpoint is rounding, i.e., do not consider abs numbers to double-check code
# n_elig_excluded_midpoint6 <- eligibility$n_elig_excluded_midpoint6
# data_processed <- eligibility$data_processed

# for feasibility check, apply the eligibility window for incident T2DM until mid2013
data_processed_feasibility <- data_processed # since input and output are both data_processed, the loop would otherwise take the output of previous loop as input
# count
n_elig_excluded_all_windows_midpoint6 <- 
  map(.x = list(0, 366, 731, 1096, 1461, 1827), # define study window (mid_years 2018 until 2013)
      .f = ~ fn_elig_criteria_midpoint6(data_processed_feasibility, study_dates, years_in_days = .x)$n_elig_excluded_midpoint6)
names(n_elig_excluded_all_windows_midpoint6) <- c("elig_mid2018_midpoint6", "elig_mid2017_midpoint6", "elig_mid2016_midpoint6", "elig_mid2015_midpoint6", "elig_mid2014_midpoint6", "elig_mid2013_midpoint6")

# apply eligibility criteria
data_processed_all_windows <-
  map(.x = list(0, 366, 731, 1096, 1461, 1827),
      .f = ~ fn_elig_criteria_midpoint6(data_processed_feasibility, study_dates, years_in_days = .x)$data_processed)
names(data_processed_all_windows) <- c("elig_mid2018", "elig_mid2017", "elig_mid2016", "elig_mid2015", "elig_mid2014", "elig_mid2013")

################################################################################
# 7 Double-check feasibility: Assign treatment/exposure and main outcome
################################################################################
# assign treatment/exposure
n_exp_severecovid_midpoint6 <- map(
  .x = data_processed_all_windows,
  .f = ~ .x %>% 
    # main exposure
    mutate(exp_bin_treat = case_when(exp_date_metfin_first <= study_dates$landmark_date 
                                     & exp_date_metfin_last >= study_dates$landmark_date - days(183) ~ 1, # 1 if started/treated/exposed and still on it
                                     TRUE ~ 0)) %>% # 0 if not started/treated/exposed until landmark or not on it anymore
    # investigate if on other antidiabetic just before landmark
    mutate(exp_bin_sulfo = case_when(cov_date_sulfo >= study_dates$landmark_date - days(183) ~ 1,
                                     TRUE ~ 0),
    exp_bin_dpp4 = case_when(cov_date_dpp4 >= study_dates$landmark_date - days(183) ~ 1,
                             TRUE ~ 0),
    exp_bin_tzd = case_when(cov_date_tzd >= study_dates$landmark_date - days(183) ~ 1,
                            TRUE ~ 0),
    exp_bin_sglt2 = case_when(cov_date_sglt2 >= study_dates$landmark_date - days(183) ~ 1,
                              TRUE ~ 0),
    exp_bin_glp1 = case_when(cov_date_glp1 >= study_dates$landmark_date - days(183) ~ 1,
                             TRUE ~ 0),
    exp_bin_megli = case_when(cov_date_megli >= study_dates$landmark_date - days(183) ~ 1,
                              TRUE ~ 0),
    exp_bin_agi = case_when(cov_date_agi >= study_dates$landmark_date - days(183) ~ 1,
                            TRUE ~ 0),
    exp_bin_insulin = case_when(cov_date_insulin >= study_dates$landmark_date - days(183) ~ 1,
                                TRUE ~ 0)) %>%
    # investigate combo antidiabetic
    mutate(exp_bin_treat_nothing = case_when(exp_bin_treat == 0 # no antidiabetic medication at all at landmark
                                                 & exp_bin_sulfo == 0
                                                 & exp_bin_dpp4 == 0
                                                 & exp_bin_tzd == 0
                                                 & exp_bin_sglt2 == 0
                                                 & exp_bin_glp1 == 0
                                                 & exp_bin_megli == 0
                                                 & exp_bin_agi == 0
                                                 & exp_bin_insulin == 0 ~ 1,
                                                 TRUE ~ 0),
    exp_bin_treat_only_metfin = case_when(exp_bin_treat == 1 # mono therapy metformin only (CAVE: double-check codelist!)
                                                 & exp_bin_sulfo == 0
                                                 & exp_bin_dpp4 == 0
                                                 & exp_bin_tzd == 0
                                                 & exp_bin_sglt2 == 0
                                                 & exp_bin_glp1 == 0
                                                 & exp_bin_megli == 0
                                                 & exp_bin_agi == 0
                                                 & exp_bin_insulin == 0 ~ 1,
                                                 TRUE ~ 0),
    exp_bin_treat_metfin_sulfo = case_when(exp_bin_treat == 1 
                                           & exp_bin_sulfo == 1 # most common, e.g. glucovance (glibenclamide + metformin)
                                           & exp_bin_dpp4 == 0
                                           & exp_bin_tzd == 0
                                           & exp_bin_sglt2 == 0
                                           & exp_bin_glp1 == 0
                                           & exp_bin_megli == 0
                                           & exp_bin_agi == 0
                                           & exp_bin_insulin == 0 ~ 1,
                                           TRUE ~ 0),
    exp_bin_treat_metfin_dpp4 = case_when(exp_bin_treat == 1 
                                          & exp_bin_sulfo == 0 
                                          & exp_bin_dpp4 == 1 # most common, e.g. janumet (sitagliptin + metformin) or kombiglyze XR (saxagliptin + metformin) or jentadueto (linagliptin + metformin)
                                          & exp_bin_tzd == 0
                                          & exp_bin_sglt2 == 0
                                          & exp_bin_glp1 == 0
                                          & exp_bin_megli == 0
                                          & exp_bin_agi == 0
                                          & exp_bin_insulin == 0 ~ 1,
                                          TRUE ~ 0),
    exp_bin_treat_metfin_tzd = case_when(exp_bin_treat == 1 
                                         & exp_bin_sulfo == 0 
                                         & exp_bin_dpp4 == 0 
                                         & exp_bin_tzd == 1 # most common, e.g. actoplusmet (pioglitazone + metformin)
                                         & exp_bin_sglt2 == 0
                                         & exp_bin_glp1 == 0
                                         & exp_bin_megli == 0
                                         & exp_bin_agi == 0
                                         & exp_bin_insulin == 0 ~ 1,
                                         TRUE ~ 0),
    exp_bin_treat_metfin_sglt2 = case_when(exp_bin_treat == 1 
                                           & exp_bin_sulfo == 0 
                                           & exp_bin_dpp4 == 0 
                                           & exp_bin_tzd == 0 
                                           & exp_bin_sglt2 == 1 # most common, e.g. invokamet (canaglifozin + metformin) or xigduo XR (dapaglifozin + metformin)
                                           & exp_bin_glp1 == 0
                                           & exp_bin_megli == 0
                                           & exp_bin_agi == 0
                                           & exp_bin_insulin == 0 ~ 1,
                                           TRUE ~ 0),
    exp_bin_treat_metfin_insulin = case_when(exp_bin_treat == 1 
                                         & exp_bin_sulfo == 0 
                                         & exp_bin_dpp4 == 0 
                                         & exp_bin_tzd == 0
                                         & exp_bin_sglt2 == 0
                                         & exp_bin_glp1 == 0
                                         & exp_bin_megli == 0
                                         & exp_bin_agi == 0
                                         & exp_bin_insulin == 1 ~ 1, # with insulin
                                         TRUE ~ 0)) %>%
    # add primary outcome
    mutate(out_bin_severecovid = case_when(out_date_covid19_severe > study_dates$landmark_date ~ 1, # 1 if severe covid outcome
                                           TRUE ~ 0)) %>%                     # 0 if no severe covid outcome
    
    # summarise everything
    summarise(
      n_exp_bin_treat_midpoint6 = fn_roundmid_any(sum(exp_bin_treat, na.rm = TRUE), threshold), 
      n_exp_bin_sulfo_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo, na.rm = TRUE), threshold),
      n_exp_bin_dpp4_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4, na.rm = TRUE), threshold),
      n_exp_bin_tzd_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd, na.rm = TRUE), threshold),
      n_exp_bin_sglt2_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2, na.rm = TRUE), threshold),
      n_exp_bin_glp1_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1, na.rm = TRUE), threshold),
      n_exp_bin_megli_midpoint6 = fn_roundmid_any(sum(exp_bin_megli, na.rm = TRUE), threshold),
      n_exp_bin_agi_midpoint6 = fn_roundmid_any(sum(exp_bin_agi, na.rm = TRUE), threshold),
      n_exp_bin_insulin_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin, na.rm = TRUE), threshold),
      n_exp_bin_treat_nothing_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_nothing, na.rm = TRUE), threshold),
      n_exp_bin_treat_only_metfin_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_only_metfin, na.rm = TRUE), threshold),
      n_exp_bin_treat_metfin_sulfo_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_metfin_sulfo, na.rm = TRUE), threshold),
      n_exp_bin_treat_metfin_dpp4_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_metfin_dpp4, na.rm = TRUE), threshold),
      n_exp_bin_treat_metfin_tzd_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_metfin_tzd, na.rm = TRUE), threshold),
      n_exp_bin_treat_metfin_sglt2_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_metfin_sglt2, na.rm = TRUE), threshold),
      n_exp_bin_treat_metfin_insulin_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_metfin_insulin, na.rm = TRUE), threshold),
      n_severeCOVID_midpoint6 = fn_roundmid_any(sum(out_bin_severecovid, na.rm = TRUE), threshold)
      )
  )
names(n_exp_severecovid_midpoint6) <- c("treat_outcome_mid2018_midpoint6", "treat_outcome_mid2017_midpoint6", "treat_outcome_mid2016_midpoint6", "treat_outcome_mid2015_midpoint6", "treat_outcome_mid2014_midpoint6", "treat_outcome_mid2013_midpoint6")

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
  .x = n_exp_severecovid_midpoint6, 
  .y = paste0(names(n_exp_severecovid_midpoint6), ".csv"),
  .f = ~ write.csv(.x, 
                   file = here::here("output", "data_properties", .y), 
                   row.names = FALSE)
)
