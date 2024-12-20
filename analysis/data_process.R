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
data_processed_dm_algo <- readRDS("data_processed_dm_algo.rds")
input_filename <- "dataset.arrow"

################################################################################
# 2 Reformat the imported data
################################################################################
data_extracted <- fn_extract_data(input_filename)

################################################################################
# 3 Process the data and apply diabetes algorithm
################################################################################
data_processed <- data_extracted %>%
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

# combine the two datasets data_processed_dm_algo and data_processed

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
n_completeness_excluded <- completeness$n_completeness_excluded
n_completeness_excluded_midpoint6 <- completeness$n_completeness_excluded_midpoint6
data_processed <- completeness$data_processed # CAVE: Being alive and registration based on mid2018, not landmark!

################################################################################
# 6 Apply the eligibility criteria
################################################################################
# Our primary eligibility window to define incident T2DM is mid2018-mid2019, but maybe we may want to extend the window until max. mid2013 later on => if so, use function with loop that can be mapped to other windows
eligibility <- fn_elig_criteria_midpoint6(data_processed, study_dates, years_in_days = 0)
n_elig_excluded <- eligibility$n_elig_excluded
n_elig_excluded_midpoint6 <- eligibility$n_elig_excluded_midpoint6
data_processed <- eligibility$data_processed

# # for feasibility check, apply the eligibility window for incident T2DM until mid2013
# data_processed_feasibility <- data_processed # since input and output are both data_processed, the loop would otherwise take the output of previous loop as input
# # count
# n_elig_excluded_all_windows <- 
#   map(.x = list(0, 366, 731, 1096, 1461, 1827), # define study window (mid_years 2018 until 2013)
#       .f = ~ fn_elig_criteria_midpoint6(data_processed_feasibility, study_dates, years_in_days = .x)$n_elig_excluded)
# names(n_elig_excluded_all_windows) <- c("elig_mid2018", "elig_mid2017", "elig_mid2016", "elig_mid2015", "elig_mid2014", "elig_mid2013")
# # count_midpoint6
# n_elig_excluded_all_windows_midpoint6 <- 
#   map(.x = list(0, 366, 731, 1096, 1461, 1827), # define study window (mid_years 2018 until 2013)
#       .f = ~ fn_elig_criteria_midpoint6(data_processed_feasibility, study_dates, years_in_days = .x)$n_elig_excluded_midpoint6)
# names(n_elig_excluded_all_windows_midpoint6) <- c("elig_mid2018_midpoint6", "elig_mid2017_midpoint6", "elig_mid2016_midpoint6", "elig_mid2015_midpoint6", "elig_mid2014_midpoint6", "elig_mid2013_midpoint6")
# 
# # apply eligibility criteria
# data_processed_all_windows <-
#   map(.x = list(0, 366, 731, 1096, 1461, 1827),
#       .f = ~ fn_elig_criteria_midpoint6(data_processed_feasibility, study_dates, years_in_days = .x)$data_processed)
# names(data_processed_all_windows) <- c("elig_mid2018", "elig_mid2017", "elig_mid2016", "elig_mid2015", "elig_mid2014", "elig_mid2013")

################################################################################
# 7 Double-check feasibility: Assign treatment/exposure and main outcome
################################################################################
# assign treatment/exposure and one outcome measure
data_processed <- data_processed %>% 
  # if on any metformin prescription in 6m prior to landmark (allow for T2DM diagnosis 6m prior)
  mutate(exp_bin_metfin_last = case_when(exp_date_metfin_last >= study_dates$landmark_date - days(183) ~ 1, 
                                   TRUE ~ 0),
         exp_date_metfin_last = case_when(exp_date_metfin_last >= study_dates$landmark_date - days(183) ~ exp_date_metfin_last, 
                                         TRUE ~ NA_Date_),
  # if started any metformin before landmark and still on it (prescription in 6m prior to landmark)
         exp_bin_metfin_last_first = case_when(exp_date_metfin_first <= study_dates$landmark_date 
                                               & exp_date_metfin_last >= study_dates$landmark_date - days(183) ~ 1, 
                                   TRUE ~ 0),
  # if started any metformin before landmark (those with exp_bin_metfin_first before mid2018 are already excluded)
  # take the difference to the two above and that's those that stopped again before landmark
         exp_bin_metfin_first = case_when(exp_date_metfin_first <= study_dates$landmark_date ~ 1, 
                                   TRUE ~ 0),
  # if on metformin monotherapy prescription in 6m prior to landmark
  # should be less than exp_bin_metfin_last
         exp_bin_metfin_mono_last = case_when(exp_date_metfin_mono_last >= study_dates$landmark_date - days(183) ~ 1, 
                                   TRUE ~ 0),
  # if started any metformin before OR after landmark (i.e. any initiation until study end date)
  # CAVE: irrespective those that stopped just before landmark
  # also calculate time between T2DM and metformin start
  # also check for metformin monotherapy
         exp_bin_metfin_any = case_when((exp_date_metfin_first <= study_dates$landmark_date)
                                          | (out_date_metfin_first > study_dates$landmark_date) ~ 1, 
                                   TRUE ~ 0),
         exp_date_metfin_any = case_when((exp_date_metfin_first <= study_dates$landmark_date)
                                          | (out_date_metfin_first > study_dates$landmark_date) ~ out_date_metfin_first, 
                                   TRUE ~ NA_Date_),
         exp_bin_metfin_mono = case_when((exp_date_metfin_mono_first <= study_dates$landmark_date)
                                  | (out_date_metfin_mono_first > study_dates$landmark_date) ~ 1, 
                                  TRUE ~ 0),
         tb_T2DMdiag_metfin_any = case_when(exp_bin_metfin_any == 1 ~ as.numeric(difftime(exp_date_metfin_any, elig_date_t2dm, units = "days")),
                                   TRUE ~ NA_real_),
         exp_bin_metfin_any_3m = case_when(tb_T2DMdiag_metfin_any <= 90 ~ 1,
                                   TRUE ~ 0),
         exp_bin_metfin_any_6m = case_when(tb_T2DMdiag_metfin_any <= 180 ~ 1,
                                   TRUE ~ 0),
         tb_T2DMdiag_metfin = case_when(exp_bin_metfin_last == 1 ~ as.numeric(difftime(exp_date_metfin_last, elig_date_t2dm, units = "days")),
                                     TRUE ~ NA_real_),
         exp_bin_metfin_3m = case_when(tb_T2DMdiag_metfin <= 90 ~ 1,
                                    TRUE ~ 0),
         exp_bin_metfin_6m = case_when(tb_T2DMdiag_metfin <= 180 ~ 1,
                                    TRUE ~ 0)) %>%
  # other antidiabetics (they include combo possibilities)
  mutate(exp_bin_sulfo = case_when(cov_date_sulfo >= study_dates$landmark_date - days(183) ~ 1,
                                   TRUE ~ 0),
         exp_bin_dpp4 = case_when(cov_date_dpp4 >= study_dates$landmark_date - days(183) ~ 1,
                                  TRUE ~ 0),
         exp_bin_dpp4_mono = case_when(cov_date_dpp4_mono >= study_dates$landmark_date - days(183) ~ 1,
                                  TRUE ~ 0),
         exp_bin_tzd = case_when(cov_date_tzd >= study_dates$landmark_date - days(183) ~ 1,
                                 TRUE ~ 0),
         exp_bin_tzd_mono = case_when(cov_date_tzd_mono >= study_dates$landmark_date - days(183) ~ 1,
                                 TRUE ~ 0),
         exp_bin_sglt2 = case_when(cov_date_sglt2 >= study_dates$landmark_date - days(183) ~ 1,
                                   TRUE ~ 0),
         exp_bin_sglt2_mono = case_when(cov_date_sglt2_mono >= study_dates$landmark_date - days(183) ~ 1,
                                   TRUE ~ 0),
         exp_bin_glp1 = case_when(cov_date_glp1 >= study_dates$landmark_date - days(183) ~ 1,
                                  TRUE ~ 0),
         exp_bin_megli = case_when(cov_date_megli >= study_dates$landmark_date - days(183) ~ 1,
                                   TRUE ~ 0),
         exp_bin_agi = case_when(cov_date_agi >= study_dates$landmark_date - days(183) ~ 1,
                                 TRUE ~ 0),
         exp_bin_insulin = case_when(cov_date_insulin >= study_dates$landmark_date - days(183) ~ 1,
                                     TRUE ~ 0)) %>%
  ## investigate treatment strategies
        # no antidiabetic medication at all at landmark
  mutate(exp_bin_treat_nothing = case_when(exp_bin_metfin_last == 0 # covers exp_bin_metfin_mono_last (sub-codelist of exp_bin_metfin_last)
                                           & exp_bin_sulfo == 0
                                           & exp_bin_dpp4 == 0
                                           & exp_bin_tzd == 0
                                           & exp_bin_sglt2 == 0
                                           & exp_bin_glp1 == 0
                                           & exp_bin_megli == 0
                                           & exp_bin_agi == 0
                                           & exp_bin_insulin == 0 ~ 1,
                                           TRUE ~ 0),
         # mono therapy metformin only (metformin combo excluded)
         exp_bin_treat_only_metfin = case_when(exp_bin_metfin_mono_last == 1 
                                               & exp_bin_sulfo == 0
                                               & exp_bin_dpp4 == 0
                                               & exp_bin_tzd == 0
                                               & exp_bin_sglt2 == 0
                                               & exp_bin_glp1 == 0
                                               & exp_bin_megli == 0
                                               & exp_bin_agi == 0
                                               & exp_bin_insulin == 0 ~ 1,
                                               TRUE ~ 0),
         # on metformin + sulfo (don't use metformin mono codelist, otherwise someone with a combo with sulfo would not be counted)  
         exp_bin_treat_metfin_sulfo = case_when(exp_bin_metfin_last == 1 
                                                & exp_bin_sulfo == 1 # e.g. glucovance (glibenclamide + metformin), but not commonly used anymore, due to increased risk of hypo (see: https://openprescribing.net/bnf/060102/) and therefore also not part of codelist
                                                & exp_bin_dpp4 == 0
                                                & exp_bin_tzd == 0
                                                & exp_bin_sglt2 == 0
                                                & exp_bin_glp1 == 0
                                                & exp_bin_megli == 0
                                                & exp_bin_agi == 0
                                                & exp_bin_insulin == 0 ~ 1,
                                                TRUE ~ 0),
         # on metformin + dpp4
         exp_bin_treat_metfin_dpp4 = case_when(exp_bin_metfin_last == 1 
                                               & exp_bin_sulfo == 0 
                                               & exp_bin_dpp4 == 1 # e.g. janumet (sitagliptin + metformin) or kombiglyze XR (saxagliptin + metformin) or jentadueto (linagliptin + metformin) or with vildagliptin or alogliptin
                                               & exp_bin_tzd == 0
                                               & exp_bin_sglt2 == 0
                                               & exp_bin_glp1 == 0
                                               & exp_bin_megli == 0
                                               & exp_bin_agi == 0
                                               & exp_bin_insulin == 0 ~ 1,
                                               TRUE ~ 0),
         # on metformin + tzd 
         exp_bin_treat_metfin_tzd = case_when(exp_bin_metfin_last == 1 
                                              & exp_bin_sulfo == 0 
                                              & exp_bin_dpp4 == 0 
                                              & exp_bin_tzd == 1 # e.g. actoplusmet (pioglitazone + metformin), maybe a few Metformin hydrochloride/rosiglitazone (0601023V0)
                                              & exp_bin_sglt2 == 0
                                              & exp_bin_glp1 == 0
                                              & exp_bin_megli == 0
                                              & exp_bin_agi == 0
                                              & exp_bin_insulin == 0 ~ 1,
                                              TRUE ~ 0),
         # on metformin + sglt2 
         exp_bin_treat_metfin_sglt2 = case_when(exp_bin_metfin_last == 1 
                                                & exp_bin_sulfo == 0 
                                                & exp_bin_dpp4 == 0 
                                                & exp_bin_tzd == 0 
                                                & exp_bin_sglt2 == 1 # e.g. invokamet (canaglifozin + metformin) or xigduo XR (dapaglifozin + metformin) or synjardy (Empagliflozin/metformin (0601023AR))
                                                & exp_bin_glp1 == 0
                                                & exp_bin_megli == 0
                                                & exp_bin_agi == 0
                                                & exp_bin_insulin == 0 ~ 1,
                                                TRUE ~ 0),
         # on metformin + insulin
         exp_bin_treat_metfin_insulin = case_when(exp_bin_metfin_last == 1 
                                                  & exp_bin_sulfo == 0 
                                                  & exp_bin_dpp4 == 0 
                                                  & exp_bin_tzd == 0
                                                  & exp_bin_sglt2 == 0
                                                  & exp_bin_glp1 == 0
                                                  & exp_bin_megli == 0
                                                  & exp_bin_agi == 0
                                                  & exp_bin_insulin == 1 ~ 1,
                                                  TRUE ~ 0),
         # SGLT2 mono
         exp_bin_treat_sglt2_mono = case_when(exp_bin_metfin_last == 0 
                                                  & exp_bin_sulfo == 0 
                                                  & exp_bin_dpp4 == 0 
                                                  & exp_bin_tzd == 0
                                                  & exp_bin_sglt2_mono == 1
                                                  & exp_bin_glp1 == 0
                                                  & exp_bin_megli == 0
                                                  & exp_bin_agi == 0
                                                  & exp_bin_insulin == 0 ~ 1,
                                                  TRUE ~ 0),
         # DPP4 mono
         exp_bin_treat_dpp4_mono = case_when(exp_bin_metfin_last == 0 
                                              & exp_bin_sulfo == 0 
                                              & exp_bin_dpp4_mono == 1 
                                              & exp_bin_tzd == 0
                                              & exp_bin_sglt2 == 0
                                              & exp_bin_glp1 == 0
                                              & exp_bin_megli == 0
                                              & exp_bin_agi == 0
                                              & exp_bin_insulin == 0 ~ 1,
                                              TRUE ~ 0),
         # TZD mono
         exp_bin_treat_tzd_mono = case_when(exp_bin_metfin_last == 0 
                                              & exp_bin_sulfo == 0 
                                              & exp_bin_dpp4 == 0 
                                              & exp_bin_tzd_mono == 1
                                              & exp_bin_sglt2 == 0
                                              & exp_bin_glp1 == 0
                                              & exp_bin_megli == 0
                                              & exp_bin_agi == 0
                                              & exp_bin_insulin == 0 ~ 1,
                                              TRUE ~ 0),
         # Sulfo mono
         exp_bin_treat_sulfo_mono = case_when(exp_bin_metfin_last == 0 
                                              & exp_bin_sulfo == 1 # does not contain any combinations with metformin (no FDCs!)
                                              & exp_bin_dpp4 == 0 
                                              & exp_bin_tzd == 0
                                              & exp_bin_sglt2 == 0
                                              & exp_bin_glp1 == 0
                                              & exp_bin_megli == 0
                                              & exp_bin_agi == 0
                                              & exp_bin_insulin == 0 ~ 1,
                                              TRUE ~ 0),
         # Insulin mono
         exp_bin_treat_insulin_mono = case_when(exp_bin_metfin_last == 0 
                                              & exp_bin_sulfo == 0 
                                              & exp_bin_dpp4 == 0 
                                              & exp_bin_tzd == 0
                                              & exp_bin_sglt2 == 0
                                              & exp_bin_glp1 == 0
                                              & exp_bin_megli == 0
                                              & exp_bin_agi == 0
                                              & exp_bin_insulin == 1 ~ 1, # does not contain any combinations with metformin (no FDCs!)
                                              TRUE ~ 0)) %>%
  ## add primary outcome
  mutate(out_bin_severecovid = case_when(out_date_covid19_severe > study_dates$landmark_date ~ 1, # severe covid outcome (hosp or death)
                                         TRUE ~ 0))

n_exp_out <- data_processed %>% 
  summarise(
    n_exp_bin_metfin_last = sum(exp_bin_metfin_last), 
    n_exp_bin_metfin_last_first = sum(exp_bin_metfin_last_first), 
    n_exp_bin_metfin_first = sum(exp_bin_metfin_first), 
    n_exp_bin_metfin_mono_last = sum(exp_bin_metfin_mono_last), 
    n_exp_bin_metfin_any = sum(exp_bin_metfin_any), 
    n_exp_bin_metfin_mono = sum(exp_bin_metfin_mono), 
    n_exp_bin_metfin_any_3m = sum(exp_bin_metfin_any_3m), 
    n_exp_bin_metfin_any_6m = sum(exp_bin_metfin_any_6m), 
    n_exp_bin_metfin_3m = sum(exp_bin_metfin_3m), 
    n_exp_bin_metfin_6m = sum(exp_bin_metfin_6m), 
    
    n_exp_bin_sulfo = sum(exp_bin_sulfo),
    n_exp_bin_dpp4 = sum(exp_bin_dpp4),
    n_exp_bin_dpp4_mono = sum(exp_bin_dpp4_mono), 
    n_exp_bin_tzd = sum(exp_bin_tzd),
    n_exp_bin_tzd_mono = sum(exp_bin_tzd_mono),
    n_exp_bin_sglt2 = sum(exp_bin_sglt2),
    n_exp_bin_sglt2_mono = sum(exp_bin_sglt2_mono),
    n_exp_bin_glp1 = sum(exp_bin_glp1),
    n_exp_bin_megli = sum(exp_bin_megli),
    n_exp_bin_agi = sum(exp_bin_agi),
    n_exp_bin_insulin = sum(exp_bin_insulin),
    
    n_exp_bin_treat_nothing = sum(exp_bin_treat_nothing),
    n_exp_bin_treat_only_metfin = sum(exp_bin_treat_only_metfin),
    n_exp_bin_treat_metfin_sulfo = sum(exp_bin_treat_metfin_sulfo),
    n_exp_bin_treat_metfin_dpp4 = sum(exp_bin_treat_metfin_dpp4),
    n_exp_bin_treat_metfin_tzd = sum(exp_bin_treat_metfin_tzd),
    n_exp_bin_treat_metfin_sglt2 = sum(exp_bin_treat_metfin_sglt2),
    n_exp_bin_treat_metfin_insulin = sum(exp_bin_treat_metfin_insulin),
    n_exp_bin_treat_sglt2_mono = sum(exp_bin_treat_sglt2_mono),
    n_exp_bin_treat_dpp4_mono = sum(exp_bin_treat_dpp4_mono),
    n_exp_bin_treat_tzd_mono = sum(exp_bin_treat_tzd_mono),
    n_exp_bin_treat_sulfo_mono = sum(exp_bin_treat_sulfo_mono),
    n_exp_bin_treat_insulin_mono = sum(exp_bin_treat_insulin_mono),
    n_out_severeCOVID = sum(out_bin_severecovid),
    
    median_tb_T2DMdiag_metfin_any = median(tb_T2DMdiag_metfin_any, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_any = quantile(tb_T2DMdiag_metfin_any, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_any = quantile(tb_T2DMdiag_metfin_any, 0.75, na.rm = TRUE),
    median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE)
    )

# midpoint6 rounded
n_exp_out_midpoint6 <- data_processed %>% 
  summarise(
    n_exp_bin_metfin_last_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_last, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_last_first_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_last_first, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_first_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_first, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_mono_last_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono_last, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_any_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_any, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono, na.rm = TRUE), threshold),
    n_exp_bin_metfin_any_3m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_any_3m, na.rm = TRUE), threshold),
    n_exp_bin_metfin_any_6m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_any_6m, na.rm = TRUE), threshold),
    n_exp_bin_metfin_3m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_3m, na.rm = TRUE), threshold),
    n_exp_bin_metfin_6m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_6m, na.rm = TRUE), threshold),
    
    n_exp_bin_sulfo_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo, na.rm = TRUE), threshold), 
    n_exp_bin_dpp4_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4, na.rm = TRUE), threshold), 
    n_exp_bin_dpp4_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4_mono, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd_mono, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2_mono, na.rm = TRUE), threshold), 
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
    n_exp_bin_treat_sglt2_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_sglt2_mono, na.rm = TRUE), threshold), 
    n_exp_bin_treat_dpp4_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_dpp4_mono, na.rm = TRUE), threshold), 
    n_exp_bin_treat_tzd_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_tzd_mono, na.rm = TRUE), threshold), 
    n_exp_bin_treat_sulfo_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_sulfo_mono, na.rm = TRUE), threshold), 
    n_exp_bin_treat_insulin_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_insulin_mono, na.rm = TRUE), threshold), 
    n_out_severeCOVID_midpoint6 = fn_roundmid_any(sum(out_bin_severecovid, na.rm = TRUE), threshold),
    
    median_tb_T2DMdiag_metfin_any = median(tb_T2DMdiag_metfin_any, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_any = quantile(tb_T2DMdiag_metfin_any, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_any = quantile(tb_T2DMdiag_metfin_any, 0.75, na.rm = TRUE),
    median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE)
  )


# n_exp_severecovid_midpoint6 <- map(
#   .x = data_processed_all_windows,
#   .f = ~ .x %>% 
#     # main exposure
#     mutate(exp_bin_treat = case_when(exp_date_metfin_first <= study_dates$landmark_date 
#                                      & exp_date_metfin_last >= study_dates$landmark_date - days(183) ~ 1, # 1 if started/treated/exposed and still on it
#                                      TRUE ~ 0)) %>% # 0 if not started/treated/exposed until landmark or not on it anymore
#     # investigate if on other antidiabetic just before landmark
#     mutate(exp_bin_sulfo = case_when(cov_date_sulfo >= study_dates$landmark_date - days(183) ~ 1,
#                                      TRUE ~ 0),
#     exp_bin_dpp4 = case_when(cov_date_dpp4 >= study_dates$landmark_date - days(183) ~ 1,
#                              TRUE ~ 0),
#     exp_bin_tzd = case_when(cov_date_tzd >= study_dates$landmark_date - days(183) ~ 1,
#                             TRUE ~ 0),
#     exp_bin_sglt2 = case_when(cov_date_sglt2 >= study_dates$landmark_date - days(183) ~ 1,
#                               TRUE ~ 0),
#     exp_bin_glp1 = case_when(cov_date_glp1 >= study_dates$landmark_date - days(183) ~ 1,
#                              TRUE ~ 0),
#     exp_bin_megli = case_when(cov_date_megli >= study_dates$landmark_date - days(183) ~ 1,
#                               TRUE ~ 0),
#     exp_bin_agi = case_when(cov_date_agi >= study_dates$landmark_date - days(183) ~ 1,
#                             TRUE ~ 0),
#     exp_bin_insulin = case_when(cov_date_insulin >= study_dates$landmark_date - days(183) ~ 1,
#                                 TRUE ~ 0)) %>%
#     # investigate combo antidiabetic
#     mutate(exp_bin_treat_nothing = case_when(exp_bin_treat == 0 # no antidiabetic medication at all at landmark
#                                                  & exp_bin_sulfo == 0
#                                                  & exp_bin_dpp4 == 0
#                                                  & exp_bin_tzd == 0
#                                                  & exp_bin_sglt2 == 0
#                                                  & exp_bin_glp1 == 0
#                                                  & exp_bin_megli == 0
#                                                  & exp_bin_agi == 0
#                                                  & exp_bin_insulin == 0 ~ 1,
#                                                  TRUE ~ 0),
#     exp_bin_treat_only_metfin = case_when(exp_bin_treat == 1 # mono therapy metformin only (CAVE: double-check codelist!)
#                                                  & exp_bin_sulfo == 0
#                                                  & exp_bin_dpp4 == 0
#                                                  & exp_bin_tzd == 0
#                                                  & exp_bin_sglt2 == 0
#                                                  & exp_bin_glp1 == 0
#                                                  & exp_bin_megli == 0
#                                                  & exp_bin_agi == 0
#                                                  & exp_bin_insulin == 0 ~ 1,
#                                                  TRUE ~ 0),
#     exp_bin_treat_metfin_sulfo = case_when(exp_bin_treat == 1 
#                                            & exp_bin_sulfo == 1 # most common, e.g. glucovance (glibenclamide + metformin)
#                                            & exp_bin_dpp4 == 0
#                                            & exp_bin_tzd == 0
#                                            & exp_bin_sglt2 == 0
#                                            & exp_bin_glp1 == 0
#                                            & exp_bin_megli == 0
#                                            & exp_bin_agi == 0
#                                            & exp_bin_insulin == 0 ~ 1,
#                                            TRUE ~ 0),
#     exp_bin_treat_metfin_dpp4 = case_when(exp_bin_treat == 1 
#                                           & exp_bin_sulfo == 0 
#                                           & exp_bin_dpp4 == 1 # most common, e.g. janumet (sitagliptin + metformin) or kombiglyze XR (saxagliptin + metformin) or jentadueto (linagliptin + metformin)
#                                           & exp_bin_tzd == 0
#                                           & exp_bin_sglt2 == 0
#                                           & exp_bin_glp1 == 0
#                                           & exp_bin_megli == 0
#                                           & exp_bin_agi == 0
#                                           & exp_bin_insulin == 0 ~ 1,
#                                           TRUE ~ 0),
#     exp_bin_treat_metfin_tzd = case_when(exp_bin_treat == 1 
#                                          & exp_bin_sulfo == 0 
#                                          & exp_bin_dpp4 == 0 
#                                          & exp_bin_tzd == 1 # most common, e.g. actoplusmet (pioglitazone + metformin)
#                                          & exp_bin_sglt2 == 0
#                                          & exp_bin_glp1 == 0
#                                          & exp_bin_megli == 0
#                                          & exp_bin_agi == 0
#                                          & exp_bin_insulin == 0 ~ 1,
#                                          TRUE ~ 0),
#     exp_bin_treat_metfin_sglt2 = case_when(exp_bin_treat == 1 
#                                            & exp_bin_sulfo == 0 
#                                            & exp_bin_dpp4 == 0 
#                                            & exp_bin_tzd == 0 
#                                            & exp_bin_sglt2 == 1 # most common, e.g. invokamet (canaglifozin + metformin) or xigduo XR (dapaglifozin + metformin)
#                                            & exp_bin_glp1 == 0
#                                            & exp_bin_megli == 0
#                                            & exp_bin_agi == 0
#                                            & exp_bin_insulin == 0 ~ 1,
#                                            TRUE ~ 0),
#     exp_bin_treat_metfin_insulin = case_when(exp_bin_treat == 1 
#                                          & exp_bin_sulfo == 0 
#                                          & exp_bin_dpp4 == 0 
#                                          & exp_bin_tzd == 0
#                                          & exp_bin_sglt2 == 0
#                                          & exp_bin_glp1 == 0
#                                          & exp_bin_megli == 0
#                                          & exp_bin_agi == 0
#                                          & exp_bin_insulin == 1 ~ 1, # with insulin
#                                          TRUE ~ 0)) %>%
#     # add primary outcome
#     mutate(out_bin_severecovid = case_when(out_date_covid19_severe > study_dates$landmark_date ~ 1, # 1 if severe covid outcome
#                                            TRUE ~ 0)) %>%                     # 0 if no severe covid outcome
#     
#     # summarise everything
#     summarise(
#       n_exp_bin_treat_midpoint6 = fn_roundmid_any(sum(exp_bin_treat, na.rm = TRUE), threshold), 
#       n_exp_bin_sulfo_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo, na.rm = TRUE), threshold),
#       n_exp_bin_dpp4_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4, na.rm = TRUE), threshold),
#       n_exp_bin_tzd_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd, na.rm = TRUE), threshold),
#       n_exp_bin_sglt2_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2, na.rm = TRUE), threshold),
#       n_exp_bin_glp1_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1, na.rm = TRUE), threshold),
#       n_exp_bin_megli_midpoint6 = fn_roundmid_any(sum(exp_bin_megli, na.rm = TRUE), threshold),
#       n_exp_bin_agi_midpoint6 = fn_roundmid_any(sum(exp_bin_agi, na.rm = TRUE), threshold),
#       n_exp_bin_insulin_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin, na.rm = TRUE), threshold),
#       n_exp_bin_treat_nothing_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_nothing, na.rm = TRUE), threshold),
#       n_exp_bin_treat_only_metfin_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_only_metfin, na.rm = TRUE), threshold),
#       n_exp_bin_treat_metfin_sulfo_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_metfin_sulfo, na.rm = TRUE), threshold),
#       n_exp_bin_treat_metfin_dpp4_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_metfin_dpp4, na.rm = TRUE), threshold),
#       n_exp_bin_treat_metfin_tzd_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_metfin_tzd, na.rm = TRUE), threshold),
#       n_exp_bin_treat_metfin_sglt2_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_metfin_sglt2, na.rm = TRUE), threshold),
#       n_exp_bin_treat_metfin_insulin_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_metfin_insulin, na.rm = TRUE), threshold),
#       n_severeCOVID_midpoint6 = fn_roundmid_any(sum(out_bin_severecovid, na.rm = TRUE), threshold)
#       )
#   )
# names(n_exp_severecovid_midpoint6) <- c("treat_outcome_mid2018_midpoint6", "treat_outcome_mid2017_midpoint6", "treat_outcome_mid2016_midpoint6", "treat_outcome_mid2015_midpoint6", "treat_outcome_mid2014_midpoint6", "treat_outcome_mid2013_midpoint6")

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
write.csv(n_elig_excluded_midpoint6, file = here::here("output", "data_properties", "n_elig_excluded_midpoint6.csv"))
write.csv(n_elig_excluded, file = here::here("output", "data_properties", "n_elig_excluded.csv"))
# descriptive exposure and 1 outcome
write.csv(n_exp_out_midpoint6, file = here::here("output", "data_properties", "n_exp_out_midpoint6.csv"))
write.csv(n_exp_out, file = here::here("output", "data_properties", "n_exp_out.csv"))

# purrr::walk2(
#   .x = n_elig_excluded_all_windows_midpoint6, 
#   .y = paste0(names(n_elig_excluded_all_windows_midpoint6), ".csv"),
#   .f = ~ write.csv(.x, 
#                    file = here::here("output", "data_properties", .y), 
#                    row.names = FALSE)
# )
# # Just to double-check re feasibility: Assign treatment/exposure and main outcome to above data frames going back 6 years
# purrr::walk2(
#   .x = n_exp_severecovid_midpoint6, 
#   .y = paste0(names(n_exp_severecovid_midpoint6), ".csv"),
#   .f = ~ write.csv(.x, 
#                    file = here::here("output", "data_properties", .y), 
#                    row.names = FALSE)
# )
