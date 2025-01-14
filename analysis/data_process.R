################################################################################
## This script does the following:
# 1. Import/extract feather dataset from OpenSAFELY 
# 2. Basic type formatting of variables -> fn_extract_data.R()
# 3. Process some covariates
# 4. Import the processed dataset with the DM variables (and ethnicity and qa_num_birth_year) and merge
# 5. Evaluate/apply the quality assurance criteria -> fn_quality_assurance_midpoint6()
# 6. Evaluate/apply the completeness criteria: -> fn_completeness_criteria_midpoint6()
# 7. Evaluate/apply the eligibility criteria: -> fn_elig_criteria_midpoint6()
# 8. Assign treatment, various treatment regimen patterns and main outcome
# 9. Output for cumulative incidence plots re treatment regimen pattern
## Save all output
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
# 1 Import the dataset definition
################################################################################
input_filename <- "dataset.arrow"

################################################################################
# 2 Reformat the imported data
################################################################################
data_extracted <- fn_extract_data(input_filename)

################################################################################
# 3 Process the data
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

################################################################################
# 4 Import the processed DM algo dataset and merge
################################################################################
data_processed_dm_algo <- readRDS(here::here("output", "data", "data_processed_dm_algo.rds"))
data_processed <- merge(data_processed, data_processed_dm_algo, 
                        by = "patient_id", 
                        all.x = TRUE)

################################################################################
# 5 Apply the quality assurance criteria
################################################################################
qa <- fn_quality_assurance_midpoint6(data_processed, study_dates, threshold)
n_qa_excluded_midpoint6 <- qa$n_qa_excluded_midpoint6
data_processed <- qa$data_processed

################################################################################
# 6 Apply the completeness criteria
################################################################################
completeness <- fn_completeness_criteria_midpoint6(data_processed, threshold)
n_completeness_excluded <- completeness$n_completeness_excluded
n_completeness_excluded_midpoint6 <- completeness$n_completeness_excluded_midpoint6
data_processed <- completeness$data_processed # CAVE: Being alive and registration based on mid2018, not landmark!

################################################################################
# 7 Apply the eligibility criteria
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
# 8 Assign treatment/exposure and main outcome
################################################################################
# assign treatment/exposure and main outcome measure
data_processed <- data_processed %>% 
  mutate(
    # started any metformin before pandemic start, among those with a T2DM diagnosis, mid2018 onwards (those with exp_bin_metfin_first before mid2018 were already excluded above via fn_elig_criteria_midpoint6)
    exp_bin_metfin_first = case_when(exp_date_metfin_first <= study_dates$pandemicstart_date ~ 1, 
                                     TRUE ~ 0),
    exp_bin_metfin_mono_first = case_when(exp_date_metfin_mono_first <= study_dates$pandemicstart_date ~ 1, 
                                     TRUE ~ 0),
    # any metformin prescription in 6m prior to pandemic start, among those that started after a T2DM diagnosis, mid2018 onwards
    # should be less than exp_bin_metfin_first; difference => those that stopped again before pandemic start 
    exp_bin_metfin_last = case_when(exp_bin_metfin_first == 1 
                                    & exp_date_metfin_last >= study_dates$pandemicstart_date - days(183) ~ 1, 
                                    TRUE ~ 0),
    # any metformin mono-therapy prescription in 6m prior to pandemic start
    # should be less than exp_bin_metfin_last
    exp_bin_metfin_mono_last = case_when(exp_bin_metfin_first == 1 # ensures initiation of ANY metformin (broad codelist, any combo)
                                         & exp_date_metfin_mono_last >= study_dates$pandemicstart_date - days(183) ~ 1, # reduces them to metformin mono
                                              TRUE ~ 0),
    # if started any metformin before OR after pandemic start (i.e. any initiation from T2DM diagnosis until study end date)
    # CAVE: irrespective those who stopped just before pandemic start
    exp_bin_metfin_anytime = case_when(exp_date_metfin_first <= study_dates$pandemicstart_date
                                   | out_date_metfin_first > study_dates$pandemicstart_date ~ 1, 
                                   TRUE ~ 0),
    exp_date_metfin_anytime = case_when(exp_bin_metfin_anytime == 1 ~ pmin(exp_date_metfin_first, out_date_metfin_first, na.rm = TRUE), 
                                    TRUE ~ as.Date(NA)),
    tb_T2DMdiag_metfin_anytime = case_when(exp_bin_metfin_anytime == 1 ~ as.numeric(difftime(exp_date_metfin_anytime, elig_date_t2dm, units = "days")),
                                       TRUE ~ NA_real_),
    exp_bin_metfin_anytime_3m = case_when(!is.na(tb_T2DMdiag_metfin_anytime) & tb_T2DMdiag_metfin_anytime <= 90 ~ 1,
                                      TRUE ~ 0),
    exp_bin_metfin_anytime_6m = case_when(!is.na(tb_T2DMdiag_metfin_anytime) & tb_T2DMdiag_metfin_anytime <= 180 ~ 1,
                                      TRUE ~ 0),
    # if started any metformin mono-therapy before OR after pandemic start (i.e. any initiation from T2DM diagnosis until study end date)
    exp_bin_metfin_mono_anytime = case_when(exp_date_metfin_mono_first <= study_dates$pandemicstart_date
                                    | out_date_metfin_mono_first > study_dates$pandemicstart_date ~ 1, 
                                    TRUE ~ 0),
    exp_date_metfin_mono_anytime = case_when(exp_bin_metfin_mono_anytime == 1 ~ pmin(exp_date_metfin_mono_first, out_date_metfin_mono_first, na.rm = TRUE), 
                                    TRUE ~ as.Date(NA)),
    
    
    ## NOW, let's investigate those who did not start any metfin combo, OVER ENTIRE STUDY PERIOD, i.e. exp_bin_metfin_anytime == 0
    # DPP4 mono (or combo with SGLT2)
    exp_bin_dpp4_mono_anytime = case_when(exp_bin_metfin_anytime == 0 # codelist entails all combo with dpp4 (e.g. Janumet)
                                 & (exp_date_dpp4_first <= study_dates$pandemicstart_date # codelist entails all combo with metformin, does not matter, since they are all also part of the metfin combo list and thus are set to exp_bin_metfin_anytime == 1
                                    | out_date_dpp4_first > study_dates$pandemicstart_date) ~ 1, 
                                 TRUE ~ 0),
    exp_date_dpp4_mono_anytime = case_when(exp_bin_dpp4_mono_anytime == 1 ~ pmin(exp_date_dpp4_first, out_date_dpp4_first, na.rm = TRUE), 
                                    TRUE ~ as.Date(NA)),
    # TZD mono
    exp_bin_tzd_mono_anytime = case_when(exp_bin_metfin_anytime == 0 # codelist entails all combo with tzd (e.g. actoplusmet, or combo with Rosiglitazone, but formally not in use anymore in NHS after 2010: https://www.gov.uk/drug-device-alerts/drug-alert-recall-of-avandia-4mg-8mg-avandamet-1mg-500mg-2mg-500mg-2mg-1000mg-4mg-1000mg)
                                          & (exp_date_tzd_first <= study_dates$pandemicstart_date # codelist entails all combo with metformin, does not matter, since they are all also part of the metfin combo list and thus are set to exp_bin_metfin_anytime == 1
                                             | out_date_tzd_first > study_dates$pandemicstart_date) ~ 1, 
                                          TRUE ~ 0),
    exp_date_tzd_mono_anytime = case_when(exp_bin_tzd_mono_anytime == 1 ~ pmin(exp_date_tzd_first, out_date_tzd_first, na.rm = TRUE), 
                                           TRUE ~ as.Date(NA)),
    # SGLT2 mono (or combo with DPP4)
    exp_bin_sglt2_mono_anytime = case_when(exp_bin_metfin_anytime == 0 # codelist entails all combo with sglt2 (e.g. synjardy)
                                         & (exp_date_sglt2_first <= study_dates$pandemicstart_date # codelist entails all combo with metformin, does not matter, since they are all also part of the metfin combo list and thus are set to exp_bin_metfin_anytime == 1
                                            | out_date_sglt2_first > study_dates$pandemicstart_date) ~ 1, 
                                         TRUE ~ 0),
    exp_date_sglt2_mono_anytime = case_when(exp_bin_sglt2_mono_anytime == 1 ~ pmin(exp_date_sglt2_first, out_date_sglt2_first, na.rm = TRUE), 
                                          TRUE ~ as.Date(NA)),
    # sulfo mono
    exp_bin_sulfo_mono_anytime = case_when(exp_bin_metfin_anytime == 0 # codelist only entails sulfo mono (glucovance/glibenclamid + metfin not in use anymore)
                                           & (exp_date_sulfo_first <= study_dates$pandemicstart_date
                                              | out_date_sulfo_first > study_dates$pandemicstart_date) ~ 1, 
                                           TRUE ~ 0),
    exp_date_sulfo_mono_anytime = case_when(exp_bin_sulfo_mono_anytime == 1 ~ pmin(exp_date_sulfo_first, out_date_sulfo_first, na.rm = TRUE), 
                                            TRUE ~ as.Date(NA)),
    # glp1 mono
    exp_bin_glp1_mono_anytime = case_when(exp_bin_metfin_anytime == 0 # codelist only entails glp1 mono (no combinations with metformin)
                                           & (exp_date_glp1_first <= study_dates$pandemicstart_date
                                              | out_date_glp1_first > study_dates$pandemicstart_date) ~ 1, 
                                           TRUE ~ 0),
    exp_date_glp1_mono_anytime = case_when(exp_bin_glp1_mono_anytime == 1 ~ pmin(exp_date_glp1_first, out_date_glp1_first, na.rm = TRUE), 
                                            TRUE ~ as.Date(NA)),
    # megli mono
    exp_bin_megli_mono_anytime = case_when(exp_bin_metfin_anytime == 0 # codelist only entails megli mono (no combinations with metformin)
                                          & (exp_date_megli_first <= study_dates$pandemicstart_date
                                             | out_date_megli_first > study_dates$pandemicstart_date) ~ 1, 
                                          TRUE ~ 0),
    exp_date_megli_mono_anytime = case_when(exp_bin_megli_mono_anytime == 1 ~ pmin(exp_date_megli_first, out_date_megli_first, na.rm = TRUE), 
                                           TRUE ~ as.Date(NA)),
    # agi mono
    exp_bin_agi_mono_anytime = case_when(exp_bin_metfin_anytime == 0 # codelist only entails megli mono (no combinations with metformin)
                                           & (exp_date_agi_first <= study_dates$pandemicstart_date
                                              | out_date_agi_first > study_dates$pandemicstart_date) ~ 1, 
                                           TRUE ~ 0),
    exp_date_agi_mono_anytime = case_when(exp_bin_agi_mono_anytime == 1 ~ pmin(exp_date_agi_first, out_date_agi_first, na.rm = TRUE), 
                                            TRUE ~ as.Date(NA)),
    # insulin mono
    exp_bin_insulin_mono_anytime = case_when(exp_bin_metfin_anytime == 0 # codelist only entails megli mono (no combinations with metformin)
                                         & (exp_date_insulin_first <= study_dates$pandemicstart_date
                                            | out_date_insulin_first > study_dates$pandemicstart_date) ~ 1, 
                                         TRUE ~ 0),
    exp_date_insulin_mono_anytime = case_when(exp_bin_insulin_mono_anytime == 1 ~ pmin(exp_date_insulin_first, out_date_insulin_first, na.rm = TRUE), 
                                          TRUE ~ as.Date(NA)),
    
    ## NOW, let's investigate those who did not start any metfin combo UNTIL PANDEMIC START, i.e. exp_bin_metfin_first == 0
    # Of course, they might initiate metfin later, after pandemic start
    # DPP4 mono (or combo with SGLT2)
    exp_bin_dpp4_mono_first = case_when(exp_date_metfin_first == 0
                                        & exp_date_dpp4_first <= study_dates$pandemicstart_date ~ 1, # combo with SGLT2 possible
                                        TRUE ~ 0),
    # TZD mono
    exp_bin_tzd_mono_first = case_when(exp_date_metfin_first == 0
                                        & exp_date_tzd_first <= study_dates$pandemicstart_date ~ 1,
                                        TRUE ~ 0),
    # SGLT2 mono (or combo with DPP4)
    exp_bin_sglt2_mono_first = case_when(exp_date_metfin_first == 0
                                       & exp_date_sglt2_first <= study_dates$pandemicstart_date ~ 1, # but combo with DPP4 possible
                                       TRUE ~ 0),
    # sulfo mono
    exp_bin_sulfo_mono_first = case_when(exp_date_metfin_first == 0
                                         & exp_date_sulfo_first <= study_dates$pandemicstart_date ~ 1,
                                         TRUE ~ 0),
    # glp1 mono
    exp_bin_glp1_mono_first = case_when(exp_date_metfin_first == 0
                                         & exp_date_glp1_first <= study_dates$pandemicstart_date ~ 1,
                                         TRUE ~ 0),
    # megli mono
    exp_bin_megli_mono_first = case_when(exp_date_metfin_first == 0
                                        & exp_date_megli_first <= study_dates$pandemicstart_date ~ 1,
                                        TRUE ~ 0),
    # agi mono
    exp_bin_agi_mono_first = case_when(exp_date_metfin_first == 0
                                         & exp_date_agi_first <= study_dates$pandemicstart_date ~ 1,
                                         TRUE ~ 0),
    # insulin mono
    exp_bin_insulin_mono_first = case_when(exp_date_metfin_first == 0
                                       & exp_date_insulin_first <= study_dates$pandemicstart_date ~ 1,
                                       TRUE ~ 0),
    ## NOW, let's see who had no prescription at all UNTIL PANDEMIC START
    exp_bin_treat_nothing_first = case_when(exp_bin_metfin_first == 0 # covers exp_bin_metfin_mono_first (sub-codelist of exp_bin_metfin_first)
                                           & exp_bin_dpp4_mono_first == 0
                                           & exp_bin_tzd_mono_first == 0 
                                           & exp_bin_sglt2_mono_first == 0 
                                           & exp_bin_sulfo_mono_first == 0
                                           & exp_bin_glp1_mono_first == 0 
                                           & exp_bin_megli_mono_first == 0
                                           & exp_bin_agi_mono_first == 0
                                           & exp_bin_insulin_mono_first == 0 ~ 1,
                                           TRUE ~ 0),
    
    
    ## NOW, let's go backwards from pandemicstart_date, with a last prescription in past 6m, among those who do not have a any metfin combo in 6m prior to pandemic start, i.e. exp_bin_metfin_last == 0
    ## these are less relevant...since treatment status will be defined at landmark = 6m after T2DM diagnosis and not at pandemic start !
    # DPP4 mono (or combo with SGLT2)
    exp_bin_dpp4_mono_last = case_when(exp_bin_metfin_last == 0
                                      & cov_date_dpp4_last >= study_dates$pandemicstart_date - days(183) ~ 1, # but combo with SGLT2 possible
                                      TRUE ~ 0),
    # TZD mono
    exp_bin_tzd_mono_last = case_when(exp_bin_metfin_last == 0
                                     & cov_date_tzd_last >= study_dates$pandemicstart_date - days(183) ~ 1,
                                     TRUE ~ 0),
    # SGLT2 mono (or combo with DPP4)
    exp_bin_sglt2_mono_last = case_when(exp_bin_metfin_last == 0
                                       & cov_date_sglt2_last >= study_dates$pandemicstart_date - days(183) ~ 1, # but combo with DPP4 possible
                                       TRUE ~ 0),
    # sulfo mono
    exp_bin_sulfo_mono_last = case_when(exp_bin_metfin_last == 0
                                       & cov_date_sulfo_last >= study_dates$pandemicstart_date - days(183) ~ 1,
                                       TRUE ~ 0),
    # glp1 mono
    exp_bin_glp1_mono_last = case_when(exp_bin_metfin_last == 0
                                      & cov_date_glp1_last >= study_dates$pandemicstart_date - days(183) ~ 1,
                                      TRUE ~ 0),
    # megli mono
    exp_bin_megli_mono_last = case_when(exp_bin_metfin_last == 0
                                       & cov_date_megli_last >= study_dates$pandemicstart_date - days(183) ~ 1,
                                       TRUE ~ 0),
    # agi mono
    exp_bin_agi_mono_last = case_when(exp_bin_metfin_last == 0
                                     & cov_date_agi_last >= study_dates$pandemicstart_date - days(183) ~ 1,
                                     TRUE ~ 0),
    # insulin mono
    exp_bin_insulin_mono_last = case_when(exp_bin_metfin_last == 0
                                         & cov_date_insulin_last >= study_dates$pandemicstart_date - days(183) ~ 1,
                                         TRUE ~ 0),
    
    ## NOW, let's see who had no prescription at all from the above, at pandemicstart_date (any of these prescription within 6m prior to landmark)
    exp_bin_treat_nothing_last = case_when(exp_bin_metfin_last == 0 # covers exp_bin_metfin_mono_last (sub-codelist of exp_bin_metfin_last)
                                               & exp_bin_dpp4_mono_last == 0
                                               & exp_bin_tzd_mono_last == 0 
                                               & exp_bin_sglt2_mono_last == 0 
                                               & exp_bin_sulfo_mono_last == 0
                                               & exp_bin_glp1_mono_last == 0 
                                               & exp_bin_megli_mono_last == 0
                                               & exp_bin_agi_mono_last == 0
                                               & exp_bin_insulin_mono_last == 0 ~ 1,
                                               TRUE ~ 0)
    
    # e.g. invokamet (canaglifozin + metformin) or xigduo XR (dapaglifozin + metformin) or synjardy (Empagliflozin/metformin (0601023AR))
    # e.g. janumet (sitagliptin + metformin) or kombiglyze XR (saxagliptin + metformin) or jentadueto (linagliptin + metformin) or with vildagliptin or alogliptin
    # e.g. actoplusmet (pioglitazone + metformin), maybe a few Metformin hydrochloride/rosiglitazone (0601023V0)
    # e.g. glucovance (glibenclamide + metformin), but not commonly used anymore, due to increased risk of hypo (see: https://openprescribing.net/bnf/060102/) and therefore also not part of codelist
  
    ) %>%
  
  ## add primary outcome
  mutate(out_bin_severecovid = case_when(out_date_covid19_severe > study_dates$pandemicstart_date ~ 1, # severe covid outcome (hosp or death)
                                         TRUE ~ 0))

n_exp_out <- data_processed %>% 
  summarise(
    n_exp_bin_metfin_first = sum(exp_bin_metfin_first), 
    n_exp_bin_metfin_last = sum(exp_bin_metfin_last), 
    n_exp_bin_metfin_mono_last = sum(exp_bin_metfin_mono_last), 
    n_exp_bin_metfin_anytime = sum(exp_bin_metfin_anytime),
    n_exp_bin_metfin_anytime_3m = sum(exp_bin_metfin_anytime_3m), 
    n_exp_bin_metfin_anytime_6m = sum(exp_bin_metfin_anytime_6m),
    n_exp_bin_metfin_mono_anytime = sum(exp_bin_metfin_mono_anytime),
    n_exp_bin_dpp4_mono_anytime = sum(exp_bin_dpp4_mono_anytime),
    n_exp_bin_tzd_mono_anytime = sum(exp_bin_tzd_mono_anytime),
    n_exp_bin_sglt2_mono_anytime = sum(exp_bin_sglt2_mono_anytime),
    n_exp_bin_sulfo_mono_anytime = sum(exp_bin_sulfo_mono_anytime),
    n_exp_bin_glp1_mono_anytime = sum(exp_bin_glp1_mono_anytime),
    n_exp_bin_megli_mono_anytime = sum(exp_bin_megli_mono_anytime),
    n_exp_bin_agi_mono_anytime = sum(exp_bin_agi_mono_anytime),
    n_exp_bin_insulin_mono_anytime = sum(exp_bin_insulin_mono_anytime),
    
    n_exp_bin_metfin_mono_first = sum(exp_bin_metfin_mono_first),
    n_exp_bin_dpp4_mono_first = sum(exp_bin_dpp4_mono_first),
    n_exp_bin_tzd_mono_first = sum(exp_bin_tzd_mono_first),
    n_exp_bin_sglt2_mono_first = sum(exp_bin_sglt2_mono_first),
    n_exp_bin_sulfo_mono_first = sum(exp_bin_sulfo_mono_first),
    n_exp_bin_glp1_mono_first = sum(exp_bin_glp1_mono_first),
    n_exp_bin_megli_mono_first = sum(exp_bin_megli_mono_first),
    n_exp_bin_agi_mono_first = sum(exp_bin_agi_mono_first),
    n_exp_bin_insulin_mono_first = sum(exp_bin_insulin_mono_first),
    n_exp_bin_treat_nothing_first = sum(exp_bin_treat_nothing_first),
    
    n_exp_bin_dpp4_mono_last = sum(exp_bin_dpp4_mono_last),
    n_exp_bin_tzd_mono_last = sum(exp_bin_tzd_mono_last),
    n_exp_bin_sglt2_mono_last = sum(exp_bin_sglt2_mono_last),
    n_exp_bin_sulfo_mono_last = sum(exp_bin_sulfo_mono_last),
    n_exp_bin_glp1_mono_last = sum(exp_bin_glp1_mono_last),
    n_exp_bin_megli_mono_last = sum(exp_bin_megli_mono_last),
    n_exp_bin_agi_mono_last = sum(exp_bin_agi_mono_last),
    n_exp_bin_insulin_mono_last = sum(exp_bin_insulin_mono_last),
    n_exp_bin_treat_nothing_last = sum(exp_bin_treat_nothing_last),
    
    n_out_severeCOVID = sum(out_bin_severecovid),
    
    median_tb_T2DMdiag_metfin_anytime = median(tb_T2DMdiag_metfin_anytime, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_anytime = quantile(tb_T2DMdiag_metfin_anytime, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_anytime = quantile(tb_T2DMdiag_metfin_anytime, 0.75, na.rm = TRUE)
    )

# midpoint6 rounded
n_exp_out_midpoint6 <- data_processed %>% 
  summarise(
    n_exp_bin_metfin_last_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_last, na.rm = TRUE), threshold), 
    
    n_exp_bin_metfin_first_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_first, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_last_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_last, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_mono_last_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono_last, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_anytime_3m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_anytime_3m, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_anytime_6m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_anytime_6m, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_dpp4_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_sulfo_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_glp1_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_megli_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_megli_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_agi_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_agi_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_insulin_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin_mono_anytime, na.rm = TRUE), threshold), 
    
    n_exp_bin_metfin_mono_first_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono_first, na.rm = TRUE), threshold),
    n_exp_bin_dpp4_mono_first_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4_mono_first, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_mono_first_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd_mono_first, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_mono_first_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2_mono_first, na.rm = TRUE), threshold), 
    n_exp_bin_sulfo_mono_first_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo_mono_first, na.rm = TRUE), threshold), 
    n_exp_bin_glp1_mono_first_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1_mono_first, na.rm = TRUE), threshold), 
    n_exp_bin_megli_mono_first_midpoint6 = fn_roundmid_any(sum(exp_bin_megli_mono_first, na.rm = TRUE), threshold), 
    n_exp_bin_agi_mono_first_midpoint6 = fn_roundmid_any(sum(exp_bin_agi_mono_first, na.rm = TRUE), threshold), 
    n_exp_bin_insulin_mono_first_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin_mono_first, na.rm = TRUE), threshold), 
    n_exp_bin_treat_nothing_first_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_nothing_first, na.rm = TRUE), threshold), 
    
    n_exp_bin_dpp4_mono_last_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4_mono_last, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_mono_last_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd_mono_last, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_mono_last_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2_mono_last, na.rm = TRUE), threshold), 
    n_exp_bin_sulfo_mono_last_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo_mono_last, na.rm = TRUE), threshold), 
    n_exp_bin_glp1_mono_last_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1_mono_last, na.rm = TRUE), threshold), 
    n_exp_bin_megli_mono_last_midpoint6 = fn_roundmid_any(sum(exp_bin_megli_mono_last, na.rm = TRUE), threshold), 
    n_exp_bin_agi_mono_last_midpoint6 = fn_roundmid_any(sum(exp_bin_agi_mono_last, na.rm = TRUE), threshold), 
    n_exp_bin_insulin_mono_last_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin_mono_last, na.rm = TRUE), threshold), 
    n_exp_bin_treat_nothing_last_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_nothing_last, na.rm = TRUE), threshold), 
    
    n_out_severeCOVID_midpoint6 = fn_roundmid_any(sum(out_bin_severecovid, na.rm = TRUE), threshold),
    
    median_tb_T2DMdiag_metfin_anytime = median(tb_T2DMdiag_metfin_anytime, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_anytime = quantile(tb_T2DMdiag_metfin_anytime, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_anytime = quantile(tb_T2DMdiag_metfin_anytime, 0.75, na.rm = TRUE)
  )

# n_exp_severecovid_midpoint6 <- map(
#   .x = data_processed_all_windows,
#   .f = ~ .x %>% 
#     # main exposure
#     mutate(exp_bin_treat = case_when(exp_date_metfin_first <= study_dates$pandemicstart_date 
#                                      & exp_date_metfin_last >= study_dates$pandemicstart_date - days(183) ~ 1, # 1 if started/treated/exposed and still on it
#                                      TRUE ~ 0)) %>% # 0 if not started/treated/exposed until landmark or not on it anymore
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
# 9 Output for cumulative incidence plots re treatment regimen pattern
################################################################################
data_plots <- data_processed %>%
  dplyr::select(patient_id, elig_date_t2dm, exp_date_metfin_anytime, exp_bin_metfin_anytime, exp_date_metfin_mono_anytime, exp_bin_metfin_mono_anytime, exp_date_dpp4_mono_anytime, exp_bin_dpp4_mono_anytime,
         exp_date_tzd_mono_anytime, exp_bin_tzd_mono_anytime, exp_date_sglt2_mono_anytime, exp_bin_sglt2_mono_anytime, exp_date_sulfo_mono_anytime, exp_bin_sulfo_mono_anytime,
         exp_date_glp1_mono_anytime, exp_bin_glp1_mono_anytime, exp_date_megli_mono_anytime, exp_bin_megli_mono_anytime, exp_date_agi_mono_anytime, exp_bin_agi_mono_anytime,
         exp_date_insulin_mono_anytime, exp_bin_insulin_mono_anytime,
         out_date_covid19_severe, out_date_dereg_any)

# # ensure no event_time is == 0
# data_plots <- data_plots %>%
#   dplyr::filter(elig_date_t2dm < exp_date_metfin_anytime | is.na(exp_date_metfin_anytime)) %>% 
#   dplyr::filter(elig_date_t2dm < out_date_covid19_severe | is.na(out_date_covid19_severe)) %>% 
#   dplyr::filter(elig_date_t2dm < as.Date("2020-02-01") | is.na(elig_date_t2dm))
  
################################################################################
# 10 Save output
################################################################################
# the full data
write_rds(data_processed, here::here("output", "data", "data_processed.rds"))
# data for cumulative incidence plots re treatment regimen pattern
write_feather(data_plots, here::here("output", "data", "data_plots.feather"))

# flow chart quality assurance
write.csv(n_qa_excluded_midpoint6, file = here::here("output", "data_properties", "n_qa_excluded_midpoint6.csv"))
# flow chart completeness criteria
write.csv(n_completeness_excluded_midpoint6, file = here::here("output", "data_properties", "n_completeness_excluded_midpoint6.csv"))
# flow chart eligibility criteria
write.csv(n_elig_excluded_midpoint6, file = here::here("output", "data_properties", "n_elig_excluded_midpoint6.csv"))
write.csv(n_elig_excluded, file = here::here("output", "data_properties", "n_elig_excluded.csv"))
# descriptive/feasibility data re treatment patterns, events between index_date and pandemic start, and main outcome
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

# 
# # Load necessary libraries
# library(survival)  # For Kaplan-Meier analysis
# library(ggplot2)   # For plotting
# library(arrow)     # For handling Feather files
# 
# # File paths
# df_input <- "output/data/data_plots.feather"
# dir_output <- "output/km_estimates/km_estimates_metfin.feather"
# 
# # Variable definitions
# exposure <- "exp_bin_metfin_anytime"
# origin_date <- "elig_date_t2dm"
# event_date <- "exp_date_metfin_anytime"
# censor_date <- "out_date_covid19_severe"
# 
# # Read input data
# data <- data_plots
# 
# # Prepare survival object
# surv_obj <- Surv(
#   time = as.numeric(difftime(data[[event_date]], data[[origin_date]], units = "days")),
#   event = !is.na(data[[event_date]]) & 
#     (is.na(data[[censor_date]]) | data[[event_date]] <= data[[censor_date]])
# )
# 
# # Kaplan-Meier fit
# km_fit <- survfit(surv_obj ~ data[[exposure]], data = data)
# 
# # Save KM estimates to Feather file
# km_estimates <- data.frame(
#   time = km_fit$time,
#   n_risk = km_fit$n.risk,
#   n_event = km_fit$n.event,
#   n_censor = km_fit$n.censor,
#   survival = km_fit$surv,
#   std_err = km_fit$std.err,
#   conf_lower = km_fit$lower,
#   conf_upper = km_fit$upper
# )
# 
# plot(
#   km_fit,
#   col = c("blue", "red"),  # Colors for different groups
#   lty = 1:2,               # Line types for different groups
#   xlab = "Days",
#   ylab = "Survival Probability",
#   main = "Kaplan-Meier Survival Curves",
#   conf.int = F          # Add confidence intervals
# )
