####
## This script does the following:
# 1. Import/extract feather dataset from OpenSAFELY
# 2. Basic type formatting of variables -> fn_extract_data.R()
# 3. Process covariates
# 4. Detour depending if run locally or on real data
# 5. Evaluate/apply the quality assurance criteria -> fn_quality_assurance_midpoint6()
# 6. Evaluate/apply the completeness criteria: -> fn_completeness_criteria_midpoint6()
# 7. Evaluate/apply the eligibility criteria: -> fn_elig_criteria_midpoint6()
# 8. Add various treatment regimen patterns, assign treatment strategies, outcomes and ICEs
# 9. Separate output for cumulative incidence plots re treatment regimen pattern
# 10. Save all output
####

####
# Import libraries and user functions ----
####
library('arrow')
library('readr')
library('here')
library('lubridate')
library('dplyr')
library('tidyr')
library('purrr')
library('forcats')
library('ggplot2')
library('jsonlite')
library('skimr') # for skim
source(here::here("analysis", "functions", "fn_extract_data.R"))
source(here::here("analysis", "functions", "utility.R"))
source(here::here("analysis", "functions", "fn_quality_assurance_midpoint6.R"))
source(here::here("analysis", "functions", "fn_completeness_criteria_midpoint6.R"))
source(here::here("analysis", "functions", "fn_elig_criteria_midpoint6.R"))

####
# Create directories for output ----
####
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_description"))

####
# Import command-line arguments and dates ----
####
# args <- commandArgs(trailingOnly=TRUE) # if needed at a later stage to define local testing
# study_dates <- fromJSON(here::here("output", "study_dates.json")) # does not work locally, use below instead, try later
source(here::here("analysis", "metadates.R"))
# Convert the meta-dates into Date objects
study_dates <- lapply(study_dates, function(x) as.Date(x))

####
# Define redaction threshold ----
####
threshold <- 6

####
# Import the dataset definition ----
####
input_filename <- "dataset.arrow"

####
# Reformat the imported data ----
####
data_extracted <- fn_extract_data(input_filename)

####
# Process the data ----
####
data_processed <- data_extracted %>%
  mutate(
    # POPULATION/DEMOGRAPHIC
    cov_cat_age = cut(
      cov_num_age,
      breaks = c(18, 40, 60, 80, 110),
      labels = c("18-39", "40-59", "60-79", "80+"),
      right = FALSE),
    
    cov_cat_sex = fn_case_when(
      cov_cat_sex == "female" ~ "Female",
      cov_cat_sex == "male" ~ "Male",
      TRUE ~ "Unknown"),
    
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
      TRUE ~ "Unknown"),
    
    cov_cat_stp = as.factor(cov_cat_stp),
    
    # Finalize smoking status
    cov_cat_smoking_status = fn_case_when(
      cov_cat_smoking_status == "S" ~ "Smoker",
      cov_cat_smoking_status == "E" ~ "Ever",
      cov_cat_smoking_status == "N" ~ "Never",
      TRUE ~ "Unknown"),
    
    # Create one history of obesity variable
    cov_bin_obesity = fn_case_when(
      cov_bin_obesity == TRUE | cov_cat_bmi_groups == "Obese (>30)" ~ "Obese (>30)",
      cov_bin_obesity == FALSE & (cov_cat_bmi_groups == "Underweight" | cov_cat_bmi_groups == "Healthy weight (18.5-24.9)" | cov_cat_bmi_groups == "Overweight (25-29.9)") ~ "Not Obese (<=30)",
      TRUE ~ "Unknown"),
    
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
    
    # TC/HDL ratio categories: https://www.urmc.rochester.edu/encyclopedia/content?ContentTypeID=167&ContentID=lipid_panel_hdl_ratio#:~:text=Most%20healthcare%20providers%20want%20the,1%20is%20considered%20very%20good.
    cov_cat_tc_hdl_ratio = cut(
      cov_num_tc_hdl_ratio,
      breaks = c(1, 3.5, 5.1, 50), # 50 is upper limit, see above -> NA
      labels = c("below 3.5:1" ,"3.5:1 to 5:1", "above 5:1"),
      right = FALSE),
    cov_cat_tc_hdl_ratio = case_when(is.na(cov_num_tc_hdl_ratio) ~ factor("Unknown", 
                                                                          levels = c("below 3.5:1", "3.5:1 to 5:1", "above 5:1", "Unknown")), TRUE ~ cov_cat_tc_hdl_ratio),
    
    # HbA1c categories: https://www.southtees.nhs.uk/resources/the-hba1c-test/
    ## remove HbA1c > 120
    ## remove HbA1c below 0
    cov_num_hba1c_mmol_mol = replace(cov_num_hba1c_mmol_mol, cov_num_hba1c_mmol_mol < 0.00 | cov_num_hba1c_mmol_mol > 120.00, NA_real_),
    cov_cat_hba1c_mmol_mol = cut(
      cov_num_hba1c_mmol_mol,
      breaks = c(0, 42, 59, 76, 120), # 120 is upper limit, above NA
      labels = c("below 42" ,"42-58", "59-75", "above 75"),
      right = FALSE),
    cov_cat_hba1c_mmol_mol = case_when(is.na(cov_cat_hba1c_mmol_mol) ~ factor("Unknown", 
                                                                              levels = c("below 42", "42-58", "59-75", "above 75", "Unknown")), TRUE ~ cov_cat_hba1c_mmol_mol),
    )

####
# Modify dummy data ----
####
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  message("Running locally, adapt dummy data")
  source("analysis/modify_dummy_data.R")
  message("Dummy data successfully modified")
}

####
# Apply the quality assurance criteria ----
####
qa <- fn_quality_assurance_midpoint6(data_processed, study_dates, threshold)
n_qa_excluded_midpoint6 <- qa$n_qa_excluded_midpoint6
data_processed <- qa$data_processed

####
# Apply the completeness criteria ----
####
completeness <- fn_completeness_criteria_midpoint6(data_processed, threshold)
n_completeness_excluded <- completeness$n_completeness_excluded
n_completeness_excluded_midpoint6 <- completeness$n_completeness_excluded_midpoint6
data_processed <- completeness$data_processed
  
####
# Apply the eligibility criteria (Real Data)
####
# Our primary eligibility window to define incident T2DM is mid2018-mid2019, but maybe we may want to extend the window until max. mid2013 later on 
# => if so, use function with loop that can be mapped to other windows
eligibility <- fn_elig_criteria_midpoint6(data_processed, study_dates, years_in_days = 0, dummydata = FALSE)
n_elig_excluded <- eligibility$n_elig_excluded
n_elig_excluded_midpoint6 <- eligibility$n_elig_excluded_midpoint6
data_processed <- eligibility$data_processed

####
# Assign treatment patterns, outcomes and treatment strategies ----
####
data_processed <- data_processed %>% 
  mutate(
    # started any metformin within 6 months after T2DM, among those with a T2DM diagnosis, mid2018 onwards (those with exp_bin_metfin_first before mid2018 were already excluded above via fn_elig_criteria_midpoint6 [and for dummy data: this was modified accordingly in modify_dummy_data])
    # combo
    exp_bin_metfin = case_when(exp_date_metfin_first <= elig_date_t2dm + days(183) ~ TRUE, 
                                     TRUE ~ FALSE),
    # mono
    exp_bin_metfin_mono = case_when(exp_date_metfin_mono_first <= elig_date_t2dm + days(183) ~ TRUE,
                                     TRUE ~ FALSE),
    # any metformin prescription (combo and mono) in 6m prior to pandemic start, among those that started after a T2DM diagnosis, mid2018 onwards
    # should be less than exp_bin_metfin; difference => those that stopped again before pandemic start 
    # combo
    exp_bin_metfin_pandemicstart = case_when(exp_bin_metfin == TRUE 
                                    & exp_date_metfin_last >= study_dates$pandemicstart_date - days(183) ~ TRUE, 
                                    TRUE ~ FALSE),
    # mono
    exp_bin_metfin_mono_pandemicstart = case_when(exp_bin_metfin == TRUE # ensures initiation of ANY metformin (broad codelist, any combo)
                                             & exp_date_metfin_mono_last >= study_dates$pandemicstart_date - days(183) ~ TRUE, # reduces them to metformin mono
                                             TRUE ~ FALSE),
    # if started any metformin (from T2DM diagnosis until study end date)
    # CAVE: irrespective those who stopped just before pandemic start
    # combo
    exp_bin_metfin_anytime = case_when(!is.na(exp_date_metfin_first) ~ TRUE, 
                                   TRUE ~ FALSE),
    exp_date_metfin_anytime = case_when(exp_bin_metfin_anytime == TRUE ~ exp_date_metfin_first, 
                                    TRUE ~ as.Date(NA)),
    tb_T2DMdiag_metfin_anytime = case_when(exp_bin_metfin_anytime == TRUE ~ as.numeric(difftime(exp_date_metfin_anytime, elig_date_t2dm, units = "days")),
                                       TRUE ~ NA_real_),
    exp_bin_metfin_anytime_3m = case_when(!is.na(tb_T2DMdiag_metfin_anytime) & tb_T2DMdiag_metfin_anytime <= 90 ~ TRUE,
                                      TRUE ~ FALSE),
    exp_bin_metfin_anytime_6m = case_when(!is.na(tb_T2DMdiag_metfin_anytime) & tb_T2DMdiag_metfin_anytime <= 183 ~ TRUE,
                                      TRUE ~ FALSE),
    # mono
    exp_bin_metfin_mono_anytime = case_when(!is.na(exp_date_metfin_mono_first) ~ TRUE, 
                                       TRUE ~ FALSE),
    exp_date_metfin_mono_anytime = case_when(exp_bin_metfin_mono_anytime == TRUE ~ exp_date_metfin_mono_first, 
                                        TRUE ~ as.Date(NA)),
    tb_T2DMdiag_metfin_mono_anytime = case_when(exp_bin_metfin_mono_anytime == TRUE ~ as.numeric(difftime(exp_date_metfin_mono_anytime, elig_date_t2dm, units = "days")),
                                           TRUE ~ NA_real_),
    exp_bin_metfin_mono_anytime_3m = case_when(!is.na(tb_T2DMdiag_metfin_mono_anytime) & tb_T2DMdiag_metfin_mono_anytime <= 90 ~ TRUE,
                                          TRUE ~ FALSE),
    exp_bin_metfin_mono_anytime_6m = case_when(!is.na(tb_T2DMdiag_metfin_mono_anytime) & tb_T2DMdiag_metfin_mono_anytime <= 183 ~ TRUE,
                                          TRUE ~ FALSE),
    
    ## Let's investigate those who did not start any metfin combo, OVER ENTIRE STUDY PERIOD, i.e. exp_bin_metfin_anytime == FALSE
    # DPP4 mono (or combo with SGLT2)
    exp_bin_dpp4_mono_anytime = case_when(exp_bin_metfin_anytime == FALSE # codelist entails all combo with dpp4 (e.g. Janumet)
                                          & !is.na(exp_date_dpp4_first) # codelist entails all combo with metformin, does not matter, since they are all also part of the metfin combo list and thus are set to exp_bin_metfin_anytime == TRUE
                                          & exp_date_dpp4_first >= elig_date_t2dm ~ TRUE, # need to add this line since this was not an exclusion criteria
                                          TRUE ~ FALSE),
    exp_date_dpp4_mono_anytime = case_when(exp_bin_dpp4_mono_anytime == TRUE ~ exp_date_dpp4_first, 
                                    TRUE ~ as.Date(NA)),
    # TZD mono
    exp_bin_tzd_mono_anytime = case_when(exp_bin_metfin_anytime == FALSE
                                         & !is.na(exp_date_tzd_first) # codelist entails all combo with metformin, does not matter, since they are all also part of the metfin combo list and thus are set to exp_bin_metfin_anytime == TRUE (e.g. actoplusmet, or combo with Rosiglitazone, but formally not in use anymore in NHS after 2010: https://www.gov.uk/drug-device-alerts/drug-alert-recall-of-avandia-4mg-8mg-avandamet-1mg-500mg-2mg-500mg-2mg-1000mg-4mg-1000mg)
                                         & exp_date_tzd_first >= elig_date_t2dm ~ TRUE, 
                                         TRUE ~ FALSE),
    exp_date_tzd_mono_anytime = case_when(exp_bin_tzd_mono_anytime == TRUE ~ exp_date_tzd_first, 
                                           TRUE ~ as.Date(NA)),
    # SGLT2 mono (or combo with DPP4)
    exp_bin_sglt2_mono_anytime = case_when(exp_bin_metfin_anytime == FALSE # codelist entails all combo with sglt2 (e.g. synjardy)
                                           & !is.na(exp_date_sglt2_first)
                                           & exp_date_sglt2_first >= elig_date_t2dm ~ TRUE, 
                                           TRUE ~ FALSE),
    exp_date_sglt2_mono_anytime = case_when(exp_bin_sglt2_mono_anytime == TRUE ~ exp_date_sglt2_first, 
                                          TRUE ~ as.Date(NA)),
    # sulfo mono
    exp_bin_sulfo_mono_anytime = case_when(exp_bin_metfin_anytime == FALSE # codelist only entails sulfo mono (glucovance/glibenclamid + metfin not in use anymore)
                                           & !is.na(exp_date_sulfo_first)
                                           & exp_date_sulfo_first >= elig_date_t2dm ~ TRUE, 
                                           TRUE ~ FALSE),
    exp_date_sulfo_mono_anytime = case_when(exp_bin_sulfo_mono_anytime == TRUE ~ exp_date_sulfo_first, 
                                            TRUE ~ as.Date(NA)),
    # glp1 mono
    exp_bin_glp1_mono_anytime = case_when(exp_bin_metfin_anytime == FALSE # codelist only entails glp1 mono (no combinations with metformin)
                                          & !is.na(exp_date_glp1_first)
                                          & exp_date_glp1_first >= elig_date_t2dm ~ TRUE, 
                                          TRUE ~ FALSE),
    exp_date_glp1_mono_anytime = case_when(exp_bin_glp1_mono_anytime == TRUE ~ exp_date_glp1_first, 
                                            TRUE ~ as.Date(NA)),
    # megli mono
    exp_bin_megli_mono_anytime = case_when(exp_bin_metfin_anytime == FALSE # codelist only entails megli mono (no combinations with metformin)
                                           & !is.na(exp_date_megli_first)
                                           & exp_date_megli_first >= elig_date_t2dm ~ TRUE, 
                                           TRUE ~ FALSE),
    exp_date_megli_mono_anytime = case_when(exp_bin_megli_mono_anytime == TRUE ~ exp_date_megli_first, 
                                           TRUE ~ as.Date(NA)),
    # agi mono
    exp_bin_agi_mono_anytime = case_when(exp_bin_metfin_anytime == FALSE # codelist only entails megli mono (no combinations with metformin)
                                         & !is.na(exp_date_agi_first)
                                         & exp_date_agi_first >= elig_date_t2dm ~ TRUE, 
                                         TRUE ~ FALSE),
    exp_date_agi_mono_anytime = case_when(exp_bin_agi_mono_anytime == TRUE ~ exp_date_agi_first, 
                                            TRUE ~ as.Date(NA)),
    # insulin mono
    exp_bin_insulin_mono_anytime = case_when(exp_bin_metfin_anytime == FALSE # codelist only entails megli mono (no combinations with metformin)
                                             & !is.na(exp_date_insulin_first)
                                             & exp_date_insulin_first >= elig_date_t2dm ~ TRUE, 
                                             TRUE ~ FALSE),
    exp_date_insulin_mono_anytime = case_when(exp_bin_insulin_mono_anytime == TRUE ~ exp_date_insulin_first, 
                                          TRUE ~ as.Date(NA)),
    ## No prescription at all OVER ENTIRE STUDY PERIOD after T2DM diagnosis
    exp_bin_treat_nothing_anytime = case_when(exp_bin_metfin_anytime == FALSE
                                            & exp_bin_dpp4_mono_anytime == FALSE
                                            & exp_bin_tzd_mono_anytime == FALSE
                                            & exp_bin_sglt2_mono_anytime == FALSE 
                                            & exp_bin_sulfo_mono_anytime == FALSE
                                            & exp_bin_glp1_mono_anytime == FALSE 
                                            & exp_bin_megli_mono_anytime == FALSE
                                            & exp_bin_agi_mono_anytime == FALSE
                                            & exp_bin_insulin_mono_anytime == FALSE ~ TRUE,
                                            TRUE ~ FALSE),
    
    ## Let's investigate those who did not start any metfin COMBO UNTIL 6M LANDMARK, i.e. exp_bin_metfin == FALSE
    # Of course, they might initiate metfin later
    # DPP4 mono (or combo with SGLT2)
    exp_bin_dpp4_mono = case_when(exp_bin_metfin == FALSE # if we use the metfin_mono then we allow to count people who initiated DPP4 + metformin
                                  & exp_date_dpp4_first <= elig_date_t2dm + days(183)
                                  & exp_date_dpp4_first >= elig_date_t2dm ~ TRUE, # need to add this line since this was not an exclusion criteria
                                  TRUE ~ FALSE),
    # TZD mono
    exp_bin_tzd_mono = case_when(exp_bin_metfin == FALSE
                                 & exp_date_tzd_first <= elig_date_t2dm + days(183)
                                 & exp_date_tzd_first >= elig_date_t2dm ~ TRUE,
                                 TRUE ~ FALSE),
    # SGLT2 mono (or combo with DPP4)
    exp_bin_sglt2_mono = case_when(exp_bin_metfin == FALSE
                                   & exp_date_sglt2_first <= elig_date_t2dm + days(183)
                                   & exp_date_sglt2_first >= elig_date_t2dm ~ TRUE,
                                   TRUE ~ FALSE),
    # sulfo mono
    exp_bin_sulfo_mono = case_when(exp_bin_metfin == FALSE
                                   & exp_date_sulfo_first <= elig_date_t2dm + days(183)
                                   & exp_date_sulfo_first >= elig_date_t2dm ~ TRUE,
                                   TRUE ~ FALSE),
    # glp1 mono
    exp_bin_glp1_mono = case_when(exp_bin_metfin == FALSE
                                  & exp_date_glp1_first <= elig_date_t2dm + days(183)
                                  & exp_date_glp1_first >= elig_date_t2dm ~ TRUE,
                                  TRUE ~ FALSE),
    # megli mono
    exp_bin_megli_mono = case_when(exp_bin_metfin == FALSE
                                   & exp_date_megli_first <= elig_date_t2dm + days(183)
                                   & exp_date_megli_first >= elig_date_t2dm ~ TRUE,
                                   TRUE ~ FALSE),
    # agi mono
    exp_bin_agi_mono = case_when(exp_bin_metfin == FALSE
                                 & exp_date_agi_first <= elig_date_t2dm + days(183)
                                 & exp_date_agi_first >= elig_date_t2dm ~ TRUE,
                                 TRUE ~ FALSE),
    # insulin mono
    exp_bin_insulin_mono = case_when(exp_bin_metfin == FALSE
                                     & exp_date_insulin_first <= elig_date_t2dm + days(183)
                                     & exp_date_insulin_first >= elig_date_t2dm ~ TRUE,
                                     TRUE ~ FALSE),
    ## Who has no prescription at all UNTIL 6M after T2DM diagnosis
    exp_bin_treat_nothing = case_when(exp_bin_metfin == FALSE # covers exp_bin_metfin_mono (sub-codelist of exp_bin_metfin)
                                           & exp_bin_dpp4_mono == FALSE
                                           & exp_bin_tzd_mono == FALSE 
                                           & exp_bin_sglt2_mono == FALSE 
                                           & exp_bin_sulfo_mono == FALSE
                                           & exp_bin_glp1_mono == FALSE 
                                           & exp_bin_megli_mono == FALSE
                                           & exp_bin_agi_mono == FALSE
                                           & exp_bin_insulin_mono == FALSE ~ TRUE,
                                           TRUE ~ FALSE),
    ## OAD prescription (except metformin combo) UNTIL 6M after T2DM diagnosis (i.e. will not count metfin combo, only oad mono stuff!)
    exp_bin_oad = case_when(exp_bin_metfin == FALSE # covers metfin_mono
                                              & (exp_bin_dpp4_mono == TRUE
                                                 | exp_bin_tzd_mono == TRUE 
                                                 | exp_bin_sglt2_mono == TRUE 
                                                 | exp_bin_sulfo_mono == TRUE
                                                 | exp_bin_glp1_mono == TRUE 
                                                 | exp_bin_megli_mono == TRUE
                                                 | exp_bin_agi_mono == TRUE
                                                 | exp_bin_insulin_mono == TRUE) ~ TRUE,
                                              TRUE ~ FALSE),
    ## OAD prescription (except metformin combo) UNTIL 6M after T2DM diagnosis (i.e. will not have metfin combo in control arm) OR nothing
    # should be more than exp_bin_treat_nothing
    exp_bin_oad_nothing = case_when(exp_bin_metfin == FALSE
                                    | exp_bin_dpp4_mono == TRUE
                                    | exp_bin_tzd_mono == TRUE
                                    | exp_bin_sglt2_mono == TRUE
                                    | exp_bin_sulfo_mono == TRUE
                                    | exp_bin_glp1_mono == TRUE
                                    | exp_bin_megli_mono == TRUE
                                    | exp_bin_agi_mono == TRUE
                                    | exp_bin_insulin_mono == TRUE ~ TRUE,
                                    TRUE ~ FALSE),
    
    ## Let's investigate those who did not start any metfin MONO UNTIL 6M LANDMARK, i.e. exp_bin_metfin_mono == FALSE
    # Of course, they might initiate metfin later
    # DPP4 (or combo with SGLT2) +/- metformin
    exp_bin_dpp4 = case_when(exp_bin_metfin_mono == FALSE # Now we may have people who initiated DPP4 + metformin in comparison arm
                                  & exp_date_dpp4_first <= elig_date_t2dm + days(183)
                                  & exp_date_dpp4_first >= elig_date_t2dm ~ TRUE, # need to add this line since this was not an exclusion criteria
                                  TRUE ~ FALSE),
    # TZD +/- metformin
    exp_bin_tzd = case_when(exp_bin_metfin_mono == FALSE
                                 & exp_date_tzd_first <= elig_date_t2dm + days(183)
                                 & exp_date_tzd_first >= elig_date_t2dm ~ TRUE,
                                 TRUE ~ FALSE),
    # SGLT2 (or combo with DPP4) +/- metformin
    exp_bin_sglt2 = case_when(exp_bin_metfin_mono == FALSE
                                   & exp_date_sglt2_first <= elig_date_t2dm + days(183)
                                   & exp_date_sglt2_first >= elig_date_t2dm ~ TRUE,
                                   TRUE ~ FALSE),
    # sulfo +/- metformin
    exp_bin_sulfo = case_when(exp_bin_metfin_mono == FALSE
                                   & exp_date_sulfo_first <= elig_date_t2dm + days(183)
                                   & exp_date_sulfo_first >= elig_date_t2dm ~ TRUE,
                                   TRUE ~ FALSE),
    # glp1 +/- metformin
    exp_bin_glp1 = case_when(exp_bin_metfin_mono == FALSE
                                  & exp_date_glp1_first <= elig_date_t2dm + days(183)
                                  & exp_date_glp1_first >= elig_date_t2dm ~ TRUE,
                                  TRUE ~ FALSE),
    # megli +/- metformin
    exp_bin_megli = case_when(exp_bin_metfin_mono == FALSE
                                   & exp_date_megli_first <= elig_date_t2dm + days(183)
                                   & exp_date_megli_first >= elig_date_t2dm ~ TRUE,
                                   TRUE ~ FALSE),
    # agi +/- metformin
    exp_bin_agi = case_when(exp_bin_metfin_mono == FALSE
                                 & exp_date_agi_first <= elig_date_t2dm + days(183)
                                 & exp_date_agi_first >= elig_date_t2dm ~ TRUE,
                                 TRUE ~ FALSE),
    # insulin +/- metformin
    exp_bin_insulin = case_when(exp_bin_metfin_mono == FALSE
                                     & exp_date_insulin_first <= elig_date_t2dm + days(183)
                                     & exp_date_insulin_first >= elig_date_t2dm ~ TRUE,
                                     TRUE ~ FALSE),
    ## No prescription at all UNTIL 6M after T2DM diagnosis | Should give the same result as exp_bin_treat_nothing, just a different distribution across the arms
    exp_bin_treat_nothing2 = case_when(exp_bin_metfin_mono == FALSE
                                      & exp_bin_dpp4 == FALSE
                                      & exp_bin_tzd == FALSE 
                                      & exp_bin_sglt2 == FALSE 
                                      & exp_bin_sulfo == FALSE
                                      & exp_bin_glp1 == FALSE 
                                      & exp_bin_megli == FALSE
                                      & exp_bin_agi == FALSE
                                      & exp_bin_insulin == FALSE ~ TRUE,
                                      TRUE ~ FALSE),
    ## OAD prescription (except metformin mono) UNTIL 6M after T2DM diagnosis (i.e. might have some metfin combo in control arm)
    # should be more than exp_bin_oad
    exp_bin_oad_metfincombo = case_when(exp_bin_metfin_mono == FALSE
                                       & (exp_bin_dpp4 == TRUE
                                       | exp_bin_tzd == TRUE 
                                       | exp_bin_sglt2 == TRUE 
                                       | exp_bin_sulfo == TRUE
                                       | exp_bin_glp1 == TRUE 
                                       | exp_bin_megli == TRUE
                                       | exp_bin_agi == TRUE
                                       | exp_bin_insulin == TRUE) ~ TRUE,
                                       TRUE ~ FALSE),
    ## OAD prescription (except metformin mono) UNTIL 6M after T2DM diagnosis (i.e. might have some metfin combo in control arm) OR nothing
    # should be more than exp_bin_oad_nothing
    exp_bin_oad_metfincombo_nothing = case_when(exp_bin_metfin_mono == FALSE
                                                | exp_bin_dpp4 == TRUE
                                                | exp_bin_tzd == TRUE 
                                                | exp_bin_sglt2 == TRUE 
                                                | exp_bin_sulfo == TRUE
                                                | exp_bin_glp1 == TRUE 
                                                | exp_bin_megli == TRUE
                                                | exp_bin_agi == TRUE
                                                | exp_bin_insulin == TRUE ~ TRUE,
                                        TRUE ~ FALSE)
    ) %>%
  
  ## add outcomes and ICEs
  mutate(out_bin_severecovid = case_when(out_date_covid19_severe > elig_date_t2dm ~ TRUE, # severe covid outcome (hosp or death)
                                         TRUE ~ FALSE),
         out_date_severecovid = case_when(out_bin_severecovid == TRUE ~ out_date_covid19_severe, 
                                             TRUE ~ as.Date(NA)),
         out_bin_covid_hosp = case_when(out_date_covid19_hosp > elig_date_t2dm ~ TRUE,
                                          TRUE ~ FALSE),
         out_date_covid_hosp = case_when(out_bin_covid_hosp == TRUE ~ out_date_covid19_hosp, 
                                           TRUE ~ as.Date(NA)),
         out_bin_covid_death = case_when(out_date_covid19_death > elig_date_t2dm ~ TRUE,
                                        TRUE ~ FALSE),
         out_date_covid_death = case_when(out_bin_covid_death == TRUE ~ out_date_covid19_death, 
                                         TRUE ~ as.Date(NA)),
         out_bin_covid = case_when(out_date_covid19 > elig_date_t2dm ~ TRUE,
                                         TRUE ~ FALSE),
         out_date_covid = case_when(out_bin_covid == TRUE ~ out_date_covid19, 
                                          TRUE ~ as.Date(NA)),
         out_bin_longcovid = case_when(out_date_long_covid19 > elig_date_t2dm ~ TRUE,
                                   TRUE ~ FALSE),
         out_date_longcovid = case_when(out_bin_longcovid == TRUE ~ out_date_long_covid19, 
                                    TRUE ~ as.Date(NA)),
         out_bin_virfat = case_when(out_date_viral_fatigue > elig_date_t2dm ~ TRUE,
                                        TRUE ~ FALSE),
         out_date_virfat = case_when(out_bin_virfat == TRUE ~ out_date_viral_fatigue, 
                                         TRUE ~ as.Date(NA)),
         out_bin_longcovid_virfat = case_when(out_date_long_fatigue > elig_date_t2dm ~ TRUE,
                                     TRUE ~ FALSE),
         out_date_longcovid_virfat = case_when(out_bin_longcovid_virfat == TRUE ~ out_date_long_fatigue, 
                                      TRUE ~ as.Date(NA)),
         # deaths between landmark and pandemic start
         out_bin_death_pandemicstart = case_when(!is.na(qa_date_of_death)
                                                 & qa_date_of_death <= study_dates$pandemicstart_date 
                                                 & qa_date_of_death > elig_date_t2dm + days(183) ~ TRUE,
                                                 TRUE ~ FALSE),
         out_date_death_pandemicstart = case_when(out_bin_death_pandemicstart == TRUE ~ qa_date_of_death, 
                                           TRUE ~ as.Date(NA)),
         # LTFU between landmark and pandemic start
         out_bin_ltfu_pandemicstart = case_when(!is.na(cens_date_dereg)
                                                & cens_date_dereg <= study_dates$pandemicstart_date
                                                & cens_date_dereg > elig_date_t2dm + days(183) ~ TRUE,
                                                TRUE ~ FALSE),
         out_date_ltfu_pandemicstart = case_when(out_bin_ltfu_pandemicstart == TRUE ~ cens_date_dereg, 
                                                  TRUE ~ as.Date(NA))
     ) %>%
  
  ## add treatment strategy variables
  mutate(exp_bin_treat = case_when(exp_bin_metfin_mono == TRUE ~ 1,
                                   exp_bin_treat_nothing == TRUE ~ 2,
                                   TRUE ~ NA_real_),
         exp_bin_treat_all = case_when(exp_bin_metfin_mono == TRUE ~ 1, 
                                   exp_bin_treat_nothing == TRUE ~ 2,
                                   exp_bin_oad == TRUE ~ 3,
                                   exp_bin_oad_metfincombo == TRUE ~ 4,
                                   TRUE ~ NA_real_),
         exp_bin_treat_3groups = case_when(exp_bin_metfin_mono == TRUE ~ 1, 
                                       exp_bin_treat_nothing == TRUE ~ 2,
                                       exp_bin_oad == TRUE | exp_bin_metfin == TRUE ~ 3,
                                       TRUE ~ NA_real_)
  )

n_exp_out <- data_processed %>% 
  summarise(
    n_exp_bin_metfin = sum(exp_bin_metfin), 
    n_exp_bin_metfin_mono = sum(exp_bin_metfin_mono), 
    n_exp_bin_metfin_pandemicstart = sum(exp_bin_metfin_pandemicstart), 
    n_exp_bin_metfin_mono_pandemicstart = sum(exp_bin_metfin_mono_pandemicstart), 
    n_exp_bin_metfin_anytime = sum(exp_bin_metfin_anytime),
    n_exp_bin_metfin_anytime_3m = sum(exp_bin_metfin_anytime_3m), 
    n_exp_bin_metfin_anytime_6m = sum(exp_bin_metfin_anytime_6m),
    n_exp_bin_metfin_mono_anytime = sum(exp_bin_metfin_mono_anytime),
    n_exp_bin_metfin_mono_anytime_3m = sum(exp_bin_metfin_mono_anytime_3m), 
    n_exp_bin_metfin_mono_anytime_6m = sum(exp_bin_metfin_mono_anytime_6m),
    
    n_exp_bin_dpp4_mono_anytime = sum(exp_bin_dpp4_mono_anytime),
    n_exp_bin_tzd_mono_anytime = sum(exp_bin_tzd_mono_anytime),
    n_exp_bin_sglt2_mono_anytime = sum(exp_bin_sglt2_mono_anytime),
    n_exp_bin_sulfo_mono_anytime = sum(exp_bin_sulfo_mono_anytime),
    n_exp_bin_glp1_mono_anytime = sum(exp_bin_glp1_mono_anytime),
    n_exp_bin_megli_mono_anytime = sum(exp_bin_megli_mono_anytime),
    n_exp_bin_agi_mono_anytime = sum(exp_bin_agi_mono_anytime),
    n_exp_bin_insulin_mono_anytime = sum(exp_bin_insulin_mono_anytime),
    n_exp_bin_treat_nothing_anytime = sum(exp_bin_treat_nothing_anytime),
    
    n_exp_bin_dpp4_mono = sum(exp_bin_dpp4_mono),
    n_exp_bin_tzd_mono = sum(exp_bin_tzd_mono),
    n_exp_bin_sglt2_mono = sum(exp_bin_sglt2_mono),
    n_exp_bin_sulfo_mono = sum(exp_bin_sulfo_mono),
    n_exp_bin_glp1_mono = sum(exp_bin_glp1_mono),
    n_exp_bin_megli_mono = sum(exp_bin_megli_mono),
    n_exp_bin_agi_mono = sum(exp_bin_agi_mono),
    n_exp_bin_insulin_mono = sum(exp_bin_insulin_mono),
    n_exp_bin_oad = sum(exp_bin_oad),
    n_exp_bin_oad_nothing = sum(exp_bin_oad_nothing),
    n_exp_bin_treat_nothing = sum(exp_bin_treat_nothing),
    
    n_exp_bin_dpp4 = sum(exp_bin_dpp4),
    n_exp_bin_tzd = sum(exp_bin_tzd),
    n_exp_bin_sglt2 = sum(exp_bin_sglt2),
    n_exp_bin_sulfo = sum(exp_bin_sulfo),
    n_exp_bin_glp1 = sum(exp_bin_glp1),
    n_exp_bin_megli = sum(exp_bin_megli),
    n_exp_bin_agi = sum(exp_bin_agi),
    n_exp_bin_insulin = sum(exp_bin_insulin),
    n_exp_bin_oad_metfincombo = sum(exp_bin_oad_metfincombo),
    n_exp_bin_oad_metfincombo_nothing = sum(exp_bin_oad_metfincombo_nothing),
    n_exp_bin_treat_nothing2 = sum(exp_bin_treat_nothing2),
    
    n_out_bin_severeCOVID = sum(out_bin_severecovid),
    n_out_bin_covid_hosp = sum(out_bin_covid_hosp),
    n_out_bin_covid_death = sum(out_bin_covid_death),
    n_out_bin_covid = sum(out_bin_covid),
    n_out_bin_longcovid = sum(out_bin_longcovid),
    n_out_bin_virfat = sum(out_bin_virfat),
    n_out_bin_longcovid_virfat = sum(out_bin_longcovid_virfat),
    n_out_bin_death_pandemicstart = sum(out_bin_death_pandemicstart),
    n_out_bin_ltfu_pandemicstart = sum(out_bin_ltfu_pandemicstart),
    
    median_tb_T2DMdiag_metfin_anytime = median(tb_T2DMdiag_metfin_anytime, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_anytime = quantile(tb_T2DMdiag_metfin_anytime, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_anytime = quantile(tb_T2DMdiag_metfin_anytime, 0.75, na.rm = TRUE),
    
    median_tb_T2DMdiag_metfin_mono_anytime = median(tb_T2DMdiag_metfin_mono_anytime, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_mono_anytime = quantile(tb_T2DMdiag_metfin_mono_anytime, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_mono_anytime = quantile(tb_T2DMdiag_metfin_mono_anytime, 0.75, na.rm = TRUE)
    
    ) %>% 
  
  # pivot (for easier data review in L4)
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Value"
  )


# midpoint6 rounded
n_exp_out_midpoint6 <- data_processed %>% 
  summarise(
    n_exp_bin_metfin_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_pandemicstart_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_pandemicstart, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_mono_pandemicstart_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono_pandemicstart, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_anytime_3m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_anytime_3m, na.rm = TRUE), threshold), 
    # n_exp_bin_metfin_anytime_6m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_anytime_6m, na.rm = TRUE), threshold), # matched with n_exp_bin_metfin_midpoint6: Correct
    n_exp_bin_metfin_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_mono_anytime_3m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono_anytime_3m, na.rm = TRUE), threshold), 
    # n_exp_bin_metfin_mono_anytime_6m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono_anytime_6m, na.rm = TRUE), threshold), # matched with n_exp_bin_metfin_mono_midpoint6: Correct
    
    n_exp_bin_dpp4_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_sulfo_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_glp1_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_megli_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_megli_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_agi_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_agi_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_insulin_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_treat_nothing_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_nothing_anytime, na.rm = TRUE), threshold), 
    
    n_exp_bin_dpp4_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4_mono, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd_mono, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2_mono, na.rm = TRUE), threshold), 
    n_exp_bin_sulfo_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo_mono, na.rm = TRUE), threshold), 
    n_exp_bin_glp1_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1_mono, na.rm = TRUE), threshold), 
    n_exp_bin_megli_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_megli_mono, na.rm = TRUE), threshold), 
    n_exp_bin_agi_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_agi_mono, na.rm = TRUE), threshold), 
    n_exp_bin_insulin_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin_mono, na.rm = TRUE), threshold), 
    n_exp_bin_oad_midpoint6 = fn_roundmid_any(sum(exp_bin_oad, na.rm = TRUE), threshold),     
    n_exp_bin_treat_nothing_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_nothing, na.rm = TRUE), threshold), 
    
    n_exp_bin_dpp4_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2, na.rm = TRUE), threshold), 
    n_exp_bin_sulfo_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo, na.rm = TRUE), threshold), 
    n_exp_bin_glp1_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1, na.rm = TRUE), threshold), 
    n_exp_bin_megli_midpoint6 = fn_roundmid_any(sum(exp_bin_megli, na.rm = TRUE), threshold), 
    n_exp_bin_agi_midpoint6 = fn_roundmid_any(sum(exp_bin_agi, na.rm = TRUE), threshold), 
    n_exp_bin_insulin_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin, na.rm = TRUE), threshold), 
    n_exp_bin_oad_metfincombo_midpoint6 = fn_roundmid_any(sum(exp_bin_oad_metfincombo, na.rm = TRUE), threshold),
    n_exp_bin_oad_metfincombo_nothing_midpoint6 = fn_roundmid_any(sum(exp_bin_oad_metfincombo_nothing, na.rm = TRUE), threshold),
    n_exp_bin_treat_nothing2_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_nothing2, na.rm = TRUE), threshold), 
    
    n_out_bin_severeCOVID_midpoint6 = fn_roundmid_any(sum(out_bin_severecovid, na.rm = TRUE), threshold), 
    n_out_bin_covid_hosp_midpoint6 = fn_roundmid_any(sum(out_bin_covid_hosp, na.rm = TRUE), threshold), 
    n_out_bin_covid_death_midpoint6 = fn_roundmid_any(sum(out_bin_covid_death, na.rm = TRUE), threshold), 
    n_out_bin_covid_midpoint6 = fn_roundmid_any(sum(out_bin_covid, na.rm = TRUE), threshold), 
    n_out_bin_longcovid_midpoint6 = fn_roundmid_any(sum(out_bin_longcovid, na.rm = TRUE), threshold), 
    n_out_bin_virfat_midpoint6 = fn_roundmid_any(sum(out_bin_virfat, na.rm = TRUE), threshold), 
    n_out_bin_longcovid_virfat_midpoint6 = fn_roundmid_any(sum(out_bin_longcovid_virfat, na.rm = TRUE), threshold),
    n_out_bin_death_pandemicstart_midpoint6 = fn_roundmid_any(sum(out_bin_death_pandemicstart, na.rm = TRUE), threshold), 
    n_out_bin_ltfu_pandemicstart_midpoint6 = fn_roundmid_any(sum(out_bin_ltfu_pandemicstart, na.rm = TRUE), threshold), 
    
    median_tb_T2DMdiag_metfin_anytime = median(tb_T2DMdiag_metfin_anytime, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_anytime = quantile(tb_T2DMdiag_metfin_anytime, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_anytime = quantile(tb_T2DMdiag_metfin_anytime, 0.75, na.rm = TRUE),
    
    median_tb_T2DMdiag_metfin_mono_anytime = median(tb_T2DMdiag_metfin_mono_anytime, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_mono_anytime = quantile(tb_T2DMdiag_metfin_mono_anytime, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_mono_anytime = quantile(tb_T2DMdiag_metfin_mono_anytime, 0.75, na.rm = TRUE)
    
  ) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Value"
    )

# Define labels
labels <- c(
  n_exp_bin_metfin_midpoint6 = "Metformin (combo) within 6m",
  n_exp_bin_metfin_mono_midpoint6 = "Metformin mono within 6m",
  n_exp_bin_metfin_pandemicstart_midpoint6 = "Metformin (combo) in 6m prior to pandemic start",
  n_exp_bin_metfin_mono_pandemicstart_midpoint6 = "Metformin mono in 6m prior to pandemic start",
  n_exp_bin_metfin_anytime_midpoint6 = "Metformin (combo) anytime",
  n_exp_bin_metfin_anytime_3m_midpoint6 = "Metformin (combo) within 3m",
  # n_exp_bin_metfin_anytime_6m_midpoint6 = "Metformin (combo) within 6m", # matched with n_exp_bin_metfin_midpoint6: Correct
  n_exp_bin_metfin_mono_anytime_midpoint6 = "Metformin mono anytime",
  n_exp_bin_metfin_mono_anytime_3m_midpoint6 = "Metformin mono within 3m",
  # n_exp_bin_metfin_mono_anytime_6m_midpoint6 = "Metformin mono within 6m", # matched with n_exp_bin_metfin_mono_midpoint6: Correct
  
  n_exp_bin_dpp4_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, DPP4",
  n_exp_bin_tzd_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, TZD",
  n_exp_bin_sglt2_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, SGLT2",
  n_exp_bin_sulfo_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, sulfonylurea",
  n_exp_bin_glp1_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, GLP1",
  n_exp_bin_megli_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, meglitinide",
  n_exp_bin_agi_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, alpha-glucosidase",
  n_exp_bin_insulin_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, insulin",
  n_exp_bin_treat_nothing_anytime_midpoint6 = "No metformin (combo) or any other antidiabetic anytime",
  
  n_exp_bin_dpp4_mono_midpoint6 = "Among those without metformin (combo) 6m, DPP4",
  n_exp_bin_tzd_mono_midpoint6 = "Among those without metformin (combo) 6m, TZD",
  n_exp_bin_sglt2_mono_midpoint6 = "Among those without metformin (combo) 6m, SGLT2",
  n_exp_bin_sulfo_mono_midpoint6 = "Among those without metformin (combo) 6m, sulfonylurea",
  n_exp_bin_glp1_mono_midpoint6 = "Among those without metformin (combo) 6m, GLP1",
  n_exp_bin_megli_mono_midpoint6 = "Among those without metformin (combo) 6m, meglitinide",
  n_exp_bin_agi_mono_midpoint6 = "Among those without metformin (combo) 6m, alpha-glucosidase",
  n_exp_bin_insulin_mono_midpoint6 = "Among those without metformin (combo) 6m, insulin",
  n_exp_bin_oad_midpoint6 = "Among those without metformin (combo) 6m, any antidiabetic (mono)",
  n_exp_bin_treat_nothing_midpoint6 = "No metformin (combo) or any other antidiabetic within 6m",
  
  n_exp_bin_dpp4_midpoint6 = "Among those without metformin mono 6m, DPP4 (+/- metformin)",
  n_exp_bin_tzd_midpoint6 = "Among those without metformin mono 6m, TZD (+/- metformin)",
  n_exp_bin_sglt2_midpoint6 = "Among those without metformin mono 6m, SGLT2 (+/- metformin)",
  n_exp_bin_sulfo_midpoint6 = "Among those without metformin mono 6m, sulfonylurea (+/- metformin)",
  n_exp_bin_glp1_midpoint6 = "Among those without metformin mono 6m, GLP1 (+/- metformin)",
  n_exp_bin_megli_midpoint6 = "Among those without metformin mono 6m, meglitinide (+/- metformin)",
  n_exp_bin_agi_midpoint6 = "Among those without metformin mono 6m, alpha-glucosidase (+/- metformin)",
  n_exp_bin_insulin_midpoint6 = "Among those without metformin mono 6m, insulin (+/- metformin)",
  n_exp_bin_oad_metfincombo_midpoint6 = "Among those without metformin mono 6m, any antidiabetic (+/- metformin)",
  n_exp_bin_oad_metfincombo_nothing_midpoint6 = "Among those without metformin mono 6m, any antidiabetic (+/- metformin) or nothing",
  n_exp_bin_treat_nothing2_midpoint6 = "No metformin mono or any other antidiabetic (+/- metformin) within 6m",

  n_out_bin_severeCOVID_midpoint6 = "COVID hosp or death (after baseline)",
  n_out_bin_covid_hosp_midpoint6 = "COVID hosp (after baseline)",
  n_out_bin_covid_death_midpoint6 = "COVID death (after baseline)",
  n_out_bin_covid_midpoint6 = "COVID diagnosis, pos test or hosp (after baseline)",
  n_out_bin_longcovid_midpoint6 = "Long COVID diagnosis (after baseline)",
  n_out_bin_virfat_midpoint6 = "Viral Fatigue (after baseline)",
  n_out_bin_longcovid_virfat_midpoint6 = "Long COVID or Viral Fatigue (after baseline)",
  n_out_bin_death_pandemicstart_midpoint6 = "Deaths between landmark and pandemic start",
  n_out_bin_ltfu_pandemicstart_midpoint6 = "LTFU between landmark and pandemic start",
  
  median_tb_T2DMdiag_metfin_anytime = "Median time from T2DM diagnosis to metformin (combo) start",
  IQR_lower_tb_T2DMdiag_metfin_anytime = "IQR lower bound: T2DM diagnosis to metformin (combo)",
  IQR_upper_tb_T2DMdiag_metfin_anytime = "IQR upper bound: T2DM diagnosis to metformin (combo)",
  
  median_tb_T2DMdiag_metfin_mono_anytime = "Median time from T2DM diagnosis to metformin mono start",
  IQR_lower_tb_T2DMdiag_metfin_mono_anytime = "IQR lower bound: T2DM diagnosis to metformin mono",
  IQR_upper_tb_T2DMdiag_metfin_mono_anytime = "IQR upper bound: T2DM diagnosis to metformin mono"
)

# Apply them
n_exp_out_midpoint6 <- n_exp_out_midpoint6 %>%
  mutate(Variable = labels[Variable])

####
# Output for cumulative incidence plots re treatment regimen pattern ----
####
# data_plots <- data_processed %>%
#   dplyr::select(patient_id, elig_date_t2dm, exp_date_metfin_anytime, exp_bin_metfin_anytime, exp_date_metfin_mono_anytime, exp_bin_metfin_mono_anytime,
#                 exp_date_dpp4_mono_anytime, exp_bin_dpp4_mono_anytime, exp_date_tzd_mono_anytime, exp_bin_tzd_mono_anytime,
#                 exp_date_sglt2_mono_anytime, exp_bin_sglt2_mono_anytime, exp_date_sulfo_mono_anytime, exp_bin_sulfo_mono_anytime,
#                 exp_date_glp1_mono_anytime, exp_bin_glp1_mono_anytime, exp_date_megli_mono_anytime, exp_bin_megli_mono_anytime,
#                 exp_date_agi_mono_anytime, exp_bin_agi_mono_anytime, exp_date_insulin_mono_anytime, exp_bin_insulin_mono_anytime,
#                 out_date_severecovid, qa_date_of_death)
# 
# data_plots <- data_plots %>% # double-check for plot to avoid event_time is == 0
#   dplyr::filter(elig_date_t2dm < exp_date_metfin_anytime | is.na(exp_date_metfin_anytime)) %>%
#   dplyr::filter(elig_date_t2dm < exp_date_metfin_mono_anytime | is.na(exp_date_metfin_mono_anytime)) %>%
#   dplyr::filter(elig_date_t2dm < exp_date_dpp4_mono_anytime | is.na(exp_date_dpp4_mono_anytime)) %>%
#   dplyr::filter(elig_date_t2dm < exp_date_tzd_mono_anytime | is.na(exp_date_tzd_mono_anytime)) %>%
#   dplyr::filter(elig_date_t2dm < exp_date_sglt2_mono_anytime | is.na(exp_date_sglt2_mono_anytime)) %>%
#   dplyr::filter(elig_date_t2dm < exp_date_sulfo_mono_anytime | is.na(exp_date_sulfo_mono_anytime)) %>%
#   dplyr::filter(elig_date_t2dm < exp_date_glp1_mono_anytime | is.na(exp_date_glp1_mono_anytime)) %>%
#   dplyr::filter(elig_date_t2dm < exp_date_megli_mono_anytime | is.na(exp_date_megli_mono_anytime)) %>%
#   dplyr::filter(elig_date_t2dm < exp_date_agi_mono_anytime | is.na(exp_date_agi_mono_anytime)) %>%
#   dplyr::filter(elig_date_t2dm < exp_date_insulin_mono_anytime | is.na(exp_date_insulin_mono_anytime)) %>%
#   dplyr::filter(elig_date_t2dm < qa_date_of_death | is.na(qa_date_of_death))

####
# Structure of entire dataset to double-check variables ----
####
variable_desc <- skim(data_processed)

####
# Restrict dataset ----
####
# Include only those fulfilling the decided treatment strategy
data_processed <- data_processed %>%
  filter(!is.na(exp_bin_treat))

# Drop unnecessary variables going forward
data_processed <- data_processed %>% 
  select(patient_id, elig_date_t2dm, qa_date_of_death,
         starts_with("exp_"), # Exposures
         starts_with("cov_"), # Covariates
         starts_with("out_"), # Outcomes
         starts_with("cens_"), # Censoring variables
         contains("_landmark") # ICEs between baseline and landmark
  )

####
# Save output ----
####
# the full data
arrow::write_feather(data_processed, here::here("output", "data", "data_processed.arrow"))
# data for cumulative incidence plots re treatment regimen pattern
# write_feather(data_plots, here::here("output", "data", "data_plots.feather"))
# description of each variable in the dataset
write.csv(variable_desc, file = here::here("output", "data_description", "variable_codebook.csv")) # for L4 reviewing only, not for release
# flow chart eligibility criteria
write.csv(n_elig_excluded_midpoint6, file = here::here("output", "data_description", "n_elig_excluded_midpoint6.csv"))
write.csv(n_elig_excluded, file = here::here("output", "data_description", "n_elig_excluded.csv"))
# descriptive data re treatment patterns, events between index_date and pandemic start, and main outcome
write.csv(n_exp_out_midpoint6, file = here::here("output", "data_description", "n_exp_out_midpoint6.csv"))
write.csv(n_exp_out, file = here::here("output", "data_description", "n_exp_out.csv"))
# flow chart quality assurance
write.csv(n_qa_excluded_midpoint6, file = here::here("output", "data_description", "n_qa_excluded_midpoint6.csv"))
# flow chart completeness criteria
write.csv(n_completeness_excluded_midpoint6, file = here::here("output", "data_description", "n_completeness_excluded_midpoint6.csv"))