####
## This script does the following:
# 1. Import/extract feather dataset from OpenSAFELY
# 2. Basic type formatting of variables -> fn_extract_data.R()
# 3. Process covariates
# 4. Modify dummy data if run locally
# 5. Evaluate/apply the quality assurance criteria -> fn_quality_assurance_midpoint6()
# 6. Evaluate/apply the completeness criteria: -> fn_completeness_criteria_midpoint6()
# 7. Evaluate/apply the eligibility criteria: -> fn_elig_criteria_midpoint6()
# 8. Add the landmark date, treatment/exposure and outcome variables
# 9. Restrict the dataset
# 10. Save the full dataset, the restricted dataset, the codebook and all flowchart tables
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
source(here::here("analysis", "functions", "fn_elig_criteria_midpoint6.R"))


# Create directories for output -------------------------------------------
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_description"))


# Import dates ------------------------------------------------------------
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


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
      cov_cat_ethnicity == "Other" | cov_cat_ethnicity == "Unknown" | cov_cat_ethnicity == "Mixed" ~ "Other", # collapsed due to too few events in "Unknown" and "Mixed", to avoid variable exclusion in cox model
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
    
    cov_bin_obesity = cov_bin_obesity == TRUE | cov_cat_bmi_groups == "Obese (>30)",
    
    # TC/HDL ratio values: remove biologically implausible values: https://doi.org/10.1093/ije/dyz099
    ## remove TC < 1.75 or > 20; remove HDL < 0.4 or > 5; remove ratios < 1 or > 50
    tmp_cov_num_cholesterol = replace(tmp_cov_num_cholesterol, tmp_cov_num_cholesterol < 1.75 | tmp_cov_num_cholesterol > 20, NA_real_),
    tmp_cov_num_hdl_cholesterol = replace(tmp_cov_num_hdl_cholesterol, tmp_cov_num_hdl_cholesterol < 0.4 | tmp_cov_num_hdl_cholesterol > 5, NA_real_),
    cov_num_tc_hdl_ratio = tmp_cov_num_cholesterol / tmp_cov_num_hdl_cholesterol,
    cov_num_tc_hdl_ratio = replace(cov_num_tc_hdl_ratio, cov_num_tc_hdl_ratio > 50 | cov_num_tc_hdl_ratio < 1, NA_real_),
    
    # TC/HDL ratio categories: https://www.urmc.rochester.edu/encyclopedia/content?ContentTypeID=167&ContentID=lipid_panel_hdl_ratio#:~:text=Most%20healthcare%20providers%20want%20the,1%20is%20considered%20very%20good.
    cov_cat_tc_hdl_ratio = cut(
      cov_num_tc_hdl_ratio,
      breaks = c(1, 3.5, 5.11, 50), # 50 is upper limit, see above -> NA
      labels = c("below 3.5:1" ,"3.5:1 to 5:1", "above 5:1"),
      right = FALSE),
    cov_cat_tc_hdl_ratio = case_when(is.na(cov_num_tc_hdl_ratio) ~ factor("Unknown", 
                                                                          levels = c("below 3.5:1", "3.5:1 to 5:1", "above 5:1", "Unknown")), TRUE ~ cov_cat_tc_hdl_ratio),
    
    # HbA1c categories: https://www.southtees.nhs.uk/resources/the-hba1c-test/
    ## remove HbA1c > 120; remove HbA1c below 0
    cov_num_hba1c_mmol_mol = replace(cov_num_hba1c_mmol_mol, cov_num_hba1c_mmol_mol < 0.00 | cov_num_hba1c_mmol_mol > 120.00, NA_real_),
    cov_cat_hba1c_mmol_mol = cut(
      cov_num_hba1c_mmol_mol,
      breaks = c(0, 42, 59, 76, 120), # 120 is upper limit, above NA
      labels = c("below 42" ,"42-58", "59-75", "above 75"),
      right = FALSE),
    cov_cat_hba1c_mmol_mol = case_when(is.na(cov_cat_hba1c_mmol_mol) ~ factor("Unknown", 
                                                                              levels = c("below 42", "42-58", "59-75", "above 75", "Unknown")), TRUE ~ cov_cat_hba1c_mmol_mol),
    elig_num_hba1c_landmark_mmol_mol = replace(elig_num_hba1c_landmark_mmol_mol, elig_num_hba1c_landmark_mmol_mol < 0.00 | elig_num_hba1c_landmark_mmol_mol > 120.00, NA_real_),
    elig_cat_hba1c_landmark_mmol_mol = cut(
      elig_num_hba1c_landmark_mmol_mol,
      breaks = c(0, 42, 59, 76, 120), # 120 is upper limit, above NA
      labels = c("below 42" ,"42-58", "59-75", "above 75"),
      right = FALSE),
    elig_cat_hba1c_landmark_mmol_mol = case_when(is.na(elig_cat_hba1c_landmark_mmol_mol) ~ factor("Unknown", 
                                                                              levels = c("below 42", "42-58", "59-75", "above 75", "Unknown")), TRUE ~ elig_cat_hba1c_landmark_mmol_mol),
    )


# Modify dummy data -------------------------------------------------------
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  message("Running locally, adapt dummy data")
  source("analysis/modify_dummy_data.R")
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
eligibility <- fn_elig_criteria_midpoint6(data_processed, study_dates, years_in_days = 0)
n_elig_excluded <- eligibility$n_elig_excluded
n_elig_excluded_midpoint6 <- eligibility$n_elig_excluded_midpoint6
data_processed <- eligibility$data_processed


# Assign the landmark date and max follow-up date -------------------------
data_processed <- data_processed %>% 
  mutate(landmark_date = elig_date_t2dm + days(183)) %>% 
  mutate(max_fup_date = landmark_date + days(730))


# Assign treatment/exposure -----------------------------------------------
data_processed <- data_processed %>% 
  mutate(
    ## (1) Starting metformin within 6 months after T2DM, among those with a T2DM diagnosis, 
    # mid2018 onwards (those with exp_bin_metfin_first before mid2018 were excluded above)
    # metformin combo
    exp_bin_metfin = !is.na(exp_date_metfin_first) & exp_date_metfin_first <= landmark_date,
    # metformin mono
    exp_bin_metfin_mono = !is.na(exp_date_metfin_mono_first) & exp_date_metfin_mono_first <= landmark_date,
                                    
    ## (2) Let's investigate those who did not start any metformin until 6 month landmark, i.e. exp_bin_metfin == FALSE
    # We include exp_bin_metfin == FALSE to avoid counting people who initiated metfin and shortly after DPP4 (clinically rare but possible)
    # We use exp_bin_metfin to cover all metfin prescriptions (combo and mono), otherwise we would allow to count people who initiated DPP4 + metformin
    # DPP4 mono (or combo with SGLT2, but no combo with metfin)
    exp_bin_dpp4_mono = exp_bin_metfin == FALSE & (!is.na(exp_date_dpp4_first) & exp_date_dpp4_first <= landmark_date), 
    # TZD mono
    exp_bin_tzd_mono = exp_bin_metfin == FALSE & (!is.na(exp_date_tzd_first) & exp_date_tzd_first <= landmark_date), 
    # SGLT2 mono (or combo with DPP4)
    exp_bin_sglt2_mono = exp_bin_metfin == FALSE & (!is.na(exp_date_sglt2_first) & exp_date_sglt2_first <= landmark_date),
    # sulfo mono
    exp_bin_sulfo_mono = exp_bin_metfin == FALSE & (!is.na(exp_date_sulfo_first) & exp_date_sulfo_first <= landmark_date),
    # glp1 mono
    exp_bin_glp1_mono = exp_bin_metfin == FALSE & (!is.na(exp_date_glp1_first) & exp_date_glp1_first <= landmark_date),
    # megli mono
    exp_bin_megli_mono = exp_bin_metfin == FALSE & (!is.na(exp_date_megli_first) & exp_date_megli_first <= landmark_date),
    # agi mono
    exp_bin_agi_mono = exp_bin_metfin == FALSE & (!is.na(exp_date_agi_first) & exp_date_agi_first <= landmark_date),
    # insulin mono
    exp_bin_insulin_mono = exp_bin_metfin == FALSE & (!is.na(exp_date_insulin_first) & exp_date_insulin_first <= landmark_date),
    
    ## (3) No prescription at all until 6 months after T2DM diagnosis (landmark)
    exp_bin_treat_nothing = case_when(exp_bin_metfin == FALSE # covers metformin combo and mono
                                      & exp_bin_dpp4_mono == FALSE
                                      & exp_bin_tzd_mono == FALSE 
                                      & exp_bin_sglt2_mono == FALSE 
                                      & exp_bin_sulfo_mono == FALSE
                                      & exp_bin_glp1_mono == FALSE 
                                      & exp_bin_megli_mono == FALSE
                                      & exp_bin_agi_mono == FALSE
                                      & exp_bin_insulin_mono == FALSE ~ TRUE,
                                      TRUE ~ FALSE)
  )


# Define final treatment variable -----------------------------------------
data_processed <- data_processed %>% 
  mutate(exp_bin_treat = case_when(exp_bin_metfin_mono == TRUE ~ 1,
                                   exp_bin_treat_nothing == TRUE ~ 0,
                                   TRUE ~ NA_real_))


# Assign outcomes ---------------------------------------------------------
data_processed <- data_processed %>% 
  mutate(
    # COVID events. These should not happen before landmark date by design - but just in case.
    out_bin_severecovid_afterlandmark = !is.na(out_date_severecovid) & out_date_severecovid > landmark_date,
    out_date_severecovid_afterlandmark = case_when(out_bin_severecovid_afterlandmark == TRUE ~ out_date_severecovid, 
                                                   TRUE ~ as.Date(NA)),
    out_bin_covid_hosp_afterlandmark = !is.na(out_date_covid_hosp) & out_date_covid_hosp > landmark_date,
    out_date_covid_hosp_afterlandmark = case_when(out_bin_covid_hosp_afterlandmark == TRUE ~ out_date_covid_hosp, 
                                                   TRUE ~ as.Date(NA)),
    out_bin_covid_death_afterlandmark = !is.na(out_date_covid_death) & out_date_covid_death > landmark_date,
    out_date_covid_death_afterlandmark = case_when(out_bin_covid_death_afterlandmark == TRUE ~ out_date_covid_death, 
                                                  TRUE ~ as.Date(NA)),
    out_bin_covid_afterlandmark = !is.na(out_date_covid) & out_date_covid > landmark_date,
    out_date_covid_afterlandmark = case_when(out_bin_covid_afterlandmark == TRUE ~ out_date_covid, 
                                                   TRUE ~ as.Date(NA)),
    out_bin_longcovid_afterlandmark = !is.na(out_date_longcovid) & out_date_longcovid > landmark_date,
    out_date_longcovid_afterlandmark = case_when(out_bin_longcovid_afterlandmark == TRUE ~ out_date_longcovid, 
                                             TRUE ~ as.Date(NA)),
    out_bin_virfat_afterlandmark = !is.na(out_date_virfat) & out_date_virfat > landmark_date,
    out_date_virfat_afterlandmark = case_when(out_bin_virfat_afterlandmark == TRUE ~ out_date_virfat, 
                                                 TRUE ~ as.Date(NA)),
    out_bin_longcovid_virfat_afterlandmark = !is.na(out_date_longcovid_virfat) & out_date_longcovid_virfat > landmark_date,
    out_date_longcovid_virfat_afterlandmark = case_when(out_bin_longcovid_virfat_afterlandmark == TRUE ~ out_date_longcovid_virfat, 
                                              TRUE ~ as.Date(NA)),
    # Other events. These may happen between eligibility and landmark date (but rare).
    out_bin_death_afterlandmark = (!is.na(qa_date_of_death) & qa_date_of_death > landmark_date) & is.na(out_date_covid_death),
    out_date_death_afterlandmark = case_when(out_bin_death_afterlandmark == TRUE ~ qa_date_of_death, 
                                             TRUE ~ as.Date(NA)),
    cens_bin_ltfu_afterlandmark = !is.na(cens_date_dereg) & cens_date_dereg > landmark_date,
    cens_date_ltfu_afterlandmark = case_when(cens_bin_ltfu_afterlandmark == TRUE ~ cens_date_dereg, 
                                            TRUE ~ as.Date(NA)),
    # In INTERVENTION: Identify all metformin prescription (combo and mono) in 6m prior to pandemic start, may be used to censor those who stopped before pandemic
    cens_bin_metfin_pandemicstart = exp_bin_metfin_mono == TRUE & !is.na(exp_date_metfin_mono_last) & exp_date_metfin_mono_last >= study_dates$pandemicstart_date - days(183),
    cens_date_metfin_pandemicstart = case_when(cens_bin_metfin_pandemicstart == TRUE ~ exp_date_metfin_mono_last, 
                                             TRUE ~ as.Date(NA)),
    # In CONTROL: Identify all who started any metformin after landmark
    cens_bin_metfin_start_cont = exp_bin_treat_nothing == TRUE & !is.na(exp_date_metfin_mono_first) & exp_date_metfin_mono_first > landmark_date,
    cens_date_metfin_start_cont = case_when(cens_bin_metfin_start_cont == TRUE ~ exp_date_metfin_mono_first, 
                                               TRUE ~ as.Date(NA))
    )

# HbA1c covariate -----------------------------------------------------------
# Above, I excluded all with HbA1c >75 mmol/mol. Unfortunately, the cox RA still "sees" this level and would exclude the variable if left as is (due to too 0 events in that subgroup)
data_processed <- data_processed %>% 
  mutate(cov_cat_hba1c_mmol_mol = fn_case_when(
    cov_cat_hba1c_mmol_mol == "below 42" ~ "below 42",
    cov_cat_hba1c_mmol_mol == "42-58" ~ "42-58",
    cov_cat_hba1c_mmol_mol == "59-75" ~ "59-75",
    cov_cat_hba1c_mmol_mol == "Unknown" ~ "Unknown",
    TRUE ~ NA_character_))

# Assign Cox variables ------------------------------------------------------
data_processed <- data_processed %>% 
  mutate(
    cox_date_severecovid = pmin(out_date_severecovid_afterlandmark, 
                                out_date_death_afterlandmark,
                                cens_date_ltfu_afterlandmark,
                                max_fup_date,
                                studyend_date,
                                na.rm = TRUE),
    cox_tt_severecovid = difftime(cox_date_severecovid,
                                  landmark_date,
                                  units = "days") %>% as.numeric(),
    cox_date_covid = pmin(out_date_covid_afterlandmark, 
                          out_date_death_afterlandmark,
                          cens_date_ltfu_afterlandmark,
                          max_fup_date,
                          studyend_date,
                          na.rm = TRUE),
    cox_tt_covid = difftime(cox_date_covid,
                            landmark_date,
                            units = "days") %>% as.numeric(),
    cox_date_longcovid_virfat = pmin(out_date_longcovid_virfat_afterlandmark, 
                                     out_date_death_afterlandmark,
                                     out_date_covid_death,
                                     cens_date_ltfu_afterlandmark,
                                     max_fup_date,
                                     studyend_date,
                                     na.rm = TRUE),
    cox_tt_longcovid_virfat = difftime(cox_date_longcovid_virfat,
                                                     landmark_date,
                                                     units = "days") %>% as.numeric(),
    # for cox reusable action: exposure start date for those who start, missing for all others 
    cox_date_metfin_start_within6m = case_when(exp_bin_metfin_mono == TRUE ~ landmark_date, 
                                               TRUE ~ as.Date(NA))
  )


# Save and inspect full processed dataset ---------------------------------
data_processed_full <- data_processed
data_processed_full_desc <- skim(data_processed_full)
write.csv(data_processed_full_desc, file = here::here("output", "data_description", "data_processed_full_desc.csv")) # for L4 reviewing only, not for release
arrow::write_feather(data_processed_full, here::here("output", "data", "data_processed_full.arrow"))


# Restrict the dataset for pipeline onwards -------------------------------
# (1) Include only those fulfilling the final treatment strategy
data_processed <- data_processed %>%
  filter(!is.na(exp_bin_treat))

# (2) Drop unnecessary variables
data_processed <- data_processed %>% 
  select(patient_id, 
         elig_date_t2dm, 
         landmark_date,
         qa_date_of_death, # To identify other deaths, e.g. between eligibility and landmark
         starts_with("exp_"), # Exposures
         starts_with("strat_"), # Stratification variable
         starts_with("cov_"), # Covariates
         starts_with("out_"), # Outcomes
         starts_with("cens_"), # Censoring variable
         starts_with("cox_"), # Cox variables
  )


# Save aggregate output and restricted processed dataset -------------------
# flow chart quality assurance
write.csv(n_qa_excluded_midpoint6, file = here::here("output", "data_description", "n_qa_excluded_midpoint6.csv"))
# flow chart completeness criteria
write.csv(n_completeness_excluded_midpoint6, file = here::here("output", "data_description", "n_completeness_excluded_midpoint6.csv"))
# flow chart eligibility criteria
write.csv(n_elig_excluded_midpoint6, file = here::here("output", "data_description", "n_elig_excluded_midpoint6.csv"))
write.csv(n_elig_excluded, file = here::here("output", "data_description", "n_elig_excluded.csv"))
# final (restricted) dataset
arrow::write_feather(data_processed, here::here("output", "data", "data_processed.arrow"))
