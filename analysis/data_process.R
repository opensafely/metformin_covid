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
print('Import libraries and functions')
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
library('ggplot2')
source(here::here("analysis", "functions", "fn_extract_data.R"))
source(here::here("analysis", "functions", "utility.R"))
source(here::here("analysis", "functions", "fn_quality_assurance_midpoint6.R"))
source(here::here("analysis", "functions", "fn_completeness_criteria_midpoint6.R"))
source(here::here("analysis", "functions", "fn_elig_criteria_midpoint6.R"))


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_description"))


# Import dates ------------------------------------------------------------
print('Import dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")
mid2018_date <- as.Date(study_dates$mid2018_date, format = "%Y-%m-%d")
mid2019_date <- as.Date(study_dates$mid2019_date, format = "%Y-%m-%d")
pandemicstart_date <- as.Date(study_dates$pandemicstart_date, format = "%Y-%m-%d")


# Define redaction threshold ----------------------------------------------
threshold <- 6


# Import the dataset and pre-process --------------------------------------
print('Import the dataset and pre-process')
input_filename <- "dataset.arrow"
data_extracted <- fn_extract_data(input_filename)


# Process the data --------------------------------------------------------
print('Process the data')
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
                                                                              levels = c("below 42", "42-58", "59-75", "above 75", "Unknown")), TRUE ~ elig_cat_hba1c_landmark_mmol_mol)
    )


# Assign calendar period of T2DM diagnosis --------------------------------
print('Assign calendar period of T2DM diagnosis')
# Sequence of dates, in period of possible eligible T2DM dates, to use as monthly break points
seq_dates_start_interval_month <- seq(
  from = mid2018_date,
  to = mid2019_date,
  by = "1 month"
)
# Create period_month column
data_processed <- data_processed %>%
  mutate(
    cov_num_period_month = cut(
      elig_date_t2dm,
      breaks = seq_dates_start_interval_month,
      include.lowest = TRUE,
      right = FALSE,
      labels = FALSE
    ),
    cov_num_period_month = as.numeric(cov_num_period_month)
  )


# Modify dummy data -------------------------------------------------------
print('Modify dummy data')
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  message("Running locally, adapt dummy data")
  source("analysis/modify_dummy_data.R")
  message("Dummy data successfully modified")
}


# Assign the landmark date and max follow-up date -------------------------
print('Assign the landmark date and max follow-up date')
data_processed <- data_processed %>% 
  mutate(landmark_date = elig_date_t2dm + days(183)) %>% 
  mutate(max_fup_date = landmark_date + days(730))


# Apply the quality assurance criteria ------------------------------------
print('Apply the quality assurance criteria')
qa <- fn_quality_assurance_midpoint6(data_processed, study_dates, threshold)
n_qa_excluded_midpoint6 <- qa$n_qa_excluded_midpoint6
data_processed <- qa$data_processed


# Apply the completeness criteria -----------------------------------------
print('Apply the completeness criteria')
completeness <- fn_completeness_criteria_midpoint6(data_processed, threshold)
n_completeness_excluded_midpoint6 <- completeness$n_completeness_excluded_midpoint6
data_processed <- completeness$data_processed


# Apply the eligibility criteria ------------------------------------------
print('Apply the eligibility criteria')
# Our primary eligibility window to define incident T2DM is mid2018-mid2019, but maybe we may want to extend the window until max. mid2013 later on 
# => use function with loop that can be mapped to other windows, depending on "years_in_days" input
eligibility <- fn_elig_criteria_midpoint6(data_processed, study_dates, years_in_days = 0)
n_elig_excluded <- eligibility$n_elig_excluded
n_elig_excluded_midpoint6 <- eligibility$n_elig_excluded_midpoint6
data_processed <- eligibility$data_processed


# Assign treatment/exposure -----------------------------------------------
print('Assign treatment/exposure')
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
print('Define final treatment variable')
data_processed <- data_processed %>% 
  mutate(exp_bin_treat = case_when(exp_bin_metfin_mono == TRUE ~ 1,
                                   exp_bin_treat_nothing == TRUE ~ 0,
                                   TRUE ~ NA_real_))


# Assign outcomes ---------------------------------------------------------
print('Assign outcomes')
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
    # Other events.
    out_bin_death_afterlandmark = out_bin_covid_death_afterlandmark | (!is.na(qa_date_of_death) & (qa_date_of_death > out_date_covid_death_afterlandmark)), # this additional step is needed to ensure order of death dates in dummy data (and no impact on real data)
    out_date_death_afterlandmark = case_when(out_bin_death_afterlandmark == TRUE ~ qa_date_of_death, 
                                             TRUE ~ as.Date(NA)),
    out_bin_noncoviddeath_afterlandmark = (!is.na(qa_date_of_death) & qa_date_of_death > landmark_date) & is.na(out_date_covid_death_afterlandmark),
    out_date_noncoviddeath_afterlandmark = case_when(out_bin_noncoviddeath_afterlandmark == TRUE ~ qa_date_of_death, 
                                                     TRUE ~ as.Date(NA)),
    cens_bin_ltfu_afterlandmark = (!is.na(cens_date_dereg) & cens_date_dereg > landmark_date) & (is.na(qa_date_of_death) | (!is.na(qa_date_of_death) & qa_date_of_death >= cens_date_dereg)),
    cens_date_ltfu_afterlandmark = case_when(cens_bin_ltfu_afterlandmark == TRUE ~ cens_date_dereg, 
                                             TRUE ~ as.Date(NA)),
    # Negative control, after landmark and after pandemic
    out_bin_fracture_afterlandmark = !is.na(out_date_fracture) & out_date_fracture > landmark_date,
    out_date_fracture_afterlandmark = case_when(out_bin_fracture_afterlandmark == TRUE ~ out_date_fracture, 
                                             TRUE ~ as.Date(NA)),
    out_bin_fracture_afterpandemic = !is.na(out_date_fracture) & out_date_fracture > pandemicstart_date,
    out_date_fracture_afterpandemic= case_when(out_bin_fracture_afterpandemic == TRUE ~ out_date_fracture, 
                                                TRUE ~ as.Date(NA)),
    # Positive control, after landmark and after pandemic
    out_bin_dm_death_afterlandmark = !is.na(out_date_dm_death) & out_date_dm_death > landmark_date,
    out_date_dm_death_afterlandmark = case_when(out_bin_dm_death_afterlandmark == TRUE ~ out_date_dm_death, 
                                                TRUE ~ as.Date(NA)),
    out_bin_dm_death_afterpandemic = !is.na(out_date_dm_death) & out_date_dm_death > pandemicstart_date,
    out_date_dm_death_afterpandemic= case_when(out_bin_dm_death_afterpandemic == TRUE ~ out_date_dm_death, 
                                               TRUE ~ as.Date(NA)),
    # In INTERVENTION: Identify all metformin prescription (combo and mono) in 6m prior to pandemic start, may be used to censor those who stopped before pandemic
    cens_bin_metfin_pandemicstart = exp_bin_metfin_mono == TRUE & !is.na(exp_date_metfin_mono_last) & exp_date_metfin_mono_last >= study_dates$pandemicstart_date - days(183),
    cens_date_metfin_pandemicstart = case_when(cens_bin_metfin_pandemicstart == TRUE ~ exp_date_metfin_mono_last, 
                                             TRUE ~ as.Date(NA)),
    # In CONTROL: Identify all who started metformin mono after landmark
    cens_bin_metfin_start_cont = exp_bin_treat_nothing == TRUE & !is.na(exp_date_metfin_mono_first) & exp_date_metfin_mono_first > landmark_date,
    cens_date_metfin_start_cont = case_when(cens_bin_metfin_start_cont == TRUE ~ exp_date_metfin_mono_first, 
                                               TRUE ~ as.Date(NA))
    )

# HbA1c covariate and HbA1c & lipids plot --------------------------------
print('HbA1c covariate and HbA1c & lipids plot')
# Above, I excluded all with HbA1c >75 mmol/mol. Unfortunately, the cox RA still "sees" this level and would exclude the variable if left as is (due to too 0 events in that subgroup)
data_processed <- data_processed %>% 
  mutate(cov_cat_hba1c_mmol_mol = fn_case_when(
    cov_cat_hba1c_mmol_mol == "below 42" ~ "below 42",
    cov_cat_hba1c_mmol_mol == "42-58" ~ "42-58",
    cov_cat_hba1c_mmol_mol == "59-75" ~ "59-75",
    cov_cat_hba1c_mmol_mol == "Unknown" ~ "Unknown",
    TRUE ~ NA_character_))

# plot the numeric lab values to check distribution (to see if there is any grouping indicating the use of comparators)
hba1c_plot <- data_processed %>% 
  drop_na(cov_num_hba1c_mmol_mol) %>% 
  ggplot(aes(x = cov_num_hba1c_mmol_mol)) +
  geom_density(fill = "blue", color = "black") +
  labs(title = "Density Plot of HbA1c",
       x = "HbA1c",
       y = "Density")
totchol_plot <- data_processed %>% 
  drop_na(tmp_cov_num_cholesterol) %>% 
  ggplot(aes(x = tmp_cov_num_cholesterol)) +
  geom_density(fill = "blue", color = "black") +
  labs(title = "Density Plot of Total Cholesterol",
       x = "Tot Cholesterol",
       y = "Density")
hdlchol_plot <- data_processed %>% 
  drop_na(tmp_cov_num_hdl_cholesterol) %>% 
  ggplot(aes(x = tmp_cov_num_hdl_cholesterol)) +
  geom_density(fill = "blue", color = "black") +
  labs(title = "Density Plot of HDL Cholesterol",
       x = "HDL Cholesterol",
       y = "Density")

# Assign Cox variables ------------------------------------------------------
print('Assign Cox variables')
data_processed <- data_processed %>% 
  mutate(
    cox_date_severecovid = pmin(out_date_severecovid_afterlandmark, 
                                out_date_noncoviddeath_afterlandmark,
                                cens_date_ltfu_afterlandmark,
                                max_fup_date,
                                studyend_date,
                                na.rm = TRUE),
    cox_tt_severecovid = difftime(cox_date_severecovid,
                                  landmark_date,
                                  units = "days") %>% as.numeric(),
    cox_date_covid = pmin(out_date_covid_afterlandmark, 
                          out_date_noncoviddeath_afterlandmark,
                          cens_date_ltfu_afterlandmark,
                          max_fup_date,
                          studyend_date,
                          na.rm = TRUE),
    cox_tt_covid = difftime(cox_date_covid,
                            landmark_date,
                            units = "days") %>% as.numeric(),
    cox_date_longcovid_virfat = pmin(out_date_longcovid_virfat_afterlandmark, 
                                     out_date_death_afterlandmark,
                                     cens_date_ltfu_afterlandmark,
                                     max_fup_date,
                                     studyend_date,
                                     na.rm = TRUE),
    cox_tt_longcovid_virfat = difftime(cox_date_longcovid_virfat,
                                       landmark_date,
                                       units = "days") %>% as.numeric(),
    # for cox reusable action: exposure start date for those who start, missing for all others 
    cox_date_metfin_start_within6m = case_when(exp_bin_metfin_mono == TRUE ~ landmark_date, 
                                               TRUE ~ as.Date(NA)),
    # for cox reusable action: define cox_stop for negative control (fracture) and positive control (diabetes deaths)
    cox_date_fracture_landmark = pmin(out_date_fracture_afterlandmark, 
                                      out_date_death_afterlandmark,
                                      cens_date_ltfu_afterlandmark,
                                      max_fup_date,
                                      studyend_date,
                                      na.rm = TRUE),
    cox_date_fracture_pandemic = pmin(out_date_fracture_afterpandemic, 
                                      out_date_death_afterlandmark,
                                      cens_date_ltfu_afterlandmark,
                                      max_fup_date,
                                      studyend_date,
                                      na.rm = TRUE),
    cox_date_dm_death_landmark = pmin(out_date_dm_death_afterlandmark, 
                                      out_date_death_afterlandmark,
                                      cens_date_ltfu_afterlandmark,
                                      max_fup_date,
                                      studyend_date,
                                      na.rm = TRUE),
    cox_date_dm_death_pandemic = pmin(out_date_dm_death_afterpandemic, 
                                      out_date_death_afterlandmark,
                                      cens_date_ltfu_afterlandmark,
                                      max_fup_date,
                                      studyend_date,
                                      na.rm = TRUE)
  )


# Save and inspect full processed dataset ---------------------------------
print('Save and inspect full processed dataset')
data_processed_full <- data_processed
data_processed_full_desc <- skim(data_processed_full)
write.csv(data_processed_full_desc, file = here::here("output", "data_description", "data_processed_full_desc.csv")) # for L4 reviewing only, not for release
arrow::write_feather(data_processed_full, here::here("output", "data", "data_processed_full.arrow"))


# Restrict the dataset for the pipeline onwards, and explore deaths/ltfu before landmark & pandemic ----
print('Restrict the dataset for the pipeline onwards, and explore deaths/ltfu before landmark & pandemic')
# (1) Explore deaths/ltfu
# died/ltfu between elig_date_t2dm and landmark -> add to data_processed to filter out
data_processed <- data_processed %>%
  mutate(death_landmark = (!is.na(qa_date_of_death) & (qa_date_of_death > elig_date_t2dm) & (qa_date_of_death <= landmark_date))) %>% 
  mutate(ltfu_landmark = (!is.na(cens_date_dereg) 
                          & (cens_date_dereg > elig_date_t2dm) & (cens_date_dereg <= landmark_date) 
                          & (is.na(qa_date_of_death) | (!is.na(qa_date_of_death) & qa_date_of_death >= cens_date_dereg)))) %>% 
  mutate(death_ltfu_landmark = death_landmark | ltfu_landmark)
# died/ltfu between landmark and pandemic start -> create and add to separate dataset for descriptive table ones
data_processed_death_ltfu <- data_processed %>%
  filter(!is.na(exp_bin_treat)) %>% # Filter out those with missing exp_bin_treat, but retain those who died/ltfu before landmark for descriptive purposes
  mutate(death_ltfu_pandemic = (!is.na(out_date_death_afterlandmark) & (out_date_death_afterlandmark < pandemicstart_date)) # do not include day of start of pandemic, since these will count (and start of pandemic is a bit arbitrary chosen)
         | (!is.na(cens_date_ltfu_afterlandmark) & (cens_date_ltfu_afterlandmark < pandemicstart_date)))
# overall flag
data_processed_death_ltfu <- data_processed_death_ltfu %>%
  mutate(
    death_ltfu_landmark_pandemic = death_ltfu_landmark | death_ltfu_pandemic,
    death_ltfu_pandemic_without_landmark = case_when(death_ltfu_pandemic & !death_ltfu_landmark ~ TRUE,
                                                     !death_ltfu_pandemic & !death_ltfu_landmark ~ FALSE,
                                                     TRUE ~ NA)
    )
# count died/ltfu between elig_date_t2dm and landmark and those fulfilling one of the two final treatment strategies
count <- data_processed %>%
  summarise(
    n_death_landmark = sum(death_landmark, na.rm = TRUE),
    n_ltfu_landmark = sum(ltfu_landmark, na.rm = TRUE),
    n_exp_bin_treat = sum(is.na(exp_bin_treat), na.rm = TRUE))
    
# (2) Filter: only keep those fulfilling one of the two final treatment strategies and still alive and in care at landmark 
data_processed <- data_processed %>%
  filter(
    !is.na(exp_bin_treat),
    !death_ltfu_landmark
    )
# tibble out incl. midpoint6 rounding (!)
n_restricted_midpoint6 <- tibble(
  n_before_exclusion_midpoint6 = fn_roundmid_any(nrow(data_processed_death_ltfu), threshold),
  n_exp_bin_treat_midpoint6 = fn_roundmid_any(count$n_exp_bin_treat, threshold),
  n_death_landmark_midpoint6 = fn_roundmid_any(count$n_death_landmark, threshold),
  n_ltfu_landmark_midpoint6 = fn_roundmid_any(count$n_ltfu_landmark, threshold),
  n_after_exclusion_midpoint6 = fn_roundmid_any(nrow(data_processed), threshold))

# (3) Filter main dataset: Only keep necessary variables
data_processed <- data_processed %>% 
  select(patient_id, 
         elig_date_t2dm, 
         landmark_date,
         max_fup_date,
         qa_date_of_death, # To identify other deaths, e.g. between eligibility and landmark
         starts_with("exp_"), # Exposures
         starts_with("strat_"), # Stratification variable
         starts_with("cov_"), # Covariates
         starts_with("out_"), # Outcomes
         starts_with("cens_"), # Censoring variable
         starts_with("cox_"), # Cox variables
  )


# Save aggregate output and restricted processed dataset -------------------
print('Save aggregate output and restricted processed dataset')
# flow chart quality assurance
write.csv(n_qa_excluded_midpoint6, file = here::here("output", "data_description", "n_qa_excluded_midpoint6.csv"))
# flow chart completeness criteria
write.csv(n_completeness_excluded_midpoint6, file = here::here("output", "data_description", "n_completeness_excluded_midpoint6.csv"))
# flow chart eligibility criteria
write.csv(n_elig_excluded_midpoint6, file = here::here("output", "data_description", "n_elig_excluded_midpoint6.csv"))
write.csv(n_elig_excluded, file = here::here("output", "data_description", "n_elig_excluded.csv"))
# flow chart eligibility criteria 2 (restrict to eligible treatment strategies and alive/still registered at landmark)
write.csv(n_restricted_midpoint6, file = here::here("output", "data_description", "n_restricted_midpoint6.csv"))
# final (restricted) dataset
arrow::write_feather(data_processed, here::here("output", "data", "data_processed.arrow"))
# for descriptive purpose, to create table ones: Incl. flags for those who died/LTFU during landmark period and between landmark and pandemic start
arrow::write_feather(data_processed_death_ltfu, here::here("output", "data", "data_processed_death_ltfu.arrow"))
# Lab value dens plots
ggsave(filename = here::here("output", "data_description", "hba1c_plot.png"), hba1c_plot, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "data_description", "totchol_plot.png"), totchol_plot, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "data_description", "hdlchol_plot.png"), hdlchol_plot, width = 20, height = 20, units = "cm")
