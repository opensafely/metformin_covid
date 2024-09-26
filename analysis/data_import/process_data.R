################################################################################
# A custom made function to process the extracted data incl. application of the diabetes algorithm
# diabetes_algorithm.R based on and Credits to https://github.com/opensafely/post-covid-diabetes/tree/main
################################################################################
# load libraries
library(forcats)
# load custom user functions
source(here::here("analysis", "data_import", "fct_case_when.R"))
source(here::here("analysis", "data_import", "diabetes_algorithm.R"))
source(here::here("analysis", "data_import", "define_status_and_fu_primary.R"))
source(here::here("analysis", "data_import", "add_period_cuts.R"))

# create function "process_data"
process_data <- function(data_extracted, study_dates, treat_window_days){ # grace period 10 days (baseline_date + 9)
  data_processed <- data_extracted %>%
    mutate(
      # POPULATION/DEMOGRAPHIC ----
      cov_cat_age = cut(
        cov_num_age,
        breaks = c(18, 40, 60, 80, Inf),
        labels = c("18-39", "40-59", "60-79", "80+"),
        right = FALSE),

      cov_cat_sex = fct_case_when(
        cov_cat_sex == "female" ~ "Female",
        cov_cat_sex == "male" ~ "Male",
        TRUE ~ NA_character_), # any other gender is again excluded at a later stage by applying the QA criteria

      cov_cat_ethnicity = fct_case_when(
        cov_cat_ethnicity == "1" ~ "White",
        cov_cat_ethnicity == "4" ~ "Black",
        cov_cat_ethnicity == "3" ~ "South Asian",
        cov_cat_ethnicity == "2" ~ "Mixed",
        cov_cat_ethnicity == "5" ~ "Other",
        cov_cat_ethnicity == "0" ~ "Unknown",
        TRUE ~ NA_character_),

      cov_cat_deprivation_5 = fct_case_when(
        cov_cat_deprivation_5 == "5 (least deprived)" ~ "5 (least deprived)",
        cov_cat_deprivation_5 == "4" ~ "4",
        cov_cat_deprivation_5 == "3" ~ "3",
        cov_cat_deprivation_5 == "2" ~ "2",
        cov_cat_deprivation_5 == "1 (most deprived)" ~ "1 (most deprived)",
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

      cov_cat_rural_urban = fct_case_when(
        cov_cat_rural_urban %in% c(1,2) ~ "Urban conurbation",
        cov_cat_rural_urban %in% c(3,4) ~ "Urban city or town",
        cov_cat_rural_urban %in% c(5,6,7,8) ~ "Rural town or village",
        TRUE ~ NA_character_),

      cov_cat_stp = as.factor(cov_cat_stp),

      
      # MAIN ELIGIBILITY - HISTORY OF T2DM ----
      ## First, define helper variables needed, esp. for step 5 in diabetes algorithm
      tmp_cov_year_latest_diabetes_diag = as.integer(format(tmp_cov_date_latest_diabetes_diag,"%Y")),
      tmp_age_1st_diag = tmp_cov_year_latest_diabetes_diag - qa_num_birth_year,
      tmp_age_1st_diag = replace(tmp_age_1st_diag, which(tmp_age_1st_diag < 0), NA),
      tmp_age_under_35_30_1st_diag = ifelse(!is.na(tmp_age_1st_diag) & 
                                              (tmp_age_1st_diag < 35 & (cov_cat_ethnicity == "White" | cov_cat_ethnicity == "Mixed" | cov_cat_ethnicity == "Other")) | 
                                              (tmp_age_1st_diag < 30), "Yes", "No"),
      # earliest HbA1c date for only those with >=47.5
      tmp_hba1c_date_step7 = as_date(case_when(tmp_cov_num_max_hba1c_mmol_mol >= 47.5 ~ pmin(tmp_cov_date_max_hba1c, na.rm = TRUE))),
      # take the first process code date in those individuals that have 5 or more process codes
      tmp_over5_pocc_step7 = as_date(case_when(tmp_cov_count_poccdm_ctv3 >= 5 ~ pmin(tmp_cov_date_poccdm, na.rm = TRUE))),
    ) %>%
      ## Second, apply the diabetes algorithm
    diabetes_algo() %>%
    mutate(
      ## Third, extract T2DM as a separate variable
      cov_bin_t2dm = case_when(cov_cat_diabetes == "T2DM" ~ TRUE, TRUE ~ FALSE),


      # TREATMENT ---- keep the structure as it is, may want to add more treatment strategies (other OADs) in future
      # Time-between positive test and day of treatment
      exp_num_tb_postest_treat =
        if_else(!is.na(exp_date_first_metfin),
                difftime(exp_date_first_metfin, baseline_date, units = "days") %>%
                  as.numeric(), NA_real_),

      # Treatment window
      exp_treat_window = baseline_date + days(treat_window_days), # e.g. grace period 10 (baseline_date + 9 days)

      # Any treatment in follow-up
      exp_any_treatment_cat =
        case_when(!is.na(exp_date_first_metfin) ~ "Treated",
                  TRUE ~ "Untreated") %>%
        factor(levels = c("Untreated", "Treated")),

      # Flag where treatment date falls in treatment assignment window, 1/0
      exp_treat_check =
        if_else(exp_date_first_metfin >= baseline_date &
                  exp_date_first_metfin <= exp_treat_window,1,0),

      # Flag where treatment date falls after treat_window, 1/0
      exp_treat_after_treat_window =
        if_else(exp_date_first_metfin >= baseline_date &
                  exp_date_first_metfin > exp_treat_window,1,0),

      # Flag where treatment date falls in treatment assignment window, and assign "Treated" vs "Untreated"
      exp_treatment =
        case_when(!is.na(exp_date_first_metfin) & exp_treat_check == 1 ~ "Treated",
                  TRUE ~ "Untreated") %>%
        factor(levels = c("Untreated", "Treated")),

      # Treatment date, any window
      exp_any_treatment_date =
        if_else(exp_any_treatment_cat == "Treated",
                exp_date_first_metfin,
                NA_Date_),

      # Treatment date, if treatment date falls in treatment assignment window
      exp_treatment_date =
        if_else(exp_treatment == "Treated",
                exp_date_first_metfin,
                NA_Date_),


      # CONFOUNDERS / COVARIATES ----
      # Remove biologically implausible numerical BMI values
      cov_num_bmi = replace(cov_num_bmi, cov_num_bmi > 70 | cov_num_bmi < 12, NA),
      # Combine BMI variables to create one history of obesity variable
      cov_bin_obesity = fct_case_when(
        cov_bin_obesity == TRUE | cov_cat_bmi_groups == "Obese (>30)" ~ "Obese (>30)",
        TRUE ~ NA_character_),

      cov_cat_smoking_status = fct_case_when(
        cov_cat_smoking_status == "S" ~ "Smoker",
        cov_cat_smoking_status == "E" ~ "Ever",
        cov_cat_smoking_status == "N" ~ "Never",
        cov_cat_smoking_status == "M" ~ "Unknown",
        TRUE ~ NA_character_),

      cov_cat_smoking_status_comb = fct_case_when(
        cov_cat_smoking_status == "S" ~ "Current",
        cov_cat_smoking_status == "E" ~ "Former",
        cov_cat_smoking_status %in% c("N", "M") ~ "Never and unknown",
        TRUE ~ NA_character_),

      # Remove implausible consultation rates to 365 max (average of one per day) / this variable is derived from count_for_patient() so does not contain any missing (only 0s).
      # cov_num_consultation_rate = replace(cov_num_consultation_rate, cov_num_consultation_rate > 365, 365),

      # TC/HDL ratio values: remove biologically implausible values: https://doi.org/10.1093/ije/dyz099
      ## remove TC < 1.75 or > 20
      ## remove HDL < 0.4 or > 5
      ## remove ratios < 1 or > 50
      tmp_cov_num_cholesterol = replace(tmp_cov_num_cholesterol, tmp_cov_num_cholesterol < 1.75 | tmp_cov_num_cholesterol > 20, NA),
      tmp_cov_num_hdl_cholesterol = replace(tmp_cov_num_hdl_cholesterol, tmp_cov_num_hdl_cholesterol < 0.4 | tmp_cov_num_hdl_cholesterol > 5, NA),
      cov_num_tc_hdl_ratio = tmp_cov_num_cholesterol / tmp_cov_num_hdl_cholesterol,
      cov_num_tc_hdl_ratio = replace(cov_num_tc_hdl_ratio, cov_num_tc_hdl_ratio > 50 | cov_num_tc_hdl_ratio < 1, NA),

    ) %>%

    mutate(
      # Outcome prep --> outcomes are added in add_status_and_fu_primary() functions below
      # Separate non-covid deaths from covid deaths, since non-covid death is a censoring event and covid death is an outcome
      study_window = baseline_date + days(28),
      out_date_dereg_28 = case_when(!is.na(out_date_dereg) & (out_date_dereg < study_window) ~ out_date_dereg,
                             TRUE ~ NA_Date_),
      out_date_covid_hosp_28 = case_when(!is.na(out_date_covid_hosp) & (out_date_covid_hosp < study_window) ~ out_date_covid_hosp,
                                    TRUE ~ NA_Date_),
      out_date_death_28 = case_when(!is.na(qa_date_of_death) & (qa_date_of_death < study_window) ~ qa_date_of_death,
                                    TRUE ~ NA_Date_),
      out_date_covid_death_28 = case_when(!is.na(qa_date_of_death) & out_bin_death_cause_covid == TRUE & (qa_date_of_death < study_window) ~ qa_date_of_death,
                                    TRUE ~ NA_Date_),
      out_date_noncovid_death_28 = case_when(!is.na(qa_date_of_death) & out_bin_death_cause_covid == FALSE & (qa_date_of_death < study_window) ~ qa_date_of_death,
                                    TRUE ~ NA_Date_),
      # irrespective of window
      out_date_covid_death = case_when(!is.na(qa_date_of_death) & out_bin_death_cause_covid == TRUE ~ qa_date_of_death,
                                       TRUE ~ NA_Date_),
      out_date_noncovid_death = case_when(!is.na(qa_date_of_death) & out_bin_death_cause_covid == FALSE ~ qa_date_of_death, # FALSE in a _bin_ (ehrQL logical) includes missing! => !is.na(qa_date_of_death) is needed
                                          TRUE ~ NA_Date_),

      # make distinction between noncovid hosp admission and covid hosp admissions?

    ) %>%

    # adds column status_primary and fu_primary
    add_status_and_fu_primary() %>%

    ## some patients might experience one of the outcomes (prim outcome) on or before day of treatment. If so, patients will be categorised as untreated!!

    add_period_cuts(study_dates = study_dates) %>%

  # drop unnecessary helper variables during the data processing (esp. above during diabetes algo)
    select(-contains("tmp"), -contains("step"))
}
