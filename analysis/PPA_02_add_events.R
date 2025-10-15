####
## This script does the following:
# 1. Import multiple-event-per-person long format dataset (df_months, which originates from data_processed.arrow)
# 2. Import the additional dataset with stable and dynamic time-updated covariates needed for the per-protocol analysis (data_processed_ppa.arrow -> df_ppa)
# 3. Assign treatment to df_months
# 4. Merge the stable time-updated (first-ever, constant) covariates from df_ppa to df_months. To test edge cases, further adapt dummy data.
# 5. Assign the stable time-updated (first-ever, constant) in df_months. 
# 6. Reformatting of the dynamic time-updated covariates in a long-format df_ppa, before merging to df_months. To test edge cases, further adapt dummy data.
# 7. Merge and assign the dynamic time-updated covariates from df_ppa to df_months.
# 8. Fill up the Swiss cheese and reassign the lipid ratio and other categorical variables
# 9. Drop unnecessary helper variables created along the way
# 10. Assign and shift outcomes, competing, and censoring events
# 11. Define eligibility
# 12. Save datasets
####


# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(zoo) # for fn_assign_stable_tu_cov.R

source(here::here("analysis", "functions", "fn_assign_treatment.R"))
source(here::here("analysis", "functions", "fn_assign_stable_tu_cov.R"))
source(here::here("analysis", "functions", "fn_assign_dynamic_tu_cov.R"))
source(here::here("analysis", "functions", "fn_add_and_shift_out_comp_cens_events.R"))
source(here::here("analysis", "functions", "fn_fill_forward_incl_baseline.R"))


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data"))


# Import dates ------------------------------------------------------------
print('Import dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# 1 + 2. Import the data ---------------------------------------------------------
print('Import the data')
df_months <- read_feather(here("output", "data", "df_months.arrow"))
df_ppa <- read_feather(here("output", "data", "data_processed_ppa.arrow")) 


# 3. Assign treatment to df_months -------------------------------------------
print('Assign treatment to df_months')
## Definition: Assign treatment status in interval k based on the latest treatment initiation recorded within that same interval. 
## Rationale: Treatment is the focal point of each interval; no lagging or leading is necessary for its definition.

### Moreover, we assume someone stays on treatment once is initiated, i.e., first-ever, time-fixed event 
### Moreover, here, we have defined treatment at baseline already (exp_bin_treat), we only assign treatment in control to flag those who initiate during follow-up (cens_date_metfin_start_cont)
### I can use my preexisting function, which works perfectly for this, based on these rules:
# a) if event date is not NA and happened before the minimum start date of all intervals of a person, then assign 1 (in corresponding flag variable) to all person-intervals
# b) if event date is not NA and happened after the maximum end date of all intervals of a person, then assign 0 (in corresponding flag variable) to all person-intervals
# c) if no event date is recorded (date variable == NA), then assign 0 (in corresponding flag variable) to all person-intervals (assuming no documentation = no event)
# d) if event date happened during follow-up, then assign 1 (in corresponding flag variable) to corresponding person-interval, and 0 to all person-intervals before, and 1 to all person-intervals after (stable/time-fixed event)
date_treat_vars <- c("cens_date_metfin_start_cont")
df_months <- fn_assign_treatment(df_months, 
                                 date_treat_vars,
                                 start_var = "start_date_month", 
                                 end_var = "end_date_month")
# I only have treatment starting date data from the control group. 
# But adding this extra step I get the overall time-varying "treat" variable for both groups (always treat = 1 in intervention)
# My PP question is: "starting metformin monotherapy (within 6 months of T2D diagnosis) vs never starting metformin" and since all those who did not start metformin are assigned to control, only the control group can deviate from protocol. This PP question/estimand also means we are not interested in what is happening after initiation (all subsequent "treat" rows are assigned 1 once initiated)
df_months <- df_months %>% 
  mutate(treat = case_when(exp_bin_treat == 0 & cens_bin_metfin_start_cont == 0 ~ 0,
                           exp_bin_treat == 0 & cens_bin_metfin_start_cont == 1 ~ 1,
                           exp_bin_treat == 1 ~ 1))

# To double-check
# df_months %>%
#   dplyr::select(patient_id, elig_date_t2dm, landmark_date, month, start_date_month, end_date_month,
#                 exp_bin_treat, cens_date_metfin_start_cont, cens_bin_metfin_start_cont, treat,
#                 out_date_severecovid_afterlandmark, out_date_death_afterlandmark,
#                 cens_date_ltfu_afterlandmark,
#                 cov_cat_sex, cov_num_age
#   ) %>%
#   # dplyr::filter(exp_bin_treat == 0) %>%
#   # dplyr::filter(cens_bin_metfin_start_cont == 1)
#   View()


# There are two types of time-updated covariates:
## a) stable time-updated (first-ever, constant) covariates: E.g. new stroke diagnosis, or new CKD diagnosis during follow-up
### These are stable over time, we only track the first ever occurring and then a person has it -> "0" in before the event, "1" after it happened
### These are simple to assign, do not require a long-format (event-level) dataset, and (in our case) make up a big share of the time-updated covariates
### Hence, deal with them first in a simple and efficient way
## b) dynamic time-update covariates: E.g. Hba1c, BMI, lipids (of course something like smoking could feature here as well, we ignore this for now)
### These are dynamic over time, may increase/decrease, in the person-interval dataset we carry over the last value until something new is happening
### These are more complex to assign, we need event-level dataset and work with them in a long format dataset until we have reduced them to 1 measurement (the correct one! see rules) per person/interval
### Once done, we merge them to the main person-month dataset with all other variables. But work on them in a separate long format dataset for as long as possible to be efficient.


# 4. Merge the stable time-updated (first-ever, constant) covariates from df_ppa to df_months --------
print('Merge the stable time-updated (first-ever, constant) covariates from df_ppa to df_months')
# all covariates from df_ppa, except hba1c/lipids/bmi
# it will result in a dataset with "cov_bin_* (comorbidities)" as baseline indicators and "cov_date_* (comorbidities)" for the stable one-time time-updated (first-ever, constant) covariate
# Ignore the "cov_bin_* (comorbidities)" from df_months since "cov_date_* (comorbidities)" from df_ppa was constructed as minimmum of first ever (see dataset_definition_ppa.py). It would have been cleaner to do this directly in df_months - if we would have know about the PPA from the start.
df_months <- df_months %>%
  left_join(df_ppa %>% 
              select(patient_id, starts_with("cov_date") & !matches("cov_date_bmi|cov_date_hba1c|cov_date_hdl_chol|cov_date_chol")),
            by = "patient_id"
            )

# Modify dummy data
## This script modifies dummy data of df_months to: 
## - randomly select 50% of cov_date_ami 
## - change them to be very close to the treatment information change, i.e., cens_date_metfin_start_cont (within 5 days before/after cens_date_metfin_start_cont)
## - this will test edge cases for fn_assign_time_fixed_cov
print('Modify dummy data')
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  message("Running locally, adapt dummy data")
  source("analysis/PPA_modify_dummy_data_df_months.R")
  message("Dummy data successfully modified")
}


# 5. Assign the stable time-updated (first-ever, constant) covariates --------
print('Assign the stable time-updated (first-ever, constant) covariates')
## Definition: Assign covariate status for interval k using latest covariate information before treatment status in interval k. 
## If no covariate status update before treatment status in interval k or no treatment status update, then carry over the latest covariate information from interval k−1.
## Rationale: This ensures that covariate information temporally precedes treatment assignment and in general reflects the latest available information at the start of interval k.

### Use a function with the following rules:
# First step:
# a) if covariate date (e.g. cov_date_ami) is not NA and happened before the minimum start date of all intervals of a person, then assign 1 (in corresponding flag variable) to all person-intervals
# b) if covariate date is not NA and happened after the maximum end date of all intervals of a person, then assign 0 (in corresponding flag variable) to all person-intervals
# c) if covariate date is NA, then assign 0 (in corresponding flag variable) to all person-intervals (assuming no documentation = no event)
# d) if covariate date is a date that happened during follow-up, then assign 1 (in flag variable) to next person-interval (using lag), 0 to all person-intervals before (incl. the current person-interval), and 1 to all person-intervals after (stable/time-fixed event)
## CAVE: This will now overwrite the preexisting cov_bin_ variables from the main dataset, which is fine (see comment above), but just to be aware of
# Second step:
# e) if covariate date is a date that happened during follow-up, check if this date is in same interval as the censoring event (cens_date_metfin_start_cont, i.e. started treatment in control group) and 
# if it happened on or before cens_date_metfin_start_cont, then assign 1 (in cov flag variable) to the current (not only next) person-interval. Otherwise, keep all as it is.

date_stable_tu_vars <- c("cov_date_ami",
                      "cov_date_hypertension",
                      "cov_date_all_stroke", "cov_date_other_arterial_embolism",
                      "cov_date_vte", "cov_date_hf", "cov_date_angina", "cov_date_dementia",
                      "cov_date_cancer", "cov_date_depression",
                      "cov_date_copd", "cov_date_liver_disease", "cov_date_chronic_kidney_disease",
                      "cov_date_pcos", "cov_date_prediabetes","cov_date_diabetescomp"
                      )
df_months <- fn_assign_stable_tu_cov(df_months, 
                                     date_stable_tu_vars,
                                     start_var = "start_date_month", 
                                     end_var = "end_date_month",
                                     cens_date_var = "cens_date_metfin_start_cont")

# To double-check
# df_months %>%
#   dplyr::select(patient_id, elig_date_t2dm, landmark_date, month, start_date_month, end_date_month,
#                 exp_bin_treat, cens_date_metfin_start_cont, cens_bin_metfin_start_cont, treat,
#                 cov_date_ami, cov_bin_ami,
#                 cov_date_hypertension, cov_bin_hypertension,
#                 cov_date_all_stroke, cov_bin_all_stroke,
#                 # cov_date, min_start, max_end, event_here, cov_flag_temp
#   ) %>%
#   dplyr::filter(!is.na(cov_date_ami)) %>%
#   dplyr::filter(!is.na(cens_date_metfin_start_cont)) %>%
#   View()


# 6. Reformatting of the dynamic time-updated covariates in a long-format df_ppa ---------
print('Reformatting of the dynamic time-updated covariates in a long-format df_ppa')
# The aim is to get a one-person-per-row dataset that we can merge to the main dataset
# So first, we work with df_ppa separately, merging to df_months only at the end; this avoids unnecessary creating/working with the main df_months as a long-format dataset (computational issues)

# Filter all dynamic time-updated covariates: hba1c/lipids/bmi
df_ppa <- df_ppa %>%
  select(patient_id, matches("_bmi_|_hba1c_|_chol_|_hdl_chol_"))
df_treat <- df_months %>% 
  filter(month == 0) %>% 
  select(patient_id, cens_date_metfin_start_cont)
# merge the cens_date_metfin_start_cont from df_months to adjust the dummy dates
# if there are patient_id that are in df_ppa but not in df_month (df_treat), then ignore these. Should not be the case in the real data, but is the case in the dummy data due to how the @table_from_file function works
df_ppa <- right_join(df_ppa, df_treat[, c("cens_date_metfin_start_cont", "patient_id")], by = "patient_id")

# Reshape to long format, easier to work with
df_ppa_long <- df_ppa %>%
  pivot_longer(
    cols = c(-patient_id, -cens_date_metfin_start_cont),
    names_to = c(".value", "variable", "instance"),
    names_pattern = "cov_(date|num)_([a-z0-9_]+)_(\\d+)",
    values_drop_na = TRUE
  )
df_ppa_long <- df_ppa_long %>% filter(!is.na(num) & !is.na(date))
df_ppa_long <- df_ppa_long %>% arrange(patient_id, variable, date)

# Modify dummy data
## This script modifies dummy data of df_ppa_long to: 
## - randomly select 5% of "date" (which is the date variable for all measurements in df_ppa_long)
## - change them to be very close to cens_date_metfin_start_cont (within 3 days before/after cens_date_metfin_start_cont)
## - this will test edge cases
print('Modify dummy data')
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  message("Running locally, adapt dummy data")
  source("analysis/PPA_modify_dummy_data_df_ppa_long.R")
  message("Dummy data successfully modified")
}
df_ppa_long <- df_ppa_long %>%
  select(!cens_date_metfin_start_cont)

# Remove illogical values, see protocol (and data_process.R, same rules/code applied there for baseline)
df_ppa_long <- df_ppa_long %>%
  mutate(
    num = case_when(
      variable == "chol" & (num > 20 | num < 1.75) ~ NA_real_,
      variable == "hdl_chol" & (num > 5  | num < 0.4) ~ NA_real_,
      variable == "hba1c" & (num > 120 | num <= 0) ~ NA_real_,
      variable == "bmi" & (num > 70 | num < 12) ~ NA_real_,
      TRUE ~ num)) %>%
  filter(!is.na(num))

# Perform the non-equi join with df_months to attach the num values from df_ppa to the correct monthly intervals from df_months, while still working in a reduced dataset
df_ppa_long <- df_ppa_long %>%
  # left_join to keep all rows=months from df_months => explicitly acknowledge many-to-many relationship to silence warning
  left_join(df_months, by = "patient_id", relationship = "many-to-many") %>%
  filter(
    date >= start_date_month, # left open, right open is correct the way the interval data is set up coming from fn_expand_intervals.R
    date <= end_date_month
  ) %>%
  select(
    patient_id,
    month,
    start_date_month,
    end_date_month,
    variable,
    date,
    num,
    cens_date_metfin_start_cont
  )
df_ppa_long <- df_ppa_long %>% arrange(patient_id, variable, date)

# Now, if there are several of the same dynamic time-updated covariate events happening in the same interval/month (e.g. HbA1c several times recorded/measured in same month):
## => carefully choose the relevant one, according to the following rules, taking our treatment information change/update variable (cens_date_metfin_start_cont) into account: 
### a) If cov_date is in same interval as cens_date_metfin_start_cont and cov_date > cens_date_metfin_start_cont, then flag the one closest to end_date_month to be moved to the next month (and discard the others that are in same interval as cens_date_metfin_start_cont & with cov_date > cens_date_metfin_start_cont)
### b) If cov_date is in same interval as cens_date_metfin_start_cont and cov_date <= cens_date_metfin_start_cont, then flag the one on or closest to cens_date_metfin_start_cont to keep (and discard the others that are also in same interval as cens_date_metfin_start_cont & cov_date <= cens_date_metfin_start_cont)
#### If at the same time, there were also cov_date > cens_date_metfin_start_cont in same interval, they will have been dealt with in rule (a) already
### c) If cov_date is not in same interval as cens_date_metfin_start_cont, then flag the one closest to end_date_month to be moved to the next month (and discard all others)
### d) Then, shift all covariate info flagged as "move" to the next month, while keeping all flagged as "keep" in the original month

##### Then, deal with edge cases: 
### (i) We may have months whereby we moved a covariate due to rule (c) or rule (a) into a month with a "keep" covariate (rule b) => >1 covariate info update in same person-interval
#### => Solution: "keep" beats "move", i.e., the one that was flagged as "keep" (because it contains time-update covariate info JUST before treatment change) is more important and wins
##### Side-note: In our case, we only have treatment update once, so a move due to rule (a) ending up in a month with a "keep" covariate (rule b) is not possible - but the function/rules work with more complex dynamic treatment updates, too!
### (ii) We may have the scenario whereby treatment update info and covariate update info change on same date (and no time-stamp to differentiate)
#### => Solution: Regard these covariate info as measured BEFORE treatment info change (see rule b above). Alternatively (sensitivity analyses), simply adapt rule a and b to regard them as measured AFTER treatment info change

##### NOTE: This function defines what to do with all OBSERVED covariate info (and to select to correct one in the correct interval). For all UNOBSERVED covariate in between we use "fill forward" (see below -> filling up the Swiss cheese), but only after merging to the existing monthly-interval dataset.
##### At the end, this results in only 1 covariate info update per person-interval for each measurement (BMI, lipids, HbA1c)
df_ppa_long <- fn_assign_dynamic_tu_cov(data = df_ppa_long,
                                         patient_id_col = patient_id,
                                         variable_col = variable,
                                         date_col = date,
                                         treat_date_col = cens_date_metfin_start_cont,
                                         month_col = month,
                                         start_date_col = start_date_month,
                                         end_date_col = end_date_month)

# To double-check
dupl_same_id_variable_month <- df_ppa_long %>%
  group_by(patient_id, variable, new_month) %>%
  add_count(name = "count") %>%
  filter(count > 1) %>%
  arrange(patient_id, variable, new_month)
if (nrow(dupl_same_id_variable_month) > 0) {
  warning("Error: Duplicate rows detected in df_ppa_long, i.e. several covariate info updates per person-interval within one variable category (BMI, lipids, or HbA1c)")
} else {
  message("Success: No duplicates found in df_ppa_long")
}

# Now, pivot it to wide format to have the individual measurements as separate columns to be able to merge with df_months
df_ppa_wide <- df_ppa_long %>%
  group_by(patient_id, new_month, variable) %>%
  mutate(instance = row_number()) %>%
  ungroup() %>%
  # The instance column ensures that each measurement gets its own row (i.e. if there are several measurements of the same measure in the same month)
  # We need to remove start_date_month, end_date_month, otherwise it creates two rows (because the dates can now be different, we rely solely on new_month!)
  pivot_wider(
    id_cols = c(patient_id, new_month, cens_date_metfin_start_cont, instance),
    names_from = variable,
    values_from = c(date, num),
    names_glue = "cov_{.value}_{variable}"
  ) %>%
  select(-instance) %>%
  arrange(patient_id, new_month)

duplicates_in_final <- df_ppa_wide %>%
  group_by(patient_id, new_month) %>%
  add_count(name = "count") %>%
  filter(count > 1) %>%
  ungroup()
if (nrow(duplicates_in_final) > 0) {
  warning("Error: Duplicate rows detected in df_ppa_wide, i.e. several covariate info updates per person-interval due to either BMI, lipids, or HbA1c")
} else {
  message("Success: No duplicates found in df_ppa_wide")
}


# 7. Merge and assign the dynamic time-updated covariates from df_ppa to df_months --------
print('Merge and assign the dynamic time-updated covariates covariates from df_ppa to df_months')
df_ppa_wide <- df_ppa_wide %>%
  rename(month = new_month)
df_ppa_wide <- df_ppa_wide %>%
  select(!cens_date_metfin_start_cont) # take out the cens_date (otherwise, duplicate variables)

# CAVE: due to new_month, it could be that we reach beyond max follow-up => ignore such rows for merging
df_months <- df_months %>%
  # Keep all df_months months; if there is a month in df_ppa_wide that is not in df_months do not include it. 
  left_join(df_ppa_wide, by = c("patient_id", "month"))


# 8. Fill up the Swiss cheese and reassign the lipid ratio --------
print('Fill up the Swiss cheese and reassign the lipid ratio')

# fill forward, always taking the latest value into account
# if the first interval (month 0) is missing, then take that value from a separate corresponding baseline variable (e.g. use the info from the time-fixed cov_num_bmi_b to fill in the time-updated cov_num_bmi at baseline and then take it from there)
tv_cov_num_cols <- c("cov_num_bmi", "cov_num_hba1c", "cov_num_chol", "cov_num_hdl_chol")
baseline_lookup <- c(
  cov_num_bmi = "cov_num_bmi_b",
  cov_num_hba1c = "cov_num_hba1c_b",
  cov_num_chol = "cov_num_chol_b",
  cov_num_hdl_chol = "cov_num_hdl_chol_b"
)

df_months <- df_months %>%
  group_by(patient_id) %>%
  arrange(start_date_month, .by_group = TRUE) %>%
  mutate(across(
    all_of(tv_cov_num_cols),
    ~ fn_fill_forward_incl_baseline(.x, baseline_lookup[[cur_column()]])
  )) %>%
  ungroup()

# Recreate the lipid ratio, HbA1c cat, and BMI cat (based on the monthly filled up values) - but without data cleaning, this happened at an earlier step.
# All data formatting was done on the underlying numeric value. We may loose some information by using the lipid ratio instead of using HDL and Tot Chol separately, but it makes clinically more sense, and when HDL is measured also Tot Chol is measured, so both should be available or not.
df_months <- df_months %>%
  mutate(
    # See data_process.R
    cov_num_tc_hdl_ratio = cov_num_chol / cov_num_hdl_chol,
    cov_cat_tc_hdl_ratio = cut(
      cov_num_tc_hdl_ratio,
      breaks = c(1, 3.5, 5.11, 50), # should not have any illogical values since underlying values have been cleaned (below 1, above 50), but in case it does, it replaces it with "Unknown"
      labels = c("below 3.5:1" ,"3.5:1 to 5:1", "above 5:1"),
      right = FALSE) %>% 
      forcats::fct_expand("Unknown"),
    cov_cat_tc_hdl_ratio = case_when(is.na(cov_cat_tc_hdl_ratio) ~ factor("Unknown", 
                                                                          levels = c("below 3.5:1", "3.5:1 to 5:1", "above 5:1", "Unknown")), TRUE ~ cov_cat_tc_hdl_ratio)
  ) %>% 
  mutate(
    # See data_process.R
    cov_cat_hba1c = cut(
      cov_num_hba1c,
      breaks = c(0, 42, 59, 76, 120), # should not have any illogical values since underlying values have been cleaned (below 0, above 120), but in case it does, it replaces it with "Unknown"
      labels = c("below 42" ,"42-58", "59-75", "above 75"),
      right = FALSE) %>% 
      forcats::fct_expand("Unknown"),
    cov_cat_hba1c = case_when(is.na(cov_cat_hba1c) ~ factor("Unknown", 
                                                            levels = c("below 42", "42-58", "59-75", "above 75", "Unknown")), TRUE ~ cov_cat_hba1c)
  ) %>% 
  mutate(
    cov_cat_bmi = cut(
      cov_num_bmi,
      breaks = c(12, 30, 70), # should not have any illogical values since underlying values have been cleaned, but in case it does (below 12, above 70), it replaces it with "Unknown"
      labels = c("Not obese", "Obese"),
      right = FALSE) %>% 
      forcats::fct_expand("Unknown"),
    cov_cat_bmi = case_when(is.na(cov_cat_bmi) ~ factor("Unknown", 
                                                        levels = c("Not obese", "Obese", "Unknown")), TRUE ~ cov_cat_bmi)
  )

# To double-check
# df_months %>%
#   dplyr::select(patient_id, elig_date_t2dm, landmark_date, month, start_date_month, end_date_month,
#                 cov_date_chol, cov_num_chol, cov_date_hdl_chol, cov_num_hdl_chol, cov_num_chol_b, cov_num_hdl_chol_b,
#                 cov_num_tc_hdl_ratio_b,
#                 cov_cat_tc_hdl_ratio_b,
#                 cov_num_tc_hdl_ratio,
#                 cov_cat_tc_hdl_ratio,
#                 cov_date_hba1c, cov_num_hba1c,
#                 cov_cat_hba1c,
#                 cov_num_hba1c_b, cov_cat_hba1c_b,
#                 cov_date_bmi, cov_num_bmi,
#                 cov_num_bmi_b,
#                 cov_cat_bmi_groups, cov_bin_obesity,
#                 cov_cat_bmi,
#                 cov_bin_ami, cov_bin_hypertension, cov_bin_all_stroke,
#                 out_date_severecovid_afterlandmark, out_date_death_afterlandmark,
#                 cens_date_ltfu_afterlandmark,
#                 cov_cat_sex, cov_num_age
#   ) %>%
#   # filter(is.na(cov_cat_tc_hdl_ratio)) %>%
#   View()


# 9. Drop unnecessary helper variables created along the way -------
print('Drop unnecessary helper variables created along the way')
# drop the baseline helper variables
df_months <- df_months %>%
  select(-c(cov_num_hdl_chol_b, cov_num_chol_b,
            cov_num_tc_hdl_ratio_b, cov_cat_tc_hdl_ratio_b,
            cov_num_hba1c_b, cov_cat_hba1c_b,
            cov_num_bmi_b, cov_cat_bmi_groups, cov_bin_obesity))
# drop all unnecessary covariate date variables used for assigning events to correct intervals
df_months <- df_months %>%
  select(-starts_with("cov_date"))
# drop all unnecessary bin flags of outcomes and censoring events
df_months <- df_months %>%
  select(-c(contains("out_bin"), contains("cens_bin")))
# drop outcomes we are not interested in the PPA
df_months <- df_months %>%        
  select(-c(contains("fracture"), contains("dm_death")))

# 10. Assign and shift outcome, censoring due to LTFU and competing events --------
print('Assign and shift outcome, censoring due to LTFU, and competing events')
# CAVE: use it only after the dataset has been interval-expanded (fn_expand_intervals.R) and after all events have been added (fn_assign* functions)

# Define columns called "outcome", censor_LTFU", "comp_event"
# Shift them as follows:

## a) outcome is happening in CURRENT (k + 1) interval:
### i) assign to CURRENT (k + 1) row/person-interval: outcome=NA, censor_LTFU=NA, comp_event=NA
### ii) assign to the PREVIOUS (k) row/person-interval: outcome=1, censor_LTFU=0, comp_event=0

## b) censor_LTFU is happening in CURRENT (k + 1) interval:
### i) assign to CURRENT (k + 1) row/person-interval: outcome=NA, censor_LTFU=NA, comp_event=NA
### ii) assign to the PREVIOUS (k) row/person-interval: outcome=NA, censor_LTFU=1, comp_event=NA # I still think it is possible to replace all NA with 0 without making any difference, since the outcome model is restricted to only rows with df_months[df_months$censor==0 & df_months$comp_event == 0,], but need to check. This is how Miguel's data set is set up, maybe it will have a reason later on, e.g. when looking at competing events. TBD

## c) comp_event is happening in CURRENT (k + 1) interval:
### i) assign to CURRENT (k + 1) row/person-interval: outcome=NA, censor_LTFU=NA, comp_event=NA
### ii) assign to the PREVIOUS (k) row/person-interval : outcome=NA, censor_LTFU=0, comp_event=1 # same comment re NA as above

## Edge case RULES
### d) if outcome and censor_LTFU have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order in function.
### e) if outcome and competing event (i.e. non-covid death) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order in function.
### f) if censoring and competing event event have the EXACT same date, then pick the competing event as the defining event. Ensured by case_when() order in function.

## Then, delete all rows/person-intervals with NAs in all three variables (outcome & censor_LTFU & comp_event), i.e., the person-interval the event happened (and by doing this, this also drops the covariate/eligibility/treatment info from that interval)
## This is the reason, why it needs a dataset per outcome, because the datasets will have different lengths. Alternatively tell the model to ignore info happening after outcome of interest (even if there is data afterwards), but I think this is cleaner, easier and safer.

# Currently, only calendar month and calendar week are implemented

# Primary outcome dataset: stops either at severe covid (includes death and hosp), non-covid death, or deregistration
df_months_severecovid <- fn_add_and_shift_out_comp_cens_events(df_months,
                                                               outcome_date_variable = "out_date_severecovid_afterlandmark",
                                                               comp_date_variable = "out_date_noncoviddeath_afterlandmark",
                                                               censor_date_variable = "cens_date_ltfu_afterlandmark",
                                                               studyend_date,
                                                               interval_type = "month")

# Secondary outcome dataset: stops either at severe covid (includes covid deaths, hosp, diagnoses, tests), non-covid death, or deregistration
df_months_covid_event <- fn_add_and_shift_out_comp_cens_events(df_months,
                                                               outcome_date_variable = "out_date_covid",
                                                               comp_date_variable = "out_date_noncoviddeath_afterlandmark",
                                                               censor_date_variable = "cens_date_ltfu_afterlandmark",
                                                               studyend_date,
                                                               interval_type = "month")

# Secondary outcome dataset: stops either at Long covid diagnosis (includes viral fatigue), any death, or deregistration
df_months_longcovid <- fn_add_and_shift_out_comp_cens_events(df_months,
                                                             outcome_date_variable = "out_date_longcovid_virfat",
                                                             comp_date_variable = "out_date_death_afterlandmark", # all-cause mortality as competing event for this analysis
                                                             censor_date_variable = "cens_date_ltfu_afterlandmark",
                                                             studyend_date,
                                                             interval_type = "month")

# To double-check
# df_months_severecovid %>%
#   dplyr::select(patient_id, elig_date_t2dm, start_date_month, end_date_month, month,
#                 outcome, comp_event, censor_LTFU, stop_date,
#                 # event_date, is_event_interval, is_outcome_event, is_comp_event, is_cens_event, row_num, first_event_row,
#                 exp_bin_treat,
#                 out_date_severecovid_afterlandmark,
#                 out_date_noncoviddeath_afterlandmark,
#                 cens_date_ltfu_afterlandmark,
#                 out_date_death_afterlandmark,
#                 out_date_covid_afterlandmark,
#                 cens_date_metfin_start_cont,
#                 out_date_longcovid_virfat_afterlandmark,
#                 cov_num_tc_hdl_ratio, cov_cat_tc_hdl_ratio,
#                 cov_num_chol, cov_num_hdl_chol,
#                 cov_num_hba1c,
#                 cov_num_bmi,
#                 cov_cat_sex, cov_num_age,
#                 cov_bin_ami, cov_bin_hypertension, cov_bin_all_stroke
#                 ) %>%
#   # dplyr::filter(!is.na(out_date_severecovid_afterlandmark)) %>%
#   # dplyr::filter(!is.na(out_date_noncoviddeath_afterlandmark)) %>%
#   # dplyr::filter(!is.na(cens_date_ltfu_afterlandmark)) %>%
#   # dplyr::filter(!is.na(cens_date_metfin_start_cont)) %>%
#   # dplyr::filter(is.na(censor_LTFU)) %>% # should be empty
#   View()


# 10. Define eligibility --------
print('Define eligibility')
# NA in our case, everyone is eligible in this dataset, single-trial approach not sequential trial


# 11. Save output -------------------------------------------------------------
print('Save output')
arrow::write_feather(df_months_severecovid, here::here("output", "data", "df_months_severecovid.arrow"))
arrow::write_feather(df_months_covid_event, here::here("output", "data", "df_months_covid_event.arrow"))
arrow::write_feather(df_months_longcovid, here::here("output", "data", "df_months_longcovid.arrow"))
