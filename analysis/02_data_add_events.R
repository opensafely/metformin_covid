####
## This script does the following:
# 1. Import multiple-event-per-person long format dataset, with time intervals as per study protocol
# 2. Import the additional time-updated covariate dataset (data_processed_ppa.arrow)
# 3. Assign treatment
# 4. Assign first-ever, time-fixed covariates (shift from k to k+1)
# 5. Assign time-updated covariates (shift from k to k+1)
# 6. Assign first-ever, time-fixed outcomes/competing/censoring events (shift from k to k-1)
# 7. Define eligibility (and censoring?)
# 8. Save datasets
####


# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(zoo) # for fn_add_firstever_events_to_intervals.R
source(here::here("analysis", "functions", "fn_assign_treatment.R"))
source(here::here("analysis", "functions", "fn_assign_time_fixed_cov.R"))
source(here::here("analysis", "functions", "fn_add_and_shift_out_comp_cens_events.R"))


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data"))


# Import dates ------------------------------------------------------------
print('Import dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# Import the data ---------------------------------------------------------
print('Import the data')
df_months <- read_feather(here("output", "data", "df_months.arrow"))
df_ppa <- read_feather(here("output", "data", "data_processed_ppa.arrow")) 


# Merge (some) data from the PPA dataset ----------------------------------
# merge all the cov_date* that occur once and then remain stable, i.e., all exceot hba1c/lipids/bmi, these will be added later
df_months <- df_months %>%
  left_join(df_ppa %>% 
              select(patient_id, starts_with("cov_date") & !matches("cov_date_bmi|cov_date_hba1c|cov_date_hdl_chol|cov_date_chol")),
            by = "patient_id"
            )


# Assign treatment --------------------------------------------------------
print('Assign treatment')
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
df_months <- fn_assign_treatment(df_months, date_treat_vars,
                                 start_var = "start_date_month", 
                                 end_var = "end_date_month")
# To double-check
# df_months %>%
#   dplyr::select(patient_id, elig_date_t2dm, landmark_date, month, start_date_month, end_date_month, 
#                 exp_bin_treat, cens_date_metfin_start_cont, cens_bin_metfin_start_cont,
#                 out_date_severecovid_afterlandmark, out_date_death_afterlandmark,
#                 cens_date_ltfu_afterlandmark, 
#                 cov_cat_sex, cov_num_age
#   ) %>%
#   View()

# Assign first-ever, time-fixed covariates -------------------------------
print('Assign first-ever, time-fixed covariates')
## Definition: Assign covariate status for interval k using covariate information from interval k−1.
## Rationale: This lag ensures that covariate values reflect the latest available information at the start of each interval - prior to eligibility or treatment assignment in interval k, avoiding looking into the future within an interval.

### First, add all first-ever time-fixed covariates, i.e., those occurring once and remaining stable (these are usually the highest number of variables in a dataset)
### I can on my preexisting function, based on these rules, but in addition implement the lag:
# a) if event date is not NA and happened before the minimum start date of all intervals of a person, then assign 1 (in corresponding flag variable) to all person-intervals
# b) if event date is not NA and happened after the maximum end date of all intervals of a person, then assign 0 (in corresponding flag variable) to all person-intervals
# c) if no event date is recorded (date variable == NA), then assign 0 (in corresponding flag variable) to all person-intervals (assuming no documentation = no event)
# d) if event date happened during follow-up, then assign 1 (in flag variable) to next person-interval (lag), and 0 to all person-intervals before, and 1 to all person-intervals after (stable/time-fixed event)
## CAVE: This will overwrite the preexisting cov_bin_ variables from the main dataset, which is fine, but just to be aware of
date_tf_cov_vars <- c("cov_date_hypertension", "cov_date_ami", 
                      "cov_date_all_stroke", "cov_date_other_arterial_embolism",
                      "cov_date_vte", "cov_date_hf", "cov_date_angina", "cov_date_dementia", 
                      "cov_date_cancer", "cov_date_depression",
                      "cov_date_copd", "cov_date_liver_disease", "cov_date_chronic_kidney_disease",
                      "cov_date_pcos", "cov_date_prediabetes","cov_date_diabetescomp")
df_months <- fn_assign_time_fixed_cov(df_months, date_tf_cov_vars,
                                      start_var = "start_date_month", 
                                      end_var = "end_date_month")
# # To double-check
# df_months %>%
#   dplyr::select(patient_id, landmark_date, month, start_date_month, end_date_month,
#                 cov_date_ami, cov_bin_ami,
#                 cov_date_hypertension, cov_bin_hypertension,
#                 cov_date_all_stroke, cov_bin_all_stroke,
#                 elig_date, event_here, 
#                 # event_cum
#   ) %>%
#   # dplyr::filter(cov_bin_hypertension == 1) %>%
#   # dplyr::filter(!is.na(cov_date_hypertension)) %>%
#   dplyr::filter(!is.na(cov_date_ami)) %>%
#   View()

date_vars <- c("out_date_severecovid_afterlandmark",
               "out_date_noncoviddeath_afterlandmark",
               "cens_date_ltfu_afterlandmark",
               "cens_date_metfin_start_cont", 
               "out_date_covid_afterlandmark",
               "out_date_longcovid_virfat_afterlandmark",
               "out_date_death_afterlandmark",
               
               "cov_date_ami", "cov_date_all_stroke", "cov_date_other_arterial_embolism",
               "cov_date_vte", "cov_date_hf", "cov_date_angina", "cov_date_dementia", 
               "cov_date_cancer", "cov_date_hypertension", "cov_date_depression",
               "cov_date_copd", "cov_date_liver_disease", "cov_date_chronic_kidney_disease",
               "cov_date_pcos", "cov_date_prediabetes","cov_date_diabetescomp")

# take baseline covariate info (binary variables) into account?!!!



# (2) Assign and shift eligibility, treatment, and censoring due to treatment ------
print('Assign and shift treatment and censoring due to treatment')
# Side note: In this simple scenario we assigned eligibility and treatment (exp_bin_treat) at baseline/landmark, no need to do anything (no sequential trial)
# We only need to assign and shift the event of censoring due to metformin start in control group, which we want to tackle in the per-protocol analysis
# cens_bin_metfin_start_cont comes from fn_add_firstever_events_to_intervals, so either is 0 or 1
df_months <- df_months %>%
  group_by(patient_id) %>%
  arrange(patient_id, month) %>%
  mutate(
    censor_ppa = as.integer(cumany(lag(cens_bin_metfin_start_cont, default = 0) == 1))
  ) %>%
  ungroup()
# Note: Once treatment is started, it stays TRUE (=1) for all future rows/person-intervals (cumulative OR)
# For the PPA model, we will then have to restrict the rows to after the first row with censor_ppa = 1 ! -> TBD

# To double-check
# df_months %>%
#   dplyr::select(patient_id, elig_date_t2dm, start_date_month, end_date_month, month, censor_ppa, cens_bin_metfin_start_cont, stop_date,
#                 out_date_severecovid_afterlandmark, out_bin_severecovid_afterlandmark,
#                 out_date_noncoviddeath_afterlandmark, out_bin_noncoviddeath_afterlandmark,
#                 cens_date_ltfu_afterlandmark, cens_bin_ltfu_afterlandmark,
#                 out_date_death_afterlandmark, out_bin_death_afterlandmark,
#                 out_date_covid_afterlandmark, out_bin_covid_afterlandmark,
#                 cens_date_metfin_start_cont, cens_bin_metfin_start_cont,
#                 out_date_longcovid_virfat_afterlandmark, out_bin_longcovid_virfat_afterlandmark,
#                 cov_cat_sex, cov_num_age) %>%
#   View()


# (3) Assign and shift outcome, censoring due to LTFU and competing events --------
print('Assign and shift outcome, censoring due to LTFU and competing events')
# Define columns called "outcome", censor_LTFU", "comp_event"
# Shift them up as discussed above and assign info to CURRENT (k + 1) and PREVIOUS (k) interval as follows:

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

## Now, delete all rows/person-intervals with NAs in all three variables (outcome & censor_LTFU & comp_event), i.e., the person-interval the event happened (and with this also drop the covariate/eligibility/treatment info from that interval)
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
#                 censor_ppa,
#                 out_date_severecovid_afterlandmark, out_bin_severecovid_afterlandmark,
#                 out_date_noncoviddeath_afterlandmark, out_bin_noncoviddeath_afterlandmark,
#                 cens_date_ltfu_afterlandmark, cens_bin_ltfu_afterlandmark,
#                 out_date_death_afterlandmark, out_bin_death_afterlandmark,
#                 out_date_covid_afterlandmark, out_bin_covid_afterlandmark,
#                 cens_date_metfin_start_cont, cens_bin_metfin_start_cont,
#                 censor_ppa,
#                 out_date_longcovid_virfat_afterlandmark, out_bin_longcovid_virfat_afterlandmark,
#                 cov_cat_sex, cov_num_age) %>%
#   # dplyr::filter(!is.na(out_date_severecovid_afterlandmark)) %>%
#   # dplyr::filter(!is.na(out_date_noncoviddeath_afterlandmark)) %>%
#   # dplyr::filter(!is.na(cens_date_ltfu_afterlandmark)) %>%
#   # dplyr::filter(!is.na(cens_date_metfin_start_cont)) %>%
#   # dplyr::filter(is.na(censor_LTFU)) %>% # should be empty
#   View()


# Save output -------------------------------------------------------------
print('Save output')
arrow::write_feather(df_months_severecovid, here::here("output", "data", "df_months_severecovid.arrow"))
arrow::write_feather(df_months_covid_event, here::here("output", "data", "df_months_covid_event.arrow"))
arrow::write_feather(df_months_longcovid, here::here("output", "data", "df_months_longcovid.arrow"))
