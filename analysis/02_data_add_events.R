####
## This script does the following:
# 1. Import multiple-event-per-person long format dataset, with time intervals as per study protocol
# 2. Assign first-ever, time-fixed events to intervals (outcomes, competing events, censoring events, eligibility events, time-fixed covariates), that occur once and then remain stable
# 3. Assign and shift eligibility, treatment, and censoring due to treatment
# 4. Assign and shift outcome, censoring due to LTFU and competing events and create outcome-specific datasets
# 5. Save datasets
####


# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(zoo) # for fn_add_firstever_events_to_intervals.R
source(here::here("analysis", "functions", "fn_add_firstever_events_to_intervals.R"))
source(here::here("analysis", "functions", "fn_add_and_shift_out_comp_cens_events.R"))


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data"))


# Import the data ---------------------------------------------------------
print('Import the data')
df_months <- read_feather(here("output", "data", "df_months.arrow"))


# Import dates ------------------------------------------------------------
print('Import dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# RULES to assign all events to the correct intervals ---------------------
## Always think: covariate info BEFORE treatment assignment BEFORE outcome
# RULE 1: Assign the time-fixed & time-varying covariate info to the current interval (these are usually the highest number of variables in a dataset)
# RULE 2: Assign the treatment and eligibility info to the next interval (using lagged info, so k - 1 info to k) 
## This ensures we are not looking into the future to define eligibility and treatment, 
## and we do ensure covariate info in same interval as treatment assignment is before treatment assignment and taken into account as such when calculating the treatment weights (model with y = treat == 1 ~ ...) and covariate info from k interval is seen as after treatment start
# RULE 3: For outcomes and competing events, shift info 1 interval up (i.e. info from k + 1 in interval k)
## This avoids using covariate/treatment info that might happen after the outcome in same interval. 
## (side note: This is not possible with outcomes/competing events such as death and we would not have to do this for these, but in many other outcome scenarios this is possible, e.g. long covid, let's adhere to it as a general rule). 
# RULE 4: For censoring events, we need to distinguish between different censoring events:
## If the censoring event is something like LTFU/deregistration, then use RULE 3 to be on the safe side (even though technically not necessarily needed)
## If the censoring event is something like starting metformin in the control group as a censoring event for the per-protocol analyses, then use RULE 2

### Does this shifting up and down interfere with each other?
# 1. Shifting the outcomes/competing events and LTFU-censoring events up 1 person-interval will eliminate the covariate info from the current person-interval. That's ok and desirable, so that the models with y = censor_ltfu == 0 ~ ... calculate on correct covariate info
# 2. Shifting the treatment and treatment-censoring events down 1 person-interval ensures that covariate info in same interval as treatment assignment is before treatment assignment and taken into account as such when calculating the treatment/censoring weights (model with y = treat == 1 ~ ... or y = censor_treat == 0 ~ ...) and covariate info from k interval is seen as after treatment start


# (1) Assign first-ever, time-fixed events to intervals -------------------
print('Assign first-ever, time-fixed events to intervals')
## These are all "first-ever" events, that happen once, and then remain TRUE during the entire follow-up
## This is true for the majority of events. 
## RULES:
# a) if event date is not NA and happened before the minimum start date of all intervals of a person, then assign 1 (in corresponding flag variable) to all person-intervals
# b) if event date is not NA and happened after the maximum end date of all intervals of a person, then assign 0 (in corresponding flag variable) to all person-intervals
# c) if no event date is recorded (date variable == NA), then assign 0 (in corresponding flag variable) to all person-intervals (assuming no documentation = no event)
# d) if event date happened during follow-up, then assign 1 (in corresponding flag variable) to corresponding person-interval, and 0 to all person-intervals before, and 1 to all person-intervals after (stable/time-fixed event)
# Side note: Deal with shifting information for outcome/censoring/competing event in a next step
# Side note: This function can also be used to add time-fixed eligibility events (e.g. T2DM diagnosis or metformin allergy), but this is only relevant for a sequential trial setup. Here, we defined the eligible population once, at baseline.
# Side note: I will also add time-fixed covariates for the per-protocol analysis to this step, but at a later stage (need to extract them first) -> make sure to take baseline covariate info (binary variables) into account when adding.
# We defined all our outcomes of interest in previous script, see 01_data_add_intervals
date_vars <- c("out_date_severecovid_afterlandmark",
               "out_date_noncoviddeath_afterlandmark",
               "cens_date_ltfu_afterlandmark",
               "cens_date_metfin_start_cont", 
               "out_date_covid_afterlandmark",
               "out_date_longcovid_virfat_afterlandmark",
               "out_date_death_afterlandmark")
df_months <- fn_add_firstever_events_to_intervals(df_months, 
                                                  date_vars,
                                                  start_var = "start_date_month", 
                                                  end_var = "end_date_month")
## NOTE: The function names the flag variables by simply replacing _date_ with _bin_
## NOTE: If a variable with the same name already existis (e.g. because defined before for the Cox model), it will overwrite it.

# To double-check
# df_months %>%
#   dplyr::select(patient_id, elig_date_t2dm, start_date_month, end_date_month, month, stop_date,
#                 out_date_severecovid_afterlandmark, out_bin_severecovid_afterlandmark,
#                 out_date_noncoviddeath_afterlandmark, out_bin_noncoviddeath_afterlandmark,
#                 cens_date_ltfu_afterlandmark, cens_bin_ltfu_afterlandmark,
#                 out_date_death_afterlandmark, out_bin_death_afterlandmark,
#                 out_date_covid_afterlandmark, out_bin_covid_afterlandmark,
#                 cens_date_metfin_start_cont, cens_bin_metfin_start_cont,
#                 out_date_longcovid_virfat_afterlandmark, out_bin_longcovid_virfat_afterlandmark,
#                 cov_cat_sex, cov_num_age) %>%
#   # dplyr::filter(!is.na(out_date_severecovid_afterlandmark)) %>%
#   # dplyr::filter(!is.na(cens_date_ltfu_afterlandmark)) %>%
#   # dplyr::filter(!is.na(out_date_covid_afterlandmark)) %>%
#   # dplyr::filter(!is.na(cens_date_metfin_start_cont)) %>%
#   # dplyr::filter(!is.na(out_date_death_afterlandmark)) %>%
#   View()


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
