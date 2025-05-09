####
## This script does the following:
# 1. Import processed data
# 2. Expand it to long format with interval data (weekly or monthly)
# 3. Add outcome, competing and censoring events and shift this info to prior person-interval as per rules outlined
# 4. Add other events, i.e., "first-ever" events
# 5. Save output
####

# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(lubridate) # for fn_expand_intervals.R
library(zoo) # for fn_add_firstever_events_to_intervals.R
library(purrr) # for fn_expand_intervals.R and fn_add_firstever_events_to_intervals.R and fn_expand_intervals.R
library(rlang) # for fn_expand_intervals.R and fn_add_firstever_events_to_intervals.R
library(stringr) # for fn_add_firstever_events_to_intervals.R
source(here::here("analysis", "functions", "fn_expand_intervals.R"))
source(here::here("analysis", "functions", "fn_add_out_comp_cens_events_to_intervals.R"))
source(here::here("analysis", "functions", "fn_add_firstever_events_to_intervals.R"))

# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_description_seqtrials"))

# Import dates ------------------------------------------------------------
print('Import the dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")
pandemicstart_date <- as.Date(study_dates$pandemicstart_date, format = "%Y-%m-%d")

# Import the data ---------------------------------------------------------
print('Import the data')
df <- read_feather(here("output", "data", "data_processed_seqtrials.arrow"))

# Expand the dataset ------------------------------------------------------
print('Expand the dataset')
# expand, using monthly intervals. My interval stopping events are: 
# out_date_covid_death (one of the outcomes), out_date_noncovid_death (the main competing event), cens_date_dereg (estimand: "assume no-one is LTFU" with IPCW; as part of ITT)
# COVID-related hospitalization per se is not an interval stopping event (but an important outcome and elig criteria), will be added later.
# Other outcomes, censoring events (e.g. starting metformin in control group) and time-updated covariates will be added later.

# Apply the function
## CAVE: stop_date_columns shall only have dates in the future of start_date_variable
## CAVE: only choose stopping events that are distinct (e.g. out_date_noncovid_death is not part of out_date_covid_death)

df_months <- fn_expand_intervals(df, 
                                 start_date_variable = pandemicstart_date,
                                 stop_date_columns = c("out_date_covid_death", "out_date_noncovid_death", "cens_date_dereg"), 
                                 studyend_date,
                                 interval_type = "month")
## NOTE: the function will add the final interval, even if the final interval is just 1 day (check stop_date == 2022-04-01)


# Add outcome, censoring and competing events the dataset ------------------
print('Add outcome, censoring and competing events the dataset')
## See detailed rules in function script, in brief:
# Assign Outcome/Censor/Competing columns:
## a) If outcome is the event: 
### i) assign to the PREVIOUS interval row/person-interval: outcome=1, censor=0, comp_event=0
### ii) assign to CURRENT interval row/person-interval: outcome=NA, censor=NA, comp_event=NA

## b) If censoring event is the event: 
### i) assign to the PREVIOUS interval row/person-interval: outcome=NA, censor=1, comp_event=NA
### ii) assign to CURRENT interval row/person-interval: outcome=NA, censor=NA, comp_event=NA

## c) If competing event is the event: 
### i) assign to the PREVIOUS interval row/person-interval: outcome=NA, censor=0, comp_event=1
### ii) assign to CURRENT interval row/person-interval: outcome=NA, censor=NA, comp_event=NA

df_months_out <- fn_add_out_comp_cens_events_to_intervals(df_months, 
                                                           outcome_date_variable = "out_date_covid_death",
                                                           comp_date_variable = "out_date_noncovid_death",
                                                           censor_date_variable = "cens_date_dereg",
                                                           studyend_date,
                                                           interval_type = "month")

# df_months_out %>%
#   dplyr::select(patient_id, elig_date_t2dm, out_date_covid_death, out_date_noncovid_death, cens_date_dereg,
#                 stop_date, start_date_month, end_date_month, month, outcome, comp_event, censor, is_event_interval, is_outcome_event, dupl_interval,
#                 cov_cat_sex, cov_num_age, everything()) %>%
#   # dplyr::filter(!is.na(out_date_covid_death)) %>%
#   dplyr::filter(!is.na(out_date_noncovid_death)) %>%
#   # dplyr::filter(!is.na(cens_date_dereg)) %>%
#   # dplyr::filter(is.na(censor)) %>% # should be empty
#   # dplyr::filter(dupl_interval == TRUE) %>% # should be empty
#   View()


# Add all other date events to correct interval with a bin flag ---------
print('Add all other date events to correct interval with a bin flag')
## remember, these are all "first ever" events, that happen once and remain 1

# a) if event is not NA and happening before min start date intervals, then assign 1 (in corresponding flag variable) to all person-intervals
# b) if event date is not NA and happening after max end date intervals, then assign 0 (in corresponding flag variable) to all person-intervals
# c) if no event date is recorded (date variable == NA), then assign 0 (in corresponding flag variable) to all person-intervals
# d) if event date is happening during follow-up, then assign 1 (in corresponding flag variable) to corresponding person-interval, 0 to all person-intervals before, and 1 to all person-intervals after

date_vars <- c("elig_date_metfin_allergy_first", "elig_date_ckd_45_first", "elig_date_liver_cirrhosis_first",
               "elig_date_metfin_interaction_first", "exp_date_metfin_first", "exp_date_metfin_mono_first", "exp_date_sulfo_first",
               "exp_date_dpp4_first", "exp_date_tzd_first", "exp_date_sglt2_first", "exp_date_glp1_first", "exp_date_megli_first", "exp_date_agi_first", 
               "exp_date_insulin_first", 
               "out_date_severecovid" # here are the covid hospitalizations "hidden", but will not use this, instead will use ELD for hosp (currently added just to see)
               )

df_months_tot <- fn_add_firstever_events_to_intervals(df_months_out, 
                                                      date_vars,
                                                      start_var = "start_date_month", 
                                                      end_var = "end_date_month")

# To double-check
# df_months_tot %>%
#   dplyr::select(patient_id, elig_date_t2dm, start_date_month, end_date_month, month, out_date_severecovid, out_bin_severecovid,
#                 elig_date_ckd_45_first, elig_bin_ckd_45_first, exp_date_sulfo_first, exp_bin_sulfo_first,
#                 exp_date_insulin_first, exp_bin_insulin_first,
#                 outcome, comp_event, censor, is_event_interval, is_outcome_event,
#                 cov_cat_sex, cov_num_age) %>%
#   View()


# Drop helper variables -------------------------------------------------
print('Drop helper variables')
df_months_tot <- df_months_tot %>% 
  select(-is_event_interval, -is_outcome_event, -is_cens_event, -is_comp_event, -start_date, -stop_date, -event_date, -dupl_interval)


# Save the dataset ------------------------------------------------------
arrow::write_feather(df_months_tot, here::here("output", "data", "df_months_tot.arrow"))