####
## This script does the following:
# 1. Import pre-processed one-person-per-row dataset
# 2. Expand it to long format, creating intervals (weekly or monthly)
# 3. Define/add non-ELD (i.e. "one-time only / first-ever") events
# 4. Save output
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
# 1) Death, recorded in: a) out_date_covid_death (one of the outcomes), b) out_date_noncovid_death (the main competing event)
# 2) De-registration, i.e. lost-to-follow-up: cens_date_dereg (estimand: "assume no-one is LTFU" with IPCW; as part of ITT)

# COVID-related hospitalization per se is not an interval stopping event (but an important outcome and elig criteria), will be added later, using ELD.
# Other outcomes (e.g. COVID diagnosis), censoring events (e.g. starting metformin in control group) and time-updated covariates will be added later, using ELD.

# Apply the function
## CAVE: stop_date_columns shall only have dates AFTER start_date_variable
## CAVE: only choose stopping events that are distinct (e.g. out_date_noncovid_death is NOT part of out_date_covid_death)
df_months <- fn_expand_intervals(df, 
                                 start_date_variable = pandemicstart_date,
                                 stop_date_columns = c("out_date_covid_death", "out_date_noncovid_death", "cens_date_dereg"), 
                                 studyend_date,
                                 interval_type = "month")
## NOTE: the function will add the final interval, even if the final interval is just 1 day (check stop_date == 2022-04-01)


# Add first-ever date events to correct interval with a bin flag ----------
print('Add first-ever date events to correct interval with a bin flag')
## remember, these are all "first ever" events, that happen once and remain 1 in the follow-up

# a) if event date is not NA and happened before min start date intervals, then assign 1 (in corresponding flag variable) to all person-intervals
# b) if event date is not NA and happened after max end date intervals, then assign 0 (in corresponding flag variable) to all person-intervals
# c) if no event date is recorded (date variable == NA), then assign 0 (in corresponding flag variable) to all person-intervals
# d) if event date happened during follow-up (i.e. existing interval date), then assign 1 (in corresponding flag variable) to corresponding person-interval, and 0 to all person-intervals before, and 1 to all person-intervals after

date_vars <- c("elig_date_metfin_allergy_first", "elig_date_ckd_45_first", "elig_date_liver_cirrhosis_first",
               "elig_date_metfin_interaction_first", 
               "exp_date_metfin_first", "exp_date_metfin_mono_first", "exp_date_sulfo_first",
               "exp_date_dpp4_first", "exp_date_tzd_first", "exp_date_sglt2_first", "exp_date_glp1_first", "exp_date_megli_first", "exp_date_agi_first", 
               "exp_date_insulin_first",
               # add outcomes, comp and censoring events, too. To use them for time-updated eligibility variable
               "out_date_severecovid", "out_date_covid", "out_date_longcovid_virfat", 
               "out_date_covid_death", "out_date_noncovid_death", "cens_date_dereg", 
               "out_date_longcovid" # this is needed separately, see elig variable
               )

df_months_tot <- fn_add_firstever_events_to_intervals(df_months, 
                                                      date_vars,
                                                      start_var = "start_date_month", 
                                                      end_var = "end_date_month")
# df_months_tot %>%
#   dplyr::select(patient_id, elig_date_t2dm, start_date_month, end_date_month, month, stop_date,
#                 out_date_covid_death, out_bin_covid_death,
#                 out_date_noncovid_death, out_bin_noncovid_death,
#                 cens_date_dereg, cens_bin_dereg,
#                 out_date_severecovid, out_bin_severecovid,
#                 elig_date_ckd_45_first, elig_bin_ckd_45_first,
#                 exp_date_sulfo_first, exp_bin_sulfo_first,
#                 exp_date_insulin_first, exp_bin_insulin_first,
#                 cov_cat_sex, cov_num_age) %>%
#   View()


# Save the dataset ------------------------------------------------------
arrow::write_feather(df_months_tot, here::here("output", "data", "df_months_tot.arrow"))