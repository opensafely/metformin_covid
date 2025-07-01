####
## This script does the following:
# 1. Import pre-processed data of a one-person-per-row dataset of eligible participants
# 2. Restructure the one-person-per-row data into a multiple-event-per-person long format dataset, with time intervals as per study protocol
# 3. Save dataset(s)
####


# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(lubridate) # for fn_expand_intervals.R
library(zoo) # for fn_add_firstever_events_to_intervals.R
library(purrr) # for fn_expand_intervals.R and fn_add_firstever_events_to_intervals.R
library(rlang) # for fn_expand_intervals.R and fn_add_firstever_events_to_intervals.R
library(stringr) # for fn_add_firstever_events_to_intervals.R
source(here::here("analysis", "functions", "fn_expand_intervals2.R"))
source(here::here("analysis", "functions", "fn_add_firstever_events_to_intervals.R"))


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data"))


# Import the data ---------------------------------------------------------
print('Import the data')
df <- read_feather(here("output", "data", "data_processed.arrow"))


# Import dates ------------------------------------------------------------
print('Import dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# Expand the dataset ------------------------------------------------------
print('Expand the dataset')
# In this case, these are monthly intervals. I expand the dataset for everyone until: 
# 1) Death, recorded in: a) out_date_covid_death (one of the outcomes), b) out_date_noncovid_death (the main competing event)
# 2) De-registration, i.e. lost-to-follow-up: cens_date_dereg (estimand: "assume no-one is LTFU" with IPCW; as part of ITT)











# # COVID-related hospitalization per se is not an interval stopping event (but an important outcome and elig criteria), will be added later, using ELD.
# # Other outcomes (e.g. COVID diagnosis), censoring events (e.g. starting metformin in control group) and time-updated covariates will be added later, using ELD.
# 
# # Apply the function
# ## CAVE: stop_date_columns shall only have dates AFTER start_date_variable
# ## CAVE: only choose stopping events that are distinct (e.g. out_date_noncovid_death is NOT part of out_date_covid_death)
# df_months <- fn_expand_intervals(df, 
#                                  start_date_variable = pandemicstart_date,
#                                  stop_date_columns = c("out_date_covid_death", "out_date_noncovid_death", "cens_date_dereg"), 
#                                  studyend_date,
#                                  interval_type = "month")
# ## NOTE: the function will add the final interval, even if the final interval is just 1 day (check stop_date == 2022-04-01)
# 
# # Expand the dataset into intervals and follow the rules see in fn_expand_intervals.R script
# start_date_variable <- "landmark_date" # start expanding at landmark date, not elig_date_t2dm (due to landmark design)
# stop_date_columns <- c("out_date_severecovid_afterlandmark", "out_date_death_afterlandmark", "cens_date_ltfu_afterlandmark")
# outcome_date_variable <- "out_date_severecovid_afterlandmark"
# comp_date_variable <- "out_date_death_afterlandmark"
# censor_date_variable <- "cens_date_ltfu_afterlandmark"
# 
# # Apply the function, currently only using months
# # Use studyend_date from metadates.R import
# df_long_months <- fn_expand_intervals(df, 
#                                       start_date_variable,
#                                       stop_date_columns, 
#                                       studyend_date,
#                                       outcome_date_variable, 
#                                       comp_date_variable,
#                                       censor_date_variable,
#                                       interval_type = "month")
# 
# ## To double-check | to check how events are treated in interval the event happend, see stop_date and start_date_month
# # df_long_months %>%
# #   dplyr::select(patient_id, elig_date_t2dm, landmark_date, out_date_severecovid_afterlandmark, out_date_death_afterlandmark,
# #                 cens_date_ltfu_afterlandmark, stop_date, start_date_month, end_date_month, month, outcome, comp_event, censor, is_event_interval, is_outcome_event, is_cens_event,
# #                 is_comp_event, followup_stop,
# #                 cov_cat_sex, cov_num_age
# #                 ) %>%
# #   # dplyr::filter(!is.na(out_date_death_afterlandmark)) %>%
# #   # dplyr::filter(!is.na(cens_date_ltfu_afterlandmark)) %>%
# #   # dplyr::filter(!is.na(out_date_severecovid_afterlandmark)) %>%
# #   # dplyr::filter(is.na(censor)) %>% # should be empty
# #   # dplyr::filter(stop_date == "2022-04-01") %>%
# #   # dplyr::filter(is_event_interval == TRUE & !is.na(out_date_severecovid_afterlandmark)) %>% # should only have 1 row per person
# #   # dplyr::filter(followup_stop == TRUE) %>% # should be empty, otherwise it means that there are two defining event intervals
# #   View()
# 
# 
# 
# 
# # Save output -------------------------------------------------------------
# print('Save output')

