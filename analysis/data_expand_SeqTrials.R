####
## This script does the following:
# 1. Import processed data
# 2. Expand it to long format with monthly interval
# 3. Saves all output
####

# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)

# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_descriptiond_seqtrials.R"))

# Import dates ------------------------------------------------------------
print('Import the dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")

# Import the data ---------------------------------------------------------
print('Import the data')
df <- read_feather(here("output", "data", "data_processed_seqtrials.arrow"))

# Expand the dataset ------------------------------------------------------
print('Expand the dataset')
# Expand the dataset into intervals and follow these rules:
# a) If outcome (out_date_severecovid_afterlandmark) is reached first, assign outcome=1, censor=0, comp_event=0 to the interval when it happened and stop expanding
# b) If competing event (out_date_death_afterlandmark) is reached first, assign outcome=NA, censor=0, comp_event=1 to the interval when it happened and stop expanding
# c) If censoring event (out_date_ltfu_afterlandmark) is reached first, assign outcome=NA, censor=1, comp_event=NA to the interval when it happened and stop expanding
# d) If studyend_date is reached first, then assign outcome=0, censor=0, comp_event=0 to the interval when it happened and stop expanding
# Use studyend date from metadates.R import
start_date_variable <- "landmark_date" # start expanding at landmark date, not elig_date_t2dm (due to landmark design)
stop_date_columns <- c("out_date_severecovid_afterlandmark", "out_date_death_afterlandmark", "cens_date_ltfu_afterlandmark")
outcome_date_variable <- "out_date_severecovid_afterlandmark"
comp_date_variable <- "out_date_death_afterlandmark"
censor_date_variable <- "cens_date_ltfu_afterlandmark"

# Apply the function, choose either weeks or months, currently only using months, but works for both.
df_long_months <- fn_expand_intervals(df, 
                                      start_date_variable,
                                      stop_date_columns, 
                                      studyend_date,
                                      outcome_date_variable, 
                                      comp_date_variable,
                                      censor_date_variable,
                                      interval_type = "month")

## To double-check | to check how events are treated in interval the event happend, see stop_date and start_date_month
# df_long_months %>%
#   dplyr::select(patient_id, elig_date_t2dm, landmark_date, out_date_severecovid_afterlandmark, out_date_death_afterlandmark,
#                 cens_date_ltfu_afterlandmark, stop_date, start_date_month, month, outcome, comp_event, censor, is_event_interval, is_outcome_event, followup_stop,
#                 cov_cat_sex, cov_num_age) %>%
#   # dplyr::filter(!is.na(out_date_death_afterlandmark)) %>%
#   # dplyr::filter(!is.na(cens_date_ltfu_afterlandmark)) %>%
#   # dplyr::filter(!is.na(out_date_severecovid_afterlandmark)) %>%
#   # dplyr::filter(is.na(censor)) %>%
#   # dplyr::filter(stop_date == "2022-04-01") %>%
#   # dplyr::filter(is_event_interval == TRUE & !is.na(out_date_severecovid_afterlandmark)) %>% # should only have 1 row per person
#   # dplyr::filter(patient_id == "1418") %>%
#   # dplyr::filter(is.na(censor)) %>% # should be empty
#   View()
