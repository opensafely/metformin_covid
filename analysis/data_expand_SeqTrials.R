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
library(lubridate)
source(here::here("analysis", "functions", "fn_expand_intervals_SeqTrials.R"))

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
## I define the following "stopping events" for dataset expansion:
# out_date_covid_death (one of the outcomes), out_date_death (the main competing event), cens_date_dereg (to build "assume no-one is LTFU" with IPCW, as part of ITT)
# COVID-related hospitalizatin per se is not be a stopping event (but an important outcome and elig criteria), only deaths and LTFU
# Other outcomes and censoring events (e.g. starting metformin in control group) will be added to the expanded dataset, their are not a stopping event for the ITT analysis

## Expand into monthly intervals and follow these rules:
# a) If outcome (out_date_covid_death) is reached first, stop expanding, and assign to the PREVIOUS interval: outcome=1, censor=0, comp_event=0
# => this becomes an indicator for the outcome event in interval𝑘+ 1!
# b) If competing event (out_date_death) is reached first, stop expanding, assign to the CURRENT interval: outcome=NA, censor=0, comp_event=1
# c) If censoring event (cens_date_dereg) is reached first, stop expanding, assign to the CURRENT interval: outcome=NA, censor=1, comp_event=NA
# d) If studyend_date is reached first, stop expanding, assign to the CURRENT interval: outcome=0, censor=0, comp_event=0

start_date_variable <- pandemicstart_date
stop_date_columns <- c("out_date_covid_death", "out_date_death", "cens_date_dereg")
outcome_date_variable <- "out_date_covid_death"
comp_date_variable <- "out_date_death"
censor_date_variable <- "cens_date_dereg"

# df %>%
#   select(patient_id, elig_date_t2dm, out_date_covid_death, out_date_death, cens_date_dereg, exp_date_metfin_mono_first,
#          cov_num_age, cov_cat_hba1c_mmol_mol_b,
#          cov_cat_sex, cov_cat_ethnicity,
#          strat_cat_region,
#          cov_bin_all_stroke) %>%
#   View()
# summary(df$out_date_death)

# Apply the function, choose either weeks or months, currently only using months, but works for both.
## Make sure that stop_date_columns only have dates in the future of start_date_variable ! Build in a check !
## Make sure that the stopping events are unique (e.g. out_date_death is not part of out_date_covid_death!)
df_long_months <- fn_expand_intervals(df, 
                                      start_date_variable,
                                      stop_date_columns, 
                                      studyend_date,
                                      outcome_date_variable, 
                                      comp_date_variable,
                                      censor_date_variable,
                                      interval_type = "month")

## To double-check | to check how events are treated in interval the event happend
# df_long_months %>%
#   dplyr::select(patient_id, elig_date_t2dm, out_date_covid_death, out_date_death, cens_date_dereg,
#                 stop_date, start_date_month, month, outcome, comp_event, censor, is_event_interval, is_outcome_event, followup_stop,
#                 cov_cat_sex, cov_num_age) %>%
#   # dplyr::filter(!is.na(out_date_death)) %>%
#   # dplyr::filter(!is.na(cens_date_dereg)) %>%
#   # dplyr::filter(!is.na(out_date_covid_death)) %>%
#   # dplyr::filter(is.na(censor)) %>% # should be empty
#   # dplyr::filter(is_event_interval == TRUE & !is.na(out_date_covid_death)) %>% # should only have 1 row per person
#   # dplyr::filter(patient_id == "1571") %>%
#   View()

# Save the dataset ------------------------------------------------------
arrow::write_feather(df_long_months, here::here("output", "data", "df_long_months.arrow"))