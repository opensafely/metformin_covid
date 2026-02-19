####
## This script does the following:
# 1. Import pre-processed data of a one-person-per-row dataset of eligible participants (data_processed.arrow)
# 2. Restructure the one-person-per-row data into a multiple-event-per-person long format dataset, with time intervals as per study protocol
# 3. Save dataset(s)
####


# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(lubridate) # for fn_expand_intervals.R
library(rlang) # for fn_expand_intervals.R
source(here::here("analysis", "functions", "fn_expand_intervals.R"))


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


# Filter main dataset for PLR pipeline: Only keep necessary variables -----
df <- df %>% 
  select(patient_id, 
         elig_date_t2dm, 
         landmark_date,
         max_fup_date,
         exp_bin_treat,
         starts_with("strat_"), # Stratification variable
         starts_with("cov_"), # Covariates
         starts_with("out_"), # Outcomes
         starts_with("cens_") # Censoring variable
  )


# Expand the dataset ------------------------------------------------------
print('Expand the dataset')
# (1) Define the interval: In this scenario, we use monthly intervals, balance between how much is happening in a given interval vs granularity vs computational
# (2) Define the dataset stopping events: We will expand the dataset until 
## a. any death (out_date_death_afterlandmark), 
## b. deregistration from TPP (cens_date_ltfu_afterlandmark), 
## c. end of study (studyend_date = 01.04.2022)
## d. maximum individual follow-up (max_fup_date = landmark_date + days(730))
# (3) Our outcomes of interest: 
## a. out_date_severecovid_afterlandmark (includes covid deaths): Primary outcome
## b. out_date_noncoviddeath_afterlandmark: Competing event. In primary analysis, we simply treat it as censoring event ("controlled direct effect")
#### Side note re competing event: We may - at a later stage - consider targeting a "total effect" (by assigning outcome == 0 for everyone who had a non-covid death and including their follow-up time after non-covid deaths, i.e. not censoring them)
## c. cens_date_ltfu_afterlandmark: Censoring event. We may want to upweigh those who were not censored using a hypothetical treatment strategy for this intercurrent event ("assume no-one is LTFU") 
## d. cens_date_metfin_start_cont: Censoring event for the per-protocol analysis ONLY. We will upweigh those who did not start metformin in control group, using a hypothetical treatment strategy for this intercurrent event ("assume no-one started metformin in control")
## e. out_date_covid_afterlandmark: Secondary outcome (includes covid tests, diagnoses, hosp & deaths). All other competing/censoring events same as above.
## f. out_date_longcovid_virfat_afterlandmark: Secondary outcome (includes viral fatigue and long covid codes). All other competing/censoring events same as above, except we use all-cause death as competing event: 
## g. out_date_death_afterlandmark: All-cause mortality, the competing event variable only for analysis of the secondary endpoint Long covid.

# Apply the interval expansion function
## CAVE: All stop_date_columns dates must be AFTER start_date_variable
## CAVE: only choose stopping events that are distinct and overarching (e.g. all-cause death instead of cause-specific deaths or overlapping outcomes)
df_months <- fn_expand_intervals(df, 
                                 start_date_variable = "landmark_date", # start expanding from landmark date, not elig_date_t2dm (due to landmark design)
                                 stop_date_columns = c("out_date_death_afterlandmark", "cens_date_ltfu_afterlandmark", "max_fup_date"),
                                 studyend_date,
                                 interval_type = "month") # currently, function takes "month" or "week"
## NOTE: the function will add the final interval, even if the final interval is just 1 day (check stop_date == 2022-04-01)

# To double-check
# df_months %>%
#   dplyr::select(patient_id, elig_date_t2dm, landmark_date, out_date_severecovid_afterlandmark, out_date_death_afterlandmark,
#                 max_fup_date,
#                 cens_date_ltfu_afterlandmark, stop_date, start_date_month, end_date_month, month,
#                 cov_cat_sex, cov_num_age
#   ) %>%
#   # dplyr::filter(!is.na(out_date_death_afterlandmark)) %>%
#   # dplyr::filter(!is.na(cens_date_ltfu_afterlandmark)) %>%
#   # dplyr::filter(stop_date == "2022-04-01") %>%
#   View()


# Save output -------------------------------------------------------------
print('Save output')
# Save long format dataset
arrow::write_feather(df_months, here::here("output", "data", "df_months.arrow"))
