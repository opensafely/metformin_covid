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
source(here::here("analysis", "functions", "fn_extract_data.R"))


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_description_seqtrials"))


# Import dates ------------------------------------------------------------
print('Import the dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# Import the data ---------------------------------------------------------
print('Import the data')
df_long_months <- read_feather(here("output", "data", "df_long_months.arrow"))


# Import ELD tables and pre-process --------------------------------------
input_filename <- "bmi.arrow"
df_bmi <- fn_extract_data(input_filename)
input_filename <- "chol.arrow"
df_chol <- fn_extract_data(input_filename)
input_filename <- "covid_ec.arrow"
df_covid_ec <- fn_extract_data(input_filename)
input_filename <- "covid_hosp.arrow"
df_covid_hosp <- fn_extract_data(input_filename)
input_filename <- "covid_pc.arrow"
df_covid_pc <- fn_extract_data(input_filename)
input_filename <- "covid_sgss.arrow"
df_covid_sgss <- fn_extract_data(input_filename)
input_filename <- "covid_vaccinations.arrow"
df_covid_vaccinations <- fn_extract_data(input_filename)
input_filename <- "hba1c.arrow"
df_hba1c <- fn_extract_data(input_filename)
input_filename <- "hdl.arrow"
df_hdl <- fn_extract_data(input_filename)
input_filename <- "metfin_interactions.arrow"
df_metfin_interactions <- fn_extract_data(input_filename)
input_filename <- "obesity_pc.arrow"
df_obesity_pc <- fn_extract_data(input_filename)
input_filename <- "obesity_sc.arrow"
df_obesity_sc <- fn_extract_data(input_filename)

# ELD table crunching ---------------------------------------------------
df_covid_ec$eld_out_bin_covid_ec <- TRUE
df_eld <- bind_rows(
  df_bmi, df_chol, df_covid_ec, df_covid_hosp, df_covid_pc, df_covid_sgss, df_covid_vaccinations, df_hba1c, df_hdl, 
  df_metfin_interactions, df_obesity_pc, df_obesity_sc)
df_eld <- df_eld %>%
  arrange(patient_id, date)


