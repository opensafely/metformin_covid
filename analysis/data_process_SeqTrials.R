####
## This script does the following:
# 1. Import/extract feather dataset from OpenSAFELY
# 2. ...
# 3. ...
# 4. ...
####


# Import libraries and user functions -------------------------------------
library('arrow')
library('readr')
library('here')
library('lubridate')
library('dplyr')
library('tidyr')
library('purrr')
library('forcats')
library('jsonlite')
library('skimr')
library('splines')
source(here::here("analysis", "functions", "fn_extract_data.R"))
source(here::here("analysis", "functions", "utility.R"))
source(here::here("analysis", "functions", "fn_quality_assurance_midpoint6.R"))
source(here::here("analysis", "functions", "fn_completeness_criteria_midpoint6.R"))
source(here::here("analysis", "functions", "fn_elig_criteria_midpoint6.R"))


# Create directories for output -------------------------------------------
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_description"))


# Import dates ------------------------------------------------------------
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# Define redaction threshold ----------------------------------------------
threshold <- 6


# Import the dataset and pre-process --------------------------------------
input_filename <- "dataset.arrow"
data_extracted <- fn_extract_data(input_filename)


# Import ELD tables and pre-process --------------------------------------
input_filename <- "clinical_events_pc.arrow"
df_clinical_events_pc <- fn_extract_data(input_filename)
unique(df_clinical_events_pc$covid) # when NA and when FALSE?
unique(df_clinical_events_pc$obesity) 
unique(df_clinical_events_pc$hba1c) 
unique(df_clinical_events_pc$chol)
unique(df_clinical_events_pc$hdl)


input_filename <- "clinical_events_num_pc.arrow"
df_clinical_events_num_pc <- fn_extract_data(input_filename)

input_filename <- "clinical_events_ctv3_pc.arrow"
df_clinical_events_ctv3_pc <- fn_extract_data(input_filename)
unique(df_clinical_events_ctv3_pc$date)

input_filename <- "apcs_covid.arrow"
df_apcs_covid <- fn_extract_data(input_filename)
input_filename <- "apcs_obesity.arrow"
df_apcs_obesity <- fn_extract_data(input_filename)
input_filename <- "covid_vaccinations.arrow"
df_covid_vaccinations<- fn_extract_data(input_filename)
input_filename <- "ec_covid.arrow"
df_ec_covid <- fn_extract_data(input_filename)
input_filename <- "sgss_covid.arrow"
df_sgss_covid <- fn_extract_data(input_filename)




