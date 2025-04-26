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
df <- fn_extract_data(input_filename)


# Import ELD tables and pre-process --------------------------------------
input_filename <- "pc_obesity.arrow"
df_pc_obesity <- fn_extract_data(input_filename)
input_filename <- "pc_covid.arrow"
df_pc_covid <- fn_extract_data(input_filename)
input_filename <- "hba1c.arrow"
df_pc_hba1c <- fn_extract_data(input_filename)
input_filename <- "chol.arrow"
df_pc_chol <- fn_extract_data(input_filename)
input_filename <- "hdl.arrow"
df_pc_hdl <- fn_extract_data(input_filename)
input_filename <- "sc_obesity.arrow"
df_sc_obesity <- fn_extract_data(input_filename)
input_filename <- "sc_covid.arrow"
df_sc_covid <- fn_extract_data(input_filename)
input_filename <- "ec_covid.arrow"
df_ec_covid <- fn_extract_data(input_filename)
input_filename <- "sgss_covid.arrow"
df_sgss_covid <- fn_extract_data(input_filename)
input_filename <- "covid_vaccinations.arrow"
df_covid_vaccinations <- fn_extract_data(input_filename)


# ELD table crunching ---------------------------------------------------
df_ec_covid$ec_covid <- TRUE
df_eld <- bind_rows(df_pc_obesity, df_pc_covid, df_pc_hba1c, df_pc_chol, df_pc_hdl, df_sc_obesity, df_sc_covid, df_ec_covid, df_sgss_covid, df_covid_vaccinations)
df_eld <- df_eld %>%
  arrange(patient_id, date)


# Add "time-fixed" main df to ELD table ----------------------------------
## df is my main df, left_join will ensure no unnecessary ELD is joined and at the same time creates ELD if they exist in df
df_merged <- df %>%
  left_join(df_eld, by = "patient_id")

# df_merged %>% 
#   select(patient_id, date, pc_obesity, pc_covid, hba1c_num, chol_num, hdl_num, sc_obesity, sc_covid, ec_covid, sgss_covid, vaccine_product, elig_date_t2dm, out_date_covid_death, exp_date_metfin_mono_first, qa_bin_was_alive, cov_num_age, cov_cat_sex, strat_cat_region) %>% 
#   View()
