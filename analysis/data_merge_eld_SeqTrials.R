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
library(data.table)
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
print('Import the main dataset')
df_long_months <- read_feather(here("output", "data", "df_long_months.arrow"))


# Import ELD tables and pre-process --------------------------------------
print('Import the ELD and combine it')
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
df_covid_ec$eld_out_bin_covid_ec <- TRUE # because this variable is defined differently in ehrQL, revisit later 
df_eld <- bind_rows(df_bmi, df_chol, df_covid_ec, df_covid_hosp, df_covid_pc, df_covid_sgss, df_covid_vaccinations, 
                    df_hba1c, df_hdl, df_metfin_interactions, df_obesity_pc, df_obesity_sc)
df_eld <- df_eld %>%
  arrange(patient_id, date)


# Merge all ELD to long format main dataset ------------------------------
# df_long_months only contains eligible individuals, i.e. alive, 18 or above, registered, with T2DM at pandemic start (max back to mid2018), and no contraindication to metfin before pandemic start
setDT(df_long_months)
setDT(df_eld)

# Filter df_eld to only patient_ids in df_long_months
df_eld <- df_eld[patient_id %in% df_long_months$patient_id]
# Filter df_long_months to reduce to only the necessary columns
df_long_months_sub <- df_long_months[, .(patient_id, start_date_month, end_date_month, month)]

# Perform the non-equi join from df_eld to df_long_months
# CAVE: this overwrites start_date_month and end_date_month => join back using month
# We ensured no overlap of start & end in fn_expand_intervals
interval_joined <- df_long_months_sub[df_eld, 
                             on = .(patient_id, start_date_month <= date, end_date_month >= date), 
                             nomatch = 0]
setDT(interval_joined)

## Important:
## If there are several different events happening in the same interval/month (e.g. eld_out_bin_covid_test == TRUE and eld_out_bin_covid_hosp == TRUE), 
## then there will be several rows for the same patient_id and same month and showing both entries across the two different columns
### -> keep both/all events!
## If there are two of the same event happening in the same interval/month (e.g. eld_out_bin_covid_test twice recorded), 
## then there will be several rows for the same patient_id and same month and showing both entries across the same column
### -> keep only 1 event, e.g. first or random

# 2. Handle duplicate rows: if there are multiple values in the same column for the same month, keep one randomly.
interval_joined_cleaned <- interval_joined %>%
  group_by(patient_id, month) %>%
  mutate(across(everything(), ~ifelse(length(unique(.)) > 1, sample(unique(.), 1), .), .names = "{.col}_cleaned")) %>%
  ungroup()

# 3. Clean up to keep only the original columns and not the cleaned ones.
interval_joined_cleaned_final <- interval_joined_cleaned %>%
  select(-contains("_cleaned"))

# Merge the unique person/month rows with ELD info back to df_long_months
df_long_months_eld <- merge(df_long_months, interval_joined_cleaned_final, 
                            by = c("patient_id", "month"), 
                            all.x = TRUE)

df_long_months_eld %>%
  # dplyr::select(patient_id, start_date_month, month, end_date_month, stop_date, outcome, comp_event, censor,
  #               elig_date_t2dm, out_date_covid_death, out_date_death, cens_date_dereg,
  #               cov_cat_sex, cov_num_age, 
  #               eld_cov_num_bmi, eld_cov_num_chol, eld_out_bin_covid_ec, eld_out_bin_covid_hosp, eld_out_bin_covid_pc, 
  #               eld_out_bin_covid_test, eld_cov_cat_vacc, eld_cov_bin_obesity_pc, eld_cov_bin_obesity_sc
  #               ) %>%
  dplyr::filter(patient_id == 12) %>% 
  View()
interval_joined %>%
  dplyr::select(patient_id, start_date_month, month, end_date_month, stop_date, outcome, comp_event, censor,
                elig_date_t2dm, out_date_covid_death, out_date_death, cens_date_dereg,
                cov_cat_sex, cov_num_age, 
                eld_cov_num_bmi, eld_cov_num_chol, eld_out_bin_covid_ec, eld_out_bin_covid_hosp, eld_out_bin_covid_pc, 
                eld_out_bin_covid_test, eld_cov_cat_vacc, eld_cov_bin_obesity_pc, eld_cov_bin_obesity_sc
  ) %>%
  dplyr::filter(patient_id == 12) %>% 
  View()
df_long_months %>%
  dplyr::select(patient_id, start_date_month, month, end_date_month) %>%
  dplyr::filter(patient_id == 12) %>% 
  View()
df_eld %>%
  dplyr::filter(patient_id == 12) %>% 
  View()

