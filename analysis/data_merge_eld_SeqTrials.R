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
library(data.table) # only for interval joining with "on"
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


# Add intervals to ELD and restructure it to 1 row per person -------------
# df_long_months only contains eligible individuals, i.e. alive, 18 or above, registered, with T2DM at pandemic start (max back to mid2018), and no contraindication to metfin before pandemic start
setDT(df_long_months)
setDT(df_eld)

### (1) First, work on ELD dataframe
# filter df_eld to only patient_ids in df_long_months
df_eld <- df_eld[patient_id %in% df_long_months$patient_id]
# reduce df_long_months to only the necessary columns to compare with df_eld
df_long_months_sub <- df_long_months[, .(patient_id, start_date_month, end_date_month, month)]
# Perform the non-equi join to get a dataset with only person-intervals with ELD events
# CAVE: this overwrites start_date_month and end_date_month => join back using month
interval_joined <- df_long_months_sub[df_eld, 
                                      on = .(patient_id, start_date_month <= date, end_date_month >= date), # We ensured no overlap of start_date_month & end_date_month in fn_expand_intervals !
                                      nomatch = 0]
# CAVE:
## Scenario 1: If there are several different events happening in the same interval/month (e.g. eld_out_bin_covid_test == TRUE and eld_out_bin_covid_hosp == TRUE), 
## then there will be several rows for the same patient_id and same month and showing both entries across the two different columns
### -> we want to keep both/all events but in 1 row only!
## Scenario 2: If there are two of the same event happening in the same interval/month (e.g. eld_out_bin_covid_test twice recorded), 
## then there will be several rows for the same patient_id and same month and showing both entries across the same column
### -> keep only 1 event, e.g. at random

## (a) Fix scenario 2 first: If there are multiple values in the same column for the same month, choose one randomly (in _cleaned) and overwrite original. I do not reduce rows yet.
# Rules:
# If there is only one unique non-NA value in the column for a given patient_id and month, it will be retained.
# If there are multiple non-NA values in the column, one will be randomly selected using sample().
# If there are only NA values in the column, the column will retain NA because unique(.) on all NAs results in a vector of length 1 with NA, so the ifelse statement will just keep that value (NA).
# the final cleaned variable will have the selected value across all the duplicate rows!
interval_joined_cleaned <- interval_joined %>%
  group_by(patient_id, month) %>%
  mutate(across(starts_with("eld_"), ~ifelse(length(unique(.)) > 1, sample(unique(.), 1), .), .names = "{.col}_cleaned")) %>%
  ungroup() %>%
  select(patient_id, month, ends_with("_cleaned")) %>%
  rename_with(~ gsub("_cleaned", "", .x))

## (b) Now work on scenario 1 and collapse rows to keep all unique info only in 1 row
# Since above code duplicated the info within each eld column, it is easy to filter with distinct()
# But to double-check, first find how many patients have multiple non-NA entries across eld_* columns for the same month
scenario_1_patients <- interval_joined_cleaned %>%
  group_by(patient_id, month) %>%
  # Create a temporary column to count the number of non-NA values across all eld_* columns for each patient-month group
  mutate(non_na_count = rowSums(!is.na(select(cur_data(), starts_with("eld_"))))) %>%
  ungroup() %>%
  # Filter for cases where there is more than 1 non-NA value across eld_* columns
  filter(non_na_count > 1) %>%
  distinct(patient_id, month)
# Collapse across person-month intervals, to keep only 1 person-month interval
interval_joined_collapsed <- interval_joined_cleaned %>%
  distinct(patient_id, month, .keep_all = TRUE)

# interval_joined %>%
#   dplyr::filter(patient_id == 894) %>% 
#   View()
# interval_joined_cleaned %>%
#   dplyr::filter(patient_id == 894) %>% 
#   View()
# interval_joined_collapsed %>%
#   dplyr::filter(patient_id == 894) %>% 
#   View()


# Merge ELD to main dataframe ---------------------------------------------
# Left join the collapsed data back to the original data; rows/obs should remain the same
df_long_months_eld <- df_long_months %>%
  left_join(interval_joined_collapsed, by = c("patient_id", "month"))

df_long_months_eld %>%
  dplyr::select(patient_id, month, contains("eld_")) %>%
  # dplyr::filter(patient_id == 12) %>%
  View()

df_long_months_eld %>%
  dplyr::select(patient_id, month, contains("eld_")) %>%
  dplyr::filter(patient_id == 5638) %>%
  View()

df_long_months_eld %>%
  dplyr::select(patient_id, month, contains("eld_")) %>%
  dplyr::filter(patient_id == 2701) %>%
  View()



# Reformat ELD events across the main dataframe ---------------------------
# Fill up within each ELD column: We only have only num and bin (only TRUE values!) across ELD columns.
# num: Keep NA until num is recorded, then fill up until next num or max follow-up
# bin: Add 1 to TRUE and 0 to all the other person-month intervals
df_long_months_eld <- df_long_months_eld %>%
  rowwise() %>% 
  group_by(patient_id) %>%
  mutate(across(starts_with("eld_") & contains("_bin_"), 
                ~case_when(.x == TRUE ~ 1, TRUE ~ 0))) %>% 
  mutate(across(starts_with("eld_") & contains("_num_"),
                ~na.locf(.x, na.rm = FALSE), .names = "{.col}")) %>%
  ungroup()

df_long_months_eld %>%
  dplyr::select(patient_id, start_date_month, end_date_month, month, outcome, comp_event, censor,
                
                cov_bin_obesity_b, cov_num_bmi_b, cov_cat_bmi_groups_b, cov_num_tc_hdl_ratio_b, cov_cat_tc_hdl_ratio_b,
                cov_num_hba1c_mmol_mol_b, cov_cat_hba1c_mmol_mol_b, out_date_severecovid, out_bin_severecovid, 
                
                elig_bin_metfin_allergy_first, elig_bin_ckd_45_first, elig_bin_liver_cirrhosis_first, elig_bin_metfin_interaction_first,
                
                exp_bin_metfin_first, exp_bin_metfin_mono_first, exp_bin_sulfo_first, exp_bin_dpp4_first, exp_bin_tzd_first,
                exp_bin_megli_first, exp_bin_agi_first, exp_bin_insulin_first,
                
                eld_cov_num_bmi, eld_cov_num_chol, eld_out_bin_covid_ec, eld_out_bin_covid_hosp, eld_out_bin_covid_pc,
                eld_out_bin_covid_test, eld_cov_bin_vacc, eld_cov_num_hba1c, eld_cov_num_hdl, eld_elig_bin_metfin_interaction, 
                eld_cov_bin_obesity_pc, eld_cov_bin_obesity_sc) %>%
  View()


# Assign time-varying treatment and eligibility variables --------------
df_long_months_eld <- df_long_months_eld %>%
  rowwise() %>% 
  group_by(patient_id) %>%
  mutate(across(starts_with("elig_") & contains("_bin_"), 
                ~case_when(.x == TRUE ~ 1, TRUE ~ 0))) %>% 
  mutate(across(starts_with("eld_") & contains("_num_"),
                ~na.locf(.x, na.rm = FALSE), .names = "{.col}")) %>%
  ungroup()


