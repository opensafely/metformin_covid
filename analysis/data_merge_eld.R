####
## This script does the following:
# 1. Import person-month dataframe and eld data
# 2. Expand it to long format with monthly interval
# 3. Saves all output
####

# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(data.table) # only for interval joining with "on"... reconsider using tidyverse (or data.table), not both
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


# Import the main person-month dataset ------------------------------------
print('Import the main person-month dataset')
df <- read_feather(here("output", "data", "df_months_tot.arrow"))


# Import ELD tables and pre-process --------------------------------------
print('Import the ELD, reformat and combine')
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

df_covid_ec$eld_out_bin_covid_ec <- TRUE # because this variable is defined differently in ehrQL, double-check later

df_eld <- bind_rows(df_bmi, df_chol, df_covid_ec, df_covid_hosp, df_covid_pc, df_covid_sgss, df_covid_vaccinations, 
                    df_hba1c, df_hdl, df_metfin_interactions, df_obesity_pc, df_obesity_sc)
df_eld <- df_eld %>%
  arrange(patient_id, date)


# Add intervals to ELD and restructure it to 1 row per person -------------
print('Add intervals to ELD and restructure it to 1 row per person')
# Note: df only contains eligible individuals, i.e. alive, 18 or above, registered, with T2DM at pandemic start (max back to mid2018), and no contraindication to metfin before pandemic start
setDT(df)
setDT(df_eld)

### (1) First, work on ELD dataframe
# filter df_eld to only patient_ids in df (i.e. only work with eligible ELD)
df_eld <- df_eld[patient_id %in% df$patient_id]
# reduce df to only the necessary columns to compare with df_eld
df_sub <- df[, .(patient_id, start_date_month, end_date_month, month)]
# Perform the non-equi join to get a dataset with only person-intervals plus the ELD events
# Even if (end) interval is only 1 day long (i.e. same start_date_month == end_date_month) it will behave correctly
interval_joined <- df_sub[df_eld, 
                          on = .(patient_id, start_date_month <= date, end_date_month >= date),
                          nomatch = 0]
# CAVE:
## Scenario 1: If there are several different events happening in the same interval/month (e.g. eld_out_bin_covid_test == TRUE and eld_out_bin_covid_hosp == TRUE), 
## then there will be several rows for the same patient_id and same month, showing both entries across the two different columns
### -> we want to keep both/all events but shift them to 1 row only!
## Scenario 2: If there are two of the same event happening in the same interval/month (e.g. eld_out_bin_covid_test twice recorded), 
## then there will be several rows for the same patient_id and same month, showing both entries across the same column
### -> keep only 1 event, at random!

## (a) Fix scenario 2 first: If there are multiple values in the same column for the same month, choose one randomly (in _cleaned) and overwrite original. I do not reduce rows yet.
# Rules:
# If there is only one unique non-NA value in the column for a given patient_id and month, it will be retained.
# If there are multiple non-NA values in the column, one will be randomly selected using sample().
# If there are only NA values in the column, the column will retain NA because unique(.) on all NAs results in a vector of length 1 with NA, so the ifelse statement will just keep that value (NA).
# the final cleaned variable will have the selected value across all the duplicate rows!

# Convert back to tibble
interval_joined <- as_tibble(interval_joined)

interval_joined_cleaned <- interval_joined %>%
  group_by(patient_id, month) %>%
  mutate(across(starts_with("eld_"), 
                ~ {
                  vals <- unique(na.omit(.))
                  if (length(vals) == 0) {
                    NA # if all values are NA: NA
                  } else if (length(vals) == 1) {
                    vals # if only one unique non-NA value: pick that one
                  } else {
                    sample(vals, 1) # if multiple non-NA value: randomly pick one
                  }
                }, 
                .names = "{.col}")) %>%
  ungroup()


## (b) Now work on scenario 1 and collapse rows to keep all unique info only in 1 row
# Since above code duplicated the info within each eld column, we can simply filter with distinct()
# Collapse across person-month intervals, to keep only 1 person-month interval
interval_joined_collapsed <- interval_joined_cleaned %>%
  distinct(patient_id, month, .keep_all = TRUE)

# To double-check: find how many patients have multiple non-NA entries across eld_* columns for the same month
# scenario_1_patients <- interval_joined_cleaned %>%
#   group_by(patient_id, month) %>%
#   # Create a temporary column to count the number of non-NA values across all eld_* columns for each patient-month group
#   mutate(non_na_count = rowSums(!is.na(select(cur_data(), starts_with("eld_"))))) %>%
#   ungroup() %>%
#   # Filter for cases where there is more than 1 non-NA value across eld_* columns
#   filter(non_na_count > 1) %>%
#   distinct(patient_id, month)

# interval_joined %>%
#   dplyr::filter(patient_id == 764) %>%
#   View()
# interval_joined_cleaned %>%
#   dplyr::filter(patient_id == 764) %>%
#   View()
# interval_joined_collapsed %>%
#   dplyr::filter(patient_id == 764) %>%
#   View()


# Merge ELD to main dataframe ---------------------------------------------
print('Merge ELD to main dataframe')
# ensure tidyverse
interval_joined_collapsed <- as_tibble(interval_joined_collapsed)
df <- as_tibble(df)

df_months_eld <- df %>%
  left_join(
    interval_joined_collapsed %>% select(-c(start_date_month, end_date_month)),
    by = c("patient_id", "month")
  )
# should have same length as initial df

# df_months_eld %>%
#   dplyr::select(patient_id, month, contains("eld_")) %>%
#   # dplyr::filter(patient_id == 764) %>%
#   View()


# Reformat ELD columns in the main dataframe ---------------------------
print('Reformat ELD columns in the main dataframe')
# We only have numeric ELD and bin/logical (only TRUE values!) ELD
# num: Keep NA until first numeric value is recorded, then fill up until next numeric value or max follow-up
# bin: Add 1 to person-intervals with TRUE and 0 to all other person-intervals
eld_num_cols <- names(df_months_eld) %>% stringr::str_subset("^eld_.*_num_")
eld_bin_cols <- names(df_months_eld) %>% stringr::str_subset("^eld_.*_bin_")

df_months_eld <- df_months_eld %>%
  mutate(across(
    all_of(eld_bin_cols),
    ~ case_when(
      is.logical(.x) & .x == TRUE ~ 1L, 
      is.logical(.x) & is.na(.x) ~ 0L,
      is.logical(.x) & .x == FALSE ~ 0L, # but should not exist
      TRUE ~ NA_integer_
    )
  ))

# df_months_eld <- df_months_eld %>%
#   group_by(patient_id) %>%
#   arrange(start_date_month, .by_group = TRUE) %>%
#   mutate(across(all_of(eld_num_cols), ~ {
#     # Identify where first non-NA occurs
#     first_non_na <- which(!is.na(.x))[1]
#     if (is.na(first_non_na)) {
#       # all NA → return as is
#       .x
#     } else {
#       # mask values before first non-NA
#       masked <- .x
#       masked[1:(first_non_na - 1)] <- NA
#       # forward fill
#       tidyr::fill(data.frame(masked), masked, .direction = "down")$masked
#     }
#   }))


# Assign time-varying treatment and eligibility variables --------------
print('Assign time-varying treatment and eligibility variables')
# create treat and elig variables 



