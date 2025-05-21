####
## This script does the following:
# 1. Import pre-processed one-person-per-row main dataset and ELD data
# 2. Process ELD data, get it into one-person-per-row format
# 3. Merge ELD to one-person-per-row main dataset
# 4. Assign/add eligibility and treatment variables
# 5. Assign/add outcome, censoring and competing event variables and shift them up 1 interval, and create different analysis sets for the different outcomes to be analysed using PLR
# 6. Save
####

# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(data.table) # only for interval joining with "on"... reconsider using tidyverse (or data.table), not both
source(here::here("analysis", "functions", "fn_extract_data.R"))
source(here::here("analysis", "functions", "fn_add_and_shift_out_comp_cens_events.R"))


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_description_seqtrials"))


# Import dates ------------------------------------------------------------
print('Import the dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# Import the main person-month dataset -----------------------------------
print('Import the main person-month dataset')
df <- read_feather(here("output", "data", "df_months_tot.arrow"))


# Import ELD tables, using the standardized import, and combine them ------
print('Import ELD tables, using the standardized import, and combine them')

input_filename <- "bmi.arrow"
df_bmi <- fn_extract_data(input_filename)
input_filename <- "chol.arrow"
df_chol <- fn_extract_data(input_filename)
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

df_eld <- bind_rows(df_bmi, df_chol, df_covid_vaccinations, 
                    df_hba1c, df_hdl, df_metfin_interactions, df_obesity_pc, df_obesity_sc)

df_eld <- df_eld %>%
  arrange(patient_id, date)


# Add intervals to ELD ---------------------------------------------------
print('Add intervals to ELD')
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

# Restrict ELD to one-person-per-row ---------------------------------------
print('Restrict ELD to one-person-per-row')
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

interval_joined <- as_tibble(interval_joined) # Convert back to tibble

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
#   dplyr::filter(patient_id == 178) %>%
#   View()
# interval_joined_cleaned %>%
#   dplyr::filter(patient_id == 178) %>%
#   View()
# interval_joined_collapsed %>%
#   dplyr::filter(patient_id == 178) %>%
#   View()


# Merge ELD to main dataframe --------------------------------------------
print('Merge ELD to main dataframe')
# first, ensure convertion back to tidyverse
interval_joined_collapsed <- as_tibble(interval_joined_collapsed)
df <- as_tibble(df)

df_months_eld <- df %>%
  left_join(interval_joined_collapsed %>% select(-c(start_date_month, end_date_month)),
            by = c("patient_id", "month"))
# should have same length as "df"


# Reformat ELD columns in the main dataframe -----------------------------
print('Reformat ELD columns in the main dataframe')
# ELD: We only have numeric and bin/logical (only TRUEs, no FALSE)
eld_num_cols <- names(df_months_eld) %>% stringr::str_subset("^eld_.*_num_")
eld_bin_cols <- names(df_months_eld) %>% stringr::str_subset("^eld_.*_bin_")

# a) bin: Add 1 to person-intervals with TRUE and 0 to all other person-intervals
df_months_eld <- df_months_eld %>%
  mutate(across(
    all_of(eld_bin_cols),
    ~ case_when(
      is.logical(.x) & .x == TRUE ~ 1L, 
      is.logical(.x) & is.na(.x) ~ 0L,
      is.logical(.x) & .x == FALSE ~ 0L, # but should not occur
      TRUE ~ NA_integer_
    )
  ))

# b) num: Keep NAs across person-intervals until first numeric value is recorded, then fill up until next numeric value or max follow-up
df_months_eld <- df_months_eld %>%
  group_by(patient_id) %>%
  arrange(start_date_month, .by_group = TRUE) %>%
  mutate(across(all_of(eld_num_cols), ~ {
    first_non_na <- which(!is.na(.x))[1]
    if (is.na(first_non_na)) {
      # if all NA: return as is (= NA)
      .x
    } else {
      # mask values before first non-NA
      masked <- .x
      masked[1:(first_non_na - 1)] <- NA
      # forward fill
      tidyr::fill(data.frame(masked), masked, .direction = "down")$masked
    }
  }))

# c) for HbA1c also create a "first-ever" high HbA1c flag for elig below
df_months_eld <- df_months_eld %>%
  group_by(patient_id) %>%
  mutate(
    elig_bin_hba1c = if_else(!is.na(eld_cov_num_hba1c) & eld_cov_num_hba1c > 75, 1L, 0L),
    elig_bin_hba1c = as.integer(cumsum(elig_bin_hba1c) > 0)
    ) %>%
  ungroup()

# d) for obesity, create combined covariate based on eld_cov_num_bmi, eld_cov_bin_obesity_pc, eld_cov_bin_obesity_sc
df_months_eld <- df_months_eld %>%
  group_by(patient_id) %>%
  mutate(eld_cov_bin_obesity = case_when(eld_cov_bin_obesity_pc == 1 | eld_cov_bin_obesity_sc == 1 | (!is.na(eld_cov_num_bmi) & eld_cov_num_bmi > 30) ~ 1L,
    TRUE ~ 0L)) %>%
  ungroup()


# df_months_eld %>%
#   dplyr::select(patient_id, elig_date_t2dm, start_date_month, end_date_month, month, stop_date,
#                 out_date_covid_death, out_bin_covid_death,
#                 out_date_noncovid_death, out_bin_noncovid_death,
#                 cens_date_dereg, cens_bin_dereg,
#                 out_date_severecovid, out_bin_severecovid,
#                 out_date_covid, out_bin_covid,
#                 out_date_longcovid_virfat, out_bin_longcovid_virfat,
#                 elig_date_ckd_45_first, elig_bin_ckd_45_first,
#                 exp_date_sulfo_first, exp_bin_sulfo_first,
#                 contains("eld_"),
#                 cov_cat_sex, cov_num_age) %>%
#   dplyr::filter(patient_id == 4005) %>%
#   View()


# Assign time-varying treatment variable ---------------------------------
print('Assign time-varying treatment variable')
## Notes:
### 1) We assume that someone who started metformin continued metformin. Therefore, someone has treat == 1 after starting for all subsequent intervals (already implemented through fn_add_firstever_events_to_intervals)
### 2) If someone started metformin in/during the interval, that interval receives treat == 1. => An interval can have treat == 1 and outcome == 1 especially due to shift of outcomes, see below! And see Miguel Hernan DataSetup.
df_months_eld <- df_months_eld %>%
  group_by(patient_id) %>%
  mutate(treat = if_else(exp_bin_metfin_mono_first == 1, 1L, 0L)) %>%
  ungroup()


# Assign time-varying eligibility variable ------------------------------
print('Assign time-varying eligibility variable')
## Notes:
### 1) We use mix of a) "one-time / first-ever" variables whereby once it happened all subsequent person-intervals have the event (e.g. metfin allergy) and b) time-varying eligibility through ELD
### 2) Make use of lag variables where appropriate: for prior history of ...
### 3) For HbA1c: Use the flag variable created above since we look back total 2 years (see elig criteria) which is exactly the entire duration of the study.
### 4) The 4 analyses (= analysis sets) do not necessarily use the exact same eligibility criteria. E.g. for the long covid analysis, I also want to exclude prior "viral fatigue" history => Re-apply elig after dataset creation.
df_months_eld <- df_months_eld %>%
  group_by(patient_id) %>%
  mutate(elig = case_when(
    # any prior medical history that is an exclusion criteria (first-time/one-time ever events), i.e. can simply look at prior interval since these events have ==1 once it occurred
    # lag default = 0 ensures that lag code works for first interval (i.e. eligible per default)
    lag(elig_bin_metfin_allergy_first, default = 0) == 1 | 
      lag(elig_bin_ckd_45_first, default = 0) == 1 | 
      lag(elig_bin_liver_cirrhosis_first, default = 0) == 1 | 
      # any documented history of coronavirus infection, i.e. in prior interval
      lag(out_bin_covid, default = 0) == 1 | # includes death, hosp, diagnoses, tests
      lag(out_bin_longcovid, default = 0) == 1 | # includes long covid codes (but not viral fatigue)
      # dead 
      lag(out_bin_covid_death, default = 0) == 1 |
      lag(out_bin_noncovid_death, default = 0) == 1 |
      # LTFU 
      lag(cens_bin_dereg, default = 0) == 1 |
      # Ensure that prior AND current treatment with another antidiabetic treatment is excluded, otherwise my starting point is not metformin vs nothing.
      exp_bin_sulfo_first == 1 | 
      exp_bin_dpp4_first == 1 | 
      exp_bin_tzd_first == 1 | 
      exp_bin_sglt2_first == 1 | 
      exp_bin_glp1_first == 1 | 
      exp_bin_megli_first == 1 | 
      exp_bin_agi_first == 1 | 
      exp_bin_insulin_first == 1 | 
      # prior treatment with metformin
      lag(treat, default = 0) == 1 |
      # prior very high HbA1c (above 75 mmol/mol) in prior 2 years (or also in CURRENT treatment interval?, i.e. lag or no lag?)
      lag(elig_bin_hba1c, default = 0) == 1 |
      # metfin interaction in prior 14-30 days ONLY => use ELD
      lag(eld_elig_bin_metfin_interaction, default = 0) == 1 # remember the difference between eld_elig_ and other elig_ is that elig_bin are filled up with ==1 after event happened until end of follow-up
    ~ 0L,
    TRUE ~ 1L)) %>% 
  ungroup()


# df_months_eld %>%
#   dplyr::select(patient_id, elig_date_t2dm, stop_date, start_date_month, end_date_month, month,
#                 elig, out_date_covid_death, out_date_severecovid, out_date_noncovid_death, cens_date_dereg,
#                 treat, 
#                 elig_bin_metfin_allergy_first, elig_bin_ckd_45_first, elig_bin_liver_cirrhosis_first, 
#                 out_bin_covid, out_bin_longcovid, out_bin_covid_death, out_bin_noncovid_death, cens_bin_dereg, exp_bin_sulfo_first, exp_bin_dpp4_first,
#                 exp_bin_tzd_first, exp_bin_sglt2_first, exp_bin_glp1_first, exp_bin_megli_first, exp_bin_agi_first, exp_bin_insulin_first, elig_bin_hba1c, eld_elig_bin_metfin_interaction,
#                 cov_cat_sex, cov_num_age, everything()) %>%
#   # dplyr::filter(!is.na(out_date_severecovid)) %>%
#   # dplyr::filter(!is.na(out_date_noncovid_death)) %>%
#   # dplyr::filter(!is.na(cens_date_dereg)) %>%
#   # dplyr::filter(is.na(censor)) %>% # should be empty
#   View()

# Assign outcome, censoring and competing events -------------------------
print('Assign outcome, censoring and competing events')
# Assign columns called "outcome", "censor", "comp_event" as following:

## a) If outcome is the row/person-interval event: 
### i) assign to the PREVIOUS row/person-interval: outcome=1, censor=0, comp_event=0
### ii) assign to CURRENT row/person-interval: outcome=NA, censor=NA, comp_event=NA

## b) If censoring event is the row/person-interval event: 
### i) assign to the PREVIOUS row/person-interval: outcome=NA, censor=1, comp_event=NA
### ii) assign to CURRENT row/person-interval: outcome=NA, censor=NA, comp_event=NA

## c) If competing event is the row/person-interval event: 
### i) assign to the PREVIOUS row/person-interval: outcome=NA, censor=0, comp_event=1
### ii) assign to CURRENT row/person-interval: outcome=NA, censor=NA, comp_event=NA

# 3) if outcome and censoring (i.e. LTFU) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order below.
# 4) if outcome and competing event (i.e. non-covid death) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order below.
# 5) if censoring and competing event event have the EXACT same date, then pick the competing event as the defining event. Ensured by case_when() order below.

# 6) Delete all rows/person-intervals that happen after the outcome

# Currently, only calendar month and calendar week are implemented

# Dataset for PLR that stops at either out_date_severecovid (i.e. either death or hosp), out_date_noncovid_death, or cens_date_dereg
df_months_severecovid <- fn_add_and_shift_out_comp_cens_events(df_months_eld,
                                                               outcome_date_variable = "out_date_severecovid",
                                                               comp_date_variable = "out_date_noncovid_death",
                                                               censor_date_variable = "cens_date_dereg",
                                                               studyend_date,
                                                               interval_type = "month")

# Dataset for PLR that stops at either out_date_covid_death, out_date_noncovid_death, or cens_date_dereg
df_months_covid_death <- fn_add_and_shift_out_comp_cens_events(df_months_eld,
                                                          outcome_date_variable = "out_date_covid_death",
                                                          comp_date_variable = "out_date_noncovid_death",
                                                          censor_date_variable = "cens_date_dereg",
                                                          studyend_date,
                                                          interval_type = "month")

# Dataset for PLR that stops at either out_date_covid (includes covid deaths, hosp, diagnoses, tests), out_date_noncovid_death, or cens_date_dereg
df_months_covid_event <- fn_add_and_shift_out_comp_cens_events(df_months_eld,
                                                               outcome_date_variable = "out_date_covid",
                                                               comp_date_variable = "out_date_noncovid_death",
                                                               censor_date_variable = "cens_date_dereg",
                                                               studyend_date,
                                                               interval_type = "month")

# Dataset for PLR that stops at either out_date_longcovid_virfat (first time), out_date_death, or cens_date_dereg
df_months_longcovid <- fn_add_and_shift_out_comp_cens_events(df_months_eld,
                                                               outcome_date_variable = "out_date_longcovid_virfat",
                                                               comp_date_variable = "out_date_death", # all-cause mortality as competing event for this analysis
                                                               censor_date_variable = "cens_date_dereg",
                                                               studyend_date,
                                                               interval_type = "month")

# df_months_longcovid %>%
#   dplyr::select(patient_id, elig_date_t2dm, stop_date, start_date_month, end_date_month, month,
#                 elig, out_date_covid, out_bin_covid, out_date_longcovid_virfat, out_bin_longcovid_virfat, 
#                 out_date_covid_death, out_date_severecovid, out_date_noncovid_death, cens_date_dereg,
#                 outcome, comp_event, censor, elig, treat,
#                 # is_event_interval, is_outcome_event, is_comp_event, is_cens_event, row_num, first_event_row,
#                 cov_cat_sex, cov_num_age, everything()) %>%
#   dplyr::filter(!is.na(out_date_longcovid_virfat)) %>%
#   # dplyr::filter(!is.na(out_date_death)) %>%
#   # dplyr::filter(!is.na(cens_date_dereg)) %>%
#   # dplyr::filter(is.na(censor)) %>% # should be empty
#   View()


# Re-apply eligibility for certain analysis sets --------------------------
print('Re-apply eligibility for certain analysis sets')

df_months_longcovid <- df_months_longcovid %>% 
  group_by(patient_id) %>%
  mutate(elig = case_when(
    lag(out_bin_longcovid_virfat, default = 0) == 1 ~ 0L, # includes in addition to the long covid codes also viral fatigue codes
    TRUE ~ elig)) %>% 
  ungroup()


# Save the datasets ------------------------------------------------------
arrow::write_feather(df_months_severecovid, here::here("output", "data", "df_months_severecovid.arrow"))
arrow::write_feather(df_months_covid_death, here::here("output", "data", "df_months_covid_death.arrow"))
arrow::write_feather(df_months_covid_event, here::here("output", "data", "df_months_covid_event.arrow"))
arrow::write_feather(df_months_longcovid, here::here("output", "data", "df_months_longcovid.arrow"))
