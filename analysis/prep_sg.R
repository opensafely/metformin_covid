####
## This script does the following:
# 1. Import the processed dataset of all eligible participants
# 2. Create subsets to run subgroup analyses using the reusable action cox-ipw
# 3. Save all subsets
####


# Import libraries and user functions -------------------------------------
print('Import libraries and functions')
library('tidyverse')
library('here')
library('arrow')
source(here::here("analysis", "functions", "utility.R")) # midpoint rounding


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data", "sg"))
fs::dir_create(here::here("output", "te", "sg_desc"))


# Import the data ---------------------------------------------------------
print('Import the data')
df <- read_feather(here("output", "data", "data_processed.arrow"))


# Import dates ------------------------------------------------------------
print('Import dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# Subset the data ---------------------------------------------------------
print('Subset the data')
df_below60 <- df %>%
  filter(cov_num_age < 60)
df_60orabove <- df %>%
  filter(cov_num_age >= 60)
df_female <- df %>%
  filter(cov_cat_sex == "Female")
df_male <- df %>%
  filter(cov_cat_sex == "Male")
df_white <- df %>%
  filter(cov_cat_ethnicity == "White")
df_nonwhite <- df %>%
  filter(cov_cat_ethnicity != "White") # includes unknown
df_imd1 <- df %>%
  filter(cov_cat_deprivation_5 == "1 (most deprived)")
df_nonimd1 <- df %>%
  filter(cov_cat_deprivation_5 != "1 (most deprived)") # no unknown by definition
df_obese <- df %>%
  filter(cov_bin_obesity == TRUE)
df_nonobese <- df %>%
  filter(cov_bin_obesity == FALSE) # includes unknown
df_HbA1c59orabove <- df %>%
  filter(cov_cat_hba1c_mmol_mol == "59-75") # everyone above 75 was part excluded
df_belowHbA1c59 <- df %>%
  filter(cov_cat_hba1c_mmol_mol != "59-75") # includes unknown - probably ok, but rediscuss


# Re-define cox_stop in each subset ---------------------------------------
print('Re-define cox_stop in each subset')
# variable already exists in dataset, will be overwritten with uptodate info
redefine_cox_stop <- function(df) {
  df %>% 
    mutate(cox_date_severecovid = pmin(out_date_severecovid_afterlandmark, 
                                       out_date_noncoviddeath_afterlandmark,
                                       cens_date_ltfu_afterlandmark,
                                       max_fup_date,
                                       studyend_date,
                                       na.rm = TRUE))
}

df_list <- list(
  df_below60 = df_below60,
  df_60orabove = df_60orabove,
  df_female = df_female,
  df_male = df_male,
  df_white = df_white,
  df_nonwhite = df_nonwhite,
  df_imd1 = df_imd1,
  df_nonimd1 = df_nonimd1,
  df_obese = df_obese,
  df_nonobese = df_nonobese,
  df_HbA1c59orabove = df_HbA1c59orabove,
  df_belowHbA1c59 = df_belowHbA1c59
)

df_list <- lapply(df_list, redefine_cox_stop)
 
# Assign back to environment (optional)
# list2env(df_list, envir = .GlobalEnv)



# Extract events by subgroup from each subset ------------------------------
print('Extract events by subgroup from each subset')



# Save subsets and descriptive (events by subgroup) -----------------------
# Subsets
lapply(names(df_list), function(name) {
  arrow::write_feather(
    df_list[[name]],
    here::here("output", "data", "sg", paste0(name, ".arrow"))
  )
})


# Events by subgroup (midpoint rounded)
# write.csv(n_events_by_sg_midpoint6, file = here::here("output", "te", "sg_desc", "n_events_by_sg_midpoint6.csv"))
