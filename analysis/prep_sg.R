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
source(here::here("analysis", "functions", "fn_process_and_save_sg.R"))


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "data", "sg"))


# Import the data ---------------------------------------------------------
print('Import the data')
df <- read_feather(here("output", "data", "data_processed.arrow"))


# Import dates ------------------------------------------------------------
print('Import dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# Define the subgroups ----------------------------------------------------
print('Define the subgroups')
subgroups <- list(
  df_below60       = quote(cov_num_age < 60),
  df_60orabove     = quote(cov_num_age >= 60),
  df_female        = quote(cov_cat_sex == "Female"),
  df_male          = quote(cov_cat_sex == "Male"),
  df_white         = quote(cov_cat_ethnicity == "White"),
  df_nonwhite      = quote(cov_cat_ethnicity != "White"),
  df_imd1          = quote(cov_cat_deprivation_5 == "1 (most deprived)"),
  df_nonimd1       = quote(cov_cat_deprivation_5 != "1 (most deprived)"),
  df_obese         = quote(cov_bin_obesity == TRUE),
  df_nonobese      = quote(cov_bin_obesity == FALSE),
  df_HbA1c59orabove = quote(cov_cat_hba1c_mmol_mol == "59-75"),
  df_belowHbA1c59  = quote(cov_cat_hba1c_mmol_mol != "59-75")
)


# Subset the data and redefine cox end date, all without saving each subset in the environment -----------
print('Run the loop & function to subset the data and redefine the cox end date in each subset')
for (name in names(subgroups)) {
  fn_process_and_save_sg(
    name = name,
    subgroup = subgroups[[name]],
    df = df,
    studyend_date = studyend_date
  )
}