####
## This script does the following:
# 1. Modify the dummy data to ensure the dummy data runs through quality assurance and completeness criteria
# 2. Modify the dummy data to ensure the covariates make sense
# 3. Modify the dummy data to ensure the order of the dates make sense
# (4. Increase proportion treated)
####

# Import libraries and functions ------------------------------------------
library(dplyr)

# Import dates ------------------------------------------------------------
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")
mid2018_date <- as.Date(study_dates$mid2018_date, format = "%Y-%m-%d")

# Set seed ----------------------------------------------------------------
set.seed(36)

# Quality assurance critera -----------------------------------------------
# let qa_bin_hrt, qa_bin_cocp and qa_bin_prostate_cancer be FALSE to avoid exclusions
data_processed <- data_processed %>%
  mutate(qa_bin_pregnancy = FALSE,
         qa_bin_hrt = FALSE,
         qa_bin_cocp = FALSE,
         qa_bin_prostate_cancer = FALSE)

# Completeness criteria ---------------------------------------------------
# Increase completeness in sex, resample with replacement
sex_categories <- c("Female", "Male", "Unknown")
sex_categories_prob <- c(0.4, 0.6, 0.0)
data_processed <- data_processed %>%
  mutate(qa_bin_is_female_or_male = TRUE) %>% 
  mutate(cov_cat_sex = sample(sex_categories, n(), replace = TRUE, prob = sex_categories_prob)) %>%
  mutate(cov_cat_sex = factor(cov_cat_sex, levels = sex_categories))

# Make all adults, and age has some odd values
data_processed <- data_processed %>%
  mutate(qa_bin_was_adult = TRUE) %>% 
  mutate(cov_num_age = round(pmax(pmin(rnorm(n(), mean = 50, sd = 20), 110), 18)))

# Complete deprivation info
cov_cat_deprivation_5_categories <- c("1 (most deprived)", "2", "3", "4", "5 (least deprived)", "unknown")
cov_cat_deprivation_5_prob <- c(0.3, 0.2, 0.2, 0.15, 0.15, 0.0) # see prelim data
data_processed <- data_processed %>%
  mutate(qa_bin_known_imd = TRUE) %>% 
  mutate(cov_cat_deprivation_5 = sample(cov_cat_deprivation_5_categories, n(), replace = TRUE, prob = cov_cat_deprivation_5_prob)) %>%
  mutate(cov_cat_deprivation_5 = factor(cov_cat_deprivation_5, levels = cov_cat_deprivation_5_categories))

# Complete region info
cov_cat_region_categories <- c("North East and Yorkshire", "North West", "Midlands", "East of England", "London", "South East", "South West")
cov_cat_region_prob <- c(0.2, 0.1, 0.2, 0.25, 0.05, 0.05, 0.15) # see prelim data
data_processed <- data_processed %>%
  mutate(cov_cat_region = sample(cov_cat_region_categories, n(), replace = TRUE, prob = cov_cat_region_prob)) %>%
  mutate(cov_cat_region = factor(cov_cat_region, levels = cov_cat_region_categories))

# All registered
data_processed <- data_processed %>%
  mutate(qa_bin_was_registered = TRUE)

# Modify covariates -------------------------------------------------------
# Define the covariates included in the PS model
# covariate_names <- names(data_processed) %>%
#   grep("^cov_", ., value = TRUE)
# print(covariate_names)

# Check for factors and logicals with only one level, for e.g. PS model on dummy data
# sapply(data_processed %>% select(all_of(covariate_names)), function(x) {
#   if (is.factor(x) || is.logical(x)) {
#     length(unique(x))
#   } else {
#     NA 
#   }
# })

# HbA1c has zero counts in all groups
HbA1c_categories <- c("below 42", "42-58", "59-75", "above 75", "Unknown")
HbA1c_prob <- c(0.05, 0.35, 0.25, 0.05, 0.3)
data_processed <- data_processed %>%
  mutate(cov_cat_hba1c_mmol_mol = sample(HbA1c_categories, n(), replace = TRUE, prob = HbA1c_prob)) %>%
  mutate(cov_cat_hba1c_mmol_mol = factor(cov_cat_hba1c_mmol_mol, levels = HbA1c_categories))

# Tot Chol/HDL ratio has zero counts in all groups
tc_hdl_ratio_categories <- c("below 3.5:1", "3.5:1 to 5:1", "above 5:1", "Unknown")
tc_hdl_ratio_prob <- c(0.1, 0.4, 0.2, 0.3)
data_processed <- data_processed %>%
  mutate(cov_cat_tc_hdl_ratio = sample(tc_hdl_ratio_categories, n(), replace = TRUE, prob = tc_hdl_ratio_prob)) %>%
  mutate(cov_cat_tc_hdl_ratio = factor(cov_cat_tc_hdl_ratio, levels = tc_hdl_ratio_categories))

# cov_bin_pcos only contains 1 TRUE
data_processed <- data_processed %>%
  mutate(cov_bin_pcos = sample(c(TRUE, FALSE), n(), replace = TRUE, prob = c(0.05, 0.95)))

# cov_bin_prediabetes only contains 1 TRUE
data_processed <- data_processed %>%
  mutate(cov_bin_prediabetes = sample(c(TRUE, FALSE), n(), replace = TRUE, prob = c(0.5, 0.5))) # see prelim data

# cov_bin_obesity only contains 1 value in 1 category
cov_bin_obesity_categories <- c("Obese (>30)", "Not Obese (<=30)", "Unknown")
cov_bin_obesity_prob <- c(0.3, 0.5, 0.2)
data_processed <- data_processed %>%
  mutate(cov_bin_obesity = sample(cov_bin_obesity_categories, n(), replace = TRUE, prob = cov_bin_obesity_prob)) %>%
  mutate(cov_bin_obesity = factor(cov_bin_obesity, levels = cov_bin_obesity_categories))

# Modify dates ------------------------------------------------------------
# 1) Ensure they all have a diabetes and elig_date_t2dm are in the window (mid2018 to mid2019)
# Define the range (mid2018 is taken from metadates.R)
end_date_elig_date_t2dm <- as.Date("2019-07-01")
# Replace missing elig_date_t2dm values with random dates
data_processed <- data_processed %>%
  mutate(elig_date_t2dm = sample(seq(mid2018_date, end_date_elig_date_t2dm, by = "day"), 
                                 n(), replace = TRUE))

# 2) Ensure all exposure and outcome dates are after baseline date (elig_date_t2dm)
## If the date is NA, leave it as NA.
## If the date is before elig_date_t2dm, replace it with a random valid date.
## If the date is valid, keep it unchanged.
adjust_dates <- function(date_var, elig_date, studyend_date) {
  ifelse(is.na(date_var), NA,
    ifelse(
      date_var <= elig_date, 
      if (as.Date(elig_date) + 1 <= as.Date(studyend_date)) {
        sample(seq(as.Date(elig_date) + 1, as.Date(studyend_date), by = "day"), 1)
      } else {
        NA # Avoid errors if range is invalid
      },
      date_var # Keep valid existing dates
    ))
  }
date_vars <- c(
  "exp_date_metfin_first", "exp_date_metfin_mono_first", "exp_date_sulfo_first",
  "exp_date_dpp4_first", "exp_date_tzd_first", "exp_date_sglt2_first",
  "exp_date_glp1_first", "exp_date_megli_first", "exp_date_agi_first",
  "exp_date_insulin_first", "out_date_covid19_hosp", "out_date_covid19_death",
  "out_date_covid19_severe", "out_date_covid19", "out_date_long_covid19",
  "out_date_viral_fatigue", "out_date_long_fatigue", "cens_date_dereg")

data_processed <- data_processed %>%
  rowwise() %>%
  mutate(across(all_of(date_vars), ~ adjust_dates(.x, elig_date_t2dm, studyend_date))) %>%
  ungroup() %>%
  mutate(across(all_of(date_vars), ~ as.Date(.x, origin = "1970-01-01")))

# Treat cens_date_dereg separately since it is completely missing; impute 5% of values
generate_cens_date_dereg <- function(elig_date, studyend_date) {
  if (runif(1) <= 0.05) {  # 5% chance to impute a random date
    if (as.Date(elig_date) + 1 <= as.Date(studyend_date)) {
      return(sample(seq(as.Date(elig_date) + 1, as.Date(studyend_date), by = "day"), 1))
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}
data_processed$cens_date_dereg <- sapply(data_processed$elig_date_t2dm, function(elig_date) {
  generate_cens_date_dereg(elig_date, studyend_date)
})
data_processed$cens_date_dereg <- as.Date(data_processed$cens_date_dereg, origin = "1970-01-01")

# Increase the number treated ---------------------------------------------
# Set to a 50:50 distribution, as per real data

# Increase sample size ----------------------------------------------------
# Sample with replacement to generate more data
# n_new_rows <- 80000
# data_processed <- data_processed[sample(1:nrow(data_processed), n_new_rows, replace = TRUE), ]