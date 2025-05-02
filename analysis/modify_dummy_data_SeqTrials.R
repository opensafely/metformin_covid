####
## This script does the following:
# 1. Modify the dummy data to ensure the dummy data runs through quality assurance and completeness criteria
# 2. Modify the dummy data to ensure the covariates make sense
# 3. Modify the dummy data to ensure the order of the dates make sense
# (4. Increase proportion treated)
####

# Import libraries and functions ------------------------------------------
library(dplyr)
source(here::here("analysis", "functions", "fn_dd_exp_out_dates.R"))
source(here::here("analysis", "functions", "fn_dd_cens_dates.R"))


# Import dates ------------------------------------------------------------
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")
mid2018_date <- as.Date(study_dates$mid2018_date, format = "%Y-%m-%d")
pandemicstart_date <- as.Date(study_dates$pandemicstart_date, format = "%Y-%m-%d")


# Set seed ----------------------------------------------------------------
set.seed(36)


# Modify quality assurance critera -----------------------------------------------
# let qa_bin_hrt, qa_bin_cocp and qa_bin_prostate_cancer have very few FALSE to avoid exclusions
data_processed <- data_processed %>%
  # Quality assurance: Pregnancy/birth codes
  mutate(qa_bin_pregnancy = rbernoulli(nrow(.), p = 0.001)) %>%
  # Quality assurance: HRT meds
  mutate(qa_bin_hrt = rbernoulli(nrow(.), p = 0.001)) %>%
  # Quality assurance: COCP meds
  mutate(qa_bin_cocp = rbernoulli(nrow(.), p = 0.001)) %>%
  # Quality assurance: Prostate cancer codes
  mutate(qa_bin_prostate_cancer = rbernoulli(nrow(.), p = 0.001))


# Modify completeness criteria ---------------------------------------------------
# Increase completeness in sex, resample with replacement
sex_categories <- c("Female", "Male", "Unknown")
sex_categories_prob <- c(0.4, 0.5, 0.0)
data_processed <- data_processed %>%
  mutate(cov_cat_sex = as.factor(sample(sex_categories, n(), replace = TRUE, prob = sex_categories_prob))) %>%
  mutate(qa_bin_is_female_or_male = cov_cat_sex %in% c("Female", "Male"))

# Make all adults, and age has some odd values
data_processed <- data_processed %>%
  mutate(cov_num_age = round(pmax(pmin(rnorm(n(), mean = 50, sd = 20), 110), 18))) %>% 
  mutate(qa_bin_was_adult = cov_num_age >= 18 & cov_num_age <= 110)

# Complete deprivation info
cov_cat_deprivation_5_categories <- c("1 (most deprived)", "2", "3", "4", "5 (least deprived)", "Unknown")
cov_cat_deprivation_5_prob <- c(0.3, 0.2, 0.2, 0.15, 0.15, 0.0) # see prelim data
data_processed <- data_processed %>%
  mutate(cov_cat_deprivation_5 = as.factor(sample(cov_cat_deprivation_5_categories, n(), replace = TRUE, prob = cov_cat_deprivation_5_prob))) %>%
  mutate(qa_bin_known_imd = cov_cat_deprivation_5 %in% c("1 (most deprived)", "2", "3", "4", "5 (least deprived)"))

# Complete region info
strat_cat_region_categories <- c("North East and Yorkshire", "North West", "Midlands", "East of England", "London", "South East", "South West", "Unknown")
strat_cat_region_prob <- c(0.2, 0.1, 0.2, 0.25, 0.05, 0.05, 0.15, 0.0) # see prelim data
data_processed <- data_processed %>%
  mutate(strat_cat_region = as.factor(sample(strat_cat_region_categories, n(), replace = TRUE, prob = strat_cat_region_prob)))

# All registered
data_processed <- data_processed %>%
  mutate(qa_bin_was_registered = TRUE)


# Modify covariates -------------------------------------------------------
# HbA1c has zero counts in all groups
HbA1c_categories <- c("below 42", "42-58", "59-75", "above 75", "Unknown")
HbA1c_prob <- c(0.05, 0.35, 0.25, 0.05, 0.3)
data_processed <- data_processed %>%
  mutate(cov_cat_hba1c_mmol_mol_b = as.factor(sample(HbA1c_categories, n(), replace = TRUE, prob = HbA1c_prob)))

# Tot Chol/HDL ratio has zero counts in all groups
tc_hdl_ratio_categories <- c("below 3.5:1", "3.5:1 to 5:1", "above 5:1", "Unknown")
tc_hdl_ratio_prob <- c(0.1, 0.4, 0.2, 0.3)
data_processed <- data_processed %>%
  mutate(cov_cat_tc_hdl_ratio_b = as.factor(sample(tc_hdl_ratio_categories, n(), replace = TRUE, prob = tc_hdl_ratio_prob)))

# cov_bin_pcos only contains 1 TRUE
data_processed <- data_processed %>%
  mutate(cov_bin_pcos = sample(c(TRUE, FALSE), n(), replace = TRUE, prob = c(0.05, 0.95)))

# cov_bin_prediabetes only contains 1 TRUE
data_processed <- data_processed %>%
  mutate(cov_bin_prediabetes = sample(c(TRUE, FALSE), n(), replace = TRUE, prob = c(0.5, 0.5))) # see prelim data

# cov_bin_obesity only contains very few TRUE
data_processed <- data_processed %>%
  mutate(cov_bin_obesity_b = sample(c(TRUE, FALSE), n(), replace = TRUE, prob = c(0.3, 0.7)))


# Modify dates ------------------------------------------------------------
# (1) Ensure more diabetes diagnoses and elig_date_t2dm are in the window (mid2018 to pandemic start)
data_processed <- data_processed %>%
  mutate(elig_date_t2dm = sample(seq(mid2018_date, pandemicstart_date, by = "day"),
                                 n(), replace = TRUE))
 
# # (2) increase the amount of completeness in the main intervention variables first
# data_processed <- data_processed %>%
#   rowwise() %>%
#   mutate(exp_date_metfin_mono_first = ifelse(
#     is.na(exp_date_metfin_mono_first) & !is.na(elig_date_t2dm) & (elig_date_t2dm < studyend_date),
#     sample(seq(elig_date_t2dm + 1, studyend_date, by = "day"), 1),
#     exp_date_metfin_mono_first
#   )) %>%
#   ungroup()
# data_processed <- data_processed %>%
#   rowwise() %>%
#   mutate(exp_date_metfin_first = ifelse(
#     is.na(exp_date_metfin_first) & !is.na(elig_date_t2dm) & (elig_date_t2dm < studyend_date),
#     sample(seq(elig_date_t2dm + 1, studyend_date, by = "day"), 1),
#     exp_date_metfin_first
#   )) %>%
#   ungroup()
# # Now, introduce ~10% missing values randomly in these
# data_processed <- data_processed %>%
#   mutate(exp_date_metfin_mono_first = ifelse(
#     runif(n()) < 0.1,  # 10% probability of NA
#     NA,
#     exp_date_metfin_mono_first
#   ))
# data_processed <- data_processed %>%
#   mutate(exp_date_metfin_first = ifelse(
#     runif(n()) < 0.1,  # 10% probability of NA
#     NA,
#     exp_date_metfin_first
#   ))
# data_processed <- data_processed %>%
#   mutate(exp_date_metfin_mono_first = as.Date(exp_date_metfin_mono_first, origin = "1970-01-01")) %>% 
#   mutate(exp_date_metfin_first = as.Date(exp_date_metfin_first, origin = "1970-01-01"))

# (3) Ensure all exposure variables are between baseline date (elig_date_t2dm) and studyend date
## If the date is NA, leave it as NA.
## If the date is before elig_date_t2dm, replace it with a random valid date.
## If the date is valid, keep it unchanged.

# exp_vars <- c(
#   "exp_date_metfin_first", "exp_date_metfin_mono_first",
#   "exp_date_sulfo_first",
#   "exp_date_dpp4_first", "exp_date_tzd_first", "exp_date_sglt2_first",
#   "exp_date_glp1_first", "exp_date_megli_first", "exp_date_agi_first",
#   "exp_date_insulin_first")
# data_processed <- data_processed %>%
#   rowwise() %>%
#   mutate(across(all_of(exp_vars), ~ fn_dd_exp_out_dates(.x, elig_date_t2dm, studyend_date))) %>%
#   ungroup() %>%
#   mutate(across(all_of(exp_vars), ~ as.Date(.x, origin = "1970-01-01")))

# (4) Ensure all outcome dates are between landmark date and studyend date
## If the date is NA, leave it as NA.
## If the date is before landmark, replace it with a random valid date.
## If the date is valid, keep it unchanged.

# out_vars <- c(
#   "out_date_covid_hosp", "out_date_covid_death",
#   "out_date_severecovid", "out_date_covid", "out_date_longcovid",
#   "out_date_virfat", "out_date_longcovid_virfat", "qa_date_of_death")
# data_processed <- data_processed %>%
#   rowwise() %>%
#   mutate(across(all_of(out_vars), ~ fn_dd_exp_out_dates(.x, elig_date_t2dm + days(183), studyend_date))) %>%
#   ungroup() %>%
#   mutate(across(all_of(out_vars), ~ as.Date(.x, origin = "1970-01-01")))

# (5) Ensure all censoring dates are between baseline date (elig_date_t2dm) and studyend date
## cens_date_dereg is completely missing in the current dummy data, so replace all
## but impute only 5% of cens_date_dereg date variables

# data_processed$cens_date_dereg <- sapply(data_processed$elig_date_t2dm, function(baseline_date) {
#   fn_dd_cens_dates(baseline_date, studyend_date)
# })
# data_processed$cens_date_dereg <- as.Date(data_processed$cens_date_dereg, origin = "1970-01-01")