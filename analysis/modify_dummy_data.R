####
## This script does the following:
# 1. Modify the dummy data to ensure the covariates make sense
# 2. Modify the dummy data to ensure the order of the dates make sense
####

# Import libraries and functions ------------------------------------------
library(dplyr)

# Set seed ----------------------------------------------------------------
set.seed(36)

# Modify covariates -------------------------------------------------------
# Define the covariates included in the PS model
covariate_names <- names(data_processed) %>%
  grep("^cov_", ., value = TRUE)
# print(covariate_names)

# age has some odd values, since no quality/completeness applied to dummy data
data_processed <- data_processed %>%
  mutate(cov_num_age = round(pmax(pmin(rnorm(n(), mean = 50, sd = 20), 110), 18)))

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
# Define probabilities for each category, more in the lower and middle
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

# cov_bin_pcos only contains FALSE
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

## Modify covariates because no quality criteria nor completeness criteria were applied to dummy data 
# cov_cat_sex has NA in dummy data, not in real data (due completeness criteria)
cov_cat_sex_categories <- c("Female", "Male")
cov_cat_sex_prob <- c(0.4, 0.6) # see prelim data
data_processed <- data_processed %>%
  mutate(cov_cat_sex = sample(cov_cat_sex_categories, n(), replace = TRUE, prob = cov_cat_sex_prob)) %>%
  mutate(cov_cat_sex = factor(cov_cat_sex, levels = cov_cat_sex_categories))

# cov_cat_deprivation_5 has NA/unknown in dummy data, not in real data (due completeness criteria)
cov_cat_deprivation_5_categories <- c("1 (most deprived)", "2", "3", "4", "5 (least deprived)")
cov_cat_deprivation_5_prob <- c(0.3, 0.2, 0.2, 0.15, 0.15) # see prelim data
data_processed <- data_processed %>%
  mutate(cov_cat_deprivation_5 = sample(cov_cat_deprivation_5_categories, n(), replace = TRUE, prob = cov_cat_deprivation_5_prob)) %>%
  mutate(cov_cat_deprivation_5 = factor(cov_cat_deprivation_5, levels = cov_cat_deprivation_5_categories))

# cov_cat_region has NA in dummy data, not in real data (due completeness criteria)
cov_cat_region_categories <- c("North East and Yorkshire", "North West", "Midlands", "East of England", "London", "South East", "South West")
cov_cat_region_prob <- c(0.2, 0.1, 0.2, 0.25, 0.05, 0.05, 0.15) # see prelim data
data_processed <- data_processed %>%
  mutate(cov_cat_region = sample(cov_cat_region_categories, n(), replace = TRUE, prob = cov_cat_region_prob)) %>%
  mutate(cov_cat_region = factor(cov_cat_region, levels = cov_cat_region_categories))

# Increase the number treated ---------------------------------------------
# Set to a 50:50 distribution, randomly set
data_processed <- data_processed %>%
  mutate(exp_bin_treat = sample(c(1, 2), nrow(data_processed), replace = TRUE, prob = c(0.5, 0.5)))

# Increase sample size ----------------------------------------------------
# Sample with replacement to generate more data
# n_new_rows <- 80000
# data_processed <- data_processed[sample(1:nrow(data_processed), n_new_rows, replace = TRUE), ]
