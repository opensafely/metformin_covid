####
## This script does the following:
# 1. Import processed data
# 2. Runs the cox model for the primary analysis
# 3. Saves all output
####

# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(lubridate)
library(splines)
library(rms) # strat
library(survival) # survival/TTE analyses
library(ggfortify) # autoplot
library(gtsummary) # tbl_regression

# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "te", "cox_scripted"))

# Import the data ---------------------------------------------------------
print('Import the data')
df <- read_feather(here("output", "data", "data_processed.arrow"))

# Import dates ------------------------------------------------------------
print('Import the dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")

# Add splines -------------------------------------------------------------
print('Add splines')
# Compute knot locations based on percentiles, according to study protocol
age_knots <- quantile(df$cov_num_age, probs = c(0.10, 0.50, 0.90))
df <- df %>%
  mutate(cov_num_age_spline = ns(cov_num_age, knots = age_knots))

# Add treatment variable --------------------------------------------------
print('Add treatment variable')
# Set format and reference to ensure proper logistic regression modeling
df$exp_bin_treat <- factor(df$exp_bin_treat, 
                           levels = c(0, 1), # Reference first
                           labels = c("nothing", "metformin"))

# Add covariates ----------------------------------------------------------
print('Add covariates')
covariate_names <- names(df) %>%
    grep("^cov_", ., value = TRUE) %>% 
    # exclude those not needed in the model:
    ## cov_bin_obesity covers for cov_num_bmi & cov_cat_bmi_groups,
    ## cov_cat_hba1c_mmol_mol covers cov_num_hba1c_mmol_mol
    ## cov_cat_tc_hdl_ratio covers cov_num_tc_hdl_ratio
    ## cov_num_age_spline covers cov_cat_age and cov_num_age
    setdiff(c("cov_num_bmi", "cov_cat_bmi_groups", "cov_num_hba1c_mmol_mol"
              , "cov_num_tc_hdl_ratio", "cov_cat_age", "cov_num_age")) 
# print(covariate_names)

# Checks before proceeding to Cox model -----------------------------------
if (any(df$cox_tt_severecovid < 0, na.rm = TRUE)) {
  stop("Error: Some values in cox_tt_severecovid are negative.")
} else {
  print("Check passed: All values in cox_tt_severecovid are non-negative.")
}

# Cox model ---------------------------------------------------------------
print('Run Cox model')
cox_formula_severecovid <- as.formula(paste("Surv(cox_tt_severecovid, out_bin_severecovid_afterlandmark) ~ exp_bin_treat +", 
                                            paste(covariate_names, collapse = " + "), "+ strat(strat_cat_region)"))
cox_model_severecovid <- coxph(cox_formula_severecovid, data = df)
cox_severecovid <- tbl_regression(cox_model_severecovid, exp = TRUE)
cox_severecovid_df <- cox_severecovid %>% 
  as_tibble()

# Save output -------------------------------------------------------------
print('Save output')
# cox model table
write.csv(cox_severecovid_df, file = here::here("output", "te", "cox_scripted", "cox_severecovid.csv"))