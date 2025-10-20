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
library(survminer) # KM curve
source(here::here("analysis", "covariates.R")) 

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

# Add covariates/confounders ----------------------------------------------
print('Add covariates/confounders')
print(covariates_names)
print(covariates_ex_strat_names)

# Checks before proceeding to Cox model -----------------------------------
if (any(df$cox_tt_severecovid < 0, na.rm = TRUE)) {
  stop("Error: Some values in cox_tt_severecovid are negative.")
} else {
  print("Check passed: All values in cox_tt_severecovid are non-negative.")
}

if (anyNA(df[, c("cox_tt_severecovid", "out_bin_severecovid_afterlandmark", covariates_names)])) {
  warning("Missing values detected in model variables")
}

# Cox model ---------------------------------------------------------------
print('Run Cox model')
cox_formula_severecovid <- as.formula(paste("Surv(cox_tt_severecovid, out_bin_severecovid_afterlandmark) ~ exp_bin_treat +", 
                                            paste(covariates_ex_strat_names, collapse = " + "), "+ strat(strat_cat_region)"))
cox_model_severecovid <- coxph(cox_formula_severecovid, data = df)
cox_severecovid <- tbl_regression(cox_model_severecovid, exp = TRUE)
cox_severecovid_df <- cox_severecovid %>% 
  as_tibble()

# Cumulative Incidence Curve ----------------------------------------------
print('Generate cumulative incidence curve')

# Fit Kaplan-Meier estimates
cuminc_fit <- survfit(Surv(cox_tt_severecovid, out_bin_severecovid_afterlandmark) ~ exp_bin_treat, data = df)

# Plot cumulative incidence (1 - survival)
cuminc_plot <- ggsurvplot(
  cuminc_fit,
  data = df,
  fun = "event",              
  risk.table = TRUE,
  pval = F,
  conf.int = TRUE,
  xlab = "Time (days)",
  ylab = "Cumulative incidence",
  ggtheme = theme_minimal(),
  legend.title = "Treatment",
  legend.labs = c("Nothing", "Metformin")
)

# Save output -------------------------------------------------------------
print('Save output')
# cox model table
write.csv(cox_severecovid_df, file = here::here("output", "te", "cox_scripted", "cox_severecovid.csv"))
# KM curve
ggsave(filename = here::here("output", "te", "cox_scripted", "cuminc_plot.png"), 
       plot = cuminc_plot$plot,
       width = 20, height = 20, units = "cm")
