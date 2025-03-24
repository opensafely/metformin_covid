####
## This script does the following:
# 1. Import processed data
# 2. Runs the cox model for the primary analysis
# 3. Saves all output
####

# Import libraries and functions ------------------------------------------
library(arrow)
library(here)
library(tidyverse)
library(lubridate)
library(splines)
library(survival) # survival/TTE analyses
library(ggfortify) # autoplot
library(gtsummary) # tbl_regression

# Create directories for output -------------------------------------------
fs::dir_create(here::here("output", "te"))

# Import the data ---------------------------------------------------------
df <- read_feather(here("output", "data", "data_processed.arrow"))

# Import dates ------------------------------------------------------------
source(here::here("analysis", "metadates.R"))
# Convert the meta-dates into Date objects
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")

# Add splines ------------------------------------------------------------- ## define in data_process
# Compute knot locations based on percentiles, according to study protocol
age_knots <- quantile(df$cov_num_age, probs = c(0.10, 0.50, 0.90))
df <- df %>%
  mutate(cov_num_age_spline = ns(cov_num_age, knots = age_knots))

# Create time-to-event variables for Cox regression ----------------------- ## define in data_process
df <- df %>% 
  filter(qa_date_of_death > elig_date_t2dm | is.na(qa_date_of_death)) # add this to data_process

df <- df %>% 
  mutate(landmark_date = elig_date_t2dm + days(183)) # add this to data_process

df <- df %>% 
  mutate(
    out_bin_severecovid_afterlandmark = case_when(!is.na(out_date_covid19_severe)
                                                  & out_date_covid19_severe > elig_date_t2dm + days(183) ~ TRUE,
                                                  TRUE ~ FALSE),
    out_date_severecovid_afterlandmark = case_when(out_bin_severecovid_afterlandmark == TRUE ~ out_date_covid19_severe, 
                                     TRUE ~ as.Date(NA)),
    out_bin_death_afterlandmark = case_when(!is.na(qa_date_of_death)
                                            & qa_date_of_death > elig_date_t2dm + days(183) ~ TRUE,
                                            TRUE ~ FALSE),
    out_date_death_afterlandmark = case_when(out_bin_death_afterlandmark == TRUE ~ qa_date_of_death, 
                                                   TRUE ~ as.Date(NA)),
    out_bin_ltfu_afterlandmark = case_when(!is.na(cens_date_dereg)
                                                  & cens_date_dereg > elig_date_t2dm + days(183) ~ TRUE,
                                                  TRUE ~ FALSE),
    out_date_ltfu_afterlandmark = case_when(out_bin_ltfu_afterlandmark == TRUE ~ cens_date_dereg, 
                                             TRUE ~ as.Date(NA))
  ) %>% 
  mutate(
    cox_date_severecovid = pmin(out_date_severecovid_afterlandmark, 
                          out_date_death_afterlandmark,
                          out_date_ltfu_afterlandmark,
                          studyend_date,
                          na.rm = TRUE),
    cox_cat_severecovid = case_when(
      # pt should not have both noncovid and covid death
      cox_date_severecovid == out_date_severecovid_afterlandmark ~ "covid_death_hosp",
      cox_date_severecovid == out_date_death_afterlandmark ~ "noncovid_death",
      cox_date_severecovid == out_date_ltfu_afterlandmark ~ "ltfu",
      TRUE ~ "none"
    ),
    cox_tt_severecovid = difftime(cox_date_severecovid,
                                  elig_date_t2dm + days(183), # count from landmark! ## define in data_process
                                  units = "days") %>% as.numeric(),
    cox_bin_severecovid = case_when(cox_cat_severecovid %in% c("noncovid_death", "ltfu", "none") ~ 0,
                                    cox_cat_severecovid == "covid_death_hosp" ~ 1,
                                    TRUE ~ NA_real_),
    cox_date_severecovid_censor = case_when(cox_bin_severecovid == 0 ~ cox_date_severecovid,
                                          TRUE ~ as.Date(NA))
  )

# Add treatment variable -------------------------------------------------- ## define in data_process
# Set format and reference to ensure proper logistic regression modeling
df$exp_bin_treat <- factor(df$exp_bin_treat, 
                           levels = c(2, 1), # Reference first
                           labels = c("nothing", "metformin"))

# Add covariates ----------------------------------------------------------
covariate_names <- names(df) %>%
    grep("^cov_", ., value = TRUE) %>% 
    # exclude those not needed in the model: 
    ## cov_cat_region covers for cov_cat_stp, 
    ## cov_bin_obesity covers for cov_num_bmi & cov_cat_bmi_groups,
    ## cov_cat_hba1c_mmol_mol covers cov_num_hba1c_mmol_mol
    ## cov_cat_tc_hdl_ratio covers cov_num_tc_hdl_ratio
    ## cov_num_age_spline covers cov_cat_age and cov_num_age
    ## take out cov_cat_region since it's used as a stratification variable instead
    setdiff(c("cov_cat_stp", "cov_num_bmi", "cov_cat_bmi_groups", "cov_num_hba1c_mmol_mol", "cov_num_tc_hdl_ratio", "cov_cat_age", "cov_num_age", "cov_cat_region")) 
print(covariate_names)

# Checks before proceeding to Cox model -----------------------------------
if (any(df$cox_tt_severecovid < 0, na.rm = TRUE)) {
  stop("Error: Some values in cox_tt_severecovid are negative.")
} else {
  print("Check passed: All values in cox_tt_severecovid are non-negative.")
}

# Cox model ---------------------------------------------------------------
cox_formula_severecovid <- as.formula(paste("Surv(cox_tt_severecovid, cox_bin_severecovid) ~ exp_bin_treat +", 
                                            paste(covariate_names, collapse = " + "), "+ strata(cov_cat_region)"))
cox_model_severecovid <- coxph(cox_formula_severecovid, data = df)
cox_severecovid <- tbl_regression(cox_model_severecovid, exp = TRUE)
cox_severecovid_df <- cox_severecovid %>% 
  as_tibble()

# Save output -------------------------------------------------------------
# cox model table
write.csv(cox_severecovid_df, file = here::here("output", "te", "cox_severecovid.csv"))