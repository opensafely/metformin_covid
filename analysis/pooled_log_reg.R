####
## This script does the following:
# 1. Import processed data
# 2. Runs the pooled logistic regression
# 3. Saves all output
####

# Import libraries and functions ------------------------------------------
library(arrow)
library(here)
library(tidyverse)
library(lubridate)
library(splines)

library(sandwich) # for robust standard errors
library(lmtest) # For hypothesis testing with robust SEs
library(car) # For deltaMethod

library(purrr) # for data wrangling
library(boot)

library(ggplot2)
# library(speedglm) # not available in OpenSAFELY

source(here::here("analysis", "functions", "fn_expand_intervals.R"))

# Create directories for output -------------------------------------------
fs::dir_create(here::here("output", "te", "pooled_log_reg"))

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

# Create the data setup needed for pooled log reg -------------------------
df <- df %>% 
  filter(qa_date_of_death > elig_date_t2dm | is.na(qa_date_of_death)) # add this to data_process

df <- df %>% 
  mutate(landmark_date = elig_date_t2dm + days(183)) # add this to data_process

df <- df %>% # add this to data_process
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

# Expand the dataset into intervals and assign primary outcome (out_date_severecovid_afterlandmark), while censoring at the other pre-defined event date
stop_date_columns <- c("out_date_severecovid_afterlandmark", "out_date_death_afterlandmark", "out_date_ltfu_afterlandmark")
outcome_date_variable <- "out_date_severecovid_afterlandmark"

# Apply the function, choose either weeks or months
df_long_weeks <- fn_expand_intervals(df, studyend_date, stop_date_columns, outcome_date_variable, interval_type = "week")
df_long_months <- fn_expand_intervals(df, studyend_date, stop_date_columns, outcome_date_variable, interval_type = "month")

# double-check
# df_long_months %>%
#   select(patient_id, elig_date_t2dm, landmark_date, out_date_severecovid_afterlandmark, out_date_death_afterlandmark,
#          out_date_ltfu_afterlandmark, stop_date, start_date_month, month, outcome,
#          qa_date_of_death, cov_cat_sex, cov_num_age, cov_cat_deprivation_5) %>%
#   View()

# Define treatment variable -------------------------------------------------- ## define in data_process
# Keep it at 1 and 0 for model below
df <- df %>% 
  mutate(exp_bin_treat = case_when(exp_bin_treat == 2 ~ 0,
                                   exp_bin_treat == 1 ~ 1))

# Define covariates ----------------------------------------------------------
covariate_names <- names(df) %>%
  grep("^cov_", ., value = TRUE) %>% 
  # exclude those not needed in the model: 
  ## cov_cat_region covers for cov_cat_stp, 
  ## cov_bin_obesity covers for cov_num_bmi & cov_cat_bmi_groups,
  ## cov_cat_hba1c_mmol_mol covers cov_num_hba1c_mmol_mol
  ## cov_cat_tc_hdl_ratio covers cov_num_tc_hdl_ratio
  ## cov_num_age_spline covers cov_cat_age and cov_num_age
  ## CAVE: Keep cov_cat_region since it's used as a stratification variable instead
  setdiff(c("cov_cat_stp", "cov_num_bmi", "cov_cat_bmi_groups", "cov_num_hba1c_mmol_mol", "cov_num_tc_hdl_ratio"
            , "cov_cat_age", "cov_num_age")) 
# print(covariate_names)

# Pooled logistic regression ----------------------------------------------
### Pooled logistic regression
## Pooled logistic regression models naturally allow for the parametric estimation of risks, and thus risk differences
## and risk ratios. Also, these models can be specified such that the effect of the treatment can vary over time, 
## as opposed to relying on a proportional hazards assumption.
## We're modeling/simulating expected (counterfactual) risks over time under the assumption that 
## individuals could have been followed until K-1. Hence everyone has risk data until until K-1.
## The aggregation at the end ensures that the true censoring distribution is respected when computing risks.
## Some general details regarding the PLR model:
## Age is included with splines
## CAVE: Make sure not to include follow-up time AFTER censoring event (here it is ok, see above, fn_expand_intervals)
## CAVE: Currently cov_cat_region is included as any other confounder - however in protocol we specified cov_cat_region as a stratification. Discuss, rethink.
## In this case our max follow-up is: K = 39 months or K = 169 weeks (from earliest possible landmark_date (01.01.2019) to study end (01.04.2022))

### Adjustment for baseline confounding
## Adjustment for baseline confounding, which can be conceptualized as an attempt to emulate randomization in an observational analysis, 
## can be accomplished using a variety of methods. 
## I will use (i) standardization and (ii) inverse probability weighting (IPW) and (iii) their combination (to obtain more precise estimates)

## (i) Standardization
## The standardized outcome among the treated and untreated groups is estimated by taking a weighted average of the conditional
## outcomes, using the prevalence of the baseline confounders in the observed study population as weights. This entails:
## (1) fitting an outcome regression model conditional on the confounders listed above and 
## (2) standardizing over the empirical distribution of the confounders to obtain marginal effect estimates. 
## We model the follow-up time in the outcome regression model using linear and quadratic terms and include product terms between the treatment group indicator and follow-up time.
## Other follow-up time modelling is possible (cubic, splines).

## (ii) IPW
## Like standardization, IPW can also be used to obtain marginal estimates of causal effects. 
## Briefly, IPW can be used to create a pseudopopulation (i.e.,a hypothetical population) in which treatment is independent of the measured confounders.
## Informally, the denominator of the inverse probability weight for each individual is the probability of receiving their observed treatment value, given their confounder history.


# Standardization ---------------------------------------------------------


# MONTH
K <- 39 # Define total follow-up (currently, all is in months)
df_long_months$monthsqr <- df_long_months$month^2
plr_formula_severecovid <- as.formula(paste("out_bin_severecovid_afterlandmark ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr) +", 
                                            paste(covariate_names, collapse = " + ")))
plr_model_severecovid <- glm(plr_formula_severecovid, 
                             family = binomial(link = 'logit'),
                             data = df_long_months) 
# summary(plr_model_severecovid)


# 95% CI using robust standard errors and the delta method for RD and RR----
## Robust SEs (sandwich estimator): Accounts for heteroskedasticity and potential model misspecifications.
## Delta method: Used for computing confidence intervals for functions of regression coefficients (RD, RR).

# Create dataset with all time points for each individual under each treatment level
df_pred <- df_long_months %>% filter(month == 0) %>% select(-month) %>% crossing(month = 0:(K-1))
df_pred$monthsqr <- df_pred$month^2

# Control group (everyone untreated)
df_pred0 <- df_pred
df_pred0$exp_bin_treat <- 0
df_pred0$p.event0 <- predict(plr_model_severecovid, df_pred0, type="response")

# Treatment group (everyone treated)
df_pred1 <- df_pred
df_pred1$exp_bin_treat <- 1
df_pred1$p.event1 <- predict(plr_model_severecovid, df_pred1, type="response")
# The above creates a person-time dataset where we have predicted discrete-time hazards

# Obtain predicted survival probabilities from discrete-time hazards
df_pred0 <- df_pred0 %>% group_by(patient_id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
df_pred1 <- df_pred1 %>% group_by(patient_id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()

# Estimate risks from survival probabilities; compute risk at month K-1
# Risk = 1 - S(t)
df_pred0$risk0 <- 1 - df_pred0$surv0
df_pred1$risk1 <- 1 - df_pred1$surv1

# Get the mean in each treatment group at each month time point
risk0 <- aggregate(df_pred0[c("exp_bin_treat", "month", "risk0")], by=list(df_pred0$month), FUN=mean)[c("exp_bin_treat", "month", "risk0")]
risk1 <- aggregate(df_pred1[c("exp_bin_treat", "month", "risk1")], by=list(df_pred1$month), FUN=mean)[c("exp_bin_treat", "month", "risk1")]

# Put all in 1 data frame (for a plot but also to extract specific time points)
graph.pred <- merge(risk0, risk1, by=c("month"))
# Edit data frame to reflect that risks are estimated at the END of each interval
graph.pred$time_0 <- graph.pred$month + 1
zero <- data.frame(cbind(0,0,0,1,0,0))
zero <- setNames(zero,names(graph.pred))
graph <- rbind(zero, graph.pred) ## can be used for the cumulative incidence plot (but without 95% CI, see below)

# Add RD and RR
graph$rd <- graph$risk1-graph$risk0
graph$rr <- graph$risk1/graph$risk0

# Extract overall risk estimate (end of follow-up)
risk0 <- graph$risk0[which(graph$month==K-1)] 
risk1 <- graph$risk1[which(graph$month==K-1)]
rd <- graph$rd[which(graph$month==K-1)]
rr <- graph$rr[which(graph$month==K-1)]

# Compute robust standard errors
vcov_robust <- vcovHC(plr_model_severecovid, type = "HC0")
robust_se <- sqrt(diag(vcov_robust))
coefs <- coef(plr_model_severecovid)

# Compute standard errors for the two arms
se_risk0 <- sqrt(vcov_robust["(Intercept)", "(Intercept)"])
se_risk1 <- sqrt(vcov_robust["exp_bin_treat", "exp_bin_treat"])

# Compute confidence intervals for predicted risks using normal approximation
ci_risk0 <- c(risk0 - 1.96 * se_risk0, risk0 + 1.96 * se_risk0)
ci_risk1 <- c(risk1 - 1.96 * se_risk1, risk1 + 1.96 * se_risk1)

# Formula for risk difference (RD)
rd_formula <- "1 / (1 + exp(-((Intercept)))) - 1 / (1 + exp(-((Intercept) + exp_bin_treat)))"

# Delta method for RD
se_rd <- deltaMethod(coefs, rd_formula, vcov_robust)$SE
ci_rd_lower <- rd - 1.96 * se_rd
ci_rd_upper <- rd + 1.96 * se_rd

# Formula for risk ratio (RR)
rr_formula <- "(1 / (1 + exp(-((Intercept) + exp_bin_treat)))) / (1 / (1 + exp(-((Intercept)))))"

# Delta method for RR
se_rr <- deltaMethod(coefs, rr_formula, vcov_robust)$SE
ci_rr_lower <- rr - 1.96 * se_rr
ci_rr_upper <- rr + 1.96 * se_rr

# Create results table
risk_estimates_from_plr_rse_tbl <- data.frame(
  Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
  Estimate = c(risk0, risk1, rd, rr),
  Lower_CI = c(ci_risk0[1], ci_risk1[1], ci_rd_lower, ci_rr_lower),
  Upper_CI = c(ci_risk0[2], ci_risk1[2], ci_rd_upper, ci_rr_upper)
)

# SUPERSEDED 95% CI using bootstrapping ----------------------------------------
# risk_estimates_from_plr_withoutCI <- function(df_long, plr_model, K, interval_type = "month") {
#   
#   # Determine column names based on interval type
#   interval_col <- ifelse(interval_type == "month", "month", "week")
#   interval_max <- K - 1  # Ensuring indexing starts from 0
#   interval_days <- ifelse(interval_type == "month", 30, 7)  # Fixed month length of 30 days, same as fn_expand_intervals
#   
#   # Expand dataset by creating time intervals 0 to K-1 for each individual
#   treat0 <- df_long %>%
#     filter(.data[[interval_col]] == 0) %>%
#     select(-all_of(interval_col)) %>%  # Remove interval column to re-add it (e.g. we used month in df_long_months, but re-create it again here)
#     crossing(!!interval_col := 0:interval_max)
#   
#   treat0 <- treat0 %>%
#     mutate(!!paste0(interval_col, "sqr") := .data[[interval_col]]^2) # re-create time squared
#   
#   # Create treatment groups (forcing treatment to 0 or 1)
#   treat0$exp_bin_treat <- 0
#   treat1 <- treat0
#   treat1$exp_bin_treat <- 1
#   
#   # Predict discrete-time hazards
#   treat0[[paste0("p.event0")]] <- predict(plr_model, treat0, type = "response")
#   treat1[[paste0("p.event1")]] <- predict(plr_model, treat1, type = "response")
#   
#   # Compute survival probabilities
#   treat0.surv <- treat0 %>%
#     group_by(patient_id) %>%
#     mutate(!!paste0("surv0") := cumprod(1 - .data[[paste0("p.event0")]])) %>%
#     ungroup()
#   
#   treat1.surv <- treat1 %>%
#     group_by(patient_id) %>%
#     mutate(!!paste0("surv1") := cumprod(1 - .data[[paste0("p.event1")]])) %>%
#     ungroup()
#   
#   # Compute risk as 1 - survival probability
#   treat0.surv <- treat0.surv %>%
#     mutate(!!paste0("risk0") := 1 - .data[[paste0("surv0")]])
#   
#   treat1.surv <- treat1.surv %>%
#     mutate(!!paste0("risk1") := 1 - .data[[paste0("surv1")]])
#   
#   # Aggregate risks by time point
#   risk0 <- aggregate(treat0.surv[c("exp_bin_treat", interval_col, "risk0")], # aggregate in addition by exp_bin_treat not needed, but does no hurt and useful for future trial emulation setup 
#                      by = list(treat0.surv[[interval_col]]), FUN = mean)[c("exp_bin_treat", interval_col, "risk0")]
#   
#   risk1 <- aggregate(treat1.surv[c("exp_bin_treat", interval_col, "risk1")],
#                      by = list(treat1.surv[[interval_col]]), FUN = mean)[c("exp_bin_treat", interval_col, "risk1")]
#   
#   # Merge risk estimates
#   graph.pred <- merge(risk0, risk1, by = interval_col)
#   
#   # Adjust time to reflect end of each interval
#   graph.pred$time_0 <- graph.pred[[interval_col]] + 1
#   zero <- data.frame(cbind(0, 0, 0, 1, 0, 0))
#   zero <- setNames(zero, names(graph.pred))
#   graph <- rbind(zero, graph.pred)
#   
#   # Compute risk differences and risk ratios
#   graph$rd <- graph$risk1 - graph$risk0
#   graph$rr <- graph$risk1 / graph$risk0
#   
#   return(graph)
# }

# Generate risk estimates for 39 months follow-up duration. nice, but not needed, since bootstrap below returns also original point estimates
# graph_months <- risk_estimates_from_plr_withoutCI(df_long_months, plr_model_severecovid, K = 39, interval_type = "month")

### Use pooled logistic regression estimates to compute causal estimates, but only point estimates
# end of follow-up estimates (i.e. "the overall risk")

# K <- 39 # total follow-up (here, months)
# graph_months$risk0[which(graph_months$month==K-1)]
# graph_months$risk1[which(graph_months$month==K-1)]
# graph_months$rd[which(graph_months$month==K-1)]
# graph_months$rr[which(graph_months$month==K-1)]


# 95% CI using bootstrapping ----------------------------------------------
# Reduce dataset as input into bootstrap?
# Adapt in a second step to dynamically include weeks instead of months

# Create input list of ids (eligible persons)
K <- 39 # Define total follow-up (currently, all is in months)
R <- 2 # Define number of bootstraps
study_ids <- data.frame(patient_id = df$patient_id)

# Create a function to obtain risks, RD, and RR from each bootstrap sample - and return 1 time point (e.g. overall risk)
te_one_timepoint_rd_rr_withCI <- function(data, indices) {
  # Select individuals into each bootstrapped sample
  ids <- data$patient_id
  boot.ids <- data.frame(patient_id = ids[indices])
  boot.ids$bid <- 1:nrow(boot.ids)
  
  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- left_join(boot.ids, df_long_months, by = "patient_id", relationship = "many-to-many")
  
  # Fit pooled logistic model to estimate discrete hazards
  plr_formula_severecovid <- as.formula(paste("out_bin_severecovid_afterlandmark ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr) +", 
                                              paste(covariate_names, collapse = " + ")))
  plr_model_severecovid <- glm(plr_formula_severecovid, 
                               family = binomial(link = 'logit'),
                               data = d)
  
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- d %>% filter(month == 0) %>% select(-month) %>% crossing(month = 0:(K-1))
  treat0$monthsqr <- treat0$month^2
  
  # In the no metformin arm ("force" everyone to be untreated)
  treat0$exp_bin_treat <- 0
  # In the metformin arm ("force" everyone to be treated)
  treat1 <- treat0
  treat1$exp_bin_treat <- 1
  
  # Extract predicted values from pooled logistic regression model for each person-time row
  # Predicted values correspond to discrete-time hazards
  treat0$p.event0 <- predict(plr_model_severecovid, treat0, type="response")
  treat1$p.event1 <- predict(plr_model_severecovid, treat1, type="response")
  # The above creates a person-time dataset where we have predicted discrete-time hazards
  # For each person-time row in the dataset
  
  # Obtain predicted survival probabilities from discrete-time hazards
  treat0.surv <- treat0 %>% group_by(bid) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
  treat1.surv <- treat1 %>% group_by(bid) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  treat0.surv$risk0 <- 1 - treat0.surv$surv0
  treat1.surv$risk1 <- 1 - treat1.surv$surv1
  
  # Get the mean in each treatment group at each month time point
  risk0 <- aggregate(treat0.surv[c("exp_bin_treat", "month", "risk0")], by=list(treat0.surv$month), FUN=mean)[c("exp_bin_treat", "month", "risk0")]
  risk1 <- aggregate(treat1.surv[c("exp_bin_treat", "month", "risk1")], by=list(treat1.surv$month), FUN=mean)[c("exp_bin_treat", "month", "risk1")]
  
  # Prepare data
  graph.pred <- merge(risk0, risk1, by=c("month"))
  # Edit data frame to reflect that risks are estimated at the END of each interval
  graph.pred$time_0 <- graph.pred$month + 1
  zero <- data.frame(cbind(0,0,0,1,0,0))
  zero <- setNames(zero,names(graph.pred))
  graph <- rbind(zero, graph.pred)
  
  graph$rd <- graph$risk1-graph$risk0
  graph$rr <- graph$risk1/graph$risk0
  return(c(graph$risk0[which(graph$month==K-1)], # adapt time point if necessary
           graph$risk1[which(graph$month==K-1)],
           graph$rd[which(graph$month==K-1)],
           graph$rr[which(graph$month==K-1)]))
}

# Run
# set.seed(423)
# te_one_timepoint_rd_rr_withCI <- boot(data = study_ids, statistic = te_one_timepoint_rd_rr_withCI, R = R)

# Function to extract bootstrapped confidence intervals - but also add a column with original point estimates (without CI)
# extract_ci_boot <- function(boot_obj, index) {
#   ci <- boot.ci(boot_obj, conf = 0.95, type = "perc", index = index)
#   if (!is.null(ci$percent)) {
#     return(c(ci$percent[4], ci$percent[5]))  # Lower and Upper CI
#   } else {
#     return(c(NA, NA))  # If CI is not available
#   }
# }
# 
# # Create results table
# risk_estimates_from_plr_boot_tbl <- data.frame(
#   Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
#   Estimate_original = te_one_timepoint_rd_rr_withCI$t0,  # Original estimates
#   Estimate_boot = colMeans(te_one_timepoint_rd_rr_withCI$t),  # Bootstrapped mean estimate
#   Lower_CI = sapply(1:4, function(i) extract_ci_boot(te_one_timepoint_rd_rr_withCI, i)[1]),
#   Upper_CI = sapply(1:4, function(i) extract_ci_boot(te_one_timepoint_rd_rr_withCI, i)[2])
# )


# Marginal parametric cumulative incidence (risk) curves ------------------
# Create input list of ids (eligible persons)
K <- 39 # Define total follow-up (currently, all is in months)
R <- 2 # Define number of bootstraps
study_ids <- data.frame(patient_id = df$patient_id)

# Bootstrap function
te_all_timepoints_withCI <- function(data, indices) {
  # Select individuals into each bootstrapped sample
  ids <- data$patient_id
  boot.ids <- data.frame(patient_id = ids[indices])
  boot.ids$bid <- 1:nrow(boot.ids)
  
  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- left_join(boot.ids, df_long_months, by = "patient_id", relationship = "many-to-many")
  
  # Fit pooled logistic model to estimate discrete hazards
  plr_formula_severecovid <- as.formula(paste("out_bin_severecovid_afterlandmark ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr) +", 
                                              paste(covariate_names, collapse = " + ")))
  plr_model_severecovid <- glm(plr_formula_severecovid, 
                               family = binomial(link = 'logit'),
                               data = d)
  
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- d %>% filter(month == 0) %>% select(-month) %>% crossing(month = 0:(K-1))
  treat0$monthsqr <- treat0$month^2
  
  # In the no metformin arm ("force" everyone to be untreated)
  treat0$exp_bin_treat <- 0
  # In the metformin arm ("force" everyone to be treated)
  treat1 <- treat0
  treat1$exp_bin_treat <- 1
  
  # Extract predicted values from pooled logistic regression model for each person-time row
  # Predicted values correspond to discrete-time hazards
  treat0$p.event0 <- predict(plr_model_severecovid, treat0, type="response")
  treat1$p.event1 <- predict(plr_model_severecovid, treat1, type="response")
  # The above creates a person-time dataset where we have predicted discrete-time hazards
  # For each person-time row in the dataset
  
  # Obtain predicted survival probabilities from discrete-time hazards
  treat0.surv <- treat0 %>% group_by(bid) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
  treat1.surv <- treat1 %>% group_by(bid) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  treat0.surv$risk0 <- 1 - treat0.surv$surv0
  treat1.surv$risk1 <- 1 - treat1.surv$surv1
  
  # Get the mean in each treatment group at each month time point
  risk0 <- aggregate(treat0.surv[c("exp_bin_treat", "month", "risk0")], by=list(treat0.surv$month), FUN=mean)[c("exp_bin_treat", "month", "risk0")]
  risk1 <- aggregate(treat1.surv[c("exp_bin_treat", "month", "risk1")], by=list(treat1.surv$month), FUN=mean)[c("exp_bin_treat", "month", "risk1")]
  
  # Prepare data
  graph.pred <- merge(risk0, risk1, by=c("month"))
  # Edit data frame to reflect that risks are estimated at the END of each interval
  graph.pred$time_0 <- graph.pred$month + 1
  zero <- data.frame(cbind(0,0,0,1,0,0))
  zero <- setNames(zero,names(graph.pred))
  graph <- rbind(zero, graph.pred)
  
  return(c(graph$risk0, graph$risk1))
}

# set.seed(423)
# te_all_timepoints_withCI <- boot(data = study_ids, statistic = te_all_timepoints_withCI, R = R)

# Check the dimensions of the bootstrapped results
# dim(te_all_timepoints_withCI$t)  # Should be (R, 2 * N=timepoints) because risk0 and risk1 are concatenated
# Inspect the first bootstrap iteration's results
# te_all_timepoints_withCI$t[1, ]

# Create an empty data frame to store the structured output for the plot
# risk_graph <- data.frame(
#   time = 0:K,
#   mean.0 = numeric(K+1),
#   ll.0 = numeric(K+1),
#   ul.0 = numeric(K+1),
#   mean.1 = numeric(K+1),
#   ll.1 = numeric(K+1),
#   ul.1 = numeric(K+1)
# )
# 
# # Calculate the mean and confidence intervals for control group (mean_risk0) and treatment group (mean_risk1)
# # We assume that `te_all_timepoints_withCI$t` holds the bootstrapped estimates
# 
# mean_risk0 <- apply(te_all_timepoints_withCI$t[, 1:40], 2, mean)  # Mean for control group
# mean_risk1 <- apply(te_all_timepoints_withCI$t[, 41:80], 2, mean) # Mean for treatment group
# K <- 39
# 
# mean_risk0 <- apply(te_all_timepoints_withCI$t[, 1:(K+1)], 2, mean)  # Mean for control group; 1:40 (first half contains risk0)
# mean_risk1 <- apply(te_all_timepoints_withCI$t[, (K+2):(2*(K+1))], 2, mean) # Mean for treatment group; 41:80 (second half contains risk1)
# 
# # Calculate 2.5th and 97.5th percentiles (CI) for control and treatment groups
# ll_risk0 <- apply(te_all_timepoints_withCI$t[, 1:(K+1)], 2, function(x) quantile(x, 0.025)) # Lower bound CI for control
# ul_risk0 <- apply(te_all_timepoints_withCI$t[, 1:(K+1)], 2, function(x) quantile(x, 0.975)) # Upper bound CI for control
# 
# ll_risk1 <- apply(te_all_timepoints_withCI$t[, (K+2):(2*(K+1))], 2, function(x) quantile(x, 0.025)) # Lower bound CI for treatment
# ul_risk1 <- apply(te_all_timepoints_withCI$t[, (K+2):(2*(K+1))], 2, function(x) quantile(x, 0.975)) # Upper bound CI for treatment
# 
# # Populate the `risk.boot.graph` data frame
# risk_graph$mean.0 <- mean_risk0
# risk_graph$ll.0 <- ll_risk0
# risk_graph$ul.0 <- ul_risk0
# risk_graph$mean.1 <- mean_risk1
# risk_graph$ll.1 <- ll_risk1
# risk_graph$ul.1 <- ul_risk1
# 
# # Create plot
# plot_cum_risk <- ggplot(risk_graph,
#                       aes(x=time)) +
#   geom_line(aes(y = mean.1, # create line for intervention group
#                 color = "Metformin"), linewidth = 1.5) +
#   geom_ribbon(aes(ymin = ll.1, ymax = ul.1, fill = "Metformin"), alpha = 0.4) +
#   geom_line(aes(y = mean.0, # create line for control group
#                 color = "No Metformin"), linewidth = 1.5) +
#   geom_ribbon(aes(ymin = ll.0, ymax = ul.0, fill = "No Metformin"), alpha=0.4) +
#   xlab("Months") +
#   # scale_x_continuous(limits = c(0, 39), # format x axis
#   #                    breaks=c(0, 6, 12, 24, 36, 39)) +
#   ylab("Cumulative Incidence (%)") + # label y axis
#   # scale_y_continuous(limits=c(0, 0.125), # format y axis
#   #                    breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125),
#   #                    labels=c("0.0%", "2.5%", "5.0%",
#   #                             "7.5%", "10.0%", "12.5%")) +
#   theme_minimal()+
#   theme(axis.text = element_text(size=14), legend.position.inside = c(0.2, 0.8),
#         axis.line = element_line(colour = "black"),
#         legend.title = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.major.y = element_blank())+
#   scale_color_manual(values=c("#E7B800", # set colors
#                               "#2E9FDF"),
#                      breaks=c('No Metformin',
#                               'Metformin')) +
#   scale_fill_manual(values=c("#E7B800", # set colors
#                              "#2E9FDF"),
#                     breaks=c('No Metformin',
#                              'Metformin'))

# Use TrialEmulation package instead -------------------------------------
# library(TrialEmulation)
# risk_estimates_from_plr_package <- trial_msm(
#   data = df_long_months, 
#   estimand_type = "ITT", 
#   outcome_cov = covariate_names,
#   model_var = "exp_bin_treat", 
#   glm_function = "parglm", 
#   use_sample_weights = FALSE, 
#   analysis_weights = "unweighted" # no artificial censoring weights
# )

# Save output -------------------------------------------------------------
# Risk estimates from plr model
write.csv(risk_estimates_from_plr_rse_tbl, file = here::here("output", "te", "pooled_log_reg", "risk_estimates_from_plr_rse_severecovid.csv"))
# write.csv(risk_estimates_from_plr_boot_tbl, file = here::here("output", "te", "pooled_log_reg", "risk_estimates_from_plr_boot_severecovid.csv"))
# Marginal parametric cumulative incidence (risk) curves from plr model, with 95% CI
# ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_cum_risk_severecovid.png"), plot_cum_risk, width = 20, height = 20, units = "cm")
