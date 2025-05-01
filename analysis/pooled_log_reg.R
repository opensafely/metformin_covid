####
## This script does the following:
# 1. Import processed data
# 2. Restructure the one-person-per-row data into a multiple-event-per-person data
# 3. Estimate adjusted but marginal risks per group, risk differences, risk ratios using (i) standardization and (ii) inverse probability weighting (IPW), 
#    with the help of pooled logistic regression, and derive 95% CI via robust standard error/delta method and bootstrapping
# 4. Create marginal parametric cumulative incidence (risk) curves incl. 95% CI (via bootstrapping)
# 5. Save all output
####


# ToDo --------------------------------------------------------------------
## a) Adapt the script to dynamically use months or weeks dynamically
## b) Adapt the script to dynamically use different outcomes
## c) Discuss re bootstrapping by arm
## d) Discuss re truncating/trimming: Use trimmed dataset from PS from the onset? Still apply truncation at 99th?
## e) Complete the marginal parametric cumulative incidence (risk) curves incl. 95% CI bootstrapping for IPTW & IPCW (currently only curves from IPTW & IPCW without 95%)
## f) Consider doing the Love/SMD plot in the separate PS branch/script only
## g) Move each function into a separate R script
## h) Consider deleting (or moving) all (i) Standardization chapters, only keep IPW approach
## i) Discuss different spline methods, ns() versus rcs(), but rcs() does not work within parglm!


# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(lubridate)
library(splines)
library(rms) # strat
library(purrr)
library(sandwich) # for robust standard errors
library(lmtest) # for hypothesis testing with robust standard errors
library(car) # for deltaMethod
library(boot)
library(ggplot2)
library(Hmisc) # for Love/SMD plot
library(parglm) # to be computationally more efficient
source(here::here("analysis", "functions", "fn_expand_intervals.R"))


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "te", "pooled_log_reg"))


# Import the data ---------------------------------------------------------
print('Import the data')
df <- read_feather(here("output", "data", "data_processed.arrow"))


# Import dates ------------------------------------------------------------
print('Import dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# Add splines -------------------------------------------------------------
print('Add splines')
# Compute knot locations based on percentiles, according to study protocol
age_knots <- quantile(df$cov_num_age, probs = c(0.10, 0.50, 0.90))
df <- df %>%
  mutate(cov_num_age_spline = ns(cov_num_age, knots = age_knots))


# Define covariates ----------------------------------------------------------
print('Define covariates')
covariate_names <- names(df) %>%
  grep("^cov_", ., value = TRUE) %>% 
  # exclude those not needed in the model: 
  ## cov_bin_obesity covers for cov_num_bmi & cov_cat_bmi_groups,
  ## cov_cat_hba1c_mmol_mol covers cov_num_hba1c_mmol_mol
  ## cov_cat_tc_hdl_ratio covers cov_num_tc_hdl_ratio
  ## cov_num_age_spline covers cov_cat_age and cov_num_age
  setdiff(c("cov_num_bmi", "cov_cat_bmi_groups", "cov_num_hba1c_mmol_mol", "cov_num_tc_hdl_ratio", "cov_cat_age", "cov_num_age")) 
print(covariate_names)

# Drop unnecessary variables to reduce dataset size -----------------------
print('Drop unnecessary variables to reduce dataset size')
df <- df %>% 
  dplyr::select(patient_id, exp_bin_treat, elig_date_t2dm, landmark_date,
                starts_with("cov_"),
                starts_with("strat_"),
                starts_with("out_") & ends_with("_afterlandmark"),
                starts_with("cens_"))


# Expand the dataset ------------------------------------------------------
print('Expand the dataset')
# Expand the dataset into intervals and follow these rules:
# a) If outcome (out_date_severecovid_afterlandmark) is reached first, assign outcome=1, censor=0, comp_event=0 to the interval when it happened and stop expanding
# b) If competing event (out_date_death_afterlandmark) is reached first, assign outcome=NA, censor=0, comp_event=1 to the interval when it happened and stop expanding
# c) If censoring event (out_date_ltfu_afterlandmark) is reached first, assign outcome=NA, censor=1, comp_event=NA to the interval when it happened and stop expanding
# d) If studyend_date is reached first, then assign outcome=0, censor=0, comp_event=0 to the interval when it happened and stop expanding
# Use studyend date from metadates.R import
start_date_variable <- "landmark_date" # start expanding at landmark date, not elig_date_t2dm (due to landmark design)
stop_date_columns <- c("out_date_severecovid_afterlandmark", "out_date_death_afterlandmark", "cens_date_ltfu_afterlandmark")
outcome_date_variable <- "out_date_severecovid_afterlandmark"
comp_date_variable <- "out_date_death_afterlandmark"
censor_date_variable <- "cens_date_ltfu_afterlandmark"

# Apply the function, choose either weeks or months, currently only using months, but works for both.
df_long_months <- fn_expand_intervals(df, 
                                      start_date_variable,
                                      stop_date_columns, 
                                      studyend_date,
                                      outcome_date_variable, 
                                      comp_date_variable,
                                      censor_date_variable,
                                      interval_type = "month")

## To double-check | to check how events are treated in interval the event happend, see stop_date and start_date_month
# df_long_months %>%
#   dplyr::select(patient_id, elig_date_t2dm, landmark_date, out_date_severecovid_afterlandmark, out_date_death_afterlandmark,
#                 cens_date_ltfu_afterlandmark, stop_date, start_date_month, month, outcome, comp_event, censor, is_event_interval, is_outcome_event, followup_stop,
#                 cov_cat_sex, cov_num_age) %>%
#   # dplyr::filter(!is.na(out_date_death_afterlandmark)) %>%
#   # dplyr::filter(!is.na(cens_date_ltfu_afterlandmark)) %>%
#   # dplyr::filter(!is.na(out_date_severecovid_afterlandmark)) %>%
#   # dplyr::filter(is.na(censor)) %>%
#   # dplyr::filter(stop_date == "2022-04-01") %>%
#   # dplyr::filter(is_event_interval == TRUE & !is.na(out_date_severecovid_afterlandmark)) %>% # should only have 1 row per person
#   # dplyr::filter(patient_id == "1418") %>%
#   # dplyr::filter(is.na(censor)) %>% # should be empty
#   View()


# Background description -----------------------------------------------------
### a) Pooled logistic regression
## Pooled logistic regression models allow for the parametric estimation of risks, and thus risk differences and risk ratios.
## Also, these models can be specified such that the effect of the treatment can vary over time, without relying on a proportional hazards assumption.
## We're modeling/simulating expected/predicted (counterfactual) risks over time under the assumption that 
## individuals could have been followed until max fup time (K-1). Hence everyone has predicted risk data until K-1.
## The aggregation at the end ensures that the true censoring distribution is respected when computing risks.
## The dataset is currently set up as in CAUSALab TTE course material (esp. see coding for "outcome", "censor" and "comp_event")

### b) Adjustment for baseline confounding
## Adjustment for baseline confounding, i.e. attempt to emulate randomization can be accomplished using a variety of methods. 
## I will use (i) standardization and (ii) inverse probability weighting (IPW), however, focus on IPW going forward

## (i) Standardization
## The standardized outcome among the treated and untreated groups is estimated by taking a weighted average of the conditional
## outcomes, using the prevalence of the baseline confounders in the observed study population as weights. This entails:
## (1) fitting an outcome regression model conditional on the confounders and 
## (2) standardizing over the empirical distribution of the confounders to obtain marginal effect estimates. 
## We model the follow-up time in the outcome regression model using linear and quadratic terms and include product terms between the treatment group indicator and follow-up time.
## Other follow-up time modelling would be possible (cubic, splines), I will stick to linear and quadratic term.

## (ii) IPW
## Like standardization, IPW can be used to obtain marginal estimates of causal effects. 
## Briefly, IPW can be used to create a pseudopopulation (i.e.,a hypothetical population) in which treatment is independent of the measured confounders.
## Informally, the denominator of the inverse probability weight for each individual is the probability of receiving their observed treatment value, given their confounder history.
## Unstabilized and stabilized weights can be used, I will focus on stabilized.
## If we want to incorporate time-varying confounding, e.g. for censoring weights or adherence weights (per protocol analysis), 
## then standardization alone is not possible anymore, we then need time-update and time-varying weights, and simply multiply all weights

## (iii) Combination of above, useful when:
## If we need to adjust for additional baseline confounding that can't be included when creating IPW (e.g. baseline calendar week/month in sequential trial setup), then we create IPW first and then standardize to the empirical distribution of baseline calendar week/month in the dataset


# Define interval data set and number of bootstraps ----------------------------
## I currently only use the monthly interval data set => max follow-up is: K = 39 months (earliest possible landmark_date [01.01.2019] to study end [01.04.2022])
## If weeks, then K = 169 weeks
K <- 39 # Total follow-up in months
df_long_months$monthsqr <- df_long_months$month^2 # add months square
R <- 10 # Total bootstraps (ideally >500)


# (i) Standardization ---------------------------------------------------------
# print('Standardization creating te_plr_stand_rse_tbl')
# ## (1) fitting an outcome regression model conditional on the confounders listed above
# df_long_months$monthsqr <- df_long_months$month^2
# plr_formula_severecovid <- as.formula(paste("outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr) +", 
#                                             paste(covariate_names, collapse = " + "), "+ strat(strat_cat_region)"))
# plr_model_severecovid <- parglm(plr_formula_severecovid, 
#                              family = binomial(link = 'logit'),
#                              data = df_long_months) 
# # summary(plr_model_severecovid)
# 
# ## (2) standardizing over the empirical distribution of the confounders to obtain marginal effect estimates. 
# # Create dataset with all time points for each individual under each treatment level
# df_pred <- df_long_months %>% filter(month == 0) %>% dplyr::select(-month) %>% crossing(month = 0:(K-1))
# df_pred$monthsqr <- df_pred$month^2
# 
# # Control group (everyone untreated)
# df_pred0 <- df_pred
# df_pred0$exp_bin_treat <- 0
# df_pred0$p.event0 <- predict(plr_model_severecovid, df_pred0, type="response")
# 
# # Treatment group (everyone treated)
# df_pred1 <- df_pred
# df_pred1$exp_bin_treat <- 1
# df_pred1$p.event1 <- predict(plr_model_severecovid, df_pred1, type="response")
# # The above creates a person-time dataset where we have predicted discrete-time hazards
# 
# # Obtain predicted survival probabilities from discrete-time hazards
# df_pred0 <- df_pred0 %>% group_by(patient_id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
# df_pred1 <- df_pred1 %>% group_by(patient_id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
# 
# # Estimate risks from survival probabilities; compute risk at month K-1
# # Risk = 1 - S(t)
# df_pred0$risk0 <- 1 - df_pred0$surv0
# df_pred1$risk1 <- 1 - df_pred1$surv1
# 
# # Get the mean in each treatment group at each month time point
# risk0 <- aggregate(df_pred0[c("exp_bin_treat", "month", "risk0")], by=list(df_pred0$month), FUN=mean)[c("exp_bin_treat", "month", "risk0")]
# risk1 <- aggregate(df_pred1[c("exp_bin_treat", "month", "risk1")], by=list(df_pred1$month), FUN=mean)[c("exp_bin_treat", "month", "risk1")]
# 
# # Put all in 1 data frame
# graph.pred <- merge(risk0, risk1, by=c("month"))
# 
# # Edit data frame to reflect that risks are estimated at the END of each interval
# graph.pred$time_0 <- graph.pred$month + 1
# zero <- data.frame(cbind(0,0,0,1,0,0))
# zero <- setNames(zero,names(graph.pred))
# graph <- rbind(zero, graph.pred) ## can be used for the cumulative incidence plot, but without 95% CI (for that, see below)
# 
# # Add RD and RR
# graph$rd <- graph$risk1-graph$risk0
# graph$rr <- graph$risk1/graph$risk0
# 
# # Extract overall risk estimate (end of follow-up K-1), but without 95% CI (for that, see below)
# risk0 <- graph$risk0[which(graph$month==K-1)] 
# risk1 <- graph$risk1[which(graph$month==K-1)]
# rd <- graph$rd[which(graph$month==K-1)]
# rr <- graph$rr[which(graph$month==K-1)]


# (i) Standardization: Add 95% CI using robust standard errors and the delta method for RD and RR | See code from Will/Fizz ----
# # Compute robust standard errors
# vcov_robust <- vcovHC(plr_model_severecovid, type = "HC0")
# robust_se <- sqrt(diag(vcov_robust))
# coefs <- coef(plr_model_severecovid)
# 
# # standard errors for the two arms
# se_risk0 <- sqrt(vcov_robust["(Intercept)", "(Intercept)"])
# se_risk1 <- sqrt(vcov_robust["exp_bin_treat", "exp_bin_treat"])
# 
# # confidence intervals for predicted risks using normal approximation
# ci_risk0 <- c(risk0 - 1.96 * se_risk0, risk0 + 1.96 * se_risk0)
# ci_risk1 <- c(risk1 - 1.96 * se_risk1, risk1 + 1.96 * se_risk1)
# 
# # Formula for risk difference (RD)
# rd_formula <- "1 / (1 + exp(-((Intercept)))) - 1 / (1 + exp(-((Intercept) + exp_bin_treat)))"
# 
# # Delta method for RD
# se_rd <- deltaMethod(coefs, rd_formula, vcov_robust)$SE
# ci_rd_lower <- rd - 1.96 * se_rd
# ci_rd_upper <- rd + 1.96 * se_rd
# 
# # Formula for risk ratio (RR)
# rr_formula <- "(1 / (1 + exp(-((Intercept) + exp_bin_treat)))) / (1 / (1 + exp(-((Intercept)))))"
# 
# # Delta method for RR
# se_rr <- deltaMethod(coefs, rr_formula, vcov_robust)$SE
# ci_rr_lower <- rr - 1.96 * se_rr
# ci_rr_upper <- rr + 1.96 * se_rr
# 
# # Create results table
# te_plr_stand_rse_tbl <- data.frame(
#   Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
#   Estimate = c(risk0, risk1, rd, rr),
#   Lower_CI = c(ci_risk0[1], ci_risk1[1], ci_rd_lower, ci_rr_lower),
#   Upper_CI = c(ci_risk0[2], ci_risk1[2], ci_rd_upper, ci_rr_upper)
# )


# (i) Standardization: Add 95% CI for risks, RD and RR using bootstrapping ---------------
# print('Standardization creating te_plr_stand_boot_tbl')
# study_ids <- data.frame(patient_id = df$patient_id)
# 
# # Create a function to obtain risks, RD, and RR from each bootstrap sample - and return 1 risk time point
# te_stand_rd_rr_withCI <- function(data, indices) {
#   # Select individuals into each bootstrapped sample
#   ids <- data$patient_id
#   boot.ids <- data.frame(patient_id = ids[indices])
#   boot.ids$bid <- 1:nrow(boot.ids)
#   
#   # Subset person-time data to individuals selected into the bootstrapped sample
#   d <- left_join(boot.ids, df_long_months, by = "patient_id", relationship = "many-to-many")
#   
#   # Fit pooled logistic model to estimate discrete hazards
#   plr_formula_severecovid <- as.formula(paste("outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr) +", 
#                                               paste(covariate_names, collapse = " + "), "+ strat(strat_cat_region)"))
#   plr_model_severecovid <- parglm(plr_formula_severecovid, 
#                                family = binomial(link = 'logit'),
#                                data = d)
#   
#   # Create dataset with all time points for each individual under each treatment level
#   treat0 <- d %>% filter(month == 0) %>% dplyr::select(-month) %>% crossing(month = 0:(K-1))
#   treat0$monthsqr <- treat0$month^2
#   
#   # In the no metformin arm ("force" everyone to be untreated)
#   treat0$exp_bin_treat <- 0
#   # In the metformin arm ("force" everyone to be treated)
#   treat1 <- treat0
#   treat1$exp_bin_treat <- 1
#   
#   # Extract predicted values from pooled logistic regression model for each person-time row
#   # Predicted values correspond to discrete-time hazards
#   treat0$p.event0 <- predict(plr_model_severecovid, treat0, type="response")
#   treat1$p.event1 <- predict(plr_model_severecovid, treat1, type="response")
#   # The above creates a person-time dataset where we have predicted discrete-time hazards
#   # For each person-time row in the dataset
#   
#   # Obtain predicted survival probabilities from discrete-time hazards
#   treat0.surv <- treat0 %>% group_by(bid) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
#   treat1.surv <- treat1 %>% group_by(bid) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
#   
#   # Estimate risks from survival probabilities
#   # Risk = 1 - S(t)
#   treat0.surv$risk0 <- 1 - treat0.surv$surv0
#   treat1.surv$risk1 <- 1 - treat1.surv$surv1
#   
#   # Get the mean in each treatment group at each month time point
#   risk0 <- aggregate(treat0.surv[c("exp_bin_treat", "month", "risk0")], by=list(treat0.surv$month), FUN=mean)[c("exp_bin_treat", "month", "risk0")]
#   risk1 <- aggregate(treat1.surv[c("exp_bin_treat", "month", "risk1")], by=list(treat1.surv$month), FUN=mean)[c("exp_bin_treat", "month", "risk1")]
#   
#   # Prepare data
#   graph.pred <- merge(risk0, risk1, by=c("month"))
#   # Edit data frame to reflect that risks are estimated at the END of each interval
#   graph.pred$time_0 <- graph.pred$month + 1
#   zero <- data.frame(cbind(0,0,0,1,0,0))
#   zero <- setNames(zero,names(graph.pred))
#   graph <- rbind(zero, graph.pred)
#   graph$rd <- graph$risk1-graph$risk0
#   graph$rr <- graph$risk1/graph$risk0
#   return(c(graph$risk0[which(graph$month==K-1)], # adapt time point if necessary
#            graph$risk1[which(graph$month==K-1)],
#            graph$rd[which(graph$month==K-1)],
#            graph$rr[which(graph$month==K-1)]))
# }
# 
# set.seed(423)
# te_stand_rd_rr_withCI_boot <- boot(data = study_ids, statistic = te_stand_rd_rr_withCI, R = R)
# # summary(te_stand_rd_rr_withCI_boot$t)
# 
# # Function to extract bootstrapped confidence intervals - and keep a column with the original point estimates (without CI)
# extract_ci_boot <- function(boot_obj, index) {
#   ci <- boot.ci(boot_obj, conf = 0.95, type = "perc", index = index)
#   if (!is.null(ci$percent)) {
#     return(c(ci$percent[4], ci$percent[5])) # Lower and Upper CI for the bootstrapped values
#   } else {
#     return(c(NA, NA)) # for the original value
#   }
# }
# 
# # Create the results table
# te_plr_stand_boot_tbl <- data.frame(
#   Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
#   Estimate_original = te_stand_rd_rr_withCI_boot$t0,  # Original estimates
#   Estimate_boot = colMeans(te_stand_rd_rr_withCI_boot$t),  # Bootstrapped mean estimate
#   Lower_CI = sapply(1:4, function(i) extract_ci_boot(te_stand_rd_rr_withCI_boot, i)[1]),
#   Upper_CI = sapply(1:4, function(i) extract_ci_boot(te_stand_rd_rr_withCI_boot, i)[2])
# )


# (i) Standardization: Marginal parametric cumulative incidence (risk) curves incl. 95% CI bootstrapping ----
# print('Standardization creating marginal parametric cumulative incidence (risk) curves incl. 95% CI bootstrapping')
# ### This chapter will be removed, instead will do the graph with IPTW & IPCW ###
# 
# study_ids <- data.frame(patient_id = df$patient_id)
# # same function as above, except that a) it returns all risk timepoints incl. 95% CI, and b) no RD and RR
# te_all_timepoints_withCI <- function(data, indices) {
#   # Select individuals into each bootstrapped sample
#   ids <- data$patient_id
#   boot.ids <- data.frame(patient_id = ids[indices])
#   boot.ids$bid <- 1:nrow(boot.ids)
#   
#   # Subset person-time data to individuals selected into the bootstrapped sample
#   d <- left_join(boot.ids, df_long_months, by = "patient_id", relationship = "many-to-many")
#   
#   # Fit pooled logistic model to estimate discrete hazards
#   plr_formula_severecovid <- as.formula(paste("outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr) +", 
#                                               paste(covariate_names, collapse = " + "), "+ strat(strat_cat_region)"))
#   plr_model_severecovid <- parglm(plr_formula_severecovid, 
#                                family = binomial(link = 'logit'),
#                                data = d)
#   
#   # Create dataset with all time points for each individual under each treatment level
#   treat0 <- d %>% filter(month == 0) %>% dplyr::select(-month) %>% crossing(month = 0:(K-1))
#   treat0$monthsqr <- treat0$month^2
#   
#   # In the no metformin arm ("force" everyone to be untreated)
#   treat0$exp_bin_treat <- 0
#   # In the metformin arm ("force" everyone to be treated)
#   treat1 <- treat0
#   treat1$exp_bin_treat <- 1
#   
#   # Extract predicted values from pooled logistic regression model for each person-time row
#   # Predicted values correspond to discrete-time hazards
#   treat0$p.event0 <- predict(plr_model_severecovid, treat0, type="response")
#   treat1$p.event1 <- predict(plr_model_severecovid, treat1, type="response")
#   # The above creates a person-time dataset where we have predicted discrete-time hazards
#   # For each person-time row in the dataset
#   
#   # Obtain predicted survival probabilities from discrete-time hazards
#   treat0.surv <- treat0 %>% group_by(bid) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
#   treat1.surv <- treat1 %>% group_by(bid) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
#   
#   # Estimate risks from survival probabilities
#   # Risk = 1 - S(t)
#   treat0.surv$risk0 <- 1 - treat0.surv$surv0
#   treat1.surv$risk1 <- 1 - treat1.surv$surv1
#   
#   # Get the mean in each treatment group at each month time point
#   risk0 <- aggregate(treat0.surv[c("exp_bin_treat", "month", "risk0")], by=list(treat0.surv$month), FUN=mean)[c("exp_bin_treat", "month", "risk0")]
#   risk1 <- aggregate(treat1.surv[c("exp_bin_treat", "month", "risk1")], by=list(treat1.surv$month), FUN=mean)[c("exp_bin_treat", "month", "risk1")]
#   
#   # Prepare data
#   graph.pred <- merge(risk0, risk1, by=c("month"))
#   # Edit data frame to reflect that risks are estimated at the END of each interval
#   graph.pred$time_0 <- graph.pred$month + 1
#   zero <- data.frame(cbind(0,0,0,1,0,0))
#   zero <- setNames(zero,names(graph.pred))
#   graph <- rbind(zero, graph.pred)
#   
#   return(c(graph$risk0, graph$risk1))
# }

# set.seed(423)
# te_all_timepoints_withCI <- boot(data = study_ids, statistic = te_all_timepoints_withCI, R = R)

# Check the dimensions of the bootstrapped results
# dim(te_all_timepoints_withCI$t) # Should be (R, 2 * N=timepoints) because risk0 and risk1 are concatenated
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
# # Calculate the mean and confidence intervals
# mean_risk0 <- apply(te_all_timepoints_withCI$t[, 1:40], 2, mean)
# mean_risk1 <- apply(te_all_timepoints_withCI$t[, 41:80], 2, mean)
# 
# mean_risk0 <- apply(te_all_timepoints_withCI$t[, 1:(K+1)], 2, mean) # Mean for control group across 1:40 (first half contains risk0)
# mean_risk1 <- apply(te_all_timepoints_withCI$t[, (K+2):(2*(K+1))], 2, mean) # Mean for treatment group across 41:80 (second half contains risk1)
# 
# ll_risk0 <- apply(te_all_timepoints_withCI$t[, 1:(K+1)], 2, function(x) quantile(x, 0.025))
# ul_risk0 <- apply(te_all_timepoints_withCI$t[, 1:(K+1)], 2, function(x) quantile(x, 0.975))
# 
# ll_risk1 <- apply(te_all_timepoints_withCI$t[, (K+2):(2*(K+1))], 2, function(x) quantile(x, 0.025))
# ul_risk1 <- apply(te_all_timepoints_withCI$t[, (K+2):(2*(K+1))], 2, function(x) quantile(x, 0.975))
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
#   theme(axis.text = element_text(size=14), legend.position = c(0.2, 0.8),
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


# (ii) IPTW ----------------------------------------------------------------
print('IPTW')
## The denominator of the IPTW for each individual is the probability of receiving their observed(!) treatment value, given their confounder history.
## Here, we deal with treatment that is defined at baseline, not time-varying uptake of treatment => one constant baseline treatment weight per individual
## And, the numerator therefore is simply the proportion of patients who were actually treated.

# Calculate denominator for "probability of receiving their observed(!) treatment value"
iptw_denom_formula <- as.formula(paste("exp_bin_treat ~ ", paste(covariate_names, collapse = " + "), "+ strat(strat_cat_region)"))
iptw.denom <- parglm(iptw_denom_formula, 
                     family = binomial(link = 'logit'),
                     data = df_long_months) 
# summary(iptw.denom)
df_long_months$iptw_denom <- predict(iptw.denom, df_long_months, type="response")

# Estimate stabilized weights
df_long_months$w_treat_stab <- ifelse(df_long_months$exp_bin_treat==1,
                                      mean(df_long_months$exp_bin_treat)/df_long_months$iptw_denom,
                                      (1-mean(df_long_months$exp_bin_treat))/(1-df_long_months$iptw_denom))
## this results in 1 weight estimate for each individual (-> constant over time)


# (ii) IPTW: Check the weights and the balance ---------------------------------------
print('IPTW: Check the weights and the balance')
summary(df_long_months$w_treat_stab)
sd(df_long_months$w_treat_stab)

# Create subsets of data, according to treat_b status
treat_b0 <- subset(df_long_months,exp_bin_treat==0)
treat_b1 <- subset(df_long_months,exp_bin_treat==1)

# List of variables to compare (include more, not only covariate_names, e.g. also the numeric variables for HbA1c and TC/HDL ratio)
varlist <- names(df) %>%
  grep("^cov", ., value = TRUE) %>% 
  # but exclude those that make no sense: 
  ## cov_num_age_spline and cov_cat_age, represented by cov_num_age
  setdiff(c("cov_num_age_spline", "cov_cat_age")) 
# print(varlist)

# Create function to take mean difference, or SMD for age
meanfctn <- function(x){
  if(x == "cov_num_age"){
    t0 <- treat_b0[[x]]
    t1 <- treat_b1[[x]]
    md <- (mean(t1) - mean(t0))/sd(t1)
  }else{
    t0 <- treat_b0[[x]]
    t1 <- treat_b1[[x]]
    md <- mean(t1) - mean(t0)}
  return(c(var = x, md = md))
}

# Calculate mean differences for covariates (SMD for age) after weighting
wmean_fctn <- function(x){
  # Convert factors to numeric
  if(is.factor(treat_b0[[x]]) || is.character(treat_b0[[x]])) {
    treat_b0[[x]] <- as.numeric(treat_b0[[x]])
  }
  if(is.factor(treat_b1[[x]]) || is.character(treat_b1[[x]])) {
    treat_b1[[x]] <- as.numeric(treat_b1[[x]])
  }
  
  # Handle missing values
  if(any(is.na(treat_b0[[x]])) || any(is.na(treat_b1[[x]]))){
    return(c(var = x, md = NA))  # Return NA if missing values exist
  }
  
  if(x == "cov_num_age"){
    md <- (weighted.mean(treat_b1[[x]], treat_b1$w_treat_stab, na.rm = TRUE) - 
             weighted.mean(treat_b0[[x]], treat_b0$w_treat_stab, na.rm = TRUE)) /
      sqrt(wtd.var(treat_b1[[x]], treat_b1$w_treat_stab, na.rm = TRUE))
  } else {
    t0 <- weighted.mean(treat_b0[[x]], treat_b0$w_treat_stab, na.rm = TRUE)
    t1 <- weighted.mean(treat_b1[[x]], treat_b1$w_treat_stab, na.rm = TRUE)
    md <- t1 - t0
  }
  
  return(c(var = x, md = md))
}

# Create the covariate plot
covplot_w <- lapply(varlist, wmean_fctn) %>% do.call(rbind,.) %>% as.data.frame()
covplot_w$md <- as.numeric(covplot_w$md)

# Plot it
covplot_weighted <- ggplot(data = covplot_w) +
  geom_point(aes(x = md, y = var), color = "steelblue") + scale_x_continuous(limits = c(-0.9, 0.9)) +
  geom_vline(xintercept = 0) +
  labs(y = "Covariates", x = "Mean Difference", title = "Covariate Balance Plot")
## However, since the numeric variables are not needed for modelling, they were not modified in the dummy data, and thus contain lots of missing


# (ii) IPTW: Truncate weights --------------------------------------------------------
print('IPTW: Truncate weights')
# Truncate stabilized weights at the 99th percentile
threshold_99 <- quantile(df_long_months$w_treat_stab, 0.99)
df_long_months$w_treat_stab_99 <- df_long_months$w_treat_stab
df_long_months$w_treat_stab_99[df_long_months$w_treat_stab_99 > threshold_99] <- threshold_99

###  Min, 25th percentile, median, mean, SD, 75th percentile, and max:
summary(df_long_months$w_treat_stab)
sd(df_long_months$w_treat_stab)
summary(df_long_months$w_treat_stab_99)
sd(df_long_months$w_treat_stab_99)


# (ii) IPTW: Fit pooled logistic regression -------------------------------------
print('IPTW: Fit pooled logistic regression and create plot_cum_risk_iptw')
# Fit pooled logistic regression, with stabilized weights (currently only IPW for treatment, i.e. baseline confounding => same weights across time)
# Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time -> implement!
# Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
# The stratification variable has been accounted for in the weights - do I still need to include it here, i.e. the have different baseline hazards by region?
plr_model_severecovid <- parglm(outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr), 
                                family = binomial(link = 'logit'),
                                data = df_long_months[df_long_months$censor==0 & df_long_months$comp_event == 0,],
                                weights = df_long_months[df_long_months$censor==0 & df_long_months$comp_event == 0,]$w_treat_stab_99)
# summary(plr_model_severecovid)

### Transform estimates to risks at each time point in each group ###

# Create dataset with all time points for each individual under each treatment level
df_pred <- df_long_months %>% filter(month == 0) %>% dplyr::select(-month) %>% crossing(month = 0:(K-1))
df_pred$monthsqr <- df_pred$month^2

# Control group (everyone untreated) with predicted discrete-time hazards
df_pred0 <- df_pred
df_pred0$exp_bin_treat <- 0
df_pred0$p.event0 <- predict(plr_model_severecovid, df_pred0, type="response")

# Treatment group (everyone treated) with predicted discrete-time hazards
df_pred1 <- df_pred
df_pred1$exp_bin_treat <- 1
df_pred1$p.event1 <- predict(plr_model_severecovid, df_pred1, type="response")

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

# Put all in 1 data frame
graph.pred <- merge(risk0, risk1, by=c("month"))

# Edit data frame to reflect that risks are estimated at the END of each interval
graph.pred$time_0 <- graph.pred$month + 1
zero <- data.frame(cbind(0,0,0,1,0,0))
zero <- setNames(zero,names(graph.pred))
graph <- rbind(zero, graph.pred) ## can be used for the cumulative incidence plot, but without 95% CI (for that, see below)

# Add RD and RR
graph$rd <- graph$risk1-graph$risk0
graph$rr <- graph$risk1/graph$risk0

# Extract overall risk estimate (end of follow-up K-1), but without 95% CI (for that, see below)
risk0 <- graph$risk0[which(graph$month==K-1)] 
risk1 <- graph$risk1[which(graph$month==K-1)]
rd <- graph$rd[which(graph$month==K-1)]
rr <- graph$rr[which(graph$month==K-1)]

### Construct marginal parametric cumulative incidence (risk) curves (without CIs) ###

# Create plot (without CIs)
plot_cum_risk_iptw <- ggplot(graph,
                  aes(x=time_0, y=risk)) + 
  geom_line(aes(y = risk1, 
                color = "Metformin"), linewidth = 1.5) +
  geom_line(aes(y = risk0,
                color = "No Metformin"), linewidth = 1.5) +
  xlab("Months") +
  # scale_x_continuous(limits = c(0, 24),
  #                    breaks=c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)) +
  ylab("Cumulative Incidence (%)") + 
  # scale_y_continuous(limits=c(0, 0.125), # format y axis
  #                    breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125),
  #                    labels=c("0.0%", "2.5%", "5.0%",
  #                             "7.5%", "10.0%", "12.5%")) +
  theme_minimal()+ # set plot theme elements
  theme(axis.text = element_text(size=14), legend.position = c(0.2, 0.8),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())+
  scale_color_manual(values=c("#E7B800","#2E9FDF"),
                     breaks=c('No Metformin', 'Metformin'))


# (ii) IPTW: Add 95% CI using robust standard errors and the delta method for RD and RR | See code from Will/Fizz ----
# print('IPTW: create te_plr_iptw_rse_tbl')
# # Compute robust standard errors
# vcov_robust <- vcovHC(plr_model_severecovid, type = "HC0")
# robust_se <- sqrt(diag(vcov_robust))
# coefs <- coef(plr_model_severecovid)
# 
# # standard errors for the two arms
# se_risk0 <- sqrt(vcov_robust["(Intercept)", "(Intercept)"])
# se_risk1 <- sqrt(vcov_robust["exp_bin_treat", "exp_bin_treat"])
# 
# # confidence intervals for predicted risks using normal approximation
# ci_risk0 <- c(risk0 - 1.96 * se_risk0, risk0 + 1.96 * se_risk0)
# ci_risk1 <- c(risk1 - 1.96 * se_risk1, risk1 + 1.96 * se_risk1)
# 
# # Formula for risk difference (RD)
# rd_formula <- "1 / (1 + exp(-((Intercept)))) - 1 / (1 + exp(-((Intercept) + exp_bin_treat)))"
# 
# # Delta method for RD
# se_rd <- deltaMethod(coefs, rd_formula, vcov_robust)$SE
# ci_rd_lower <- rd - 1.96 * se_rd
# ci_rd_upper <- rd + 1.96 * se_rd
# 
# # Formula for risk ratio (RR)
# rr_formula <- "(1 / (1 + exp(-((Intercept) + exp_bin_treat)))) / (1 / (1 + exp(-((Intercept)))))"
# 
# # Delta method for RR
# se_rr <- deltaMethod(coefs, rr_formula, vcov_robust)$SE
# ci_rr_lower <- rr - 1.96 * se_rr
# ci_rr_upper <- rr + 1.96 * se_rr
# 
# # Create results table
# te_plr_iptw_rse_tbl <- data.frame(
#   Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
#   Estimate = c(risk0, risk1, rd, rr),
#   Lower_CI = c(ci_risk0[1], ci_risk1[1], ci_rd_lower, ci_rr_lower),
#   Upper_CI = c(ci_risk0[2], ci_risk1[2], ci_rd_upper, ci_rr_upper)
# )


# (ii) IPTW & IPCW: Add censoring event weights -----------------------------------
print('IPTW & IPCW: Add censoring event weights')
## In our case we want to censor for LTFU, i.e. "the effect had no one been lost to follow-up"
## Informally, an uncensored individual’s inverse probability of censoring weight is the inverse of their probability of remaining uncensored given their treatment and covariate history. 
## In the context of loss to follow-up, each individual who was not lost to follow-up receives a weight that is proportional to the inverse of the probability of not being lost to follow-up, given their specific history.
## We will consider stabilized IP weights for censoring, in which the numerator of the weights is the probability of remaining uncensored given an individual’s treatment.

## Once an individual is censored, they receive a censoring weight of 0 from that point onward. 
## In this simplified example, we assume censoring due to loss to follow-up depends only on the baseline treatment and baseline covariates.
## Once we can incorporate time-varying covariates (from OpenSAFELY), we will adapt, including treat_1, i.e. first time treated.

# Indicator for ever being censored due to loss to follow-up
df_long_months <- df_long_months %>%
  group_by(patient_id) %>%
  mutate(censor_any = ifelse(is.na(censor), NA, max(censor, na.rm = T))) %>%
  ungroup()

## Fit a pooled logistic regression model for the denominator of the IP weights for censoring 
## This model should predict the probability of remaining uncensored (not being lost to follow-up) at each timepoint. 
## Continuous time, here months (modeled using linear and quadratic terms), will serve as the time scale for this model. 
## We do not include any product terms in the model.
ipcw_denom_formula <- as.formula(paste("censor == 0 ~ exp_bin_treat + month + monthsqr +", paste(covariate_names, collapse = " + "), "+ strat(strat_cat_region)"))
ipcw.denom <- parglm(ipcw_denom_formula, 
                  family = binomial(link = 'logit'),
                  data = df_long_months) 

# Obtain predicted probabilities of being uncensored for denominator
df_long_months$ipcw_denom <- predict(ipcw.denom, df_long_months, type="response")

## For the numerator we now fit a model to include the time-varying aspect of the weights - but without any covariates
# strat(strat_cat_region)?
ipcw_num_formula <- as.formula(paste("censor == 0 ~ exp_bin_treat + month + monthsqr"))
ipcw.num <- glm(ipcw_num_formula, 
                family = binomial(link = 'logit'),
                data = df_long_months) 

# Obtain predicted probabilities of being uncensored for numerator
df_long_months$ipcw_num <- predict(ipcw.num, df_long_months, type="response")

### Estimate stabilized inverse probability weights for censoring ###
# Take cumulative products starting at baseline
df_long_months <- df_long_months %>%
  group_by(patient_id) %>%
  mutate(w_cens_stab = cumprod(ipcw_num)/cumprod(ipcw_denom)) %>%
  ungroup() %>%
  # ensure that individuals who are not censored are assigned a weight of 1 - but does not exist in our case due to the data setup
  mutate(w_cens_stab = ifelse(is.na(censor), 1, w_cens_stab)) 

# Take product of weight for censoring and treatment for each individual (take the untruncated w_treat_stab from above)
df_long_months$w_treat_cens_stab <- df_long_months$w_treat_stab * df_long_months$w_cens_stab

### Truncate final stabilized weight at the 99th percentile ###
threshold_99 <- quantile(df_long_months$w_treat_cens_stab, 0.99)
df_long_months$w_treat_cens_stab_99 <- df_long_months$w_treat_cens_stab
df_long_months$w_treat_cens_stab_99[df_long_months$w_treat_cens_stab_99 > threshold_99] <- threshold_99

###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: truncated weights ###
summary(df_long_months$w_treat_cens_stab_99)
sd(df_long_months$w_treat_cens_stab_99)
# df_long_months %>% 
#   filter(is.na(df_long_months$w_treat_cens_stab_99)) %>% 
#   View()

# (ii) IPTW & IPCW: Fit pooled logistic regression -------------------------
print('IPTW & IPCW: Fit pooled logistic regression and create plot_cum_risk_iptw_ipcw')
# Fit pooled logistic regression, with stabilized weights for IPTW and IPCW
# Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time. -> implement!
# Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
# The stratification variable has been accounted for in the weights - do I still need to include it here, i.e. the have different baseline hazards by region?
plr_model_severecovid <- parglm(outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr), 
                                family = binomial(link = 'logit'),
                                data = df_long_months[df_long_months$censor==0 & df_long_months$comp_event == 0,],
                                weights = df_long_months[df_long_months$censor==0 & df_long_months$comp_event == 0,]$w_treat_cens_stab_99)
# summary(plr_model_severecovid)

### Transform estimates to risks at each time point in each group ###

# Create dataset with all time points for each individual under each treatment level
df_pred <- df_long_months %>% filter(month == 0) %>% dplyr::select(-month) %>% crossing(month = 0:(K-1))
df_pred$monthsqr <- df_pred$month^2

# Control group (everyone untreated) with predicted discrete-time hazards
df_pred0 <- df_pred
df_pred0$exp_bin_treat <- 0
df_pred0$p.event0 <- predict(plr_model_severecovid, df_pred0, type="response")

# Treatment group (everyone treated) with predicted discrete-time hazards
df_pred1 <- df_pred
df_pred1$exp_bin_treat <- 1
df_pred1$p.event1 <- predict(plr_model_severecovid, df_pred1, type="response")

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

# Put all in 1 data frame
graph.pred <- merge(risk0, risk1, by=c("month"))

# Edit data frame to reflect that risks are estimated at the END of each interval
graph.pred$time_0 <- graph.pred$month + 1
zero <- data.frame(cbind(0,0,0,1,0,0))
zero <- setNames(zero,names(graph.pred))
graph <- rbind(zero, graph.pred) ## can be used for the cumulative incidence plot, but without 95% CI (for that, see below)

# Add RD and RR
graph$rd <- graph$risk1-graph$risk0
graph$rr <- graph$risk1/graph$risk0

# Extract overall risk estimate (end of follow-up K-1), but without 95% CI (for that, see below)
risk0 <- graph$risk0[which(graph$month==K-1)] 
risk1 <- graph$risk1[which(graph$month==K-1)]
rd <- graph$rd[which(graph$month==K-1)]
rr <- graph$rr[which(graph$month==K-1)]

### Construct marginal parametric cumulative incidence (risk) curves (without CIs) ###

# Create plot (without CIs)
plot_cum_risk_iptw_ipcw <- ggplot(graph,
                            aes(x=time_0, y=risk)) + 
  geom_line(aes(y = risk1, 
                color = "Metformin"), linewidth = 1.5) +
  geom_line(aes(y = risk0,
                color = "No Metformin"), linewidth = 1.5) +
  xlab("Months") +
  # scale_x_continuous(limits = c(0, 24),
  #                    breaks=c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)) +
  ylab("Cumulative Incidence (%)") + 
  # scale_y_continuous(limits=c(0, 0.125), # format y axis
  #                    breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125),
  #                    labels=c("0.0%", "2.5%", "5.0%",
  #                             "7.5%", "10.0%", "12.5%")) +
  theme_minimal()+ # set plot theme elements
  theme(axis.text = element_text(size=14), legend.position = c(0.2, 0.8),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())+
  scale_color_manual(values=c("#E7B800","#2E9FDF"),
                     breaks=c('No Metformin', 'Metformin'))


# (ii) IPTW & IPCW: Add 95% CI using robust standard errors and the delta method for RD and RR | See code from Will/Fizz ----
# print('IPTW & IPCW: create te_plr_iptw_ipcw_boot_tbl')
# # Compute robust standard errors
# vcov_robust <- vcovHC(plr_model_severecovid, type = "HC0")
# robust_se <- sqrt(diag(vcov_robust))
# coefs <- coef(plr_model_severecovid)
# 
# # standard errors for the two arms
# se_risk0 <- sqrt(vcov_robust["(Intercept)", "(Intercept)"])
# se_risk1 <- sqrt(vcov_robust["exp_bin_treat", "exp_bin_treat"])
# 
# # confidence intervals for predicted risks using normal approximation
# ci_risk0 <- c(risk0 - 1.96 * se_risk0, risk0 + 1.96 * se_risk0)
# ci_risk1 <- c(risk1 - 1.96 * se_risk1, risk1 + 1.96 * se_risk1)
# 
# # Formula for risk difference (RD)
# rd_formula <- "1 / (1 + exp(-((Intercept)))) - 1 / (1 + exp(-((Intercept) + exp_bin_treat)))"
# 
# # Delta method for RD
# se_rd <- deltaMethod(coefs, rd_formula, vcov_robust)$SE
# ci_rd_lower <- rd - 1.96 * se_rd
# ci_rd_upper <- rd + 1.96 * se_rd
# 
# # Formula for risk ratio (RR)
# rr_formula <- "(1 / (1 + exp(-((Intercept) + exp_bin_treat)))) / (1 / (1 + exp(-((Intercept)))))"
# 
# # Delta method for RR
# se_rr <- deltaMethod(coefs, rr_formula, vcov_robust)$SE
# ci_rr_lower <- rr - 1.96 * se_rr
# ci_rr_upper <- rr + 1.96 * se_rr
# 
# # Create results table
# te_plr_iptw_ipcw_rse_tbl <- data.frame(
#   Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
#   Estimate = c(risk0, risk1, rd, rr),
#   Lower_CI = c(ci_risk0[1], ci_risk1[1], ci_rd_lower, ci_rr_lower),
#   Upper_CI = c(ci_risk0[2], ci_risk1[2], ci_rd_upper, ci_rr_upper)
# )


# (ii) IPTW & IPCW: Add 95% CI for risks, RD and RR using bootstrapping ----

## add more printouts and checks to help debugging when running it on real data!

# Create a function to obtain risks, RD, and RR from each bootstrap sample - and return 1 risk time point
study_ids <- data.frame(patient_id = unique(df_long_months$patient_id))

te_iptw_ipcw_rd_rr_withCI <- function(data, indices) {
  
  # capture and count unsuccessful bootstraps
  tryCatch({
  
  # Resample the dataset
  boot.ids <- data[indices, , drop = FALSE] # random resampling based on data = study_ids, & force the subset to stay a data frame instead of accidentally becoming a vector
  d <- left_join(boot.ids, df_long_months, by = "patient_id", relationship = "many-to-many")
  # d <- df_long_months %>% filter(patient_id %in% boot.ids$patient_id) # Instead of left_join, filter the df_long_months, probably faster, but does not respect "bootstrapping with replacement"
  
  ### (1) Calculate IPTW
  iptw_denom_formula <- as.formula(paste("exp_bin_treat ~ ", paste(covariate_names, collapse = " + "), "+ strat(strat_cat_region)"))
  iptw.denom <- parglm(iptw_denom_formula, 
                       family = binomial(link = 'logit'),
                       data = d) 
  d$iptw_denom <- predict(iptw.denom, d, type="response")
  
  # check that no NA in predictions that will mess up calculations later down the road
  if (any(is.na(d$iptw_denom))) {
    stop("NA in IPTW denominator predictions")
  }
  
  # Estimate stabilized weights
  df_long_months$w_treat_stab <- ifelse(d$exp_bin_treat==1,
                                        mean(d$exp_bin_treat)/d$iptw_denom,
                                        (1-mean(d$exp_bin_treat))/(1-d$iptw_denom))
  ## this results in 1 weight estimate for each individual (-> constant over time)
  
  ### (2) Calculate IPCW
  ipcw_denom_formula <- as.formula(paste("censor == 0 ~ exp_bin_treat + month + monthsqr +", paste(covariate_names, collapse = " + "), "+ strat(strat_cat_region)"))
  ipcw.denom <- parglm(ipcw_denom_formula, 
                       family = binomial(link = 'logit'),
                       data = d) 
  
  # Obtain predicted probabilities of being uncensored for denominator
  d$ipcw_denom <- predict(ipcw.denom, d, type="response")
  if (any(is.na(d$ipcw_denom))) {
    stop("NA in IPCW denominator predictions")
  }
  
  ## For the numerator we now fit a model to include the time-varying aspect of the weights
  ipcw_num_formula <- as.formula(paste("censor == 0 ~ exp_bin_treat + month + monthsqr"))
  ipcw.num <- parglm(ipcw_num_formula, 
                  family = binomial(link = 'logit'),
                  data = d) 
  
  # Obtain predicted probabilities of being uncensored for numerator
  d$ipcw_num <- predict(ipcw.num, d, type="response")
  if (any(is.na(d$ipcw_num))) {
    stop("NA in IPCW numerator predictions")
  }
  
  # Take cumulative products starting at baseline
  d <- d %>%
    group_by(patient_id) %>%
    mutate(w_cens_stab = cumprod(ipcw_num)/cumprod(ipcw_denom)) %>%
    ungroup()
  
  ### (3) Multiply IPTW * IPCW
  d$w_treat_cens_stab <- d$w_treat_stab * d$w_cens_stab
  if (any(is.na(d$w_treat_cens_stab))) {
    cat("NA in final weight w_treat_cens_stab:", sum(is.na(d$w_treat_cens_stab)), "\n")
  }
  # d <- d[!is.na(d$w_treat_cens_stab), ] # Drop NA weights if there are just a few
  
  ### (4) Truncate final stabilized weight at the 99th percentile
  threshold_99 <- quantile(d$w_treat_cens_stab, 0.99)
  d$w_treat_cens_stab_99 <- d$w_treat_cens_stab
  d$w_treat_cens_stab_99[d$w_treat_cens_stab_99 > threshold_99] <- threshold_99
  
  ### (5) Fit pooled logistic regression, with stabilized weights for IPTW and IPCW
  # Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time. -> implement!
  # Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
  # The stratification variable has been accounted for in the weights - do I still need to include it here, i.e. the have different baseline hazards by region?
  plr_model_severecovid <- parglm(outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr), 
                                  family = binomial(link = 'logit'),
                                  data = d[d$censor==0 & d$comp_event == 0,],
                                  weights = d[d$censor==0 & d$comp_event == 0,]$w_treat_cens_stab_99)
  
  ### (6) Transform estimates to risks at each time point in each group
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- d %>% filter(month == 0) %>% dplyr::select(-month) %>% crossing(month = 0:(K-1))
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
  treat0.surv <- treat0 %>% group_by(patient_id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
  treat1.surv <- treat1 %>% group_by(patient_id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
  
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
  
  }, error = function(e) {
    # If any error occurs, count as a failed bootstrap
    cat("Bootstrap sample failed due to:", conditionMessage(e), "\n")
    boot_failures <<- boot_failures + 1
    return(rep(NA, 4)) # Return NA values to avoid breaking the boot()
  })
  
}

set.seed(423)
boot_failures <- 0 # Reset failure counter

# run it
te_iptw_ipcw_rd_rr_withCI_boot <- boot(data = study_ids, statistic = te_iptw_ipcw_rd_rr_withCI, R = R)
cat("Number of failed bootstrap samples:", boot_failures, "\n")

# Function to extract bootstrapped confidence intervals - and keep a column with the original point estimates (without CI)
extract_ci_boot <- function(boot_obj, index) {
  ci <- boot.ci(boot_obj, conf = 0.95, type = "perc", index = index)
  if (!is.null(ci$percent)) {
    return(c(ci$percent[4], ci$percent[5])) # Lower and Upper CI for the bootstrapped values
  } else {
    return(c(NA, NA)) # for the original value
  }
}

# Create the results table
te_plr_iptw_ipcw_boot_tbl <- data.frame(
  Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
  Estimate_original = te_iptw_ipcw_rd_rr_withCI_boot$t0,  # Original estimates
  Estimate_boot = colMeans(te_iptw_ipcw_rd_rr_withCI_boot$t),  # Bootstrapped mean estimate
  Lower_CI = sapply(1:4, function(i) extract_ci_boot(te_iptw_ipcw_rd_rr_withCI_boot, i)[1]),
  Upper_CI = sapply(1:4, function(i) extract_ci_boot(te_iptw_ipcw_rd_rr_withCI_boot, i)[2])
)


# Save output -------------------------------------------------------------
print('Save output')
# Treatment effect (risk) estimates
# write.csv(te_plr_stand_rse_tbl, file = here::here("output", "te", "pooled_log_reg", "te_plr_stand_rse_severecovid.csv"))
# write.csv(te_plr_stand_boot_tbl, file = here::here("output", "te", "pooled_log_reg", "te_plr_stand_boot_severecovid.csv"))
# write.csv(te_plr_iptw_rse_tbl, file = here::here("output", "te", "pooled_log_reg", "te_plr_iptw_rse_severecovid.csv"))
# write.csv(te_plr_iptw_ipcw_rse_tbl, file = here::here("output", "te", "pooled_log_reg", "te_plr_iptw_ipcw_rse_severecovid.csv"))
write.csv(te_plr_iptw_ipcw_boot_tbl, file = here::here("output", "te", "pooled_log_reg", "te_plr_iptw_ipcw_boot_severecovid.csv"))

# Marginal parametric cumulative incidence (risk) curves from plr model, with 95% CI
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_cum_risk_iptw_severecovid.png"), plot_cum_risk_iptw, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_cum_risk_iptw_ipcw_severecovid.png"), plot_cum_risk_iptw_ipcw, width = 20, height = 20, units = "cm")

# Love/SMD plot
ggsave(filename = here::here("output", "te", "pooled_log_reg", "covplot_weighted.png"), covplot_weighted, width = 20, height = 20, units = "cm")
