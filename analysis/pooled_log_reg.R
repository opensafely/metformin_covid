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
## h) Discuss different spline methods, ns() versus rcs(), but rcs() does not work within parglm


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
source(here::here("analysis", "functions", "fn_smd_before.R"))
source(here::here("analysis", "functions", "fn_smd_after.R"))


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
# Expand the dataset into intervals and follow the rules see in fn_expand_intervals.R script
start_date_variable <- "landmark_date" # start expanding at landmark date, not elig_date_t2dm (due to landmark design)
stop_date_columns <- c("out_date_severecovid_afterlandmark", "out_date_death_afterlandmark", "cens_date_ltfu_afterlandmark")
outcome_date_variable <- "out_date_severecovid_afterlandmark"
comp_date_variable <- "out_date_death_afterlandmark"
censor_date_variable <- "cens_date_ltfu_afterlandmark"

# Apply the function, currently only using months
# Use studyend_date from metadates.R import
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
#                 cens_date_ltfu_afterlandmark, stop_date, start_date_month, end_date_month, month, outcome, comp_event, censor, is_event_interval, is_outcome_event, is_cens_event,
#                 is_comp_event, followup_stop,
#                 cov_cat_sex, cov_num_age
#                 ) %>%
#   # dplyr::filter(!is.na(out_date_death_afterlandmark)) %>%
#   # dplyr::filter(!is.na(cens_date_ltfu_afterlandmark)) %>%
#   # dplyr::filter(!is.na(out_date_severecovid_afterlandmark)) %>%
#   # dplyr::filter(is.na(censor)) %>% # should be empty
#   # dplyr::filter(stop_date == "2022-04-01") %>%
#   # dplyr::filter(is_event_interval == TRUE & !is.na(out_date_severecovid_afterlandmark)) %>% # should only have 1 row per person
#   # dplyr::filter(followup_stop == TRUE) %>% # should be empty, otherwise it means that there are two defining event intervals
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
df_long_months$monthsqr <- df_long_months$month^2 # add months square to model time in PLR
R <- 10 # Total bootstraps (ideally >500)


# IPTW: Fit treatment model and truncate the weights ----------------------
print('IPTW')
## The denominator of the IPTW for each individual is the probability of receiving their observed(!) treatment value, given their confounder history.
## Here, we deal with treatment that is defined at baseline, not time-varying uptake of treatment => one constant baseline treatment weight per individual
## And, the numerator therefore is simply the proportion of patients who were actually treated.

# Calculate denominator for "probability of receiving their observed(!) treatment value"
iptw_denom_formula <- as.formula(paste("exp_bin_treat ~ ", paste(covariate_names, collapse = " + "), "+ strat(strat_cat_region)"))
iptw.denom <- parglm(iptw_denom_formula, 
                     family = binomial(link = 'logit'),
                     data = df_long_months[df_long_months$month == 0, ]) 
# summary(iptw.denom)
df_long_months$iptw_denom <- predict(iptw.denom, df_long_months, type="response")

# Estimate stabilized weights
df_long_months$w_treat_stab <- ifelse(df_long_months$exp_bin_treat==1,
                                      mean(df_long_months$exp_bin_treat)/df_long_months$iptw_denom,
                                      (1-mean(df_long_months$exp_bin_treat))/(1-df_long_months$iptw_denom))
## this results in 1 weight estimate for each individual (-> constant over time)

print('IPTW: Truncate weights at 1st and 99th percentile')
lower_threshold <- quantile(df_long_months$w_treat_stab, 0.01, na.rm = TRUE)
upper_threshold <- quantile(df_long_months$w_treat_stab, 0.99, na.rm = TRUE)
df_long_months$w_treat_stab_trunc <- df_long_months$w_treat_stab
df_long_months$w_treat_stab_trunc[df_long_months$w_treat_stab_trunc < lower_threshold] <- lower_threshold
df_long_months$w_treat_stab_trunc[df_long_months$w_treat_stab_trunc > upper_threshold] <- upper_threshold

## Min, 25th percentile, median, mean, SD, 75th percentile, and max:
summary(df_long_months$w_treat_stab)
sd(df_long_months$w_treat_stab)
summary(df_long_months$w_treat_stab_trunc)
sd(df_long_months$w_treat_stab_trunc)


# IPTW: Fit pooled logistic regression -------------------------------------
print('IPTW: Fit pooled logistic regression and create plot_cum_risk_iptw')
# Fit pooled logistic regression, with stabilized weights (currently only IPW for treatment, i.e. baseline confounding => same weights across time)
# Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time -> implement!
# Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
# The stratification variable has been accounted for in the weights - do I still need to include it here, i.e. the have different baseline hazards by region?
plr_model_severecovid <- parglm(outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr), 
                                family = binomial(link = 'logit'),
                                data = df_long_months[df_long_months$censor==0 & df_long_months$comp_event == 0,],
                                weights = df_long_months[df_long_months$censor==0 & df_long_months$comp_event == 0,]$w_treat_stab_trunc)
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
  ylab("Cumulative Incidence (%)") + 
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


# IPTW & IPCW: Add censoring event weights ----------------------------------
print('IPTW & IPCW: Add censoring event weights')
## In our case we want to censor for LTFU, i.e. "the effect had no one been lost to follow-up"
## Informally, an uncensored individual’s inverse probability of censoring weight is the inverse of their probability of remaining uncensored given their treatment and covariate history. 
## In the context of loss to follow-up, each individual who was not lost to follow-up receives a weight that is proportional to the inverse of the probability of not being lost to follow-up, given their specific history.
## We will consider stabilized IP weights for censoring, in which the numerator of the weights is the probability of remaining uncensored given an individual’s treatment.

## Once an individual is censored, they receive a censoring weight of 0 from that point onward. 
## In this simplified example, we assume censoring due to loss to follow-up depends only on the baseline treatment and baseline covariates.
## Once we can incorporate time-varying covariates (from OpenSAFELY), we will adapt, including treat_1, i.e. first time treated.

# Indicator for ever being censored due to loss to follow-up (just for descriptive)
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
ipcw_num_formula <- as.formula(paste("censor == 0 ~ month + monthsqr"))
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
  ungroup()

# Take product of weight for censoring and treatment for each individual (take the untruncated w_treat_stab from above)
df_long_months$w_treat_cens_stab <- df_long_months$w_treat_stab * df_long_months$w_cens_stab

### Truncate final stabilized weight at the 1st and 99th percentile ###
lower_threshold <- quantile(df_long_months$w_treat_cens_stab, 0.01, na.rm = TRUE)
upper_threshold <- quantile(df_long_months$w_treat_cens_stab, 0.99, na.rm = TRUE)
df_long_months$w_treat_cens_stab_trunc <- df_long_months$w_treat_cens_stab
df_long_months$w_treat_cens_stab_trunc[df_long_months$w_treat_cens_stab_trunc < lower_threshold] <- lower_threshold
df_long_months$w_treat_cens_stab_trunc[df_long_months$w_treat_cens_stab_trunc > upper_threshold] <- upper_threshold

###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: truncated weights ###
summary(df_long_months$w_treat_cens_stab_trunc)
sd(df_long_months$w_treat_cens_stab_trunc)


# IPTW & IPCW: Check the summary weights and the balance ---------------------------------------
print('IPTW: Check the weights and the balance')

# Create subsets of data, according to treat_b status
treat_b0 <- subset(df_long_months, exp_bin_treat==0)
treat_b1 <- subset(df_long_months, exp_bin_treat==1)

# List of variables to compare (include more, not only covariate_names, e.g. also the numeric variables for HbA1c and TC/HDL ratio)
varlist <- names(df) %>%
  grep("^cov", ., value = TRUE) %>% 
  # but exclude those that make no sense: 
  ## cov_num_age_spline and cov_cat_age, represented by cov_num_age
  setdiff(c("cov_num_age_spline", "cov_cat_age")) 
# print(varlist)

### Check standardized mean difference BEFORE weighting
smd_before <- lapply(varlist, fn_smd_before)
smd_before_df <- as.data.frame(do.call(rbind, smd_before), stringsAsFactors = FALSE)
smd_before_df <- transform(smd_before_df,
                           var = as.character(var),
                           type = as.character(type),
                           smd = as.numeric(smd))
covplot_unweighted <- ggplot(smd_before_df, aes(x = smd, y = reorder(var, smd))) +
  geom_point(color = "steelblue", size = 3) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "red") +
  labs(title = "Love Plot: Standardized Mean Differences before weighting",
       x = "Standardized Mean Difference (SMD)",
       y = "Covariate") +
  theme_minimal()

### Check standardized mean difference AFTER weighting
w_treat_0 <- treat_b0$w_treat_cens_stab_trunc
w_treat_1 <- treat_b1$w_treat_cens_stab_trunc
smd_after <- lapply(varlist, 
                    fn_smd_after, 
                    w_treat_0 = w_treat_0, 
                    w_treat_1 = w_treat_1)
smd_after_df <- as.data.frame(do.call(rbind, smd_after), stringsAsFactors = FALSE)
smd_after_df <- transform(smd_after_df,
                          var = as.character(var),
                          type = as.character(type),
                          smd = as.numeric(smd))
covplot_weighted <- ggplot(smd_after_df, aes(x = smd, y = reorder(var, smd))) +
  geom_point(color = "darkorange", size = 3) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "red") +
  labs(title = "Love Plot: Standardized Mean Differences after weighting",
       x = "Standardized Mean Difference (SMD)",
       y = "Covariate") +
  theme_minimal()
## Note re dummy data: cov_num_tc_hdl_ratio, _hba1c, _bmi, were not modified in dummy data adaptation since not used for modelling (only their categorical counterparts) => missing


# IPTW & IPCW: Fit pooled logistic regression -------------------------
print('IPTW & IPCW: Fit pooled logistic regression and create plot_cum_risk_iptw_ipcw')
# Fit pooled logistic regression, with stabilized weights for IPTW and IPCW
# Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time. -> implement!
# Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
# The stratification variable has been accounted for in the weights - do I still need to include it here, i.e. the have different baseline hazards by region?
plr_model_severecovid <- parglm(outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr), 
                                family = binomial(link = 'logit'),
                                data = df_long_months[df_long_months$censor==0 & df_long_months$comp_event == 0,],
                                weights = df_long_months[df_long_months$censor==0 & df_long_months$comp_event == 0,]$w_treat_cens_stab_trunc)
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
  ylab("Cumulative Incidence (%)") + 
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


# IPTW & IPCW: Add 95% CI using robust standard errors and the delta method for RD and RR | See code from Will/Fizz ----
# print('IPTW & IPCW: create te_plr_iptw_ipcw_boot_tbl')


# IPTW & IPCW: Add 95% CI using bootstrapping ----
study_ids <- data.frame(patient_id = unique(df_long_months$patient_id))

te_iptw_ipcw_rd_rr_withCI <- function(data, indices) {
  
  # capture and count unsuccessful bootstraps
  tryCatch({
  
  # Resample patient IDs (with replacement)
  boot.ids <- data[indices, , drop = FALSE] %>%
    mutate(bootstrap_draw = row_number())  # Assign a unique row for each draw
  
  # Join back the full time records **for each draw** (i.e. same patient is allowed to be sampled several times)
  d <- boot.ids %>%
    left_join(df_long_months, by = "patient_id", relationship = "many-to-many") %>%
    mutate(boot_patient_id = paste0(patient_id, "_", bootstrap_draw))
  
  ### (1) Calculate IPTW
  iptw_denom_formula <- as.formula(paste("exp_bin_treat ~ ", paste(covariate_names, collapse = " + "), "+ strat(strat_cat_region)"))
  iptw.denom <- parglm(iptw_denom_formula, 
                       family = binomial(link = 'logit'),
                       data = d[d$month == 0, ]) # restrict to baseline person-month for model fitting
  d$iptw_denom <- predict(iptw.denom, d, type="response")
  
  # check that no NA in predictions that will mess up calculations later down the road
  if (any(is.na(d$iptw_denom))) {
    stop("NA in IPTW denominator predictions")
  }
  
  # Estimate stabilized weights
  d$w_treat_stab <- ifelse(d$exp_bin_treat==1,
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
  ipcw_num_formula <- as.formula(paste("censor == 0 ~ month + monthsqr")) # TBD if exp_bin_treat included or not
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
    group_by(boot_patient_id) %>%
    mutate(w_cens_stab = cumprod(ipcw_num)/cumprod(ipcw_denom)) %>%
    ungroup()
  
  ### (3) Multiply IPTW * IPCW
  d$w_treat_cens_stab <- d$w_treat_stab * d$w_cens_stab
  if (any(is.na(d$w_treat_cens_stab))) {
    cat("NA in final weight w_treat_cens_stab:", sum(is.na(d$w_treat_cens_stab)), "\n")
  }
  # d <- d[!is.na(d$w_treat_cens_stab), ] # Drop NA weights if there are just a few
  
  ### (4) Truncate final stabilized weight at the 1st and 99th percentile
  lower_threshold <- quantile(d$w_treat_cens_stab, 0.01)
  upper_threshold <- quantile(d$w_treat_cens_stab, 0.99)
  d$w_treat_cens_stab_trunc <- d$w_treat_cens_stab
  d$w_treat_cens_stab_trunc[d$w_treat_cens_stab_trunc < lower_threshold] <- lower_threshold
  d$w_treat_cens_stab_trunc[d$w_treat_cens_stab_trunc > upper_threshold] <- upper_threshold
  
  ### (5) Fit pooled logistic regression, with stabilized weights for IPTW and IPCW
  # Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time. -> implement!
  # Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
  # The stratification variable has been accounted for in the weights - do I still need to include it here, i.e. the have different baseline hazards by region?
  # If "I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr)" is excluded, then it would mirror the cox model
  plr_model_severecovid <- parglm(outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr), 
                                  family = binomial(link = 'logit'),
                                  data = d[d$censor==0 & d$comp_event == 0,],
                                  weights = d[d$censor==0 & d$comp_event == 0,]$w_treat_cens_stab_trunc)
  
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
  treat0.surv <- treat0 %>% group_by(boot_patient_id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
  treat1.surv <- treat1 %>% group_by(boot_patient_id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
  
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
  return(c(graph$risk0, graph$risk1, graph$rd, graph$rr))
  
  }, error = function(e) {
    # If any error occurs, count as a failed bootstrap
    cat("Bootstrap sample failed due to:", conditionMessage(e), "\n")
    boot_failures <<- boot_failures + 1
    return(rep(NA, length(c(graph$risk0, graph$risk1, graph$rd, graph$rr)))) # Return NA values to avoid breaking the boot()
  })
  
}

set.seed(423)
boot_failures <- 0 # Reset failure counter

# RUN it
te_iptw_ipcw_rd_rr_withCI_boot <- boot(data = study_ids, statistic = te_iptw_ipcw_rd_rr_withCI, R = R)
cat("Number of failed bootstrap samples:", boot_failures, "\n")

# Function to extract bootstrapped confidence intervals (only at final time point) - and keep a column with the original point estimates (without CI)
extract_ci_boot <- function(boot_obj, index) {
  ci <- boot.ci(boot_obj, conf = 0.95, type = "perc", index = index)
  if (!is.null(ci$percent)) {
    return(c(ci$percent[4], ci$percent[5])) # Lower and Upper CI for the bootstrapped values
  } else {
    return(c(NA, NA)) # for the original value
  }
}

# Identify indices for final time point. The final time point is index K, since R indexes from 1. (even though we shifted month to K - 1 above to start at month == 0)
final_indices <- c(K, 2*K, 3*K, 4*K)

# Create the results table for final month only
te_plr_iptw_ipcw_boot_tbl <- data.frame(
  Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
  Estimate_original = te_iptw_ipcw_rd_rr_withCI_boot$t0[final_indices],
  Estimate_boot = colMeans(te_iptw_ipcw_rd_rr_withCI_boot$t[, final_indices], na.rm = TRUE),
  Lower_CI = sapply(final_indices, function(i) extract_ci_boot(te_iptw_ipcw_rd_rr_withCI_boot, i)[1]),
  Upper_CI = sapply(final_indices, function(i) extract_ci_boot(te_iptw_ipcw_rd_rr_withCI_boot, i)[2])
)

# Create an empty data frame to store the structured output for the plot
risk_graph <- data.frame(
  time = 0:K,
  mean.0 = numeric(K+1),
  ll.0 = numeric(K+1),
  ul.0 = numeric(K+1),
  mean.1 = numeric(K+1),
  ll.1 = numeric(K+1),
  ul.1 = numeric(K+1)
)

# Calculate the mean and confidence intervals
mean_risk0 <- apply(te_iptw_ipcw_rd_rr_withCI_boot$t[, 1:(K+1)], 2, mean) # Mean for control group (first half contains risk0)
mean_risk1 <- apply(te_iptw_ipcw_rd_rr_withCI_boot$t[, (K+2):(2*(K+1))], 2, mean) # Mean for treatment group (second half contains risk1)

ll_risk0 <- apply(te_iptw_ipcw_rd_rr_withCI_boot$t[, 1:(K+1)], 2, function(x) quantile(x, 0.025))
ul_risk0 <- apply(te_iptw_ipcw_rd_rr_withCI_boot$t[, 1:(K+1)], 2, function(x) quantile(x, 0.975))

ll_risk1 <- apply(te_iptw_ipcw_rd_rr_withCI_boot$t[, (K+2):(2*(K+1))], 2, function(x) quantile(x, 0.025))
ul_risk1 <- apply(te_iptw_ipcw_rd_rr_withCI_boot$t[, (K+2):(2*(K+1))], 2, function(x) quantile(x, 0.975))

# Populate the `risk.boot.graph` data frame
risk_graph$mean.0 <- mean_risk0
risk_graph$ll.0 <- ll_risk0
risk_graph$ul.0 <- ul_risk0
risk_graph$mean.1 <- mean_risk1
risk_graph$ll.1 <- ll_risk1
risk_graph$ul.1 <- ul_risk1

# Create plot
plot_cum_risk_iptw_ipcw_ci <- ggplot(risk_graph,
                        aes(x=time)) +
  geom_line(aes(y = mean.1, # create line for intervention group
                color = "Metformin"), linewidth = 1.5) +
  geom_ribbon(aes(ymin = ll.1, ymax = ul.1, fill = "Metformin"), alpha = 0.4) +
  geom_line(aes(y = mean.0, # create line for control group
                color = "No Metformin"), linewidth = 1.5) +
  geom_ribbon(aes(ymin = ll.0, ymax = ul.0, fill = "No Metformin"), alpha=0.4) +
  xlab("Months") +
  # scale_x_continuous(limits = c(0, 39), # format x axis
  #                    breaks=c(0, 6, 12, 24, 36, 39)) +
  ylab("Cumulative Incidence (%)") + # label y axis
  # scale_y_continuous(limits=c(0, 0.125), # format y axis
  #                    breaks=c(0, 0.025, 0.05, 0.075, 0.1, 0.125),
  #                    labels=c("0.0%", "2.5%", "5.0%",
  #                             "7.5%", "10.0%", "12.5%")) +
  theme_minimal()+
  theme(axis.text = element_text(size=14), legend.position = c(0.2, 0.8),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())+
  scale_color_manual(values=c("#E7B800", # set colors
                              "#2E9FDF"),
                     breaks=c('No Metformin',
                              'Metformin')) +
  scale_fill_manual(values=c("#E7B800", # set colors
                             "#2E9FDF"),
                    breaks=c('No Metformin',
                             'Metformin'))


# Save output -------------------------------------------------------------
print('Save output')
# Treatment effect (risk) estimates
write.csv(te_plr_iptw_ipcw_boot_tbl, file = here::here("output", "te", "pooled_log_reg", "te_plr_iptw_ipcw_boot_severecovid.csv"))

# Marginal parametric cumulative incidence (risk) curves from plr model, with 95% CI
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_cum_risk_iptw_severecovid.png"), plot_cum_risk_iptw, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_cum_risk_iptw_ipcw_severecovid.png"), plot_cum_risk_iptw_ipcw, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_cum_risk_iptw_ipcw_ci_severecovid.png"), plot_cum_risk_iptw_ipcw_ci, width = 20, height = 20, units = "cm")

# Love/SMD plot
ggsave(filename = here::here("output", "te", "pooled_log_reg", "covplot_unweighted.png"), covplot_unweighted, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", "covplot_weighted.png"), covplot_weighted, width = 20, height = 20, units = "cm")
