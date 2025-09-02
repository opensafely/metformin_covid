####
## This script does the following:
# 1. Import processed data
# 2. Estimate adjusted but marginal risks per group, risk differences, risk ratios using (i) standardization and (ii) inverse probability weighting (IPW), 
#    with the help of pooled logistic regression, and derive 95% CI via bootstrapping
# 3. Create marginal parametric cumulative incidence (risk) curves incl. 95% CI (via bootstrapping)
# 4. Save all output
####


# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(parglm) # to be computationally more efficient
library(splines)
library(lubridate)
library(boot)
library(Hmisc) # for Love/SMD plot
source(here::here("analysis", "functions", "fn_smd_before.R"))
source(here::here("analysis", "functions", "fn_smd_after.R"))


# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "te", "pooled_log_reg"))


# Import the data ---------------------------------------------------------
print('Import the data')
df_months_severecovid <- read_feather(here("output", "data", "df_months_severecovid.arrow"))
# df_months_covid_event <- read_feather(here("output", "data", "df_months_covid_event.arrow"))
# df_months_longcovid <- read_feather(here("output", "data", "df_months_longcovid.arrow"))


# Import dates ------------------------------------------------------------
print('Import dates')
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# Add splines -------------------------------------------------------------
print('Add splines')
# Compute knot locations based on percentiles, according to study protocol
age_knots <- quantile(df_months_severecovid$cov_num_age, probs = c(0.10, 0.50, 0.90))
df_months_severecovid <- df_months_severecovid %>%
  mutate(cov_num_age_spline = ns(cov_num_age, knots = age_knots))


# Define covariates ----------------------------------------------------------
print('Define covariates')
covariate_names <- names(df_months_severecovid) %>%
  grep("^(cov_num|cov_cat|cov_bin|strat_)", ., value = TRUE) %>% 
  # exclude those not needed in the model: 
  ## cov_cat_bmi covers cov_num_bmi,
  ## cov_cat_hba1c_mmol_mol covers cov_num_hba1c & cov_num_hba1c_mmol_mol
  ## cov_cat_tc_hdl_ratio covers cov_num_tc_hdl_ratio & cov_num_chol & cov_num_hdl_chol
  ## cov_num_age_spline covers cov_cat_age and cov_num_age
  setdiff(c("cov_num_bmi",
            "cov_num_hba1c", "cov_num_hba1c_mmol_mol", 
            "cov_num_tc_hdl_ratio", "cov_num_chol", "cov_num_hdl_chol",
            "cov_cat_age", "cov_num_age"
            )) 
print(covariate_names)


# Background -----------------------------------------------------
### a) Pooled logistic regression
## Pooled logistic regression models allow for the parametric estimation of risks, and thus risk differences and risk ratios.
## These models can be specified such that the effect of the treatment can vary over time, without relying on a proportional hazards assumption.
## We're modeling/simulating expected/predicted (counterfactual) risks over time under the assumption that 
## individuals could have been followed until max fup time (K-1). Hence everyone has predicted risk data until K-1.
## The aggregation at the end ensures that the true censoring distribution is respected when computing risks.

### b) Adjustment for baseline confounding
## Adjustment for baseline confounding, i.e. attempt to emulate randomization, can be accomplished using a variety of methods. 
## We will use (i) inverse probability weighting (IPW) and then (ii) risk estimation via standardisation (parametric g-formula) to get marginal risks (not conditional, but adjusted for confounders)

## IPW 
## IPW for treatment is used to create a pseudopopulation (i.e.,a hypothetical population) in which treatment is independent of the measured baseline confounders.
## Informally, the denominator of the inverse probability weight for each individual is the probability of receiving their observed treatment value, given their confounder history.
## Unstabilized and stabilized weights can be used, we will focus on stabilized.
## If we want to incorporate time-varying confounding, e.g. for censoring weights or adherence weights (per protocol analysis), 
## then we then need time-updated / time-varying weights (and info), but can simply multiply all weights
## We model the follow-up time in the outcome regression model using linear and quadratic terms and include product terms between the treatment group indicator and follow-up time.
## Other follow-up time modelling would be possible (cubic, splines), we will stick to linear and quadratic term.

## Standardization
## standardizing over the empirical distribution of the confounders to obtain marginal effect estimates (weighted average of the conditional outcomes)

## We will run through the first dataset/outcome/PLR (severe covid) step by step. For the other outcomes, we use the function at the end to apply everything at once.


# Define interval and number of bootstraps -----------------------
print('Define interval data set and number of bootstraps')
## We currently only use the monthly interval data set => max follow-up is: K = 39 months (earliest possible landmark_date [01.01.2019] to study end [01.04.2022])
## If weeks, then K = 169 weeks
K <- 39 # Total follow-up in months; if we use 39, we allow for the final 1-day fup (01.04.2022-01.04.2022); should not make much difference
df_months_severecovid$monthsqr <- df_months_severecovid$month^2 # add months square to model time in PLR
R <- 10 # Total bootstraps (ideally >500)


# IPTW: Fit treatment model and truncate the weights ----------------------
print('IPTW: Fit treatment model and truncate the weights')
## The denominator of the IPTW for each individual is the probability of receiving their observed(!) treatment value, given their confounder history.
## Here, we deal with treatment that is defined at baseline, not time-varying uptake of treatment => one constant baseline treatment weight per individual
## And, the numerator therefore is simply the proportion of patients who were actually treated.

# Calculate denominator for "probability of receiving their observed(!) treatment value"
iptw_denom_formula <- as.formula(paste("I(exp_bin_treat == 1) ~ ", paste(covariate_names, collapse = " + ")))
iptw.denom <- parglm(iptw_denom_formula, 
                     family = binomial(link = 'logit'),
                     data = df_months_severecovid[df_months_severecovid$month == 0, ]) 
# summary(iptw.denom)
# summary(iptw.denom)$coefficients
coef_summary <- summary(iptw.denom)$coefficients
# Check if any estimates are NA
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some iptw.denom coefficient estimates are NA!")
} else {
  message("All iptw.denom coefficient estimates are valid (no NAs).")
}

# Get predictions (but only use baseline rows, stable weights over time!)
df_months_severecovid <- df_months_severecovid %>%
  group_by(patient_id) %>%
  mutate(
    iptw_denom = predict(
      iptw.denom,
      newdata = cur_data()[month == 0, ],
      type = "response"
    ) %>% first()  # take the first (baseline) prediction and assign to the remaining months
  ) %>%
  ungroup()

# Estimate stabilized weights
p_treat <- mean(df_months_severecovid$exp_bin_treat[df_months_severecovid$month == 0]) # prop treated, based on baseline month only, should not matter in this case, but safer
df_months_severecovid <- df_months_severecovid %>%
  mutate(w_treat_stab = ifelse(exp_bin_treat == 1,
                               p_treat / iptw_denom,
                               (1 - p_treat) / (1 - iptw_denom)))
# truncate the weights
lower_threshold <- quantile(df_months_severecovid$w_treat_stab, 0.01, na.rm = TRUE)
upper_threshold <- quantile(df_months_severecovid$w_treat_stab, 0.99, na.rm = TRUE)
df_months_severecovid$w_treat_stab_trunc <- df_months_severecovid$w_treat_stab
df_months_severecovid$w_treat_stab_trunc[df_months_severecovid$w_treat_stab_trunc < lower_threshold] <- lower_threshold
df_months_severecovid$w_treat_stab_trunc[df_months_severecovid$w_treat_stab_trunc > upper_threshold] <- upper_threshold

## Min, 25th percentile, median, mean, SD, 75th percentile, and max:
summary(df_months_severecovid$w_treat_stab)
sd(df_months_severecovid$w_treat_stab)
summary(df_months_severecovid$w_treat_stab_trunc)
sd(df_months_severecovid$w_treat_stab_trunc)


# IPTW: Check the weights and the balance ---------------------------------------
print('IPTW: Check the weights and the balance')

# Restrict to baseline rows
df_baseline <- df_months_severecovid %>% filter(month == 0)

# List of variables to compare (include more, not only covariate_names, e.g. also the numeric variables for HbA1c and TC/HDL ratio)
varlist <- names(df_months_severecovid) %>%
  grep("^(cov_num|cov_cat|cov_bin|strat_)", ., value = TRUE) %>% 
  # but exclude those that make no sense: 
  ## cov_num_age_spline and cov_cat_age, represented by cov_num_age
  setdiff(c("cov_num_age_spline", "cov_cat_age")) 
print(varlist)

# Create subsets of data, according to treat_b status, at month 0
treat_b0 <- df_baseline %>% filter(exp_bin_treat==0)
treat_b1 <- df_baseline %>% filter(exp_bin_treat==1)
w_treat_0 <- treat_b0$w_treat_stab_trunc
w_treat_1 <- treat_b1$w_treat_stab_trunc

### Check standardized mean difference BEFORE weighting
smd_before <- lapply(varlist, fn_smd_before)
smd_before_df <- as.data.frame(do.call(rbind, smd_before), stringsAsFactors = FALSE)
smd_before_df <- transform(smd_before_df,
                           var = as.character(var),
                           type = as.character(type),
                           smd = as.numeric(smd))
covplot_unweighted_severecovid <- ggplot(smd_before_df, aes(x = smd, y = reorder(var, smd))) +
  geom_point(color = "steelblue", size = 3) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "red") +
  labs(title = "Love Plot: Standardized Mean Differences before weighting",
       x = "Standardized Mean Difference (SMD)",
       y = "Covariate") +
  theme_minimal()

### Check standardized mean difference AFTER weighting
smd_after <- lapply(varlist, 
                    fn_smd_after, 
                    w_treat_0 = w_treat_0, 
                    w_treat_1 = w_treat_1)
smd_after_df <- as.data.frame(do.call(rbind, smd_after), stringsAsFactors = FALSE)
smd_after_df <- transform(smd_after_df,
                          var = as.character(var),
                          type = as.character(type),
                          smd = as.numeric(smd))
covplot_weighted_severecovid <- ggplot(smd_after_df, aes(x = smd, y = reorder(var, smd))) +
  geom_point(color = "darkorange", size = 3) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "red") +
  labs(title = "Love Plot: Standardized Mean Differences after weighting",
       x = "Standardized Mean Difference (SMD)",
       y = "Covariate") +
  theme_minimal()


# IPTW: Fit pooled logistic regression -------------------------------------
print('IPTW: Fit pooled logistic regression and create plot_cum_risk_iptw')
# Fit pooled logistic regression, with stabilized weights (currently only IPW for treatment, i.e. baseline confounding => same weights across time)
# Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time
# Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
plr_model_severecovid <- parglm(outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr), 
                                family = binomial(link = 'logit'),
                                data = df_months_severecovid[df_months_severecovid$censor_LTFU==0 & df_months_severecovid$comp_event == 0,],
                                weights = df_months_severecovid[df_months_severecovid$censor_LTFU==0 & df_months_severecovid$comp_event == 0,]$w_treat_stab_trunc)
# summary(plr_model_severecovid)
coef_summary <- summary(plr_model_severecovid)$coefficients
# Check if any estimates are NA
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some plr_model_severecovid coefficient estimates are NA!")
} else {
  message("All plr_model_severecovid coefficient estimates are valid (no NAs).")
}

### STANDARDIZATION: Transform estimates to risks at each time point in each group ###

# Create dataset with all time points for each individual under each treatment level
df_pred <- df_months_severecovid %>% filter(month == 0) %>% dplyr::select(-month) %>% crossing(month = 0:(K-1))
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
plot_cum_risk_iptw_severecovid <- ggplot(graph,
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


# IPTW & IPCW for PPA: Add censoring event weights for the per-protocol analysis ----------------------------------
print('IPTW & IPCW for PPA: Add censoring event weights for the per-protocol analysis')
## We want to use IPW for starting treatment in the control group (PPA): cens_bin_metfin_start_cont | “What would the effect have been if everyone had adhered to their assigned strategy?”
## An uncensored individual’s IPCW is the inverse of their probability of remaining uncensored given their treatment and covariate history. 
## In the context of the PPA, each individual who continued adhering to their treatment strategy receives a weight that is proportional to the inverse of the probability of adhering to the treatment strategy, given their specific history.
## Once an individual is censored, they receive a censoring weight of 0 from that point onward. 

## Fit a pooled logistic regression model for the denominator of IPCW
## This model should predict the probability of remaining uncensored (not deviating from their treatment strategy) at each timepoint. 
## Continuous time, here months (modeled using linear and quadratic terms), will serve as the time scale for this model. 
## We do not include any product terms in the model

# Denominator model: P(uncensored | covariates + time); do not include treatment history since our PPA only applies to the control group!
ipcw_denom_formula <- as.formula(paste("I(cens_bin_metfin_start_cont == 0) ~ month + monthsqr +", paste(covariate_names, collapse = " + ")))
ipcw.denom <- parglm(ipcw_denom_formula, 
                     family = binomial(link = 'logit'),
                     data = df_months_severecovid) 
# summary(ipcw.denom)
# summary(ipcw.denom)$coefficients
coef_summary <- summary(ipcw.denom)$coefficients
# Check if any estimates are NA
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some ipcw.denom coefficient estimates are NA!")
} else {
  message("All ipcw.denom coefficient estimates are valid (no NAs).")
}

# Obtain predicted probabilities of being uncensored
df_months_severecovid$ipcw_denom <- predict(ipcw.denom, df_months_severecovid, type="response")

## For the numerator we now fit a model to include the time-varying aspect of the weights - but without any covariates
# Numerator model: P(uncensored | time only)
ipcw_num_formula <- as.formula(paste("I(cens_bin_metfin_start_cont == 0) ~ month + monthsqr"))
ipcw.num <- parglm(ipcw_num_formula, 
                   family = binomial(link = 'logit'),
                   data = df_months_severecovid) 

# Obtain predicted probabilities of being uncensored
df_months_severecovid$ipcw_num <- predict(ipcw.num, df_months_severecovid, type="response")

### Estimate stabilized inverse probability weights for censoring ###
# Take cumulative products starting at baseline: cumulative product of num/denom
df_months_severecovid <- df_months_severecovid %>%
  group_by(patient_id) %>%
  mutate(
    # w_cens_stab = cumprod(ipcw_num)/cumprod(ipcw_denom),
    w_cens_stab = cumprod(ipcw_num/ipcw_denom),
    # zero out weights after censoring (cens_bin_metfin_start_cont == 1 means censored)
    w_cens_stab = ifelse(cummax(cens_bin_metfin_start_cont) == 1, 0, w_cens_stab)
  ) %>%
  ungroup()

# Take product of weight for censoring and treatment for each individual (take the untruncated w_treat_stab from above, and only then truncate both combined)
df_months_severecovid$w_treat_cens_stab <- df_months_severecovid$w_treat_stab * df_months_severecovid$w_cens_stab

### Truncate final stabilized weight at the 1st and 99th percentile ###
lower_threshold <- quantile(df_months_severecovid$w_treat_cens_stab, 0.01, na.rm = TRUE)
upper_threshold <- quantile(df_months_severecovid$w_treat_cens_stab, 0.99, na.rm = TRUE)
df_months_severecovid$w_treat_cens_stab_trunc <- df_months_severecovid$w_treat_cens_stab
df_months_severecovid$w_treat_cens_stab_trunc[df_months_severecovid$w_treat_cens_stab_trunc < lower_threshold] <- lower_threshold
df_months_severecovid$w_treat_cens_stab_trunc[df_months_severecovid$w_treat_cens_stab_trunc > upper_threshold] <- upper_threshold

###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: truncated weights ###
summary(df_months_severecovid$w_treat_cens_stab_trunc)
sd(df_months_severecovid$w_treat_cens_stab_trunc)

## quite variable weights, investigate a bit more, rare censoring events in some months?
censoring_rates <- df_months_severecovid %>%
  group_by(month) %>%
  summarise(censor_rate = mean(cens_bin_metfin_start_cont, na.rm = TRUE), .groups="drop") %>%
  mutate(censor_rate_pct = round(censor_rate * 100, 2))
plot_censoring_ppa <- ggplot(censoring_rates, aes(x = month, y = censor_rate_pct)) +
  geom_line() +
  geom_point() +
  ylab("Censoring for starting metformin per month (%)") +
  xlab("Month") +
  theme_minimal()
censoring_weights <- df_months_severecovid %>%
  group_by(month) %>%
  summarise(
    mean_w = mean(w_cens_stab, na.rm = TRUE),
    sd_w = sd(w_cens_stab, na.rm = TRUE),
    max_w = max(w_cens_stab, na.rm = TRUE)
  )
print(censoring_weights, n = Inf)

# df_months_severecovid %>%
#   select(patient_id, month, start_date_month, end_date_month, exp_bin_treat, outcome, cens_bin_metfin_start_cont, censor_LTFU, comp_event, w_treat_cens_stab_trunc, w_treat_stab, w_cens_stab, ipcw_num, ipcw_denom) %>%
#   filter(month == 38) %>%
#   # filter(w_treat_cens_stab_trunc > 3000) %>%
#   View()


# IPTW & IPCW for PPA: Fit pooled logistic regression -------------------------
print('IPTW & IPCW for PPA: Fit pooled logistic regression and create plot_cum_risk_iptw_ipcw_severecovid')
# Fit pooled logistic regression, with stabilized weights for IPTW and IPCW
# Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time
# Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
plr_model_severecovid <- parglm(outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr), 
                                family = binomial(link = 'logit'),
                                data = df_months_severecovid[df_months_severecovid$censor_LTFU==0 & df_months_severecovid$comp_event == 0,],
                                weights = df_months_severecovid[df_months_severecovid$censor_LTFU==0 & df_months_severecovid$comp_event == 0,]$w_treat_cens_stab_trunc)
# summary(plr_model_severecovid)
coef_summary <- summary(plr_model_severecovid)$coefficients
# Check if any estimates are NA
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some plr_model_severecovid with PPA censoring coefficient estimates are NA!")
} else {
  message("All plr_model_severecovid with PPA censoring coefficient estimates are valid (no NAs).")
}

### Transform estimates to risks at each time point in each group ###

# Create dataset with all time points for each individual under each treatment level
df_pred <- df_months_severecovid %>% filter(month == 0) %>% dplyr::select(-month) %>% crossing(month = 0:(K-1))
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
plot_cum_risk_iptw_ipcw_severecovid <- ggplot(graph,
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


# IPTW & IPCW for PPA & IPCW for LTFU: Add censoring event weights for LTFU ----------
print('IPTW & IPCW for PPA & IPCW for LTFU')

## We want to use IPW for two censoring events: a) LTFU, i.e. "the effect had no one been lost to follow-up", and b) starting treatment in the control group (PPA)
## Informally, an uncensored individual’s inverse probability of censoring weight is the inverse of their probability of remaining uncensored given their treatment and covariate history. 
## In the context of loss to follow-up, each individual who was not lost to follow-up receives a weight that is proportional to the inverse of the probability of not being lost to follow-up, given their specific history.
## We will consider stabilized IP weights for censoring, in which the numerator of the weights is the probability of remaining uncensored given an individual’s history
## Once an individual is censored, they receive a censoring weight of 0 from that point onward. 

## In this simplified example, we assume censoring due to loss to follow-up depends only on the baseline treatment and baseline covariates. 
## In a second step, we will incorporate time-updated covariate info, too.

## Fit a pooled logistic regression model for the denominator of the IP weights for censoring 
## This model should predict the probability of remaining uncensored (not being lost to follow-up or not having started treatment) at each timepoint. 
## Continuous time, here months (modeled using linear and quadratic terms), will serve as the time scale for this model. 
## We do not include any product terms in the model.
## Include exp_bin_treat here (different to PPA!)

ipcw_ltfu_denom_formula <- as.formula(paste("I(censor_LTFU == 0) ~ exp_bin_treat + month + monthsqr +", paste(covariate_names, collapse = " + ")))
ipcw_ltfu.denom <- parglm(ipcw_ltfu_denom_formula, 
                  family = binomial(link = 'logit'),
                  data = df_months_severecovid) 
coef_summary <- summary(ipcw_ltfu.denom)$coefficients
# Check if any estimates are NA
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some ipcw_ltfu.denom coefficient estimates are NA!")
} else {
  message("All ipcw_ltfu.denom coefficient estimates are valid (no NAs).")
}

# Obtain predicted probabilities of being uncensored for denominator
df_months_severecovid$ipcw_ltfu_denom <- predict(ipcw_ltfu.denom, df_months_severecovid, type="response")

## For the numerator we now fit a model to include the time-varying aspect of the weights - but without any covariates
ipcw_ltfu_num_formula <- as.formula(paste("I(censor_LTFU == 0) ~ month + monthsqr"))
ipcw_ltfu.num <- parglm(ipcw_ltfu_num_formula, 
                   family = binomial(link = 'logit'),
                   data = df_months_severecovid) 
coef_summary <- summary(ipcw_ltfu.num)$coefficients
# Check if any estimates are NA
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some ipcw_ltfu.num coefficient estimates are NA!")
} else {
  message("All ipcw_ltfu.num coefficient estimates are valid (no NAs).")
}

# Obtain predicted probabilities of being uncensored for numerator
df_months_severecovid$ipcw_ltfu_num <- predict(ipcw_ltfu.num, df_months_severecovid, type="response")

### Estimate stabilized inverse probability weights for censoring ###
# Take cumulative products starting at baseline
df_months_severecovid <- df_months_severecovid %>%
  group_by(patient_id) %>%
  mutate(
    w_cens_ltfu_stab = cumprod(ipcw_ltfu_num/ipcw_ltfu_denom),
    # zero out weights after censoring (cens_bin_metfin_start_cont == 1 means censored)
    w_cens_ltfu_stab = ifelse(cummax(censor_LTFU) == 1, 0, w_cens_ltfu_stab)
  ) %>% 
  ungroup()

# Take product of weight for censoring and treatment for each individual (take the untruncated w_treat_stab from above)
df_months_severecovid$w_treat_cens_ltfu_stab <- df_months_severecovid$w_treat_stab * df_months_severecovid$w_cens_stab * df_months_severecovid$w_cens_ltfu_stab

### Truncate final stabilized weight at the 1st and 99th percentile ###
lower_threshold <- quantile(df_months_severecovid$w_treat_cens_ltfu_stab, 0.01, na.rm = TRUE)
upper_threshold <- quantile(df_months_severecovid$w_treat_cens_ltfu_stab, 0.99, na.rm = TRUE)
df_months_severecovid$w_treat_cens_ltfu_stab_trunc <- df_months_severecovid$w_treat_cens_ltfu_stab
df_months_severecovid$w_treat_cens_ltfu_stab_trunc[df_months_severecovid$w_treat_cens_ltfu_stab_trunc < lower_threshold] <- lower_threshold
df_months_severecovid$w_treat_cens_ltfu_stab_trunc[df_months_severecovid$w_treat_cens_ltfu_stab_trunc > upper_threshold] <- upper_threshold

###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: truncated weights ###
summary(df_months_severecovid$w_treat_cens_ltfu_stab_trunc)
sd(df_months_severecovid$w_treat_cens_ltfu_stab_trunc)

# df_months_severecovid %>%
#   select(patient_id, outcome, exp_bin_treat, censor_LTFU, cens_bin_metfin_start_cont, month, w_treat_cens_ltfu_stab_trunc, w_cens_ltfu_stab, ipcw_ltfu_num, ipcw_ltfu_denom, w_treat_cens_stab_trunc, w_treat_stab, w_cens_stab, ipcw_num, ipcw_denom) %>%
#   # filter(cens_bin_metfin_start_cont == 1) %>%
#   # filter(w_treat_cens_stab_trunc > 3000) %>%
#   View()

censoring_rates <- df_months_severecovid %>%
  group_by(month) %>%
  summarise(censor_rate = mean(censor_LTFU, na.rm = TRUE), .groups="drop") %>%
  mutate(censor_rate_pct = round(censor_rate * 100, 2))
plot_censoring_ltfu <- ggplot(censoring_rates, aes(x = month, y = censor_rate_pct)) +
  geom_line() +
  geom_point() +
  ylab("Censoring for LTFU per month (%)") +
  xlab("Month") +
  theme_minimal()
censoring_weights <- df_months_severecovid %>%
  group_by(month) %>%
  summarise(
    mean_w = mean(w_cens_ltfu_stab, na.rm = TRUE),
    sd_w = sd(w_cens_ltfu_stab, na.rm = TRUE),
    max_w = max(w_cens_ltfu_stab, na.rm = TRUE)
  )
print(censoring_weights, n = Inf)


# IPTW & IPCW for PPA & IPCW for LTFU: Fit pooled logistic regression -------------------------
print('IPTW & IPCW for PPA & IPCW for LTFU: Fit pooled logistic regression and create plot_cum_risk_iptw_ipcw_ltfu_severecovid')
# Fit pooled logistic regression, with stabilized weights for IPTW and IPCW (both for PPA and LTFU)
# Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time
# Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
plr_model_severecovid <- parglm(outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr), 
                                family = binomial(link = 'logit'),
                                data = df_months_severecovid[df_months_severecovid$censor_LTFU==0 & df_months_severecovid$comp_event == 0,],
                                weights = df_months_severecovid[df_months_severecovid$censor_LTFU==0 & df_months_severecovid$comp_event == 0,]$w_treat_cens_ltfu_stab_trunc)
# summary(plr_model_severecovid)
coef_summary <- summary(plr_model_severecovid)$coefficients
# Check if any estimates are NA
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some plr_model_severecovid incl. PPA and LTFU censoring coefficient estimates are NA!")
} else {
  message("All plr_model_severecovid incl. PPA and LTFU censoring coefficient estimates are valid (no NAs).")
}

### Transform estimates to risks at each time point in each group ###

# Create dataset with all time points for each individual under each treatment level
df_pred <- df_months_severecovid %>% filter(month == 0) %>% dplyr::select(-month) %>% crossing(month = 0:(K-1))
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
plot_cum_risk_iptw_ipcw_ltfu_severecovid <- ggplot(graph,
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


# IPTW & IPCW for PPA & IPCW for LTFU: Do all steps together incl. 95% CI using bootstrapping -------------------------
print('IPTW & IPCW for PPA & IPCW for LTFU: Do all steps together incl. 95% CI using bootstrapping and create plot_cum_risk_iptw_ipcw_ltfu_severecovid')

te_iptw_ipcw_rd_rr_withCI <- function(data, indices, data_full, covariate_names) {

  # capture and count unsuccessful bootstraps
  tryCatch({

  # Resample patient IDs (with replacement)
  boot.ids <- data[indices, , drop = FALSE] %>%
    mutate(bootstrap_draw = row_number())  # Assign a unique row for each draw

  # Join back the full time records **for each draw** (i.e. same patient is allowed to be sampled several times, i.e. with replacement)
  d <- boot.ids %>%
    left_join(data_full, by = "patient_id", relationship = "many-to-many") %>%
    mutate(boot_patient_id = paste0(patient_id, "_", bootstrap_draw))

  ### (1) Calculate IPTW
  iptw_denom_formula <- as.formula(paste("I(exp_bin_treat == 1) ~ ", paste(covariate_names, collapse = " + ")))
  iptw.denom <- parglm(iptw_denom_formula,
                       family = binomial(link = 'logit'),
                       data = d[d$month == 0, ]) # restrict to baseline person-month for model fitting
  coef_summary <- summary(iptw.denom)$coefficients
  if (anyNA(coef_summary[, "Estimate"])) {
    warning("Some iptw.denom coefficient estimates are NA!")
  }
  
  # baseline predictions, one row per patient, based on baseline month covariates only!
  d <- d %>%
    group_by(boot_patient_id) %>%
    mutate(
      iptw_denom = predict(
        iptw.denom,
        newdata = cur_data()[month == 0, ],
        type = "response"
      ) %>% first()  # take the first (baseline) prediction and assign to the remaining months
    ) %>%
    ungroup()
  
  # d$iptw_denom <- predict(iptw.denom, d, type="response") # works but wrong if time-varying covariates

  # check that no NA in predictions that will mess up calculations later down the road
  if (any(is.na(d$iptw_denom) | is.infinite(d$iptw_denom))) {
    stop("iptw_denom has invalid stabilized weights (NA or Inf)")
  }

  # Estimate stabilized weights
  p_treat <- mean(d$exp_bin_treat[d$month == 0]) # prop treated, based on baseline month only, should not matter in this case, but safer
  d <- d %>%
    mutate(w_treat_stab = ifelse(exp_bin_treat == 1,
                                 p_treat / iptw_denom,
                                 (1 - p_treat) / (1 - iptw_denom)))

  if (any(is.na(d$w_treat_stab) | is.infinite(d$w_treat_stab))) {
    stop("w_treat_stab has invalid stabilized weights (NA or Inf)")
  }
  

  ### (2) Calculate IPCW for PPA
  ipcw_denom_formula <- as.formula(paste("I(cens_bin_metfin_start_cont == 0) ~ month + monthsqr +", paste(covariate_names, collapse = " + ")))
  ipcw.denom <- parglm(ipcw_denom_formula,
                       family = binomial(link = 'logit'),
                       data = d)
  coef_summary <- summary(ipcw.denom)$coefficients
  if (anyNA(coef_summary[, "Estimate"])) {
    warning("Some ipcw.denom coefficient estimates are NA!")
  }

  # Obtain predicted probabilities of being uncensored for denominator
  d$ipcw_denom <- predict(ipcw.denom, d, type="response")
  if (any(is.na(d$ipcw_denom))) {
    stop("NA in IPCW for PPA denominator predictions")
  }
  if(any(d$ipcw_denom <= 0 | d$ipcw_denom >= 1, na.rm=TRUE)) {
    cat("Problematic ipcw_denom values:", sum(d$ipcw_denom <= 0 | d$ipcw_denom >= 1, na.rm=TRUE), "\n")
  }

  ## For the numerator we now fit a model to include the time-varying aspect of the weights
  ipcw_num_formula <- as.formula(paste("I(cens_bin_metfin_start_cont == 0) ~ month + monthsqr"))
  ipcw.num <- parglm(ipcw_num_formula,
                  family = binomial(link = 'logit'),
                  data = d)

  # Obtain predicted probabilities of being uncensored for numerator
  d$ipcw_num <- predict(ipcw.num, d, type="response")

  if (any(is.na(d$ipcw_num))) {
    stop("NA in IPCW for PPA numerator predictions")
  }

  # Take the cumulative product
  d <- d %>%
    group_by(boot_patient_id) %>%
    mutate(
      w_cens_stab = cumprod(ipcw_num/ipcw_denom),
      w_cens_stab = ifelse(cummax(cens_bin_metfin_start_cont) == 1, 0, w_cens_stab)
    ) %>%
    ungroup()
  
  
  ### (3) Calculate IPCW for LTFU
  ipcw_ltfu_denom_formula <- as.formula(paste("I(censor_LTFU == 0) ~ exp_bin_treat + month + monthsqr +", paste(covariate_names, collapse = " + ")))
  ipcw_ltfu.denom <- parglm(ipcw_ltfu_denom_formula,
                       family = binomial(link = 'logit'),
                       data = d)
  coef_summary <- summary(ipcw_ltfu.denom)$coefficients
  if (anyNA(coef_summary[, "Estimate"])) {
    warning("Some ipcw_ltfu.denom coefficient estimates are NA!")
  }
  
  # Obtain predicted probabilities of being uncensored for denominator
  d$ipcw_ltfu_denom <- predict(ipcw_ltfu.denom, d, type="response")
  if (any(is.na(d$ipcw_ltfu_denom))) {
    stop("NA in IPCW for LTFU denominator predictions")
  }
  if(any(d$ipcw_ltfu_denom <= 0 | d$ipcw_ltfu_denom >= 1, na.rm=TRUE)) {
    cat("Problematic ipcw_ltfu_denom values:", sum(d$ipcw_ltfu_denom <= 0 | d$ipcw_ltfu_denom >= 1, na.rm=TRUE), "\n")
  }
  
  ## For the numerator we now fit a model to include the time-varying aspect of the weights
  ipcw_ltfu_num_formula <- as.formula(paste("I(censor_LTFU == 0) ~ month + monthsqr"))
  ipcw_ltfu.num <- parglm(ipcw_ltfu_num_formula,
                     family = binomial(link = 'logit'),
                     data = d)
  coef_summary <- summary(ipcw_ltfu.denom)$coefficients
  if (anyNA(coef_summary[, "Estimate"])) {
    warning("Some ipcw_ltfu.num coefficient estimates are NA!")
  }
  
  # Obtain predicted probabilities of being uncensored for numerator
  d$ipcw_ltfu_num <- predict(ipcw_ltfu.num, d, type="response")
  if (any(is.na(d$ipcw_ltfu_num))) {
    stop("NA in IPCW for PPA numerator predictions")
  }
  
  # Take the cumulative product
  d <- d %>%
    group_by(boot_patient_id) %>%
    mutate(
      w_cens_ltfu_stab = cumprod(ipcw_ltfu_num/ipcw_ltfu_denom),
      w_cens_ltfu_stab = ifelse(cummax(censor_LTFU) == 1, 0, w_cens_ltfu_stab)
    ) %>% 
    ungroup()
  

  ### (4) Multiply IPTW * IPCW for PPA * IPCW for LTFU
  d$w_treat_cens_stab <- d$w_treat_stab * d$w_cens_ltfu_stab * d$w_cens_stab 
  
  if (any(is.na(d$w_treat_cens_stab))) {
    cat("NA in final weight w_treat_cens_stab:", sum(is.na(d$w_treat_cens_stab)), "\n")
  }

  
  ### (5) Truncate final stabilized weight at the 1st and 99th percentile
  lower_threshold <- quantile(d$w_treat_cens_stab, 0.01)
  upper_threshold <- quantile(d$w_treat_cens_stab, 0.99)
  d$w_treat_cens_stab_trunc <- d$w_treat_cens_stab
  d$w_treat_cens_stab_trunc[d$w_treat_cens_stab_trunc < lower_threshold] <- lower_threshold
  d$w_treat_cens_stab_trunc[d$w_treat_cens_stab_trunc > upper_threshold] <- upper_threshold

  ### (6) Fit pooled logistic regression, with stabilized weights for IPTW and IPCW
  # Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time
  # Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
  # If "I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr)" is excluded, then it would mirror the cox model
  plr_model_severecovid <- parglm(outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr),
                                  family = binomial(link = 'logit'),
                                  data = d[d$censor_LTFU==0 & d$comp_event == 0,],
                                  weights = d[d$censor_LTFU==0 & d$comp_event == 0,]$w_treat_cens_stab_trunc)
  coef_summary <- summary(plr_model_severecovid)$coefficients
  if (anyNA(coef_summary[, "Estimate"])) {
    warning("Some plr_model_severecovid coefficient estimates are NA!")
  }

  ### (7) Transform estimates to risks at each time point in each group
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

### RUN it, for each outcome/dataset
set.seed(443)
boot_failures <- 0 # Reset failure counter

## Severe COVID / primary outcome
# extract unique identifiers, as a data frame
study_ids <- data.frame(patient_id = unique(df_months_severecovid$patient_id))
# add dataset, covariate list (see above), and boostrap runs (see above)
te_iptw_ipcw_rd_rr_withCI_severecovid_boot <- boot(
  data = study_ids,
  statistic = function(data, indices) {
    te_iptw_ipcw_rd_rr_withCI(data, indices, data_full = df_months_severecovid, covariate_names = covariate_names)
  },
  R = R
)
cat("Number of failed bootstrap samples:", boot_failures, "\n")

### risk tables
# Function to extract bootstrapped confidence intervals for a time point - and keep a column with the original point estimates (without CI)
extract_ci_boot <- function(boot_obj, index) {
  ci <- boot.ci(boot_obj, conf = 0.95, type = "perc", index = index)
  if (!is.null(ci$percent)) {
    return(c(ci$percent[4], ci$percent[5])) # Lower and Upper CI for the bootstrapped values
  } else {
    return(c(NA, NA)) # for the original value
  }
}

# Identify target months (indices) for when to extract risk estimates.
# The final time point is index K, since R indexes from 1 (even though we shifted month to K - 1 above to start at month == 0)
k_39_index <- K   # final time point (K = 39 if months 0 to 38)
k_24_index <- 25  # corresponds to month 24, because R is 1-indexed => corresponds to Cox 730 days timepoint
indices_39 <- c(k_39_index, 2*k_39_index, 3*k_39_index, 4*k_39_index)
indices_24 <- c(k_24_index, 2*k_24_index, 3*k_24_index, 4*k_24_index)

te_plr_iptw_ipcw_severecovid_39_boot_tbl <- data.frame(
  Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
  Estimate_original = te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t0[indices_39],
  Estimate_boot = colMeans(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, indices_39], na.rm = TRUE),
  Lower_CI = sapply(indices_39, function(i) extract_ci_boot(te_iptw_ipcw_rd_rr_withCI_severecovid_boot, i)[1]),
  Upper_CI = sapply(indices_39, function(i) extract_ci_boot(te_iptw_ipcw_rd_rr_withCI_severecovid_boot, i)[2])
)
te_plr_iptw_ipcw_severecovid_24_boot_tbl <- data.frame(
  Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
  Estimate_original = te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t0[indices_24],
  Estimate_boot = colMeans(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, indices_24], na.rm = TRUE),
  Lower_CI = sapply(indices_24, function(i) extract_ci_boot(te_iptw_ipcw_rd_rr_withCI_severecovid_boot, i)[1]),
  Upper_CI = sapply(indices_24, function(i) extract_ci_boot(te_iptw_ipcw_rd_rr_withCI_severecovid_boot, i)[2])
)


### risk plots
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
mean_risk0 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, 1:(K+1)], 2, mean) # Mean for control group (first half contains risk0)
mean_risk1 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, (K+2):(2*(K+1))], 2, mean) # Mean for treatment group (second half contains risk1)

ll_risk0 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, 1:(K+1)], 2, function(x) quantile(x, 0.025))
ul_risk0 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, 1:(K+1)], 2, function(x) quantile(x, 0.975))

ll_risk1 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, (K+2):(2*(K+1))], 2, function(x) quantile(x, 0.025))
ul_risk1 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, (K+2):(2*(K+1))], 2, function(x) quantile(x, 0.975))

# Populate the `risk.boot.graph` data frame
risk_graph$mean.0 <- mean_risk0
risk_graph$ll.0 <- ll_risk0
risk_graph$ul.0 <- ul_risk0
risk_graph$mean.1 <- mean_risk1
risk_graph$ll.1 <- ll_risk1
risk_graph$ul.1 <- ul_risk1

# Create plot
plot_cum_risk_iptw_ipcw_ltfu_ci_severecovid <- ggplot(risk_graph,
                        aes(x=time)) +
  geom_line(aes(y = mean.1, # create line for intervention group
                color = "Metformin"), linewidth = 1.5) +
  geom_ribbon(aes(ymin = ll.1, ymax = ul.1, fill = "Metformin"), alpha = 0.4) +
  geom_line(aes(y = mean.0, # create line for control group
                color = "No Metformin"), linewidth = 1.5) +
  geom_ribbon(aes(ymin = ll.0, ymax = ul.0, fill = "No Metformin"), alpha=0.4) +
  xlab("Months") +
  ylab("Cumulative Incidence (%)") +
  theme_minimal()+
  theme(axis.text = element_text(size=14), legend.position = c(0.2, 0.8),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())+
  scale_color_manual(values=c("#E7B800",
                              "#2E9FDF"),
                     breaks=c('No Metformin',
                              'Metformin')) +
  scale_fill_manual(values=c("#E7B800",
                             "#2E9FDF"),
                    breaks=c('No Metformin',
                             'Metformin'))


# Save output -------------------------------------------------------------
print('Save output')
## primary outcome: Severe COVID
# Love/SMD plot
ggsave(filename = here::here("output", "te", "pooled_log_reg", "covplot_unweighted_severecovid.png"), covplot_unweighted_severecovid, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", "covplot_weighted_severecovid.png"), covplot_weighted_severecovid, width = 20, height = 20, units = "cm")
# Censoring plots
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_censoring_ppa.png"), plot_censoring_ppa, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_censoring_ltfu.png"), plot_censoring_ltfu, width = 20, height = 20, units = "cm")
# Treatment effect (risk) estimates at max follow-up (39 months)
write.csv(te_plr_iptw_ipcw_severecovid_39_boot_tbl, file = here::here("output", "te", "pooled_log_reg", "te_plr_iptw_ipcw_severecovid_39_boot_tbl.csv"))
# Treatment effect (risk) estimates at 24 months (corresponding to Cox model)
write.csv(te_plr_iptw_ipcw_severecovid_24_boot_tbl, file = here::here("output", "te", "pooled_log_reg", "te_plr_iptw_ipcw_severecovid_24_boot_tbl.csv"))
# Marginal parametric cumulative incidence (risk) curves from plr model, last one including 95% CI
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_cum_risk_iptw_severecovid.png"), plot_cum_risk_iptw_severecovid, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_cum_risk_iptw_ipcw_severecovid.png"), plot_cum_risk_iptw_ipcw_severecovid, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_cum_risk_iptw_ipcw_ltfu_severecovid.png"), plot_cum_risk_iptw_ipcw_ltfu_severecovid, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_cum_risk_iptw_ipcw_ltfu_ci_severecovid.png"), plot_cum_risk_iptw_ipcw_ltfu_ci_severecovid, width = 20, height = 20, units = "cm")
