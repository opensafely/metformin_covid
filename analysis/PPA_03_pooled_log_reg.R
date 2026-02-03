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
library(rlang) # for fn_standardize_risks.R
source(here::here("analysis", "functions", "fn_truncate_by_percentile.R"))
source(here::here("analysis", "functions", "fn_standardize_risks.R"))
source(here::here("analysis", "functions", "fn_extract_ci_boot.R"))
source(here::here("analysis", "covariates.R"))


# Parse command line arguments --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Default values
include_interaction <- TRUE
R_arg <- NULL  # store the command-line R value if provided

# Parse arguments
if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "no_interaction") {
      include_interaction <- FALSE
    } else if (grepl("^R=", args[i])) {
      R_arg <- as.integer(sub("^R=", "", args[i]))
    }
  }
}

# Set R based on environment and arguments
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  # Local/test environment: always use 5 bootstraps only !!
  R <- 5
  if (!is.null(R_arg)) {
    message("Note: Running on local/test data - ignoring R argument and using R=5 by default")
  }
} else {
  # Real data environment: use argument if provided, otherwise default to 100
  R <- ifelse(!is.null(R_arg), R_arg, 100)
}

# Set output suffix based on interaction
if (include_interaction) {
  print('Running analysis WITH interaction terms')
  output_suffix <- ""
} else {
  print('Running analysis WITHOUT interaction terms')
  output_suffix <- "_no_interaction"
}

message("Number of bootstraps: ", R)


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


# Define covariates ---------------------------------------------------------
print('Define covariates')
print(covariates_tv_names)


# Define the primary time point of risk estimation, and time^2 --------------
print('Define the primary time point of risk estimation, and time^2')
# We currently use monthly intervals => 2-year risk estimation time point is K = 24 months
K <- 24
df_months_severecovid$monthsqr <- df_months_severecovid$month^2 # add months square to model time in PLR


# Define the treatment history variable --------------------------------------
print('Define the treatment history variable')
## treat is my time-varying treatment variable, see PPA_02_add_events
## exp_bin_treat is my group indicator (at baseline, constant throughout all months)
## let's build treat_lag to capture treatment history
## And by default at month 0 treat_lag shall be 0, since everyone is by design untreated before baseline -> "initiating metformin"
df_months_severecovid <- df_months_severecovid %>%
  group_by(patient_id)  %>%
  arrange(month)  %>%
  mutate(treat_lag = lag(treat, default = 0)) %>% 
  ungroup()


# Adjust for baseline confounding and time-varying protocol deviation using one time-varying IPTW model -------

## Build the weights -------
print('Adjust for baseline confounding and time-varying protocol deviation using one time-varying IPTW model: Build the weights')
## I(treat == 1) variable is the time-varying treatment indicator.
## Denominator model: probability of treatment given past treatment and time-varying covariates
treat_denom_formula <- as.formula(paste("I(treat == 1) ~ treat_lag + month + monthsqr +", paste(c(covariates_tv_names), collapse = " + ")))
treat.denom <- parglm(treat_denom_formula,
                       family = binomial(link = 'logit'),
                       data = df_months_severecovid)
coef_summary <- summary(treat.denom)$coefficients
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some coefficient estimates from treat.denom are NA!")
}
df_months_severecovid$treat_denom <- predict(treat.denom, df_months_severecovid, type = "response")

## Numerator model: probability of treatment given past treatment history only (and time modelling element)
treat_num_formula <- as.formula(paste("I(treat == 1) ~ treat_lag + month + monthsqr"))
treat.num <- parglm(treat_num_formula,
                     family = binomial(link = 'logit'),
                     data = df_months_severecovid)
coef_summary <- summary(treat.num)$coefficients
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some coefficient estimates from treat.num are NA!")
}
df_months_severecovid$treat_num <- predict(treat.num, df_months_severecovid, type = "response")

## Calculate the cumulative time-varying stabilized weight
# The treat_denom model predicts the probability of a person being in the treatment group (I(treat == 1)) at any given time, conditioned on their covariates.
# => the weight calculation cumprod(treat_num / treat_denom) is only correct for the treatment group. For the control group, the probability of their observed treatment status (treat == 0) is 1 - treat_denom.
df_months_severecovid <- df_months_severecovid %>%
  group_by(patient_id) %>%
  arrange(month, .by_group = TRUE) %>%
  mutate(
    w_treat_stab = ifelse(
      exp_bin_treat == 1, cumprod(treat_num / treat_denom),
      cumprod((1 - treat_num)/(1 - treat_denom))
    ),
    # Zero weights after censoring in control
    w_treat_stab = ifelse(exp_bin_treat == 0 & cummax(treat) == 1, 0, w_treat_stab)
  ) %>%
  ungroup()

## Truncate final stabilized weight at the 1st and 99th percentile
df_months_severecovid <- fn_truncate_by_percentile(df = df_months_severecovid, col_name = "w_treat_stab")

## Min, 25th percentile, median, mean, SD, 75th percentile, and max
summary(df_months_severecovid$w_treat_stab)
summary(df_months_severecovid$w_treat_stab_trunc)


## Fit the outcome model -------
print('Adjust for baseline confounding and time-varying protocol deviation using one time-varying IPTW model: Fit the outcome model')
outcome_formula <- if (include_interaction) {
  outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr)
} else {
  outcome ~ exp_bin_treat + month + monthsqr
}
treat.severecovid <- parglm(outcome_formula,
                            family = binomial(link = 'logit'),
                            data = df_months_severecovid[df_months_severecovid$censor_LTFU==0 & df_months_severecovid$comp_event == 0,],
                            weights = df_months_severecovid[df_months_severecovid$censor_LTFU==0 & df_months_severecovid$comp_event == 0,]$w_treat_stab_trunc)
coef_summary <- summary(treat.severecovid)$coefficients
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some coefficient estimates from the treat.severecovid with w_treat_stab_trunc are NA!")
}


## Standardization -------
print('Adjust for baseline confounding and time-varying protocol deviation using one time-varying IPTW model: Standardization')
# Transform estimates to risks at each time point in each group
cum_risk_treat_severecovid <- fn_standardize_risks(
  df = df_months_severecovid,
  K = K, # CAVE: data frame is modified to reflect that risks are estimated at the END of each interval, i.e., start at month == 0, initial zero row added => Final estimate in K-1
  model = treat.severecovid,
  group_col = "exp_bin_treat",
  time_col = "month",
  patient_id_col = "patient_id"
)


## Construct marginal parametric cumulative incidence (risk) curves (without CIs) -------
print('Adjust for baseline confounding and time-varying protocol deviation using one time-varying IPTW model: Plot')
plot_cum_risk_treat_severecovid <- ggplot(cum_risk_treat_severecovid$graph_data,
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


# Adjust additionally for time-varying LTFU ----------

## Build the weights -------
print('Adjust additionally for time-varying LTFU: Build the weights')
## We want to use IPW also for LTFU, i.e. "the effect had no one been lost to follow-up"
## Informally, an uncensored individuals' inverse probability of censoring weight is the inverse of their probability of remaining uncensored given their treatment and covariate history. 
## Once an individual is censored, they receive a censoring weight of 0 from that point onward. 
ltfu_denom_formula <- as.formula(paste("I(censor_LTFU == 0) ~ treat_lag + month + monthsqr +", paste(covariates_tv_names, collapse = " + ")))
ltfu.denom <- parglm(ltfu_denom_formula, 
                          family = binomial(link = 'logit'),
                          data = df_months_severecovid) 
coef_summary <- summary(ltfu.denom)$coefficients
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some ltfu.denom coefficient estimates are NA!")
}
df_months_severecovid$ltfu_denom <- predict(ltfu.denom, df_months_severecovid, type="response")

## For the numerator we now fit a model to include the time-varying aspect of the weights - but without any covariates
ltfu_num_formula <- as.formula(paste("I(censor_LTFU == 0) ~ treat_lag + month + monthsqr"))
ltfu.num <- parglm(ltfu_num_formula, 
                        family = binomial(link = 'logit'),
                        data = df_months_severecovid) 
coef_summary <- summary(ltfu.num)$coefficients
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some ltfu.num coefficient estimates are NA!")
}
df_months_severecovid$ltfu_num <- predict(ltfu.num, df_months_severecovid, type="response")

## Estimate stabilized inverse probability weights for censoring
## Take cumulative products starting at baseline
df_months_severecovid <- df_months_severecovid %>%
  group_by(patient_id) %>%
  arrange(month, .by_group = TRUE) %>%
  mutate(
    w_ltfu_stab = cumprod(ltfu_num/ltfu_denom),
    w_ltfu_stab = ifelse(cummax(censor_LTFU) == 1, 0, w_ltfu_stab)
  ) %>% 
  ungroup()
summary(df_months_severecovid$w_ltfu_stab)

## Investigate the censoring-for-LTFU a bit more
censoring_ltfu_rates <- df_months_severecovid %>%
  group_by(month) %>%
  summarise(censor_ltfu_rate = mean(censor_LTFU, na.rm = TRUE), .groups="drop") %>%
  mutate(
    censor_ltfu_rate_pct = round(censor_ltfu_rate * 100, 2),
    cum_censor_ltfu_rate = cumsum(censor_ltfu_rate),
    cum_censor_ltfu_rate_pct = round(cum_censor_ltfu_rate * 100, 2)
  )
plot_censoring_ltfu <- ggplot(censoring_ltfu_rates, aes(x = month)) +
  geom_line(aes(y = censor_ltfu_rate_pct), color = "blue") +
  geom_point(aes(y = censor_ltfu_rate_pct), color = "blue") +
  geom_line(aes(y = cum_censor_ltfu_rate_pct), color = "red") +
  geom_point(aes(y = cum_censor_ltfu_rate_pct), color = "red") +
  ylab("Censoring for LTFU (%)") +
  xlab("Month") +
  theme_minimal() +
  ggtitle("Monthly (blue) and cumulative (red) LTFU censoring rate")
censoring_ltfu_weights <- df_months_severecovid %>%
  group_by(month) %>%
  summarise(
    mean_w = mean(w_ltfu_stab, na.rm = TRUE),
    sd_w = sd(w_ltfu_stab, na.rm = TRUE),
    max_w = max(w_ltfu_stab, na.rm = TRUE)
  )


# Adjust additionally for time-varying competing risk (non-covid deaths) ----------

## Build the weights -------
print('Adjust additionally for time-varying competing risk (non-covid deaths)')
## We want to use IPW also for the competing risk non-covid deaths, i.e. "the effect had no one died to other causes than covid-19"
## Informally, an uncensored individuals' inverse probability of censoring weight is the inverse of their probability of remaining uncensored given their treatment and covariate history. 
## Once an individual is censored, they receive a censoring weight of 0 from that point onward. 
comp_denom_formula <- as.formula(paste("I(comp_event == 0) ~ treat_lag + month + monthsqr +", paste(covariates_tv_names, collapse = " + ")))
comp.denom <- parglm(comp_denom_formula, 
                     family = binomial(link = 'logit'),
                     data = df_months_severecovid) 
coef_summary <- summary(comp.denom)$coefficients
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some comp.denom coefficient estimates are NA!")
}
df_months_severecovid$comp_denom <- predict(comp.denom, df_months_severecovid, type="response")

## For the numerator we now fit a model to include the time-varying aspect of the weights - but without any covariates
comp_num_formula <- as.formula(paste("I(comp_event == 0) ~ treat_lag + month + monthsqr"))
comp.num <- parglm(comp_num_formula, 
                   family = binomial(link = 'logit'),
                   data = df_months_severecovid) 
coef_summary <- summary(comp.num)$coefficients
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some comp.num coefficient estimates are NA!")
}
df_months_severecovid$comp_num <- predict(comp.num, df_months_severecovid, type="response")

## Estimate stabilized inverse probability weights for censoring
## Take cumulative products starting at baseline
df_months_severecovid <- df_months_severecovid %>%
  group_by(patient_id) %>%
  arrange(month, .by_group = TRUE) %>%
  mutate(
    w_comp_stab = cumprod(comp_num/comp_denom),
    w_comp_stab = ifelse(cummax(comp_event) == 1, 0, w_comp_stab)
  ) %>% 
  ungroup()
summary(df_months_severecovid$w_comp_stab)

## Investigate the censoring-for-competing risk a bit more
censoring_comp_rates <- df_months_severecovid %>%
  group_by(month) %>%
  summarise(censor_comp_rate = mean(comp_event, na.rm = TRUE), .groups="drop") %>%
  mutate(
    censor_comp_rate_pct = round(censor_comp_rate * 100, 2),
    cum_censor_comp_rate = cumsum(censor_comp_rate),
    cum_censor_comp_rate_pct = round(cum_censor_comp_rate * 100, 2)
  )
plot_censoring_comp <- ggplot(censoring_comp_rates, aes(x = month)) +
  geom_line(aes(y = censor_comp_rate_pct), color = "blue") +
  geom_point(aes(y = censor_comp_rate_pct), color = "blue") +
  geom_line(aes(y = cum_censor_comp_rate_pct), color = "red") +
  geom_point(aes(y = cum_censor_comp_rate_pct), color = "red") +
  ylab("Censoring for competing risk (%)") +
  xlab("Month") +
  theme_minimal() +
  ggtitle("Monthly (blue) and cumulative (red) competing risk censoring rate")
censoring_comp_weights <- df_months_severecovid %>%
  group_by(month) %>%
  summarise(
    mean_w = mean(w_comp_stab, na.rm = TRUE),
    sd_w = sd(w_comp_stab, na.rm = TRUE),
    max_w = max(w_comp_stab, na.rm = TRUE)
  )



# Take the product of all 3 (untruncated) weights -----------
print('Take the product of all (untruncated) weights')
df_months_severecovid$w_treat_ltfu_comp_stab <- df_months_severecovid$w_treat_stab * df_months_severecovid$w_ltfu_stab * df_months_severecovid$w_comp_stab

## Truncate the final stabilized weight at the 1st and 99th percentile ------
print('Truncate the final stabilized weight at the 1st and 99th percentile')
df_months_severecovid <- fn_truncate_by_percentile(df = df_months_severecovid, col_name = "w_treat_ltfu_comp_stab")

## Min, 25th percentile, median, mean, SD, 75th percentile, and max
summary(df_months_severecovid$w_treat_ltfu_comp_stab)
summary(df_months_severecovid$w_treat_ltfu_comp_stab_trunc)



# Fit the outcome model -----------
print('Adjust additionally for time-varying LTFU: Fit the outcome model')
# Fit pooled logistic regression, with stabilized weights for IPTW and IPCW (both for PPA and LTFU)
# Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
outcome_formula <- if (include_interaction) {
  outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr)
} else {
  outcome ~ exp_bin_treat + month + monthsqr
}
severecovid.treat.ltfu.comp <- parglm(outcome_formula, 
                                family = binomial(link = 'logit'),
                                data = df_months_severecovid[df_months_severecovid$censor_LTFU==0 & df_months_severecovid$comp_event == 0,],
                                weights = df_months_severecovid[df_months_severecovid$censor_LTFU==0 & df_months_severecovid$comp_event == 0,]$w_treat_ltfu_comp_stab_trunc)
coef_summary <- summary(severecovid.treat.ltfu.comp)$coefficients
if (anyNA(coef_summary[, "Estimate"])) {
  warning("Some severecovid.treat.ltfu.comp (incl. PPA and LTFU censoring and Comp censoring) coefficient estimates are NA!")
}

## Standardization -------
print('Adjust additionally for time-varying LTFU: Standardization')
cum_risk_treat_ltfu_comp_severecovid <- fn_standardize_risks(
  df = df_months_severecovid,
  K = K, # CAVE: data frame is modified to reflect that risks are estimated at the END of each interval, initial zero row added => Final estimate in K-1
  model = severecovid.treat.ltfu.comp,
  group_col = "exp_bin_treat",
  time_col = "month",
  patient_id_col = "patient_id"
)

## Construct marginal parametric cumulative incidence (risk) curves (without CIs) -------
print('Adjust additionally for time-varying LTFU and competing risk: Plot')
plot_cum_risk_treat_ltfu_severecovid <- ggplot(cum_risk_treat_ltfu_comp_severecovid$graph_data,
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



# All steps together incl. 95% CI using bootstrapping -------------------------
print('All steps together incl. 95% CI using bootstrapping')

te_iptw_ipcw_rd_rr_withCI <- function(data, indices, data_full, covariates_tv_names, K, include_interaction = TRUE) {
  
  # capture and count unsuccessful bootstraps
  tryCatch({
    
    # Re-sample patient IDs (with replacement)
    boot.ids <- data[indices, , drop = FALSE] %>%
      mutate(bootstrap_draw = row_number())  # Assign a unique row for each draw
    
    # Join back the full time records **for each draw** (i.e. same patient is allowed to be sampled several times, i.e. with replacement)
    d <- boot.ids %>%
      left_join(data_full, by = "patient_id", relationship = "many-to-many") %>%
      mutate(boot_patient_id = paste0(patient_id, "_", bootstrap_draw))
    
    
    ### (1) Calculate IPTW for baseline confounding and time-varying protocol deviation
    # Denominator model
    treat_denom_formula <- as.formula(paste("I(treat == 1) ~ treat_lag + month + monthsqr +", paste(c(covariates_tv_names), collapse = " + ")))
    treat.denom <- parglm(treat_denom_formula,
                          family = binomial(link = 'logit'),
                          data = d)
    coef_summary <- summary(treat.denom)$coefficients
    if (anyNA(coef_summary[, "Estimate"])) {
      warning("Some treat.denom coefficient estimates are NA!")
    }
    d$treat_denom <- predict(treat.denom, d, type = "response")
    
    # Numerator model
    treat_num_formula <- as.formula(paste("I(treat == 1) ~ treat_lag + month + monthsqr"))
    treat.num <- parglm(treat_num_formula,
                        family = binomial(link = 'logit'),
                        data = d)
    coef_summary <- summary(treat.num)$coefficients
    if (anyNA(coef_summary[, "Estimate"])) {
      warning("Some coefficient estimates from treat.num are NA!")
    }
    d$treat_num <- predict(treat.num, d, type = "response")
    
    # Cumulative product and ensure correct inverse and censoring behaviour
    d <- d %>%
      group_by(boot_patient_id) %>%
      arrange(month, .by_group = TRUE) %>%
      mutate(
        w_treat_stab = ifelse(
          exp_bin_treat == 1, cumprod(treat_num / treat_denom),
          cumprod((1 - treat_num)/(1 - treat_denom))
        ),
        w_treat_stab = ifelse(exp_bin_treat == 0 & cummax(treat) == 1, 0, w_treat_stab)
      ) %>%
      ungroup()
    
    # check that no NA in predictions that will mess up calculations later down the road
    if (any(is.na(d$w_treat_stab) | is.infinite(d$w_treat_stab))) {
      stop("w_treat_stab has invalid stabilized weights (NA or Inf)")
    }
    
    
    ### (2) Calculate IPCW for LTFU
    # Denominator model
    ltfu_denom_formula <- as.formula(paste("I(censor_LTFU == 0) ~ treat_lag + month + monthsqr +", paste(covariates_tv_names, collapse = " + ")))
    ltfu.denom <- parglm(ltfu_denom_formula,
                              family = binomial(link = 'logit'),
                              data = d)
    coef_summary <- summary(ltfu.denom)$coefficients
    if (anyNA(coef_summary[, "Estimate"])) {
      warning("Some ltfu.denom coefficient estimates are NA!")
    }
    d$ltfu_denom <- predict(ltfu.denom, d, type="response")
    
    # Numerator model
    ltfu_num_formula <- as.formula(paste("I(censor_LTFU == 0) ~ treat_lag + month + monthsqr"))
    ltfu.num <- parglm(ltfu_num_formula,
                            family = binomial(link = 'logit'),
                            data = d)
    coef_summary <- summary(ltfu.num)$coefficients
    if (anyNA(coef_summary[, "Estimate"])) {
      warning("Some ltfu.num coefficient estimates are NA!")
    }
    d$ltfu_num <- predict(ltfu.num, d, type="response")
    
    # Cumulative product and ensure correct inverse and censoring behaviour
    d <- d %>%
      group_by(boot_patient_id) %>%
      arrange(month, .by_group = TRUE) %>%
      mutate(
        w_ltfu_stab = cumprod(ltfu_num/ltfu_denom),
        w_ltfu_stab = ifelse(cummax(censor_LTFU) == 1, 0, w_ltfu_stab)
      ) %>% 
      ungroup()
    
    
    ### (3) Calculate IPCW for the competing risk (non-covid death)
    # Denominator model
    comp_denom_formula <- as.formula(paste("I(comp_event == 0) ~ treat_lag + month + monthsqr +", paste(covariates_tv_names, collapse = " + ")))
    comp.denom <- parglm(comp_denom_formula,
                         family = binomial(link = 'logit'),
                         data = d)
    coef_summary <- summary(comp.denom)$coefficients
    if (anyNA(coef_summary[, "Estimate"])) {
      warning("Some comp.denom coefficient estimates are NA!")
    }
    d$comp_denom <- predict(comp.denom, d, type="response")
    
    # Numerator model
    comp_num_formula <- as.formula(paste("I(comp_event == 0) ~ treat_lag + month + monthsqr"))
    comp.num <- parglm(comp_num_formula,
                       family = binomial(link = 'logit'),
                       data = d)
    coef_summary <- summary(comp.num)$coefficients
    if (anyNA(coef_summary[, "Estimate"])) {
      warning("Some comp.num coefficient estimates are NA!")
    }
    d$comp_num <- predict(comp.num, d, type="response")
    
    # Cumulative product and ensure correct inverse and censoring behaviour
    d <- d %>%
      group_by(boot_patient_id) %>%
      arrange(month, .by_group = TRUE) %>%
      mutate(
        w_comp_stab = cumprod(comp_num/comp_denom),
        w_comp_stab = ifelse(cummax(comp_event) == 1, 0, w_comp_stab)
      ) %>% 
      ungroup()
    
    
    ### (4) Multiply the weights
    d$w_treat_ltfu_comp_stab <- d$w_treat_stab * d$w_ltfu_stab * d$w_comp_stab
    if (any(is.na(d$w_treat_ltfu_comp_stab))) {
      cat("NA in final weight w_treat_ltfu_comp_stab:", sum(is.na(d$w_treat_ltfu_comp_stab)), "\n")
    }
    
    
    ### (5) Truncate final stabilized weight at the 1st and 99th percentile
    d <- fn_truncate_by_percentile(
      df = d,
      col_name = "w_treat_ltfu_comp_stab"
    )
    
    
    ### (6) Fit the outcome model
    # Include the treatment group indicator, the follow-up time (linear and quadratic terms), and product terms between the treatment group indicator and follow-up time
    # Train the model only on individuals who were still at risk of the outcome, i.e. only include individuals who are uncensored and alive
    # If "I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr)" is excluded, then it mirrors the cox model
    outcome_formula <- if (include_interaction) {
      outcome ~ exp_bin_treat + month + monthsqr + I(exp_bin_treat*month) + I(exp_bin_treat*monthsqr)
    } else {
      outcome ~ exp_bin_treat + month + monthsqr
    }
    
    severecovid.treat.ltfu.comp <- parglm(outcome_formula,
                                          family = binomial(link = 'logit'),
                                          data = d[d$censor_LTFU==0 & d$comp_event == 0,],
                                          weights = d[d$censor_LTFU==0 & d$comp_event == 0,]$w_treat_ltfu_comp_stab_trunc)
    coef_summary <- summary(severecovid.treat.ltfu.comp)$coefficients
    if (anyNA(coef_summary[, "Estimate"])) {
      warning("Some severecovid.treat.ltfu.comp coefficient estimates are NA!")
    }
    
    
    ### (7) Standardization/Marginalization
    results_std <- fn_standardize_risks(
      df = d,
      K = K, # CAVE: data frame is modified to reflect that risks are estimated at the END of each interval, i.e., start at month == 0, initial zero row added => Final estimate in K-1
      model = severecovid.treat.ltfu.comp,
      group_col = "exp_bin_treat",
      time_col = "month",
      patient_id_col = "boot_patient_id"
    )
    
    return(c(results_std$graph_data$risk0, results_std$graph_data$risk1, results_std$graph_data$rd, results_std$graph_data$rr))

  }, error = function(e) {
    # If any error occurs, count as a failed bootstrap
    cat("Bootstrap sample failed due to:", conditionMessage(e), "\n")
    boot_failures <<- boot_failures + 1
    # Check if results_std$graph_data is available before attempting to get its length
    if(exists("results_std") && !is.null(results_std$graph_data)) {
      return(rep(NA, length(c(results_std$graph_data$risk0, results_std$graph_data$risk1, results_std$graph_data$rd, results_std$graph_data$rr))))
    } else {
      # If not, return NA values of a default length
      return(rep(NA, 4 * K)) # 4 metrics (risk0, risk1, rd, rr) for K time points
    }
  })
}

### RUN it, for each outcome/dataset (we run the per protocol estimate only on the primary outcome, but would also work for the other outcomes)
set.seed(443)
boot_failures <- 0 # Reset failure counter
study_ids <- data.frame(patient_id = unique(df_months_severecovid$patient_id))

te_iptw_ipcw_rd_rr_withCI_severecovid_boot <- boot(
  data = study_ids,
  statistic = function(data, indices) {
    te_iptw_ipcw_rd_rr_withCI(data, indices, data_full = df_months_severecovid, 
                              covariates_tv_names = covariates_tv_names, K = K,
                              include_interaction = include_interaction)
  },
  R = R
)
cat("Number of failed bootstrap samples:", boot_failures, "\n")

### My bootstrap object contains:
cat(sprintf("   Number of statistics per sample (risk0, risk1, RR, RD): %d\n", 
            ncol(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t)))
# 100 values per bootstrap sample
# Indices: risk0 (1-25), risk1 (26-50), rd (51-75), rr (76-100)
# the final month (23) is at position 24 in the vectors
# Row 1: month = 0 (zero row, time_0 = 0)
# Row 2: month = 0 (time_0 = 1)
# Row 3: month = 1 (time_0 = 2)

### Extract the final time point only
# Indices
indices <- c(
  K, # risk0 at month 23 (index 24)
  K + (K + 1), # risk1 at month 23 (index 49)
  K + 2 * (K + 1), # rd at month 23 (index 74)
  K + 3 * (K + 1) # rr at month 23 (index 99)
)
te_plr_iptw_ipcw_severecovid_boot_tbl <- data.frame(
  Measure = c("Risk Control", "Risk Treatment", "Risk Difference", "Risk Ratio"),
  Estimate_original = te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t0[indices],
  Estimate_boot = colMeans(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, indices], na.rm = TRUE),
  Lower_CI = sapply(indices, function(i) fn_extract_ci_boot(te_iptw_ipcw_rd_rr_withCI_severecovid_boot, i)[1]),
  Upper_CI = sapply(indices, function(i) fn_extract_ci_boot(te_iptw_ipcw_rd_rr_withCI_severecovid_boot, i)[2])
)

### Extract all time points to build the risk table for the cum risk graph
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
mean_risk0 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, 1:(K+1)], 2, mean) # Mean for control group (first half contains risk0) 1:25 correctly extracts all 25 values (months 0-23) for control
mean_risk1 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, (K+2):(2*(K+1))], 2, mean) # Mean for treatment group (second half contains risk1) 26:50 correctly extracts all 25 values (months 0-23) for intervention

ll_risk0 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, 1:(K+1)], 2, function(x) quantile(x, 0.025))
ul_risk0 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, 1:(K+1)], 2, function(x) quantile(x, 0.975))

ll_risk1 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, (K+2):(2*(K+1))], 2, function(x) quantile(x, 0.025))
ul_risk1 <- apply(te_iptw_ipcw_rd_rr_withCI_severecovid_boot$t[, (K+2):(2*(K+1))], 2, function(x) quantile(x, 0.975))

risk_graph$mean.0 <- mean_risk0
risk_graph$ll.0 <- ll_risk0
risk_graph$ul.0 <- ul_risk0
risk_graph$mean.1 <- mean_risk1
risk_graph$ll.1 <- ll_risk1
risk_graph$ul.1 <- ul_risk1

# Create plot
plot_cum_risk_treat_ltfu_comp_ci_severecovid <- ggplot(risk_graph,
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
  scale_color_manual(
    name = "",
    values = c("No Metformin" = "#E7B800",
               "Metformin" = "#2E9FDF")
  ) +
  scale_fill_manual(
    name = "",
    values = c("No Metformin" = "#E7B800",
               "Metformin" = "#2E9FDF")
  )


# Save output -------------------------------------------------------------
print('Save output')
# Treatment effect (risk) estimates at 24 months, once with interaction and once without (to mirror Cox), depending on terminal args
write.csv(te_plr_iptw_ipcw_severecovid_boot_tbl, 
          file = here::here("output", "te", "pooled_log_reg", 
                            paste0("te_plr_iptw_ipcw_severecovid_boot_tbl", output_suffix, ".csv")))
# Marginal parametric cumulative incidence (risk) curves from plr models, first two only with point estimates, last one including 95% CI (from bootstrap)
ggsave(filename = here::here("output", "te", "pooled_log_reg", 
                             paste0("plot_cum_risk_treat_severecovid", output_suffix, ".png")), 
       plot_cum_risk_treat_severecovid, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", 
                             paste0("plot_cum_risk_treat_ltfu_severecovid", output_suffix, ".png")), 
       plot_cum_risk_treat_ltfu_severecovid, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "te", "pooled_log_reg", 
                             paste0("plot_cum_risk_treat_ltfu_comp_ci_severecovid", output_suffix, ".png")), 
       plot_cum_risk_treat_ltfu_comp_ci_severecovid, width = 20, height = 20, units = "cm")
# Censoring plots and underlying data
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_censoring_ltfu.png"), plot_censoring_ltfu, width = 20, height = 20, units = "cm")
write.csv(censoring_ltfu_rates, file = here::here("output", "te", "pooled_log_reg", "censoring_ltfu_rates.csv"))
write.csv(censoring_ltfu_weights, file = here::here("output", "te", "pooled_log_reg", "censoring_ltfu_weights.csv"))
ggsave(filename = here::here("output", "te", "pooled_log_reg", "plot_censoring_comp.png"), plot_censoring_comp, width = 20, height = 20, units = "cm")
write.csv(censoring_comp_rates, file = here::here("output", "te", "pooled_log_reg", "censoring_comp_rates.csv"))
write.csv(censoring_comp_weights, file = here::here("output", "te", "pooled_log_reg", "censoring_comp_weights.csv"))
