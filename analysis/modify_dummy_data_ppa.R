#### 
## This script modifies dummy data to: 
## - Add realistic numeric values to our dynamic time-updated covariates: HbA1c, Tot Cholesterol, HDL, and BMI (currently, nearly all empty in dummy data)
## - Introduce a few outliers to test filtering rules for HbA1c, Tot Chol and HDL Chol
## - Generate sequential values with sequential dates for these dynamic time-updated covariates (e.g. 12 measurement dates over the course of 2 years)
## - Add some covariate dates between baseline date (landmark_date) and studyend date for fixed time-updated covariates
####


# Import libraries and functions ------------------------------------------
library(dplyr)
source(here::here("analysis", "functions", "fn_dd_cens_dates.R"))


# Import dates ------------------------------------------------------------
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))
studyend_date <- as.Date(study_dates$studyend_date, format = "%Y-%m-%d")


# Set seed for reproducibility --------------------------------------------
set.seed(36)


# Function to generate lab values with some outliers -----------------------
generate_lab_values <- function(n, mean, sd, min_val, max_val, outlier_prob = 0.05, outlier_range = c(-50, 200)) {
  vals <- rnorm(n, mean = mean, sd = sd)
  
  # Introduce random outliers
  is_outlier <- runif(n) < outlier_prob
  vals[is_outlier] <- runif(sum(is_outlier), min = outlier_range[1], max = outlier_range[2])
  
  # Introduce some NAs randomly
  vals[sample(1:n, size = round(0.1*n))] <- NA
  
  return(vals)
}


# Function to generate sequential dates over 2 years -------
# 12 measurements over ~2 years (20–120 days apart)
generate_dates <- function(start_date, n = 12) {
  as.Date(
    cumsum(c(0, sample(20:120, n - 1, replace = TRUE))) + as.numeric(start_date),
    origin = "1970-01-01"
  )
}

# Function to generate measurement dates and values for multiple variables (and overwrite old one) ------
generate_measurements <- function(df, vars, n = 12) {
  for (v in vars) {
    var_name <- v$name
    
    # ---- Generate dates ----
    var_dates <- t(sapply(df$landmark_date, v$date_fun, n = n))
    var_dates <- as.data.frame(var_dates)
    names(var_dates) <- paste0("cov_date_", var_name, "_", 1:n)
    var_dates[] <- lapply(var_dates, as.Date, origin = "1970-01-01")  # ensure Date
    
    # ---- Generate values ----
    var_values <- t(sapply(1:nrow(df), function(i) v$value_fun(n)))
    var_values <- as.data.frame(var_values)
    names(var_values) <- paste0("cov_num_", var_name, "_", 1:n)
    
    # ---- Overwrite existing columns or create new ----
    for (i in 1:n) {
      df[[paste0("cov_date_", var_name, "_", i)]] <- var_dates[[i]]
      df[[paste0("cov_num_", var_name, "_", i)]]  <- var_values[[i]]
    }
  }
  
  return(df)
}


vars <- list(
  list(
    name = "bmi",
    date_fun = function(landmark_date, n) generate_dates(landmark_date, n),
    value_fun = function(n) pmin(pmax(rnorm(n, mean = 27, sd = 4), 16), 50),
    type = "numeric"
  ),
  list(
    name = "hba1c",
    date_fun = function(landmark_date, n) generate_dates(landmark_date, n),
    value_fun = function(n) pmin(pmax(rnorm(n, mean = 50, sd = 15), 0), 120),
    type = "numeric"
  ),
  list(
    name = "chol",
    date_fun = function(landmark_date, n) generate_dates(landmark_date, n),
    value_fun = function(n) pmin(pmax(rnorm(n, mean = 5, sd = 1.2), 1.75), 20),
    type = "numeric"
  ),
  list(
    name = "hdl_chol",
    date_fun = function(landmark_date, n) generate_dates(landmark_date, n),
    value_fun = function(n) pmin(pmax(rnorm(n, mean = 1.3, sd = 0.3), 0.4), 5),
    type = "numeric"
  )
)

# list(
#   name = "smoking",
#   date_fun = function(landmark_date, n) generate_dates(landmark_date, n),
#   value_fun = function(n) sample(c("never", "former", "current"), n, replace = TRUE),
#   type = "categorical"
# )

# Apply to df
df <- generate_measurements(df, vars, n = 12)


# Now, modify dummy data for HbA1c, Total Cholesterol, HDL Cholesterol to create outliers -------
for(i in 1:12) {
  df[[paste0("cov_num_hba1c_", i)]] <- generate_lab_values(
    n = nrow(df), mean = 50, sd = 15, min_val = 0, max_val = 120,
    outlier_prob = 0.1, outlier_range = c(-20, 200)
  )
  
  df[[paste0("cov_num_chol_", i)]] <- generate_lab_values(
    n = nrow(df), mean = 5, sd = 1.2, min_val = 1.75, max_val = 20,
    outlier_prob = 0.1, outlier_range = c(-5, 30)
  )
  
  df[[paste0("cov_num_hdl_chol_", i)]] <- generate_lab_values(
    n = nrow(df), mean = 1.3, sd = 0.3, min_val = 0.4, max_val = 5,
    outlier_prob = 0.1, outlier_range = c(-2, 10)
  )
}


# Add some covariate dates between baseline date (landmark_date) and studyend date for fixed time-updated covariates -----
## cov_date_hypertension and cov_date_ami are completely missing in the current dummy data, so replace all
## but impute only 5%
df$cov_date_hypertension <- sapply(df$landmark_date, function(baseline_date) {
  fn_dd_cens_dates(baseline_date, studyend_date)
})
df$cov_date_hypertension <- as.Date(df$cov_date_hypertension, origin = "1970-01-01")

df$cov_date_ami <- sapply(df$landmark_date, function(baseline_date) {
  fn_dd_cens_dates(baseline_date, studyend_date)
})
df$cov_date_ami <- as.Date(df$cov_date_ami, origin = "1970-01-01")

df$cov_date_all_stroke <- sapply(df$landmark_date, function(baseline_date) {
  fn_dd_cens_dates(baseline_date, studyend_date)
})
df$cov_date_all_stroke <- as.Date(df$cov_date_all_stroke, origin = "1970-01-01")

df$cov_date_dementia <- sapply(df$landmark_date, function(baseline_date) {
  fn_dd_cens_dates(baseline_date, studyend_date)
})
df$cov_date_dementia <- as.Date(df$cov_date_dementia, origin = "1970-01-01")