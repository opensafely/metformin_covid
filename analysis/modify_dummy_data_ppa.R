#### 
## This script modifies dummy data to: 
## - Add HbA1c, Tot Cholesterol, HDL, and BMI realistic numeric values
## - Introduce a few outliers to test filtering rules for HbA1c, Tot Chol and HDL Chol
## - Generate BMI data, not only 8 values, but also 8 sequential dates
## - Generate a few random cov_date_hypertension dates AFTER baseline to test the pipeline further down the line (currently completely empty)
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


# Generate dummy data for HbA1c, Total Cholesterol, HDL Cholesterol -------

for(i in 1:8) {
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


## BMI without outliers (8 sequential measurements) -----------------------

# Function to generate sequential dates over 2 years
generate_dates <- function(start_date) {
  # 8 measurements over ~2 years (20–120 days apart)
  as.Date(cumsum(c(0, sample(20:120, 7, replace = TRUE))) + as.numeric(start_date),
          origin = "1970-01-01")
}

# Generate 8 BMI measurement dates per patient - after landmark_date
bmi_dates <- t(sapply(df$landmark_date, generate_dates))
bmi_dates <- as.data.frame(bmi_dates)
names(bmi_dates) <- paste0("cov_date_bmi_", 1:8)
bmi_dates[] <- lapply(bmi_dates, as.Date, origin = "1970-01-01")

# Generate 8 BMI numeric values per patient (no outliers)
bmi_values <- t(sapply(1:nrow(df), function(i) {
  # BMI ~ Normal(27, 4), truncated 16–50
  pmin(pmax(rnorm(8, mean = 27, sd = 4), 16), 50)
}))
bmi_values <- as.data.frame(bmi_values)
names(bmi_values) <- paste0("cov_num_bmi_", 1:8)

# Add BMI dates and values back to df by overwriting the original values
for (i in 1:8) {
  df[[paste0("cov_date_bmi_", i)]] <- bmi_dates[[i]]
}
for (i in 1:8) {
  df[[paste0("cov_num_bmi_", i)]] <- bmi_values[[i]]
}


# Add some covariate dates between baseline date (landmark_date) and studyend date -----
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
