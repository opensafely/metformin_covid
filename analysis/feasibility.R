################################################################################
# This script does the following:
# 1. Assess key feasibility indicators re DM diagnoses, Metformin prescription and COVID diagnoses
################################################################################

################################################################################
# 0.0 Import libraries, functions and design aspects
################################################################################
library('arrow')
library('readr')
library('here')
library('lubridate')
library('dplyr')
library('tidyr')
library('purrr')
library('forcats')
library('ggplot2')
library('gt')

## Import custom user functions
source(here::here("analysis", "functions", "fn_extract_data.R"))
source(here::here("analysis", "functions", "utility.R"))
source(here::here("analysis", "functions", "fn_diabetes_algorithm.R"))

## Redaction threshold
threshold <- 6

################################################################################
# 0.1 Create directories for output
################################################################################
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_properties"))

################################################################################
# 1 Import data
################################################################################
data_processed <- read_rds(here::here("output", "data", "data_processed.rds"))

################################################################################
# 2 Check feasibility
################################################################################
## 1. How many with a T2DM diagnosis start metformin within 5d, 7d, 10d, 30d, and 90d after (the first) pos. SARS-CoV-2 test
# T2DM variable from diabetes algo
grace_periods <- c(5, 7, 10, 30, 90) # grace periods in days (30 and 90 are probably not useful as grace period, but still of interest)
fn_t2dm_covid_metfin_start <- function(grace_period) {
  data_processed %>%
    filter(!is.na(cov_date_t2dm)) %>%
    filter(!is.na(exp_date_metfin_first)) %>%
    filter(exp_date_metfin_first > cov_date_t2dm) %>%
    filter(!is.na(cov_date_covid19_first)) %>%
    filter(exp_date_metfin_first <= (cov_date_covid19_first + days(grace_period))) %>%
    
    mutate(tb_T2DMdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_t2dm, units = "days"))) %>%
    mutate(tb_COVIDdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_covid19_first, units = "days"))) %>%
    
    summarise(
      n_start_metfin_aCOVID_midpoint6 = fn_roundmid_any(n(), threshold), # Perform redaction
      median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
      IQR_lower = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),
      IQR_upper = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE)
    ) %>%
    
    mutate("grace period in days" = grace_period) # for identification
}
n_t2dm_covid_metfin_start_midpoint6 <- map_dfr(grace_periods, fn_t2dm_covid_metfin_start)

# gt table
# n_t2dm_covid_metfin_start_midpoint6 %>%
#   gt() %>%
#   tab_header(title = "Metformin start within 10d after COVID diagnosis among T2DM") %>%
#   cols_label(
#     n_start_metfin_aCOVID_midpoint6 = "N starting metformin 10d after COVID diagnosis",
#     median_tb_T2DMdiag_metfin = "Median time to metformin start after T2DM",
#     IQR_lower = "IQR Lower",
#     IQR_upper = "IQR Upper")


## 2. How many with a T2DM diagnosis start metformin thereafter (in monthly periods til max 6m) and median time to start?
# T2DM variable from diabetes algo
periods <- c(30, 60, 90, 120, 150, 180) # in days - but now, we are interested in e.g. monthly until max. 6-monthly intervals
fn_t2dm_metfin_start <- function(periods){
  data_processed %>%
    filter(!is.na(cov_date_t2dm)) %>%
    filter(!is.na(exp_date_metfin_first)) %>%
    filter(exp_date_metfin_first > cov_date_t2dm) %>%
    filter(exp_date_metfin_first <= (cov_date_t2dm + days(periods))) %>%
    
    mutate(tb_T2DMdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_t2dm, units = "days"))) %>%
    
    summarise(
      n_start_metfin_aT2DM_midpoint6 = fn_roundmid_any(n(), threshold), # Perform redaction
      median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
      IQR_lower = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),  
      IQR_upper = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE)
      ) %>% 
    mutate("period in days" = periods)
}
n_t2dm_metfin_start_midpoint6 <- map_dfr(periods, fn_t2dm_metfin_start)

# gt table
# n_t2dm_metfin_start_midpoint6 %>%
#   gt() %>%
#   tab_header(title = "Metformin start after T2DM diagnosis") %>%
#   cols_label(
#     n_start_metfin_aT2DM_midpoint6 = "N starting metformin after T2DM diagnosis",
#     median_tb_T2DMdiag_metfin = "Median time to metformin start (days)",
#     IQR_lower = "IQR Lower",
#     IQR_upper = "IQR Upper")

################################################################################
# 3 Save output
################################################################################
write.csv(n_t2dm_covid_metfin_start_midpoint6, file = here::here("output", "data_properties", "n_t2dm_covid_metfin_start_midpoint6.csv"))
write.csv(n_t2dm_metfin_start_midpoint6, file = here::here("output", "data_properties", "n_t2dm_metfin_start_midpoint6.csv"))
