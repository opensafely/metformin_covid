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

# landmark date
landmark_date <- as.Date("2020-02-01")

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
## How many with a T2DM diagnosis (in the eligible period, i.e. 07/2018 - 07/2019) have started metformin after T2DM diagnosis but before/on landmark
## How many with a T2DM diagnosis (in the eligible period, i.e. 07/2018 - 07/2019) have NOT started metformin after T2DM diagnosis but before/on landmark
# Use T2DM variable from diabetes algo
# Incl. outcomes
mid2018_date <- as.Date("2018-07-01")
years_in_days <- c(0, 366, 731, 1096, 1461, 1827) # Calculated from 1 year = 365.25 days, taking into account leap years. Until mid2013.

# redacted
fn_t2dm_midpoint6 <- function(years_in_days) {
  data_processed %>%
    # T2DM diagnosis 6m prior landmark, but after the looped mid-year date
    filter(cov_date_t2dm < landmark_date - days(183) &
             cov_date_t2dm > mid2018_date - days(years_in_days)) %>% 
    # did not start metformin or if then only after T2DM diagnosis
    filter(exp_date_metfin_first > cov_date_t2dm | is.na(exp_date_metfin_first)) %>% 
    # date of death is missing, i.e., alive, or only dead after landmark (start of fup)
    filter(is.na(qa_date_of_death) | qa_date_of_death > landmark_date) %>% 
    
    # among the above, define treatment status
    mutate(exp_bin_treat = case_when(exp_date_metfin_first <= landmark_date ~ 1, # 1 if started
                                     TRUE ~ NA_real_)) %>% 
    mutate(tb_T2DMdiag_metfin = case_when(exp_bin_treat == 1 ~ as.numeric(difftime(exp_date_metfin_first, cov_date_t2dm, units = "days")),
                                            TRUE ~ NA_real_)) %>% # time between T2DM and metformin start
    # among the above, define outcome status, binary and date
    mutate(out_bin_severecovid = case_when(out_date_severe_covid > landmark_date ~ 1,
                                             TRUE ~ NA_real_)) %>% # covid-related death or hospitalisation thereafter
    mutate(out_bin_longcovid = case_when(out_date_long_fatigue > landmark_date ~ 1,
                                                 TRUE ~ NA_real_)) %>% # long covid or viral fatigue code thereafter
    mutate(out_date_severecovid = case_when(out_bin_severecovid == 1 ~ out_date_severe_covid,
      TRUE ~ NA_Date_)) %>%
    mutate(out_date_longcovid = case_when(out_bin_longcovid == 1 ~ out_date_long_fatigue,
      TRUE ~ NA_Date_)) %>%
    # among those with an outcome, check if they also have a covid diagnosis (defined highly sensitive, incl. primary/secondary care) and after landmark date (just to be sure)
    mutate(out_bin_severecovid_diagnosed = case_when(out_bin_severecovid == 1 & cov_date_covid19_first < out_date_severecovid & cov_date_covid19_first > landmark_date ~ 1, 
                                               TRUE ~ NA_real_)) %>% 
    mutate(out_bin_longcovid_diagnosed = case_when(out_bin_longcovid == 1 & cov_date_covid19_first > out_date_longcovid & cov_date_covid19_first > landmark_date ~ 1, 
                                                   TRUE ~ NA_real_)) %>% 
    # time between landmark and outcomes
    mutate(tb_landmark_longcovid = case_when(
             out_bin_longcovid == 1 ~ as.numeric(difftime(out_date_longcovid, landmark_date, units = "days")),
             TRUE ~ NA_real_),
           tb_landmark_severecovid = case_when(
             out_bin_severecovid == 1 ~ as.numeric(difftime(out_date_severecovid, landmark_date, units = "days")),
             TRUE ~ NA_real_)) %>%
    
    summarise(
      n_t2dm_midpoint6 = fn_roundmid_any(n(), threshold), # Perform redaction
      
      n_metfin_by_landmark_midpoint6 = fn_roundmid_any(sum(exp_bin_treat, na.rm = TRUE), threshold), # Perform redaction
      median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
      IQR_lower_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),
      IQR_upper_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE),
      
      n_severeCOVID_midpoint6 = fn_roundmid_any(sum(out_bin_severecovid, na.rm = TRUE), threshold), # Perform redaction
      median_tb_landmark_severecovid = median(tb_landmark_severecovid, na.rm = TRUE),
      IQR_lower_tb_landmark_severecovid = quantile(tb_landmark_severecovid, 0.25, na.rm = TRUE),
      IQR_upper_tb_landmark_severecovid = quantile(tb_landmark_severecovid, 0.75, na.rm = TRUE),
      
      n_LongCOVID_midpoint6 = fn_roundmid_any(sum(out_bin_longcovid, na.rm = TRUE), threshold), # Perform redaction
      median_tb_landmark_longcovid = median(tb_landmark_longcovid, na.rm = TRUE),
      IQR_lower_tb_landmark_longcovid = quantile(tb_landmark_longcovid, 0.25, na.rm = TRUE),
      IQR_upper_tb_landmark_longcovid = quantile(tb_landmark_longcovid, 0.75, na.rm = TRUE),
      
      n_severeCOVID_COVIDdiag_midpoint6 = fn_roundmid_any(sum(out_bin_severecovid_diagnosed, na.rm = TRUE), threshold), # Perform redaction
      n_LongCOVID_COVIDdiag_midpoint6 = fn_roundmid_any(sum(out_bin_longcovid_diagnosed, na.rm = TRUE), threshold) # Perform redaction
    ) %>%
    
    mutate("years in days" = years_in_days) # for identification
}

t2dm_midpoint6 <- map_dfr(years_in_days, fn_t2dm_midpoint6)



# ## 1. How many with a T2DM diagnosis start metformin within 5d, 7d, 10d, 30d, and 90d after (the first) pos. SARS-CoV-2 test?
# # Use T2DM variable from diabetes algo
# # Add outcomes
# grace_periods <- c(5, 7, 10, 30, 90) # grace periods in days (30 and 90 are probably not useful as grace period, but still of interest)
# 
# ## function with midpoint6 rounding
# fn_t2dm_covid_metfin_start_midpoint6 <- function(grace_periods) {
#   data_processed %>%
#     filter(!is.na(cov_date_t2dm)) %>% # has a T2DM diagnosis
#     filter(exp_date_metfin_first > cov_date_t2dm) %>% # metformin AFTER T2DM diagnosis
#     filter(!is.na(cov_date_covid19_first)) %>% # has COVID-19
#     filter(cov_date_covid19_first < exp_date_metfin_first &
#              exp_date_metfin_first <= cov_date_covid19_first + days(grace_periods)) %>% # received metformin within grace period after COVID diagnosis and no metformin before (because it's _metformin_first)
#     mutate(out_bin_covid_outcome = case_when(out_date_death_covid > exp_date_metfin_first |
#                                                out_date_covid19_hes_first > exp_date_metfin_first |
#                                                out_date_covid19_hes_last > exp_date_metfin_first |
#                                                out_date_covid19_emergency_first > exp_date_metfin_first |
#                                                out_date_covid19_emergency_last > exp_date_metfin_first ~ 1,
#                                              TRUE ~ 0)) %>% # among the above, search for those with a covid-related death or hospitalisation thereafter
#     mutate(out_bin_longcovid_outcome = case_when(out_date_long_fatigue_first > exp_date_metfin_first ~ 1,
#                                              TRUE ~ 0)) %>% # among the above, search for those with a long covid or viral fatigue code thereafter
#     mutate(out_date_covid_outcome = case_when(
#       out_bin_covid_outcome == 1 ~ pmin(out_date_death_covid,
#                                         out_date_covid19_hes_first,
#                                         out_date_covid19_hes_last,
#                                         out_date_covid19_emergency_first,
#                                         out_date_covid19_emergency_last,
#                                         na.rm = TRUE), # choose the first event
#       TRUE ~ NA_Date_)) %>%
#     mutate(tb_T2DMdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_t2dm, units = "days")),
#            tb_COVIDdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_covid19_first, units = "days")),
#            tb_metfin_covid_outcome = case_when(
#              out_bin_covid_outcome == 1 ~ as.numeric(difftime(exp_date_metfin_first, out_date_covid_outcome, units = "days")),
#              TRUE ~ NA_real_),
#            tb_metfin_longcovid_outcome = case_when(
#              out_bin_longcovid_outcome == 1 ~ as.numeric(difftime(exp_date_metfin_first, out_date_long_fatigue_first, units = "days")),
#              TRUE ~ NA_real_)
#            ) %>%
#     summarise(
#       n_start_metfin_aCOVID_midpoint6 = fn_roundmid_any(n(), threshold), # Perform redaction
#       median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
#       IQR_lower_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),
#       IQR_upper_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE),
#       #median_tb_COVIDdiag_metfin = median(tb_COVIDdiag_metfin, na.rm = TRUE),
#       #IQR_lower_tb_COVIDdiag_metfin = quantile(tb_COVIDdiag_metfin, 0.25, na.rm = TRUE),
#       #IQR_upper_tb_COVIDdiag_metfin = quantile(tb_COVIDdiag_metfin, 0.75, na.rm = TRUE),
#       n_COVIDoutcome_midpoint6 = fn_roundmid_any(sum(out_bin_covid_outcome, na.rm = TRUE), threshold), # Perform redaction
#       median_tb_metfin_COVIDoutcome = median(tb_metfin_covid_outcome, na.rm = TRUE),
#       IQR_lower_tb_metfin_COVIDoutcome = quantile(tb_metfin_covid_outcome, 0.25, na.rm = TRUE),
#       IQR_upper_tb_metfin_COVIDoutcome = quantile(tb_metfin_covid_outcome, 0.75, na.rm = TRUE),
#       n_LongCOVIDoutcome_midpoint6 = fn_roundmid_any(sum(out_bin_longcovid_outcome, na.rm = TRUE), threshold), # Perform redaction
#       median_tb_metfin_LongCOVIDoutcome = median(tb_metfin_longcovid_outcome, na.rm = TRUE),
#       IQR_lower_tb_metfin_LongCOVIDoutcome = quantile(tb_metfin_longcovid_outcome, 0.25, na.rm = TRUE),
#       IQR_upper_tb_metfin_LongCOVIDoutcome = quantile(tb_metfin_longcovid_outcome, 0.75, na.rm = TRUE)
#     ) %>%
# 
#     mutate("grace period in days" = grace_periods) # for identification
# }
# n_t2dm_covid_metfin_start_midpoint6 <- map_dfr(grace_periods, fn_t2dm_covid_metfin_start_midpoint6)

# ## function without midpoint6 rounding
# fn_t2dm_covid_metfin_start <- function(grace_periods) {
#   data_processed %>%
#     filter(!is.na(cov_date_t2dm)) %>% # has a T2DM diagnosis
#     filter(exp_date_metfin_first > cov_date_t2dm) %>% # metformin AFTER T2DM diagnosis
#     filter(!is.na(cov_date_covid19_first)) %>% # has COVID-19
#     filter(cov_date_covid19_first < exp_date_metfin_first &
#              exp_date_metfin_first <= cov_date_covid19_first + days(grace_periods)) %>% # received metformin within grace period after COVID diagnosis and no metformin before (because it's _metformin_first)
#     mutate(out_bin_covid_outcome = case_when(out_date_death_covid > exp_date_metfin_first | 
#                                                out_date_covid19_hes_first > exp_date_metfin_first | 
#                                                out_date_covid19_hes_last > exp_date_metfin_first |
#                                                out_date_covid19_emergency_first > exp_date_metfin_first |
#                                                out_date_covid19_emergency_last > exp_date_metfin_first ~ 1,
#                                              TRUE ~ 0)) %>% # among the above, search for those with a covid-related death or hospitalisation thereafter
#     mutate(out_bin_longcovid_outcome = case_when(out_date_long_fatigue_first > exp_date_metfin_first ~ 1,
#                                                  TRUE ~ 0)) %>% # among the above, search for those with a long covid or viral fatigue code thereafter
#     mutate(out_date_covid_outcome = case_when(
#       out_bin_covid_outcome == 1 ~ pmin(out_date_death_covid, 
#                                         out_date_covid19_hes_first, 
#                                         out_date_covid19_hes_last, 
#                                         out_date_covid19_emergency_first, 
#                                         out_date_covid19_emergency_last, 
#                                         na.rm = TRUE), # choose the first event
#       TRUE ~ NA_Date_)) %>%
#     mutate(tb_T2DMdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_t2dm, units = "days")),
#            tb_COVIDdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_covid19_first, units = "days")),
#            tb_metfin_covid_outcome = case_when(
#              out_bin_covid_outcome == 1 ~ as.numeric(difftime(exp_date_metfin_first, out_date_covid_outcome, units = "days")),
#              TRUE ~ NA_real_),
#            tb_metfin_longcovid_outcome = case_when(
#              out_bin_longcovid_outcome == 1 ~ as.numeric(difftime(exp_date_metfin_first, out_date_long_fatigue_first, units = "days")),
#              TRUE ~ NA_real_)
#     ) %>%
#     summarise(
#       n_start_metfin_aCOVID = n(),
#       median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
#       IQR_lower_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),
#       IQR_upper_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE),
#       #median_tb_COVIDdiag_metfin = median(tb_COVIDdiag_metfin, na.rm = TRUE),
#       #IQR_lower_tb_COVIDdiag_metfin = quantile(tb_COVIDdiag_metfin, 0.25, na.rm = TRUE),
#       #IQR_upper_tb_COVIDdiag_metfin = quantile(tb_COVIDdiag_metfin, 0.75, na.rm = TRUE),
#       n_COVIDoutcome = sum(out_bin_covid_outcome, na.rm = TRUE), 
#       median_tb_metfin_COVIDoutcome = median(tb_metfin_covid_outcome, na.rm = TRUE),
#       IQR_lower_tb_metfin_COVIDoutcome = quantile(tb_metfin_covid_outcome, 0.25, na.rm = TRUE),
#       IQR_upper_tb_metfin_COVIDoutcome = quantile(tb_metfin_covid_outcome, 0.75, na.rm = TRUE),
#       n_LongCOVIDoutcome = sum(out_bin_longcovid_outcome, na.rm = TRUE), 
#       median_tb_metfin_LongCOVIDoutcome = median(tb_metfin_longcovid_outcome, na.rm = TRUE),
#       IQR_lower_tb_metfin_LongCOVIDoutcome = quantile(tb_metfin_longcovid_outcome, 0.25, na.rm = TRUE),
#       IQR_upper_tb_metfin_LongCOVIDoutcome = quantile(tb_metfin_longcovid_outcome, 0.75, na.rm = TRUE)
#     ) %>%
#     
#     mutate("grace period in days" = grace_periods) # for identification
# }
# n_t2dm_covid_metfin_start <- map_dfr(grace_periods, fn_t2dm_covid_metfin_start)
# 
# # gt table
# # n_t2dm_covid_metfin_start_midpoint6 %>%
# #   gt() %>%
# #   tab_header(title = "Metformin start within 10d after COVID diagnosis among T2DM") %>%
# #   cols_label(
# #     n_start_metfin_aCOVID_midpoint6 = "N starting metformin 10d after COVID diagnosis",
# #     median_tb_T2DMdiag_metfin = "Median time to metformin start after T2DM",
# #     IQR_lower = "IQR Lower",
# #     IQR_upper = "IQR Upper")
# 
# 
# ## 2. How many with a T2DM diagnosis start metformin thereafter (in monthly periods til max 6m)?
# # Use T2DM variable from diabetes algo
# # Add outcomes
# periods <- c(30, 60, 90, 120, 150, 180) # in days - but now, we are interested in e.g. monthly until max. 6-monthly intervals
# 
# ## function with midpoint6 rounding
# fn_t2dm_metfin_start_midpoint6 <- function(periods){
#   data_processed %>%
#     filter(!is.na(cov_date_t2dm)) %>%
#     filter(!is.na(exp_date_metfin_first)) %>%
#     filter(exp_date_metfin_first > cov_date_t2dm) %>%
#     filter(exp_date_metfin_first <= (cov_date_t2dm + days(periods))) %>%
#     
#     mutate(out_bin_covid_outcome = case_when(out_date_death_covid > exp_date_metfin_first | 
#                                                out_date_covid19_hes_first > exp_date_metfin_first | 
#                                                out_date_covid19_hes_last > exp_date_metfin_first |
#                                                out_date_covid19_emergency_first > exp_date_metfin_first |
#                                                out_date_covid19_emergency_last > exp_date_metfin_first ~ 1,
#                                              TRUE ~ 0)) %>% # among the above, search for those with a covid-related death or hospitalisation thereafter
#     mutate(out_bin_longcovid_outcome = case_when(out_date_long_fatigue_first > exp_date_metfin_first ~ 1,
#                                                  TRUE ~ 0)) %>% # among the above, search for those with a long covid or viral fatigue code thereafter
#     mutate(out_bin_covid_diagnosis = case_when(out_bin_covid_outcome == 1 & cov_date_covid19_first > exp_date_metfin_first ~ 1, 
#                                                TRUE ~ 0)) %>% # among those with a covid-related death or hosp, how many have a covid diagnosis before (but after metformin start)?
#     mutate(out_bin_longcovid_diagnosis = case_when(out_bin_longcovid_outcome == 1 & cov_date_covid19_first > exp_date_metfin_first ~ 1, 
#                                                TRUE ~ 0)) %>% # among those with long covid, how many have a covid diagnosis before (but after metformin start)?
#     mutate(out_date_covid_outcome = case_when(
#       out_bin_covid_outcome == 1 ~ pmin(out_date_death_covid, 
#                                         out_date_covid19_hes_first, 
#                                         out_date_covid19_hes_last, 
#                                         out_date_covid19_emergency_first, 
#                                         out_date_covid19_emergency_last, 
#                                         na.rm = TRUE), # choose the first event
#       TRUE ~ NA_Date_)) %>%
#     mutate(tb_T2DMdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_t2dm, units = "days")),
#            tb_metfin_covid_outcome = case_when(
#              out_bin_covid_outcome == 1 ~ as.numeric(difftime(exp_date_metfin_first, out_date_covid_outcome, units = "days")),
#              TRUE ~ NA_real_),
#            tb_metfin_longcovid_outcome = case_when(
#              out_bin_longcovid_outcome == 1 ~ as.numeric(difftime(exp_date_metfin_first, out_date_long_fatigue_first, units = "days")),
#              TRUE ~ NA_real_),
#            tb_metfin_covid_diagnosis = case_when(
#              out_bin_longcovid_outcome == 1 | out_bin_covid_outcome == 1 ~ as.numeric(difftime(exp_date_metfin_first, cov_date_covid19_first, units = "days")),
#              TRUE ~ NA_real_),
#            
#     ) %>%
#     summarise(
#       n_start_metfin_aT2DM_midpoint6 = fn_roundmid_any(n(), threshold), # Perform redaction
#       median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
#       IQR_lower_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),
#       IQR_upper_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE),
#       n_COVIDoutcome_midpoint6 = fn_roundmid_any(sum(out_bin_covid_outcome, na.rm = TRUE), threshold), # Perform redaction
#       median_tb_metfin_COVIDoutcome = median(tb_metfin_covid_outcome, na.rm = TRUE),
#       IQR_lower_tb_metfin_COVIDoutcome = quantile(tb_metfin_covid_outcome, 0.25, na.rm = TRUE),
#       IQR_upper_tb_metfin_COVIDoutcome = quantile(tb_metfin_covid_outcome, 0.75, na.rm = TRUE),
#       n_LongCOVIDoutcome_midpoint6 = fn_roundmid_any(sum(out_bin_longcovid_outcome, na.rm = TRUE), threshold), # Perform redaction
#       median_tb_metfin_LongCOVIDoutcome = median(tb_metfin_longcovid_outcome, na.rm = TRUE),
#       IQR_lower_tb_metfin_LongCOVIDoutcome = quantile(tb_metfin_longcovid_outcome, 0.25, na.rm = TRUE),
#       IQR_upper_tb_metfin_LongCOVIDoutcome = quantile(tb_metfin_longcovid_outcome, 0.75, na.rm = TRUE),
#       n_COVIDoutcome_COVIDdiag_midpoint6 = fn_roundmid_any(sum(out_bin_covid_diagnosis, na.rm = TRUE), threshold), # Perform redaction
#       n_LongCOVIDoutcome_COVIDdiag_midpoint6 = fn_roundmid_any(sum(out_bin_longcovid_diagnosis, na.rm = TRUE), threshold), # Perform redaction
#       median_tb_metfin_covid_diagnosis = median(tb_metfin_covid_diagnosis, na.rm = TRUE),
#       IQR_lower_tb_metfin_covid_diagnosis = quantile(tb_metfin_covid_diagnosis, 0.25, na.rm = TRUE),
#       IQR_upper_tb_metfin_covid_diagnosis = quantile(tb_metfin_covid_diagnosis, 0.75, na.rm = TRUE)
#     ) %>%
# 
#     mutate("period in days" = periods)
# }
# n_t2dm_metfin_start_midpoint6 <- map_dfr(periods, fn_t2dm_metfin_start_midpoint6)
# 
# ## function without midpoint6 rounding
# fn_t2dm_metfin_start <- function(periods){
#   data_processed %>%
#     filter(!is.na(cov_date_t2dm)) %>%
#     filter(!is.na(exp_date_metfin_first)) %>%
#     filter(exp_date_metfin_first > cov_date_t2dm) %>%
#     filter(exp_date_metfin_first <= (cov_date_t2dm + days(periods))) %>%
#     
#     mutate(out_bin_covid_outcome = case_when(out_date_death_covid > exp_date_metfin_first | 
#                                                out_date_covid19_hes_first > exp_date_metfin_first | 
#                                                out_date_covid19_hes_last > exp_date_metfin_first |
#                                                out_date_covid19_emergency_first > exp_date_metfin_first |
#                                                out_date_covid19_emergency_last > exp_date_metfin_first ~ 1,
#                                              TRUE ~ 0)) %>% # among the above, search for those with a covid-related death or hospitalisation thereafter
#     mutate(out_bin_longcovid_outcome = case_when(out_date_long_fatigue_first > exp_date_metfin_first ~ 1,
#                                                  TRUE ~ 0)) %>% # among the above, search for those with a long covid or viral fatigue code thereafter
#     mutate(out_bin_covid_diagnosis = case_when(out_bin_covid_outcome == 1 & cov_date_covid19_first > exp_date_metfin_first ~ 1, 
#                                                TRUE ~ 0)) %>% # among those with a covid-related death or hosp, how many have a covid diagnosis before (but after metformin start)?
#     mutate(out_bin_longcovid_diagnosis = case_when(out_bin_longcovid_outcome == 1 & cov_date_covid19_first > exp_date_metfin_first ~ 1, 
#                                                    TRUE ~ 0)) %>% # among those with long covid, how many have a covid diagnosis before (but after metformin start)?
#     mutate(out_date_covid_outcome = case_when(
#       out_bin_covid_outcome == 1 ~ pmin(out_date_death_covid, 
#                                         out_date_covid19_hes_first, 
#                                         out_date_covid19_hes_last, 
#                                         out_date_covid19_emergency_first, 
#                                         out_date_covid19_emergency_last, 
#                                         na.rm = TRUE), # choose the first event
#       TRUE ~ NA_Date_)) %>%
#     mutate(tb_T2DMdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_t2dm, units = "days")),
#            tb_metfin_covid_outcome = case_when(
#              out_bin_covid_outcome == 1 ~ as.numeric(difftime(exp_date_metfin_first, out_date_covid_outcome, units = "days")),
#              TRUE ~ NA_real_),
#            tb_metfin_longcovid_outcome = case_when(
#              out_bin_longcovid_outcome == 1 ~ as.numeric(difftime(exp_date_metfin_first, out_date_long_fatigue_first, units = "days")),
#              TRUE ~ NA_real_),
#            tb_metfin_covid_diagnosis = case_when(
#              out_bin_longcovid_outcome == 1 | out_bin_covid_outcome == 1 ~ as.numeric(difftime(exp_date_metfin_first, cov_date_covid19_first, units = "days")),
#              TRUE ~ NA_real_),
#            
#     ) %>%
#     summarise(
#       n_start_metfin_aT2DM = n(), 
#       median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
#       IQR_lower_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),
#       IQR_upper_tb_T2DMdiag_metfin = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE),
#       n_COVIDoutcome = sum(out_bin_covid_outcome, na.rm = TRUE), 
#       median_tb_metfin_COVIDoutcome = median(tb_metfin_covid_outcome, na.rm = TRUE),
#       IQR_lower_tb_metfin_COVIDoutcome = quantile(tb_metfin_covid_outcome, 0.25, na.rm = TRUE),
#       IQR_upper_tb_metfin_COVIDoutcome = quantile(tb_metfin_covid_outcome, 0.75, na.rm = TRUE),
#       n_LongCOVIDoutcome = sum(out_bin_longcovid_outcome, na.rm = TRUE), 
#       median_tb_metfin_LongCOVIDoutcome = median(tb_metfin_longcovid_outcome, na.rm = TRUE),
#       IQR_lower_tb_metfin_LongCOVIDoutcome = quantile(tb_metfin_longcovid_outcome, 0.25, na.rm = TRUE),
#       IQR_upper_tb_metfin_LongCOVIDoutcome = quantile(tb_metfin_longcovid_outcome, 0.75, na.rm = TRUE),
#       n_COVIDoutcome_COVIDdiag = sum(out_bin_covid_diagnosis, na.rm = TRUE), 
#       n_LongCOVIDoutcome_COVIDdiag = sum(out_bin_longcovid_diagnosis, na.rm = TRUE),
#       median_tb_metfin_covid_diagnosis = median(tb_metfin_covid_diagnosis, na.rm = TRUE),
#       IQR_lower_tb_metfin_covid_diagnosis = quantile(tb_metfin_covid_diagnosis, 0.25, na.rm = TRUE),
#       IQR_upper_tb_metfin_covid_diagnosis = quantile(tb_metfin_covid_diagnosis, 0.75, na.rm = TRUE)
#     ) %>%
#     
#     mutate("period in days" = periods)
# }
# n_t2dm_metfin_start <- map_dfr(periods, fn_t2dm_metfin_start)

################################################################################
# 3 Save output
################################################################################
# write.csv(n_t2dm_covid_metfin_start_midpoint6, file = here::here("output", "data_properties", "n_t2dm_covid_metfin_start_midpoint6.csv"))
# write.csv(n_t2dm_covid_metfin_start, file = here::here("output", "data_properties", "n_t2dm_covid_metfin_start.csv"))
# write.csv(n_t2dm_metfin_start_midpoint6, file = here::here("output", "data_properties", "n_t2dm_metfin_start_midpoint6.csv"))
# write.csv(n_t2dm_metfin_start, file = here::here("output", "data_properties", "n_t2dm_metfin_start.csv"))

write.csv(t2dm_midpoint6, file = here::here("output", "data_properties", "t2dm_midpoint6.csv"))
