################################################################################
# This script does the following:
# 1. Import/extract feather dataset from OpenSAFELY including basic type formatting of variables (extract_data())
# 2. Process the data (process_data()): a) combine some categories, 
#                                       b) apply the diabetes algorithm (diabetes_algo()), 
################################################################################

################################################################################
# 0.0 Import libraries + functions
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
source(here::here("analysis", "functions", "fn_case_when.R"))
source(here::here("analysis", "functions", "fn_diabetes_algorithm.R"))

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
## 1. Raw incidence curves for the different diabetes diagnoses and PCOS
# Using the diabetes algo variables for T2DM, T1DM and GestDM
data_processed <- data_processed %>% 
  mutate(
    year_t2dm = year(ymd(cov_date_t2dm)),  
    year_t1dm = year(ymd(cov_date_t1dm)),
    year_gestationaldm = year(ymd(cov_date_gestationaldm)),
    year_pcos = year(ymd(cov_date_pcos)),
    year_predm = year(ymd(cov_date_prediabetes))
  ) %>%
  pivot_longer(
    cols = starts_with("year_"),       
    names_to = "diagnosis",                  
    values_to = "year",              
    names_prefix = "year_"
  ) %>%
  filter(!is.na(year))                
# Summarize incidence by year and diagnosis
incidence_data <- data_processed %>%
  group_by(diagnosis, year) %>%
  summarise(count = n(), .groups = 'drop')
# Plot the incidence counts
plot_incidence <- ggplot(incidence_data, aes(x = year, y = count, color = diagnosis)) +
  geom_line(size = 1.2) +                   
  geom_point(size = 3) +  
  labs(title = "Raw incidence counts of diagnoses over time",
       x = "Year",
       y = "Number of Cases") +
  scale_color_manual(values = c("t2dm" = "blue", "t1dm" = "red", "gestationaldm" = "green", "pcos" = "yellow", "predm" = "orange")) + 
  theme_minimal()

## 2. How many with a T2DM diagnosis start metformin and median time to start ?
# T2DM variable from diabetes algo
n_t2dm_metfin_start <- data_processed %>%
  filter(!is.na(cov_date_t2dm)) %>%
  filter(!is.na(exp_date_metfin_first)) %>%
  filter(exp_date_metfin_first > cov_date_t2dm) %>%
  mutate(tb_T2DMdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_t2dm, units = "days"))) %>%
  summarise(
    n_start_metfin_aT2DM = n(),
    median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
    IQR_lower = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),  
    IQR_upper = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE))
# gt table
tbl_n_t2dm_metfin_start <- n_t2dm_metfin_start %>%
  gt() %>%
  tab_header(title = "Metformin start after T2DM diagnosis") %>%
  cols_label(
    n_start_metfin_aT2DM = "N starting metformin after T2DM diagnosis",
    median_tb_T2DMdiag_metfin = "Median time to metformin start (days)",
    IQR_lower = "IQR Lower",
    IQR_upper = "IQR Upper")

## 2. How many with a T2DM diagnosis start metformin within 10d after (the first) pos. SARS-CoV-2 test
# T2DM variable from diabetes algo
n_t2dm_covid_metfin_start <- data_processed %>%
  filter(!is.na(cov_date_t2dm)) %>%
  filter(!is.na(exp_date_metfin_first)) %>%
  filter(exp_date_metfin_first > cov_date_t2dm) %>%
  filter(!is.na(cov_date_covid19_first)) %>%
  filter(exp_date_metfin_first <= (cov_date_covid19_first + days(10))) %>% # grace period 10 days
  mutate(tb_T2DMdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_t2dm, units = "days"))) %>%
  mutate(tb_COVIDdiag_metfin = as.numeric(difftime(exp_date_metfin_first, cov_date_covid19_first, units = "days"))) %>%
  summarise(
    n_start_metfin_aCOVID = n(),
    median_tb_T2DMdiag_metfin = median(tb_T2DMdiag_metfin, na.rm = TRUE),
    IQR_lower = quantile(tb_T2DMdiag_metfin, 0.25, na.rm = TRUE),  
    IQR_upper = quantile(tb_T2DMdiag_metfin, 0.75, na.rm = TRUE))
# gt table
tbl_n_t2dm_covid_metfin_start <- n_t2dm_covid_metfin_start %>%
  gt() %>%
  tab_header(title = "Metformin start within 10d after COVID diagnosis among T2DM") %>%
  cols_label(
    n_start_metfin_aCOVID = "N starting metformin 10d after COVID diagnosis",
    median_tb_T2DMdiag_metfin = "Median time to metformin start after T2DM",
    IQR_lower = "IQR Lower",
    IQR_upper = "IQR Upper")

################################################################################
# 3 Save output
################################################################################
ggsave(plot_incidence, filename = here::here("output", "data_properties", "plot_incidence.png"))
gtsave(tbl_n_t2dm_metfin_start, filename = here::here("output", "data_properties", "tbl_n_t2dm_metfin_start.html"))
gtsave(tbl_n_t2dm_covid_metfin_start, filename = here::here("output", "data_properties", "tbl_n_t2dm_covid_metfin_start.html"))