################################################################################
# This script does the following:
# 1. Import/extract feather dataset from OpenSAFELY including basic type formatting of variables (extract_data())
# 2. Process the data (process_data()): a) combine some categories, 
#                                       b) apply the diabetes algorithm (diabetes_algo()), 
#                                       c) apply primary endpoint window (28 days) and create final outcome and exposure/treatment variables (add_status_and_fu_primary()),
#                                       d) add a treatment window period (grace period) that can be changed flexibly
#                                       e) define patients in periods according to their baseline_date (add_period_cuts()). Not needed for now, needed for SeqTrial if feasible 
# 3. Apply the data quality assurance and data completeness criteria (quality_assurance())
# 4. Apply the eligibility criteria (calc_n_excluded())
# 5. Summarize how many initiated treatment within applied grace period (calc_n_treated())
#
# Structure of this script based on https://github.com/opensafely/pax-non-users/tree/2dbf044472efdcfeb86f8fc2c8eea222e7eefe32/analysis
# Special credit to https://github.com/LindaNab
# Bristol diabetes algorithm based on and Credits to: https://github.com/opensafely/post-covid-diabetes/tree/main 
################################################################################

################################################################################
# 0.0 Import libraries + functions
################################################################################
library('feather')
library('readr')
library('here')
library('lubridate')
library('dplyr')
library('tidyr')
library('purrr')
## Import custom user functions
source(here::here("analysis", "data_import", "extract_data.R"))
source(here::here("analysis", "data_import", "process_data.R"))
source(here::here("analysis", "data_import", "quality_assurance.R"))
source(here::here("analysis", "data_import", "calc_n_excluded.R"))
source(here::here("analysis", "data_import", "calc_n_treated.R"))

################################################################################
# 0.1 Create directories for output
################################################################################
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_properties"))

################################################################################
# 0.2 Import command-line arguments
################################################################################
args <- commandArgs(trailingOnly=TRUE)
study_dates <-
    jsonlite::read_json(path = here::here("output", "study_dates.json")) %>%
    map(as.Date)

################################################################################
# 1 Import data
################################################################################
input_filename <- "dataset.arrow"
data_extracted <- extract_data(input_filename)

################################################################################
# 2 Process data
################################################################################
data_processed_g10 <- process_data(data_extracted, study_dates, treat_window_days = 9) # grace period 10

# add additional shorter grace periods besides the primary grace period 10 days (baseline_date + 9), 
# just to explore and for later (if SeqTrial applied, by running a trial each day within grace period)
data_processed <-
  map(.x = list(6, 7, 8, 9), 
      .f = ~ process_data(data_extracted, study_dates, treat_window_days = .x))
names(data_processed) <- c("grace7", "grace8", "grace9", "grace10")

################################################################################
# 3 Apply quality assurance criteria
################################################################################
n_qa_excluded <- quality_assurance(data_processed$grace10)

data_processed_g10 <- data_processed_g10 %>%
  filter(!is.na(qa_num_birth_year)) %>%
  filter(is.na(qa_date_of_death) | (qa_num_birth_year <= year(qa_date_of_death))) %>%
  filter(qa_num_birth_year >= 1793 & qa_num_birth_year <= year(Sys.Date())) %>%
  filter((qa_date_of_death > as.Date("1900-01-01")) | (qa_date_of_death < Sys.Date()) | is.na(qa_date_of_death)) %>%
  filter((cov_cat_sex == "Female" | is.na(cov_cat_sex)) | (cov_cat_sex == "Male" & (qa_bin_pregnancy == FALSE))) %>% # FALSE includes missing in a ehrQL logical
  filter((cov_cat_sex == "Female" | is.na(cov_cat_sex)) | (cov_cat_sex == "Male" & (qa_bin_hrt == FALSE)) | (cov_cat_sex == "Male" & (qa_bin_cocp == FALSE))) %>%
  filter((cov_cat_sex == "Male" | is.na(cov_cat_sex)) | (cov_cat_sex == "Female" & (qa_bin_prostate_cancer == FALSE)))

data_processed <-
  map(.x = data_processed,
      .f = ~ .x %>%
        filter(!is.na(qa_num_birth_year)) %>%
        filter(is.na(qa_date_of_death) | (qa_num_birth_year <= year(qa_date_of_death))) %>%
        filter(qa_num_birth_year >= 1793 & qa_num_birth_year <= year(Sys.Date())) %>%
        filter((qa_date_of_death > as.Date("1900-01-01")) | (qa_date_of_death < Sys.Date()) | is.na(qa_date_of_death)) %>%
        filter((cov_cat_sex == "Female" | is.na(cov_cat_sex)) | (cov_cat_sex == "Male" & (qa_bin_pregnancy == FALSE))) %>% # FALSE includes missing in a ehrQL logical
        filter((cov_cat_sex == "Female" | is.na(cov_cat_sex)) | (cov_cat_sex == "Male" & (qa_bin_hrt == FALSE)) | (cov_cat_sex == "Male" & (qa_bin_cocp == FALSE))) %>%
        filter((cov_cat_sex == "Male" | is.na(cov_cat_sex)) | (cov_cat_sex == "Female" & (qa_bin_prostate_cancer == FALSE)))
  )

################################################################################
# 4 Apply eligibility criteria
################################################################################
n_excluded <- calc_n_excluded(data_processed$grace10)

data_processed_g10 <- data_processed_g10 %>%
  # completeness criteria
  filter(qa_bin_was_alive == TRUE & (qa_date_of_death > baseline_date | is.na(qa_date_of_death))) %>% # additional condition since "qa_bin_was_alive == TRUE" may not cover all (e.g. pos test came out after death)
  filter(qa_bin_is_female_or_male == TRUE) %>%
  filter(qa_bin_known_imd == TRUE) %>%
  filter(!is.na(cov_cat_region)) %>%
  filter(qa_bin_was_registered == TRUE) %>%
  # inclusion criteria
  filter(qa_bin_was_adult == TRUE) %>%
  filter(cov_bin_t2dm == TRUE) %>%
  filter(!is.na(baseline_date)) %>%
  # exclusion criteria
  filter(cov_bin_hosp_baseline == FALSE) %>% # FALSE includes missing in a ehrQL logical
  filter(cov_bin_metfin_before_baseline == FALSE) %>%
  filter(cov_bin_metfin_allergy == FALSE) %>%
  filter(cov_bin_ckd_45 == FALSE) %>%
  filter(cov_bin_liver_cirrhosis == FALSE) %>%
  filter(cov_bin_metfin_interaction == FALSE) %>%
  filter(cov_bin_long_covid == FALSE)

data_processed <-
  map(.x = data_processed,
      .f = ~ .x %>%
        # completeness criteria
        filter(qa_bin_was_alive == TRUE & (qa_date_of_death > baseline_date | is.na(qa_date_of_death))) %>% # additional condition since "qa_bin_was_alive == TRUE" may not cover all (e.g. pos test came out after death)
        filter(qa_bin_is_female_or_male == TRUE) %>%
        filter(qa_bin_known_imd == TRUE) %>%
        filter(!is.na(cov_cat_region)) %>%
        filter(qa_bin_was_registered == TRUE) %>%
        # inclusion criteria
        filter(qa_bin_was_adult == TRUE) %>%
        filter(cov_bin_t2dm == TRUE) %>%
        filter(!is.na(baseline_date)) %>%
        # exclusion criteria
        filter(cov_bin_hosp_baseline == FALSE) %>% # FALSE includes missing in a ehrQL logical
        filter(cov_bin_metfin_before_baseline == FALSE) %>%
        filter(cov_bin_metfin_allergy == FALSE) %>%
        filter(cov_bin_ckd_45 == FALSE) %>%
        filter(cov_bin_liver_cirrhosis == FALSE) %>%
        filter(cov_bin_metfin_interaction == FALSE) %>%
        filter(cov_bin_long_covid == FALSE)
  )

################################################################################
# 5 Apply treatment window
################################################################################
n_treated <- calc_n_treated(data_processed$grace10)

################################################################################
# 6 Save data and steps/numbers of excluded participants
################################################################################
iwalk(.x = data_processed,
      .f = ~ write_rds(.x,
                       here::here("output", "data",
                                  paste0("data_processed", "_"[!.y == "grace10"],
                                    .y[!.y == "grace10"], ".rds"))))
write_rds(n_qa_excluded,
          here::here("output", "data_properties", "n_qa_excluded.rds"))
write_rds(n_excluded,
          here::here("output", "data_properties", "n_excluded.rds"))
write_rds(n_excluded,
          here::here("output", "data_properties", "n_treated.rds"))
