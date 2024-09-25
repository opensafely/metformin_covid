################################################################################
#
# Processing data
#
# This script can be run via an action in project.yaml
#
# The output of this script is:
# -./output/data/data_processed.rds
# - ./output/data_properties/n_excluded.rds
#
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
source(here::here("analysis", "data_import", "process_data.R"))
source(here::here("analysis", "data_import", "extract_data.R"))
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
    jsonlite::read_json(path = here::here("lib", "design", "study-dates.json")) %>%
    map(as.Date)

################################################################################
# 1 Import data
################################################################################
input_filename <- "dataset.arrow"
data_extracted <- extract_data(input_filename)

## dummy data issues?
# data_extracted %>% # why are all the deaths only covid? Why no noncovid deaths?
#   select(qa_date_of_death, out_bin_death_cause_covid) %>%
#   View()
# table(data_extracted$out_date_dereg) # why are there no dereg dates available?

# # change data if run using dummy data
# if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
#   data_extracted <-
#     data_extracted %>%
#     mutate(died_ons_covid_any_date =
#              if_else(!is.na(death_date), death_date, died_ons_covid_any_date),
#            death_date =
#              if_else(!is.na(died_ons_covid_any_date), died_ons_covid_any_date, death_date),
#            date_treated = if_else(!is.na(date_treated),
#                                   covid_test_positive_date + runif(nrow(data_extracted), 0, 4) %>% round(),
#                                   NA_Date_),
#            paxlovid_covid_therapeutics = if_else(!is.na(paxlovid_covid_therapeutics),
#                                                  date_treated,
#                                                  NA_Date_),
#            sotrovimab_covid_therapeutics = if_else(!is.na(sotrovimab_covid_therapeutics),
#                                                    date_treated,
#                                                    NA_Date_),
#            molnupiravir_covid_therapeutics = if_else(!is.na(molnupiravir_covid_therapeutics),
#                                                      date_treated,
#                                                      NA_Date_)
#     )
# }

################################################################################
# 2 Process data
################################################################################
data_processed_g10 <- process_data(data_extracted, study_dates, treat_window_days = 9) # grace period 10 dataset only

data_processed <-
  map(.x = list(6, 7, 8, 9), # add additional longer grace periods besides the primary grace period 10 days (baseline_date + 9)
      .f = ~ process_data(data_extracted, study_dates, treat_window_days = .x))
names(data_processed) <- c("grace7", "grace8", "grace9", "grace10")

# currently all deaths are covid-related (see above) and deregistration date not available -> modify dummy data?
# unique(data_processed_g10$out_date_death_28)
# unique(data_processed_g10$out_date_noncovid_death) # no noncovid deaths
# unique(data_processed_g10$out_date_covid_death)

# data_processed_g10 %>% # why do some have a date of death before baseline dates and marked qa_bin_was_alive == TRUE? dummy data?
#   select(patient_id, baseline_date, period_week, period_month, period_2month, period_3month, status_primary, fu_primary, qa_date_of_death, qa_bin_was_alive) %>%
#   View()

# # change data if run using dummy data
# if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
#   data_processed <-
#     map(.x = data_processed,
#         .f = ~ .x %>% group_by(patient_id) %>%
#           mutate(period_month = runif(1, 0, 12) %>% ceiling(),
#                  period_2month = runif(1, 0, 6) %>% ceiling(),
#                  period_3month = runif(1, 0, 4) %>% ceiling(),
#                  period_week = runif(1, 0, 52) %>% ceiling(),
#                  covid_test_positive_date = sample(seq(ymd("20220210"), ymd("20230209"), by = 1), 1),
#                  tb_postest_vacc_cat = sample(c(">= 84 days or unknown", "< 7 days", "7-27 days", "28-83 days"), 1) %>%
#                    factor(levels = c(">= 84 days or unknown", "< 7 days", "7-27 days", "28-83 days"))) %>%
#           ungroup() %>%
#           mutate(ageband = if_else(is.na(ageband), "18-39", ageband %>% as.character()) %>%
#                    factor(levels = c("18-39", "40-59", "60-79", "80+"))))
# }

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

# why are some periods (at the end) NA? Double-check!
# data_processed$grace10 %>%
#   select(patient_id, baseline_date, period_week, period_month, period_2month, period_3month, status_primary, fu_primary, qa_date_of_death, qa_bin_was_alive) %>%
#   View()

# # contraindications
# n_excluded_contraindicated <- calc_n_excluded_contraindicated(data_processed$grace10)
# data_processed <-
#   map(.x = data_processed,
#       .f = ~ .x %>% add_contraindicated_indicator())
# data_processed_excl_contraindicated  <-
#   map(.x = data_processed,
#       .f = ~ .x %>% filter(contraindicated == FALSE))

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
