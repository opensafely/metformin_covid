################################################################################
# SEQTRIALS PIPELINE - apply the T2DM window, data quality assurance (QA) and
# completeness criteria, and record a flowchart of the exclusions
#
# 1. Import dates, the diabetes algorithm output (for the starting counts) and the
#    QA/completeness dataset (people with a T2DM diagnosis in the T2DM window)
# 2. Modify the dummy data (only when run locally)
# 3. Flowchart steps, each person is counted once at the first criterion they fail:
#    a) Study population: all people -> T2DM diagnosis -> T2DM diagnosis in the window
#    b) QA criteria 1-7
#    c) Completeness criteria 1-3
# 4. Save the remaining people (input for dataset 3) and the flowchart
################################################################################

################################################################################
# 0.0 Import libraries and functions
################################################################################
library('arrow')
library('readr')
library('here')
library('lubridate')
library('dplyr')
library('jsonlite')
source(here::here("analysis", "functions", "fn_extract_data.R"))
source(here::here("analysis", "functions", "utility.R"))
source(here::here("analysis", "seqtrials", "functions", "fn_flowchart.R"))
source(here::here("analysis", "seqtrials", "functions", "fn_modify_dummy_data_qa.R"))

################################################################################
# 0.1 Create directories for output
################################################################################
fs::dir_create(here::here("output", "seqtrials", "data"))
fs::dir_create(here::here("output", "seqtrials", "data_description"))

################################################################################
# 0.2 Import dates and define the T2DM window
################################################################################
study_dates <- read_json(here::here("output", "seqtrials", "study_dates.json"))
studystart_date <- as.Date(study_dates$studystart_date)
studyend_date <- as.Date(study_dates$studyend_date)
# Trials start at diagnosis + 30*k days (k = 0, ..., 6), i.e. up to 180 days after
# diagnosis, so the window starts 180 days before study start
# (identical to the population definition in dataset_definition_t2dm_qa.py)
t2dm_window_start <- studystart_date - days(180)

# Redaction threshold
threshold <- 6

################################################################################
# 1.0 Starting counts from the diabetes algorithm output (all people in TPP)
################################################################################
print("Starting counts from the diabetes algorithm output")
data_dm_algo <- read_csv(
  here::here("output", "seqtrials", "data_processed_dm_algo.csv.gz"),
  col_select = c(patient_id, t2dm_date),
  col_types = cols(patient_id = col_double(), t2dm_date = col_date())
)
n_all <- nrow(data_dm_algo)
n_t2dm <- sum(!is.na(data_dm_algo$t2dm_date))
n_t2dm_window <- sum(!is.na(data_dm_algo$t2dm_date) &
                       data_dm_algo$t2dm_date >= t2dm_window_start &
                       data_dm_algo$t2dm_date <= studyend_date)
rm(data_dm_algo)

################################################################################
# 1.1 Import the QA/completeness dataset (people with T2DM in the window)
################################################################################
print("Import the QA/completeness dataset")
data <- fn_extract_data(file.path("seqtrials", "dataset_t2dm_qa.arrow"))

# Check: the dataset population must match the T2DM window count above
if (nrow(data) != n_t2dm_window) {
  stop("Number of people in dataset_t2dm_qa.arrow (", nrow(data), ") differs from the ",
       "number with a T2DM diagnosis in the window in the diabetes algorithm output (",
       n_t2dm_window, ").")
}

################################################################################
# 2.0 Modify dummy data (only when run locally)
################################################################################
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  print("Running locally: modify dummy data")
  data <- fn_modify_dummy_data_qa(data, studyend_date)
}

################################################################################
# 3.0 Flowchart a) Study population
################################################################################
flowchart_population <- tibble(
  group = "Study population",
  criterion = c(
    "People in TPP (diabetes algorithm dataset)",
    "No T2DM diagnosis",
    paste0("T2DM diagnosis outside the T2DM window (", t2dm_window_start, " to ", studyend_date, ")")
  ),
  n_excluded = c(NA, n_all - n_t2dm, n_t2dm - n_t2dm_window),
  n_remaining = c(n_all, n_t2dm, n_t2dm_window)
)

################################################################################
# 3.1 Flowchart b) Data quality assurance criteria
################################################################################
print("Apply the quality assurance criteria")
qa_exclusions <- rlang::exprs(
  "QA 1: Year of birth missing" =
    is.na(qa_num_birth_year),
  "QA 2: Year of birth after year of death" =
    !is.na(qa_date_of_death) & qa_num_birth_year > year(qa_date_of_death),
  "QA 3: Year of birth before 1838 or after study end" =
    qa_num_birth_year < 1838 | qa_num_birth_year > year(studyend_date),
  "QA 4: Date of death on or before 1900-01-01 or after today" =
    !is.na(qa_date_of_death) & (qa_date_of_death <= as.Date("1900-01-01") | qa_date_of_death > Sys.Date()),
  "QA 5: Men with pregnancy codes" =
    cov_cat_sex == "male" & qa_bin_pregnancy,
  "QA 6: Men with HRT or combined oral contraceptive pill" =
    cov_cat_sex == "male" & (qa_bin_hrt | qa_bin_cocp),
  "QA 7: Women with prostate cancer codes" =
    cov_cat_sex == "female" & qa_bin_prostate_cancer
)
qa <- fn_apply_exclusions(data, qa_exclusions, group = "Quality assurance")
data <- qa$data

################################################################################
# 3.2 Flowchart c) Data completeness criteria
################################################################################
print("Apply the completeness criteria")
completeness_exclusions <- rlang::exprs(
  "Completeness 1: Sex not female or male" =
    !(cov_cat_sex %in% c("female", "male")),
  "Completeness 2: Deprivation (IMD) unknown at T2DM diagnosis" =
    !qa_bin_known_imd,
  "Completeness 3: STP unknown at T2DM diagnosis" =
    !qa_bin_known_stp
)
completeness <- fn_apply_exclusions(data, completeness_exclusions, group = "Completeness")
data <- completeness$data

################################################################################
# 4.0 Save output
################################################################################
print("Save output")
flowchart <- bind_rows(flowchart_population, qa$flowchart, completeness$flowchart)
flowchart_midpoint6 <- fn_round_flowchart(flowchart, threshold) %>%
  select(group, criterion, n_excluded_midpoint6, n_remaining_midpoint6)

# Flowchart, midpoint-6 rounded (for release)
write_csv(flowchart_midpoint6,
          here::here("output", "seqtrials", "data_description", "flowchart_qa_midpoint6.csv"))
# People passing all criteria (input for dataset 3, via table_from_file): only the
# patient_id and T2DM diagnosis date, all other variables are re-defined in dataset 3
write_feather(select(data, patient_id, elig_date_t2dm),
              here::here("output", "seqtrials", "data", "data_t2dm_qa.arrow"))

print(flowchart_midpoint6, n = Inf, width = Inf)
