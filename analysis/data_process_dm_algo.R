################################################################################
## This script does the following:
# 1. Import feather dataset from OpenSAFELY 
# 2. Basic formatting of variables -> fn_extract_data.R()
# 3. Apply the diabetes algorithm -> fn_diabetes_algorithm()
# 4. Save the output: data_processed
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

## Import custom user functions
source(here::here("analysis", "functions", "fn_extract_data.R"))
source(here::here("analysis", "functions", "utility.R"))
source(here::here("analysis", "functions", "fn_diabetes_algorithm.R"))

################################################################################
# 0.1 Create directories for output
################################################################################
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_properties"))

################################################################################
# 0.2 Import command-line arguments and dates # to be adapted at a later stage
################################################################################
args <- commandArgs(trailingOnly=TRUE)
# study_dates <-
#    jsonlite::read_json(path = here::here("output", "study_dates.json")) %>%
#    map(as.Date)
source(here::here("analysis", "metadates.R"))
# Convert the meta-dates into Date objects
study_dates <- lapply(study_dates, function(x) as.Date(x))

################################################################################
# 0.3 Define redaction threshold
################################################################################
threshold <- 6

################################################################################
# 1 Import data
################################################################################
input_filename <- "dataset_dm_algo.arrow"

################################################################################
# 2 Reformat the imported data
################################################################################
data_extracted <- fn_extract_data(input_filename)

################################################################################
# 3 Apply diabetes algorithm
################################################################################
# apply diabetes algorithm and delete all helper variables (tmp & step) at the end
data_processed_dm_algo <- fn_diabetes_algorithm(data_extracted)

################################################################################
# 4 Save output
################################################################################
# the data
write_csv(data_processed_dm_algo, here::here("output", "data", "data_processed_dm_algo.csv.gz"), na = "") # important to specify (na = ""), otherwise @table_from_file("output/data/data_processed_dm_algo.csv.gz") doesn't work (accept empty or a date, but not NA)