####
## This script does the following:
# 1. Import/extract the PPA feather dataset from OpenSAFELY
# 2. Basic type formatting of variables -> fn_extract_data.R()
# 3. Process covariates
# 4. Modify dummy data if run locally
# 5. Save the processed PPA dataset
####


# Import libraries and user functions -------------------------------------
library('arrow')
library('readr')
library('here')
library('dplyr')
library('tidyr')
library('purrr')
library('skimr')
source(here::here("analysis", "functions", "fn_extract_data.R"))


# Create directories for output -------------------------------------------
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "data_description"))


# Import the dataset and pre-process --------------------------------------
input_filename <- "dataset_ppa.arrow"
df <- fn_extract_data(input_filename)


# Modify dummy data -------------------------------------------------------
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  message("Running locally, adapt dummy data")
  source("analysis/modify_dummy_data_ppa.R")
  message("Dummy data successfully modified")
}


# Process the data --------------------------------------------------------
# Define thresholds for each type of biomarker (see data_process)
hba1c_min <- 0
hba1c_max <- 120

chol_min <- 1.75
chol_max <- 20

hdl_min <- 0.4
hdl_max <- 5

df <- df %>%
  mutate(
    across(
      matches("^cov_num_hba1c_\\d+$"),
      ~ ifelse(.x < hba1c_min | .x > hba1c_max, NA_real_, .x)
    ),
    across(
      matches("^cov_num_chol_\\d+$"),
      ~ ifelse(.x < chol_min | .x > chol_max, NA_real_, .x)
    ),
    across(
      matches("^cov_num_hdl_chol_\\d+$"),
      ~ ifelse(.x < hdl_min | .x > hdl_max, NA_real_, .x)
    )
  )


# Filter main dataset: Only keep necessary variables ----------------------
df <- df %>% 
  select(!landmark_date)


# Save and inspect full processed dataset ---------------------------------
data_processed_ppa_desc <- skim(df)
write.csv(data_processed_ppa_desc, file = here::here("output", "data_description", "data_processed_ppa_desc.csv")) # for L4 reviewing only, not for release
arrow::write_feather(df, here::here("output", "data", "data_processed_ppa.arrow"))
