################################################################################
## This script does the following:
# 1. Import processed data
# 2. Create table 1
# 3. Save all output
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
library('jsonlite')
library('skimr') # for skim
## Import custom user functions and meta-dates
source(here::here("analysis", "functions", "utility.R"))

################################################################################
# 0.1 Create directories for output
################################################################################
fs::dir_create(here::here("output", "data_properties"))

################################################################################
# 0.3 Define redaction threshold
################################################################################
threshold <- 6

################################################################################
# 1 Import the data
################################################################################
data_processed <- readRDS(here("output", "data", "data_processed.rds"))

################################################################################
# 2 Label the data
################################################################################
var_labels <- list(
  N  ~ "Total N",
  exp_bin_treat ~ "Treatment",
  
  cov_num_age ~ "Age",
  cov_cat_age ~ "Age groups",
  cov_cat_sex ~ "Sex",
  cov_cat_ethnicity ~ "Ethnicity",
  cov_cat_deprivation_5 ~ "Deprivation",
  cov_cat_region ~ "Region",
  cov_cat_stp ~ "STP",
  cov_cat_rural_urban ~ "Rural/urban",
  cov_cat_smoking_status ~ "Smoking status",
  cov_bin_carehome_status ~ "Care/nursing home resident",
  cov_bin_obesity ~ "Body Mass Index > 40 kg/m^2",
  cov_num_bmi ~ "Body Mass Index",
  cov_cat_bmi_groups ~ "Body Mass Index, groups",
  cov_bin_ami ~ "History of acute myocardial infarct",
  cov_bin_all_stroke  ~ "History of stroke",
  cov_bin_other_arterial_embolism ~ "History of other arterial embolism",
  cov_bin_vte ~ "History of venous thromboembolism",
  cov_bin_hf ~ "History of heart failure",
  cov_bin_angina ~ "History of angina pectoris",
  cov_bin_dementia ~ "History of dementia",
  cov_bin_cancer ~ "History of cancer",
  cov_bin_hypertension ~ "History of arterial hypertension",
  cov_bin_depression ~ "History of depression",
  cov_bin_copd ~ "History of COPD",
  cov_bin_liver_disease ~ "History of liver disease",
  cov_bin_chronic_kidney_disease ~ "History of CKD",
  cov_bin_pcos ~ "History of PCOS",
  cov_bin_prediabetes ~ "History of prediabetes",
  cov_bin_diabetescomp ~ "Diabetes complication",
  cov_num_hba1c_mmol_mol ~ "HbA1c in mmol/mol",
  cov_num_tc_hdl_ratio ~ "TC/Chol ratio",
  cov_bin_healthcare_worker ~ "Healthcare worker",
  
  exp_bin_metfin_pandemicstart ~ "Any metformin prescription within 6m prior to pandemic start",
  exp_bin_metfin_anytime ~ "Starting metformin anytime in control and anytime after landmark for intervention",
  out_bin_severecovid ~ "COVID hosp or death",
  out_bin_death_pandemicstart ~ "Death between landmark and pandemic start",
  out_bin_ltfu_pandemicstart ~ "LTFU between landmark and pandemic start"
) %>%
  set_names(., map_chr(., all.vars))



################################################################################
# 11 Save output
################################################################################
# the full data
write_rds(data_processed, here::here("output", "data", "data_processed-test.rds"))
