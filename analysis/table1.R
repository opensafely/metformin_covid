####
## This script does the following:
# 1. Import processed data
# 2. Create table 1
# 3. Save all output
####

# Import libraries and user functions -------------------------------------
library('arrow')
library('here')
library('tidyverse')
library('gtsummary')
source(here::here("analysis", "functions", "utility.R")) # fn_roundmid_any

# Create directories for output -------------------------------------------
fs::dir_create(here::here("output", "data_description"))

# Define redaction threshold ----------------------------------------------
threshold <- 6

# Import the data ---------------------------------------------------------
df <- read_feather(here("output", "data", "data_processed.arrow"))
df_death_ltfu <- read_feather(here("output", "data", "data_processed_death_ltfu.arrow"))

# Label the data ---------------------------------------------------------
var_labels_main <- list(
  N  = "Total N",
  exp_bin_treat = "Treatment",
  
  cov_num_age = "Age",
  cov_cat_age = "Age groups",
  cov_cat_sex = "Sex",
  cov_cat_ethnicity = "Ethnicity",
  cov_cat_deprivation_5 = "Deprivation",
  strat_cat_region = "Region",
  cov_cat_rural_urban = "Rural/urban",
  cov_cat_smoking_status = "Smoking status",
  cov_bin_healthcare_worker = "Healthcare worker",
  cov_num_consrate = "Consultation rate in previous year",
  cov_cat_bmi_groups = "Body Mass Index (BMI) categories",
  cov_cat_hba1c_b = "HbA1c categories in mmol/mol",
  cov_cat_tc_hdl_ratio_b = "TC/HDL ratio categories",
  cov_bin_ami = "History of acute myocardial infarct",
  cov_bin_all_stroke  = "History of stroke",
  cov_bin_other_arterial_embolism = "History of other arterial embolism",
  cov_bin_vte = "History of venous thromboembolism",
  cov_bin_hf = "History of heart failure",
  cov_bin_angina = "History of angina pectoris",
  cov_bin_dementia = "History of dementia",
  cov_bin_cancer = "History of cancer",
  cov_bin_hypertension = "History of arterial hypertension",
  cov_bin_depression = "History of depression",
  cov_bin_copd = "History of COPD",
  cov_bin_liver_disease = "History of liver disease",
  cov_bin_chronic_kidney_disease = "History of CKD",
  cov_bin_pcos = "History of PCOS",
  cov_bin_prediabetes = "History of prediabetes",
  cov_bin_diabetescomp = "Diabetes complication",
  cov_num_hba1c_b = "HbA1c in mmol/mol",
  cov_num_tc_hdl_ratio_b = "TC/Chol ratio",
  cov_num_counthba1c = "Count HbA1c measurements, previous 2y",
  cov_num_countlifestyle = "Count lifestyle advices, previous 2y",
  
  cov_num_period_month = "Calendar period of T2DM diagnosis",
  out_bin_severecovid_afterlandmark = "COVID hosp or death",
  out_bin_covid_hosp_afterlandmark = "COVID hosp",
  out_bin_covid_death_afterlandmark = "COVID death",
  out_bin_covid_afterlandmark = "Any covid diagnosis, pos test or hosp",
  out_bin_longcovid_afterlandmark = "Any Long COVID diagnosis",
  out_bin_virfat_afterlandmark = "Any Viral Fatigue diagnosis",
  out_bin_longcovid_virfat_afterlandmark = "Any Long COVID or Viral Fatigue diagnosis",
  out_bin_death_afterlandmark = "Any death after landmark",
  cens_bin_ltfu_afterlandmark = "Any LTFU after landmark",
  cens_bin_metfin_pandemicstart = "INT: Any metformin prescription within 6m prior to pandemic start",
  cens_bin_metfin_start_cont = "CONT: Any metformin start in control"
) 

var_labels_death_ltfu1 <- list(
  N  = "Total N",
  death_ltfu_landmark = "Died or LTFU until landmark",
  exp_bin_treat = "Metformin treatment",
  
  cov_num_age = "Age",
  cov_cat_age = "Age groups",
  cov_cat_sex = "Sex",
  cov_cat_ethnicity = "Ethnicity",
  cov_cat_deprivation_5 = "Deprivation",
  strat_cat_region = "Region",
  cov_cat_rural_urban = "Rural/urban",
  cov_cat_smoking_status = "Smoking status",
  cov_bin_healthcare_worker = "Healthcare worker",
  cov_num_consrate = "Consultation rate in previous year",
  cov_cat_bmi_groups = "Body Mass Index (BMI) categories",
  cov_cat_hba1c_b = "HbA1c categories in mmol/mol",
  cov_cat_tc_hdl_ratio_b = "TC/HDL ratio categories",
  cov_bin_ami = "History of acute myocardial infarct",
  cov_bin_all_stroke  = "History of stroke",
  cov_bin_other_arterial_embolism = "History of other arterial embolism",
  cov_bin_vte = "History of venous thromboembolism",
  cov_bin_hf = "History of heart failure",
  cov_bin_angina = "History of angina pectoris",
  cov_bin_dementia = "History of dementia",
  cov_bin_cancer = "History of cancer",
  cov_bin_hypertension = "History of arterial hypertension",
  cov_bin_depression = "History of depression",
  cov_bin_copd = "History of COPD",
  cov_bin_liver_disease = "History of liver disease",
  cov_bin_chronic_kidney_disease = "History of CKD",
  cov_bin_pcos = "History of PCOS",
  cov_bin_prediabetes = "History of prediabetes",
  cov_bin_diabetescomp = "Diabetes complication",
  cov_num_hba1c_b = "HbA1c in mmol/mol",
  cov_num_tc_hdl_ratio_b = "TC/Chol ratio",
  cov_num_counthba1c = "Count HbA1c measurements, previous 2y",
  cov_num_countlifestyle = "Count lifestyle advices, previous 2y",
  
  cov_num_period_month = "Calendar period of T2DM diagnosis",
  out_bin_severecovid_afterlandmark = "COVID hosp or death",
  out_bin_covid_hosp_afterlandmark = "COVID hosp",
  out_bin_covid_death_afterlandmark = "COVID death",
  out_bin_covid_afterlandmark = "Any covid diagnosis, pos test or hosp",
  out_bin_longcovid_afterlandmark = "Any Long COVID diagnosis",
  out_bin_virfat_afterlandmark = "Any Viral Fatigue diagnosis",
  out_bin_longcovid_virfat_afterlandmark = "Any Long COVID or Viral Fatigue diagnosis",
  out_bin_death_afterlandmark = "Any death after landmark",
  cens_bin_ltfu_afterlandmark = "Any LTFU after landmark",
  cens_bin_metfin_pandemicstart = "INT: Any metformin prescription within 6m prior to pandemic start",
  cens_bin_metfin_start_cont = "CONT: Any metformin start in control"
)

var_labels_death_ltfu2 <- list(
  N  = "Total N; without deaths/LTFU until landmark",
  death_ltfu_pandemic_without_landmark = "Died/LTFU between landmark and pandemic start",
  exp_bin_treat = "Metformin treatment",
  
  cov_num_age = "Age",
  cov_cat_age = "Age groups",
  cov_cat_sex = "Sex",
  cov_cat_ethnicity = "Ethnicity",
  cov_cat_deprivation_5 = "Deprivation",
  strat_cat_region = "Region",
  cov_cat_rural_urban = "Rural/urban",
  cov_cat_smoking_status = "Smoking status",
  cov_bin_healthcare_worker = "Healthcare worker",
  cov_num_consrate = "Consultation rate in previous year",
  cov_cat_bmi_groups = "Body Mass Index (BMI) categories",
  cov_cat_hba1c_b = "HbA1c categories in mmol/mol",
  cov_cat_tc_hdl_ratio_b = "TC/HDL ratio categories",
  cov_bin_ami = "History of acute myocardial infarct",
  cov_bin_all_stroke  = "History of stroke",
  cov_bin_other_arterial_embolism = "History of other arterial embolism",
  cov_bin_vte = "History of venous thromboembolism",
  cov_bin_hf = "History of heart failure",
  cov_bin_angina = "History of angina pectoris",
  cov_bin_dementia = "History of dementia",
  cov_bin_cancer = "History of cancer",
  cov_bin_hypertension = "History of arterial hypertension",
  cov_bin_depression = "History of depression",
  cov_bin_copd = "History of COPD",
  cov_bin_liver_disease = "History of liver disease",
  cov_bin_chronic_kidney_disease = "History of CKD",
  cov_bin_pcos = "History of PCOS",
  cov_bin_prediabetes = "History of prediabetes",
  cov_bin_diabetescomp = "Diabetes complication",
  cov_num_hba1c_b = "HbA1c in mmol/mol",
  cov_num_tc_hdl_ratio_b = "TC/Chol ratio",
  cov_num_counthba1c = "Count HbA1c measurements, previous 2y",
  cov_num_countlifestyle = "Count lifestyle advices, previous 2y",
  
  cov_num_period_month = "Calendar period of T2DM diagnosis",
  out_bin_severecovid_afterlandmark = "COVID hosp or death",
  out_bin_covid_hosp_afterlandmark = "COVID hosp",
  out_bin_covid_death_afterlandmark = "COVID death",
  out_bin_covid_afterlandmark = "Any covid diagnosis, pos test or hosp",
  out_bin_longcovid_afterlandmark = "Any Long COVID diagnosis",
  out_bin_virfat_afterlandmark = "Any Viral Fatigue diagnosis",
  out_bin_longcovid_virfat_afterlandmark = "Any Long COVID or Viral Fatigue diagnosis",
  out_bin_death_afterlandmark = "Any death after landmark",
  cens_bin_ltfu_afterlandmark = "Any LTFU after landmark",
  cens_bin_metfin_pandemicstart = "INT: Any metformin prescription within 6m prior to pandemic start",
  cens_bin_metfin_start_cont = "CONT: Any metformin start in control"
)

# Create the main table 1 ---------------------------------------------------------
table_1_main <-
  df %>%
  mutate(
    N=1L,
    exp_bin_treat = factor(exp_bin_treat,
                          levels = c(1,0),
                          labels = c("Metformin mono", "Nothing"))
  ) %>%
  select(names(var_labels_main)
  ) %>%
  tbl_summary(
    by = exp_bin_treat,
    label = var_labels_main,
    statistic = list(
      N ~ "{N}",
      all_continuous() ~ "{median} ({p25}, {p75});  {mean} ({sd})"
    ),
  )
table_1_main <- as_tibble(table_1_main$table_body)
table_1_main_redacted <- table_1_main %>%
  mutate(
    across(
      starts_with("stat_"),
      ~ stringr::str_replace_all(
        ., "\\d+",
        function(x) as.character(fn_roundmid_any(as.numeric(x), threshold))
      )
    ),
    var_label = factor(
      var_label,
      levels = map_chr(var_labels_main[-c(1, 2)], ~ last(as.character(.)))
    ),
    label = replace_na(as.character(label), "")
  )
table_1_main_redacted <- table_1_main_redacted %>% 
  rename(metformin = stat_1,
         control = stat_2)

# Create the death/ltfu1 table ---------------------------------------------------------
df_death_ltfu <- df_death_ltfu %>% 
  mutate(exp_bin_treat = factor(exp_bin_treat, 
                                levels = c(1,0), 
                                labels = c("Metformin mono", "Nothing")))
table_1_death_ltfu1 <-
  df_death_ltfu %>%
  mutate(
    N=1L,
    death_ltfu_landmark = factor(death_ltfu_landmark,
                           labels = c("Alive and in care at landmark", "Died or LTFU until landmark"))
  ) %>%
  select(names(var_labels_death_ltfu1)
  ) %>%
  tbl_summary(
    by = death_ltfu_landmark,
    label = var_labels_death_ltfu1[setdiff(names(.), "death_ltfu_landmark")],
    statistic = list(
      N ~ "{N}",
      all_continuous() ~ "{median} ({p25}, {p75});  {mean} ({sd})"
    ),
  )

table_1_death_ltfu1 <- as_tibble(table_1_death_ltfu1$table_body)
table_1_death_ltfu1_redacted <- table_1_death_ltfu1 %>%
  mutate(
    across(
      starts_with("stat_"),
      ~ stringr::str_replace_all(
        ., "\\d+",
        function(x) as.character(fn_roundmid_any(as.numeric(x), threshold))
      )
    ),
    var_label = factor(
      var_label,
      levels = map_chr(var_labels_main[-c(1, 2)], ~ last(as.character(.)))
    ),
    label = replace_na(as.character(label), "")
  )
table_1_death_ltfu1_redacted <- table_1_death_ltfu1_redacted %>% 
  rename(metformin = stat_1,
         control = stat_2)

# Create the death/ltfu2 table ---------------------------------------------------------
table_1_death_ltfu2 <-
  df_death_ltfu %>%
  mutate(
    N=1L,
    death_ltfu_pandemic_without_landmark = factor(death_ltfu_pandemic_without_landmark, 
                           labels = c("Alive and in care at pandemic start", "Died or LTFU between landmark and pandemic start"))
  ) %>%
  select(names(var_labels_death_ltfu2)
  ) %>%
  tbl_summary(
    by = death_ltfu_pandemic_without_landmark,
    label = var_labels_death_ltfu2[setdiff(names(.), "death_ltfu_pandemic_without_landmark")],
    statistic = list(
      N ~ "{N}",
      all_continuous() ~ "{median} ({p25}, {p75});  {mean} ({sd})"
    ),
  )

table_1_death_ltfu2 <- as_tibble(table_1_death_ltfu2$table_body)
table_1_death_ltfu2_redacted <- table_1_death_ltfu2 %>%
  mutate(
    across(
      starts_with("stat_"),
      ~ stringr::str_replace_all(
        ., "\\d+",
        function(x) as.character(fn_roundmid_any(as.numeric(x), threshold))
      )
    ),
    var_label = factor(
      var_label,
      levels = map_chr(var_labels_main[-c(1, 2)], ~ last(as.character(.)))
    ),
    label = replace_na(as.character(label), "")
  )
table_1_death_ltfu2_redacted <- table_1_death_ltfu2_redacted %>% 
  rename(metformin = stat_1,
         control = stat_2)

# Save output -------------------------------------------------------------
write_csv(table_1_main, fs::path("output", "data_description", "table1_main.csv"))
write_csv(table_1_main_redacted, fs::path("output", "data_description", "table1_main_midpoint6.csv"))
write_csv(table_1_death_ltfu1_redacted, fs::path("output", "data_description", "table1_death_ltfu1_midpoint6.csv"))
write_csv(table_1_death_ltfu2_redacted, fs::path("output", "data_description", "table1_death_ltfu2_midpoint6.csv"))
