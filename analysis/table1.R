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

# Label the data ---------------------------------------------------------
var_labels <- list(
  N  ~ "Total N",
  exp_bin_treat ~ "Treatment",
  
  cov_num_age ~ "Age",
  cov_cat_age ~ "Age groups",
  cov_cat_sex ~ "Sex",
  cov_cat_ethnicity ~ "Ethnicity",
  cov_cat_deprivation_5 ~ "Deprivation",
  strat_cat_region ~ "Region",
  cov_cat_rural_urban ~ "Rural/urban",
  cov_cat_smoking_status ~ "Smoking status",
  cov_bin_carehome_status ~ "Care/nursing home resident",
  cov_bin_healthcare_worker ~ "Healthcare worker",
  cov_num_consrate ~ "Consultation rate in previous year",
  cov_bin_obesity ~ "Body Mass Index > 30 kg/m^2",
  cov_cat_hba1c_mmol_mol ~ "HbA1c categories in mmol/mol",
  cov_cat_tc_hdl_ratio ~ "TC/HDL ratio categories",
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
  
  out_bin_severecovid_afterlandmark ~ "COVID hosp or death",
  out_bin_covid_hosp_afterlandmark ~ "COVID hosp",
  out_bin_covid_death_afterlandmark ~ "COVID death",
  out_bin_covid_afterlandmark ~ "Any covid diagnosis, pos test or hosp",
  out_bin_longcovid_afterlandmark ~ "Any Long COVID diagnosis",
  out_bin_virfat_afterlandmark ~ "Any Viral Fatigue diagnosis",
  out_bin_longcovid_virfat_afterlandmark ~ "Any Long COVID or Viral Fatigue diagnosis",
  out_bin_death_afterlandmark ~ "Any death after landmark",
  cens_bin_ltfu_afterlandmark ~ "Any LTFU after landmark",
  cens_bin_metfin_pandemicstart ~ "INT: Any metformin prescription within 6m prior to pandemic start",
  cens_bin_metfin_start_cont ~ "CONT: Any metformin start in control"
) %>%
  set_names(., map_chr(., all.vars))

# Create the table ---------------------------------------------------------
table_1 <-
  df %>%
  mutate(
    N=1L,
    exp_bin_treat = factor(exp_bin_treat, 
                          levels = c(1,0), 
                          labels = c("Metformin mono", "Nothing"))
  ) %>%
  select(
    exp_bin_treat,
    all_of(names(var_labels)),
  ) %>%
  tbl_summary(
    by = exp_bin_treat,
    label = unname(var_labels[names(.)]),
    statistic = list(
      N ~ "{N}",
      all_continuous() ~ "{median} ({p25}, {p75});  {mean} ({sd})"
    ),
  )

raw_stats <- table_1$meta_data %>%
  select(var_label, df_stats) %>%
  unnest(df_stats)

raw_stats_redacted <- raw_stats %>%
  mutate(
    n = fn_roundmid_any(n, threshold),
    N = fn_roundmid_any(N, threshold),
    p = n / N,
    N_miss = fn_roundmid_any(N_miss, threshold),
    N_obs = fn_roundmid_any(N_obs, threshold),
    p_miss = N_miss / N_obs,
    N_nonmiss = fn_roundmid_any(N_nonmiss, threshold),
    p_nonmiss = N_nonmiss / N_obs,
    var_label = factor(var_label, levels = map_chr(var_labels[-c(1, 2)], ~ last(as.character(.)))),
    variable_levels = replace_na(as.character(variable_levels), "")
  )


# Save output -------------------------------------------------------------
write_csv(raw_stats_redacted, fs::path("output", "data_description", "table1_midpoint6.csv"))
