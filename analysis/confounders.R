####
## This script defines the final covariates/confounders included in the model
####
confounder_names <- c("cov_cat_sex", "cov_num_age_spline",
                       
                       "cov_cat_ethnicity", "cov_cat_deprivation_5", "cov_cat_rural_urban",
                       "cov_cat_smoking_status", "cov_bin_healthcare_worker", "cov_num_consrate", 
                       
                       "cov_bin_ami", "cov_bin_all_stroke", "cov_bin_other_arterial_embolism", "cov_bin_vte", 
                       "cov_bin_hf", "cov_bin_angina", "cov_bin_dementia", "cov_bin_cancer", "cov_bin_hypertension",
                       "cov_bin_depression", "cov_bin_copd", "cov_bin_liver_disease", "cov_bin_chronic_kidney_disease",
                       "cov_bin_pcos", "cov_bin_prediabetes", "cov_bin_diabetescomp", 
                       
                       "cov_cat_bmi_groups", "cov_cat_tc_hdl_ratio_b", "cov_cat_hba1c_b", "cov_num_period_month") 