################################################################################
# A custom made function to extract the data from a feather file
################################################################################
extract_data <- function(input_filename) {
  data_extract <- arrow::read_feather(here::here("output", input_filename))

  data_extract <- data_extract %>%
    mutate(across(c(contains("_date")),
                  ~ floor_date(as.Date(., format="%Y-%m-%d"), unit = "days")), # rounding down the date to the nearest day
           across(contains('_birth_year'),
                   ~ format(as.Date(.), "%Y")), # specifically for birth_year, then pass it onto _num to create cov_num_birth_year
           across(contains('_num'), ~ as.numeric(.)),
           across(contains('_cat'), ~ as.factor(.)),
           across(contains('_bin'), ~ as.logical(.))
           )
}


################ If import of a csv, then apply col_types instead #################
#     col_types = cols_only(
#
#       # Identifier
#       patient_id = col_integer(),
#
#       # QUALITY ASSURANCE ----
#       qa_bin_is_female_or_male = col_logical(),
#       qa_bin_was_adult = col_logical(),
#       qa_bin_was_alive = col_logical(),
#       qa_bin_known_imd = col_logical(),
#       qa_bin_was_registered = col_logical(),
#       qa_birth_year = col_integer(),
#       qa_date_of_death = col_date(format = "%Y-%m-%d"),
#       qa_bin_pregnancy = col_logical(),
#       qa_bin_cocp = col_logical(),
#       qa_bin_hrt = col_logical(),
#       qa_bin_prostate_cancer = col_logical(),
#
#       # POPULATION/DEMOGRAPHIC ----
#       cov_num_age = col_integer(),
#       cov_cat_sex = col_character(),
#       cov_cat_ethnicity = col_character(),
#       cov_cat_deprivation_5 = col_character(),
#       cov_cat_deprivation_10 = col_character(),
#       cov_cat_region = col_character(),
#       cov_cat_rural_urban = col_character(),
#       cov_cat_stp = col_character(),
#
#       # MAIN ELIGIBILITY 1 - FIRST POSITIVE SARS-CoV-2 TEST IN PERIOD ----
#       baseline_date = col_date(format = "%Y-%m-%d"),
#       cov_bin_pos_covid = col_logical(),
#
#       # MAIN ELIGIBILITY 2 - HISTORY OF T2DM ----
#       cov_date_t2dm = col_date(format = "%Y-%m-%d"), # will be modified through algorithm
#       anxcillary/temporary variables to feed diabetes algorithm and define T2DM
#.      cov_date_otherdm = col_date(format = "%Y-%m-%d"), # change to tmp?
#.      cov_date_gestationaldm = col_date(format = "%Y-%m-%d"), # change to tmp?
#       cov_date_poccdm = col_date(format = "%Y-%m-%d"), # change to tmp?
#       tmp_cov_date_t1dm_ctv3 = col_date(format = "%Y-%m-%d"),
#       tmp_cov_count_t1dm_ctv3 = col_integer(),
#       tmp_cov_count_t1dm_hes = col_integer(),
#       tmp_cov_count_t1dm = col_integer(),
#       tmp_cov_date_t2dm_ctv3 = col_date(format = "%Y-%m-%d"),
#       tmp_cov_count_t2dm_ctv3 = col_integer(),
#       tmp_cov_count_t2dm_hes = col_integer(),
#       tmp_cov_count_t2dm = col_integer(),
#       tmp_cov_count_otherdm = col_integer(),
#       tmp_cov_count_poccdm_ctv3 = col_integer(),
#       tmp_cov_num_max_hba1c_mmol_mol = col_numeric(),
#       tmp_cov_date_max_hba1c = col_date(format = "%Y-%m-%d"),
#       tmp_cov_date_insulin_snomed = col_date(format = "%Y-%m-%d"),
#       tmp_cov_date_antidiabetic_drugs_snomed = col_date(format = "%Y-%m-%d"),
#       tmp_cov_date_nonmetform_drugs_snomed = col_date(format = "%Y-%m-%d"),
#       tmp_cov_date_diabetes_medication = col_date(format = "%Y-%m-%d"),
#       tmp_cov_date_latest_diabetes_diag = col_date(format = "%Y-%m-%d"),
#
#       # OTHER ELIGIBILITY ----
#       cov_bin_prediabetes = col_logical(), # still to decide if eligibility  - otherwise confounder
#       cov_date_prediabetes = col_date(format = "%Y-%m-%d"), # still to decide if eligibility  - otherwise confounder
#       cov_bin_hosp_baseline = col_logical(),
#       cov_bin_metfin_before_baseline = col_logical(),
#       cov_date_metfin_before_baseline = col_date(format = "%Y-%m-%d"),
#       cov_bin_metfin_allergy = col_logical(),
#       cov_bin_ckd_45 = col_logical(),
#       cov_bin_liver_cirrhosis = col_logical(),
#       cov_bin_metfin_interaction = col_logical(),
#       cov_date_metfin_interaction = col_date(format = "%Y-%m-%d"),
#       cov_bin_long_covid = col_logical(),
#       cov_date_long_covid = col_date(format = "%Y-%m-%d"),
#
#       # TREATMENT - METFORMIN ----
#       exp_date_first_metfin = col_date(format = "%Y-%m-%d"),
#.      exp_count_metfin = col_integer(),
#       # date_treated = col_date(format = "%Y-%m-%d"), # we only have 1 treatment, variable not needed (unless adding other OADs?)
#       exp_bin_7d_metfin = col_logical(),
#
#       # CONFOUNDERS ----
#       cov_cat_smoking_status = col_character(),
#       cov_bin_carehome_status = col_logical(),
#       cov_bin_obesity = col_logical(),
#       cov_bin_ami = col_logical(),
#       cov_bin_all_stroke = col_logical(),
#       cov_bin_other_arterial_embolism = col_logical(),
#       cov_bin_vte = col_logical(),
#       cov_bin_hf = col_logical(),
#       cov_bin_angina = col_logical(),
#       cov_bin_dementia = col_logical(),
#       cov_bin_cancer = col_logical(),
#       cov_bin_hypertension = col_logical(),
#       cov_bin_depression = col_logical(),
#       cov_bin_copd = col_logical(),
#       cov_bin_liver_disease = col_logical(),
#       cov_bin_chronic_kidney_disease = col_logical(),
#       cov_bin_gestationaldm = col_logical(),
#       cov_bin_pcos = col_logical(),
#       cov_bin_t1dm = col_logical(),
#       cov_bin_diabetescomp = col_logical(),
#       cov_bin_hba1c_measurement = col_logical(),
#       cov_bin_ogtt_measurement = col_logical(),
#       cov_count_covid_vaccines = col_integer(),
#       cov_cat_bmi_groups = col_character(),
#       cov_num_bmi = col_numeric(),
#       cov_num_hba1c_mmol_mol = col_numeric(),
#       cov_bin_healthcare_worker = col_logical(),
#       cov_num_consultation_rate = col_integer(),
#       # cov_num_tc_hdl_ratio # derived later
#       tmp_cov_num_cholesterol = col_numeric(),
#       tmp_cov_num_hdl_cholesterol = col_numeric(),
#
#       # OUTCOMES ----
#       out_date_covid19 = col_date(format = "%Y-%m-%d"),
#       out_date_covid19_emergency = col_date(format = "%Y-%m-%d"),
#       out_bin_long_covid = col_logical(),
#       out_date_long_covid_first = col_date(format = "%Y-%m-%d"),
#       out_bin_viral_fatigue = col_logical(),
#       out_date_viral_fatigue_first = col_date(format = "%Y-%m-%d"),
#       out_bin_death_cause_covid = col_logical(),
#       out_date_death = col_date(format = "%Y-%m-%d"),
#
#       # VACCINATION ----
#       # vaccination_status = col_character(),
#       # date_most_recent_cov_vac = col_date(format = "%Y-%m-%d"),
#       # pfizer_most_recent_cov_vac = col_logical(),
#       # az_most_recent_cov_vac = col_logical(),
#       # moderna_most_recent_cov_vac = col_logical(),
#
#       # OUTCOMES ----
#       # death_date = col_date(format = "%Y-%m-%d"),
#       # died_ons_covid_any_date = col_date(format = "%Y-%m-%d"),
#       # died_ons_covid_date = col_date(format = "%Y-%m-%d"),
#       # dereg_date = col_date(format = "%Y-%m-%d"),
#       # # hosp
#       # # covid specific
#       # covid_hosp_admission_date0 = col_date(format = "%Y-%m-%d"),
#       # covid_hosp_admission_date1 = col_date(format = "%Y-%m-%d"),
#       # covid_hosp_admission_date2 = col_date(format = "%Y-%m-%d"),
#       # covid_hosp_admission_date3 = col_date(format = "%Y-%m-%d"),
#       # covid_hosp_admission_date4 = col_date(format = "%Y-%m-%d"),
#       # covid_hosp_admission_date5 = col_date(format = "%Y-%m-%d"),
#       # covid_hosp_admission_date6 = col_date(format = "%Y-%m-%d"),
#       # covid_hosp_admission_first_date7_28 = col_date(format = "%Y-%m-%d"),
#       # covid_hosp_discharge_date = col_date(format = "%Y-%m-%d"),
#       # covid_hosp_date_mabs_procedure = col_date(format = "%Y-%m-%d"),
#       # # all cause
#       # allcause_hosp_admission_date0 = col_date(format = "%Y-%m-%d"),
#       # allcause_hosp_admission_date1 = col_date(format = "%Y-%m-%d"),
#       # allcause_hosp_admission_date2 = col_date(format = "%Y-%m-%d"),
#       # allcause_hosp_admission_date3 = col_date(format = "%Y-%m-%d"),
#       # allcause_hosp_admission_date4 = col_date(format = "%Y-%m-%d"),
#       # allcause_hosp_admission_date5 = col_date(format = "%Y-%m-%d"),
#       # allcause_hosp_admission_date6 = col_date(format = "%Y-%m-%d"),
#       # allcause_hosp_admission_first_date7_28 = col_date(format = "%Y-%m-%d"),
#       # allcause_hosp_discharge_date = col_date(format = "%Y-%m-%d"),
#       # allcause_hosp_date_mabs_procedure = col_date(format = "%Y-%m-%d"),
#       # # cause/diagnosis
#       # death_cause = col_character(),
#       # allcause_hosp_admission_diagnosis0 = col_character(),
#       # allcause_hosp_admission_diagnosis1 = col_character(),
#       # allcause_hosp_admission_diagnosis2 = col_character(),
#       # allcause_hosp_admission_diagnosis4 = col_character(),
#       # allcause_hosp_admission_diagnosis3 = col_character(),
#       # allcause_hosp_admission_diagnosis5 = col_character(),
#       # allcause_hosp_admission_diagnosis6 = col_character(),
#       # allcause_hosp_admission_first_diagnosis7_28 = col_character()
#     ),
#   )
# }
