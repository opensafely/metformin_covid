#######################################################################################
# IMPORT
#######################################################################################
## ehrQL functions
from ehrql import (
    case,
    create_dataset,
    when,
    minimum_of,
    maximum_of,
    days
)

## TPP tables
from ehrql.tables.tpp import (
    clinical_events,
    medications,
    patients,
    practice_registrations,
    ons_deaths,
    sgss_covid_all_tests,
    addresses,
    occupation_on_covid_vaccine_record
)

## All codelists from codelists.py
from codelists import *

## variable helper functions 
from variable_helper_functions import *

## json (for the dates)
import json

# random seed (ideally use numpy, but currently not working on my local environment)
#import numpy as np 
#np.random.seed(19283) # random seed
import random
random.seed(19283) # random seed

#######################################################################################
# DEFINE the feasibility study end date
#######################################################################################
with open("output/study_dates.json") as f:
  study_dates = json.load(f)
studyend_date = study_dates["studyend_date"]
pandemicstart_date = study_dates["pandemicstart_date"]
mid2018_date = study_dates["mid2018_date"]

#######################################################################################
# INITIALISE the dataset and set the dummy dataset size
#######################################################################################
dataset = create_dataset()
dataset.configure_dummy_data(population_size=5000)
dataset.define_population(patients.exists_for_patient())

#######################################################################################
# Table 2) QUALITY ASSURANCES and completeness criteria
#######################################################################################
# population variables for dataset definition 
dataset.qa_bin_is_female_or_male = patients.sex.is_in(["female", "male"]) 
dataset.qa_bin_was_adult = (patients.age_on(pandemicstart_date) >= 18) & (patients.age_on(pandemicstart_date) <= 110) 
dataset.qa_bin_was_alive = patients.is_alive_on(pandemicstart_date)
dataset.qa_bin_was_alive_mid2018 = patients.is_alive_on(mid2018_date)
dataset.qa_bin_known_imd = addresses.for_patient_on(pandemicstart_date).exists_for_patient() # known deprivation
dataset.qa_bin_was_registered = practice_registrations.spanning(pandemicstart_date - days(366), pandemicstart_date).exists_for_patient() # see https://docs.opensafely.org/ehrql/reference/schemas/tpp/#practice_registrations.spanning. Calculated from 1 year = 365.25 days, taking into account leap year.
dataset.qa_bin_was_registered_mid2018 = practice_registrations.spanning(mid2018_date - days(366), mid2018_date).exists_for_patient() 

## Date of death
dataset.qa_date_of_death = ons_deaths.date

## Pregnancy (over entire study period -> don't have helper function for this)
dataset.qa_bin_pregnancy = clinical_events.where(clinical_events.snomedct_code.is_in(pregnancy_snomed_clinical)).exists_for_patient()

## Combined oral contraceptive pill (over entire study period -> don't have helper function for this)
dataset.qa_bin_cocp = medications.where(medications.dmd_code.is_in(cocp_dmd)).exists_for_patient()

## Hormone replacement therapy (over entire study period -> don't have helper function for this)
dataset.qa_bin_hrt = medications.where(medications.dmd_code.is_in(hrt_dmd)).exists_for_patient()

## Prostate cancer (over entire study period -> don't have helper function for this)
### Primary care
prostate_cancer_snomed = clinical_events.where(clinical_events.snomedct_code.is_in(prostate_cancer_snomed_clinical)).exists_for_patient()

### HES APC
prostate_cancer_hes = apcs.where(apcs.all_diagnoses.contains_any_of(prostate_cancer_icd10)).exists_for_patient()

### ONS (stated anywhere on death certificate)
prostate_cancer_death = cause_of_death_matches(prostate_cancer_icd10)
# Combined: Any prostate cancer diagnosis
dataset.qa_bin_prostate_cancer = case(
    when(prostate_cancer_snomed).then(True),
    when(prostate_cancer_hes).then(True),
    when(prostate_cancer_death).then(True),
    otherwise=False
)


#######################################################################################
# Table 3) ELIGIBILITY criteria
#######################################################################################

###
# diabetes variables defined in previous separate action/dataset definition
###

## Known hypersensitivity / intolerance to metformin, on or before landmark
dataset.elig_date_metfin_allergy_last = last_matching_event_clinical_snomed_before(metformin_allergy_snomed_clinical, pandemicstart_date).date
dataset.elig_date_metfin_allergy_first = first_matching_event_clinical_snomed_before(metformin_allergy_snomed_clinical, pandemicstart_date).date

## Moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; stage 4/5), on or before landmark
dataset.elig_date_ckd_45_last = maximum_of(
    last_matching_event_clinical_snomed_before(ckd_snomed_clinical_45, pandemicstart_date).date,
    last_matching_event_apc_before(ckd_stage4_icd10 + ckd_stage5_icd10, pandemicstart_date).admission_date
)
dataset.elig_date_ckd_45_first = minimum_of(
    first_matching_event_clinical_snomed_before(ckd_snomed_clinical_45, pandemicstart_date).date,
    first_matching_event_apc_before(ckd_stage4_icd10 + ckd_stage5_icd10, pandemicstart_date).admission_date
)

## Advance decompensated liver cirrhosis, on or before landmark
dataset.elig_date_liver_cirrhosis_last = maximum_of(
    last_matching_event_clinical_snomed_before(advanced_decompensated_cirrhosis_snomed_clinical + ascitic_drainage_snomed_clinical, pandemicstart_date).date,
    last_matching_event_apc_before(advanced_decompensated_cirrhosis_icd10, pandemicstart_date).admission_date
)
dataset.elig_date_liver_cirrhosis_first = minimum_of(
    first_matching_event_clinical_snomed_before(advanced_decompensated_cirrhosis_snomed_clinical + ascitic_drainage_snomed_clinical, pandemicstart_date).date,
    first_matching_event_apc_before(advanced_decompensated_cirrhosis_icd10, pandemicstart_date).admission_date
)

## Use of the following medications in the last 14 days (drug-drug interaction with metformin)
dataset.elig_date_metfin_interaction_last = last_matching_med_dmd_before(metformin_interaction_dmd, pandemicstart_date).date
dataset.elig_date_metfin_interaction_first = first_matching_med_dmd_before(metformin_interaction_dmd, pandemicstart_date).date

#######################################################################################
# Table 4) INTERVENTION/EXPOSURE variables
#######################################################################################
## Last metformin prescription before landmark: Use to establish "prescription xxx months prior to landmark"
dataset.exp_date_metfin_last = last_matching_med_dmd_before(metformin_dmd, pandemicstart_date).date
dataset.exp_date_metfin_mono_last = last_matching_med_dmd_before(metformin_mono_dmd, pandemicstart_date).date 
## First metformin prescription before landmark: Use to establish initiation of metformin
dataset.exp_date_metfin_first = first_matching_med_dmd_before(metformin_dmd, pandemicstart_date).date
dataset.exp_date_metfin_mono_first = first_matching_med_dmd_before(metformin_mono_dmd, pandemicstart_date).date

## Any other antidiabetic drug exposure before pandemicstart_date
dataset.exp_date_sulfo_first = first_matching_med_dmd_before(sulfonylurea_dmd, pandemicstart_date).date 
dataset.exp_date_dpp4_first = first_matching_med_dmd_before(dpp4_dmd, pandemicstart_date).date 
dataset.exp_date_dpp4_mono_first = first_matching_med_dmd_before(dpp4_mono_dmd, pandemicstart_date).date 
dataset.exp_date_tzd_first = first_matching_med_dmd_before(tzd_dmd, pandemicstart_date).date 
dataset.exp_date_tzd_mono_first = first_matching_med_dmd_before(tzd_mono_dmd, pandemicstart_date).date 
dataset.exp_date_sglt2_first = first_matching_med_dmd_before(sglt2_dmd, pandemicstart_date).date 
dataset.exp_date_sglt2_mono_first = first_matching_med_dmd_before(sglt2_mono_dmd, pandemicstart_date).date 
dataset.exp_date_glp1_first = first_matching_med_dmd_before(glp1_dmd, pandemicstart_date).date 
dataset.exp_date_megli_first = first_matching_med_dmd_before(meglitinides_dmd, pandemicstart_date).date 
dataset.exp_date_agi_first = first_matching_med_dmd_before(agi_dmd, pandemicstart_date).date 
dataset.exp_date_insulin_first = first_matching_med_dmd_before(insulin_dmd, pandemicstart_date).date


#######################################################################################
# Table 5) Demographics, covariates and potential confounders
#######################################################################################

## Sex
dataset.cov_cat_sex = patients.sex

## Age at pandemicstart_date
dataset.cov_num_age = patients.age_on(pandemicstart_date)

## Index of Multiple Deprevation Rank (rounded down to nearest 100). 5 categories.
imd_rounded = addresses.for_patient_on(pandemicstart_date).imd_rounded
dataset.cov_cat_deprivation_5 = case(
    when((imd_rounded >=0) & (imd_rounded < int(32844 * 1 / 5))).then("1 (most deprived)"),
    when(imd_rounded < int(32844 * 2 / 5)).then("2"),
    when(imd_rounded < int(32844 * 3 / 5)).then("3"),
    when(imd_rounded < int(32844 * 4 / 5)).then("4"),
    when(imd_rounded < int(32844 * 5 / 5)).then("5 (least deprived)"),
    otherwise="unknown"
)

## Practice registration info at pandemicstart_date
registered = practice_registrations.for_patient_on(pandemicstart_date)
registered_mid2018 = practice_registrations.for_patient_on(mid2018_date)
dataset.cov_cat_region = registered.practice_nuts1_region_name ## Region
dataset.cov_cat_region_mid2018 = registered_mid2018.practice_nuts1_region_name ## Region
dataset.cov_cat_stp = registered.practice_stp ## Practice
dataset.cov_cat_stp_mid2018 = registered_mid2018.practice_stp ## Practice
dataset.cov_cat_rural_urban = addresses.for_patient_on(pandemicstart_date).rural_urban_classification ## Rurality
dataset.cov_cat_rural_urban_mid2018 = addresses.for_patient_on(mid2018_date).rural_urban_classification ## Rurality

## Smoking status at pandemicstart_date
tmp_most_recent_smoking_cat = last_matching_event_clinical_ctv3_before(smoking_clear, pandemicstart_date).ctv3_code.to_category(smoking_clear)
tmp_ever_smoked = last_matching_event_clinical_ctv3_before(ever_smoking, pandemicstart_date).exists_for_patient() # uses a different codelist with ONLY smoking codes
dataset.cov_cat_smoking_status = case(
    when(tmp_most_recent_smoking_cat == "S").then("S"),
    when(tmp_most_recent_smoking_cat == "E").then("E"),
    when((tmp_most_recent_smoking_cat == "N") & (tmp_ever_smoked == True)).then("E"),
    when(tmp_most_recent_smoking_cat == "N").then("N"),
    when((tmp_most_recent_smoking_cat == "M") & (tmp_ever_smoked == True)).then("E"),
    when(tmp_most_recent_smoking_cat == "M").then("M"),
    otherwise = "M"
)

## Care home resident at pandemicstart_date, see https://github.com/opensafely/opioids-covid-research/blob/main/analysis/define_dataset_table.py
# Flag care home based on primis (patients in long-stay nursing and residential care)
tmp_care_home_code = last_matching_event_clinical_snomed_before(carehome, pandemicstart_date).exists_for_patient()
# Flag care home based on TPP
tmp_care_home_tpp1 = addresses.for_patient_on(pandemicstart_date).care_home_is_potential_match
tmp_care_home_tpp2 = addresses.for_patient_on(pandemicstart_date).care_home_requires_nursing
tmp_care_home_tpp3 = addresses.for_patient_on(pandemicstart_date).care_home_does_not_require_nursing
# combine
dataset.cov_bin_carehome_status = case(
    when(tmp_care_home_code).then(True),
    when(tmp_care_home_tpp1).then(True),
    when(tmp_care_home_tpp2).then(True),
    when(tmp_care_home_tpp3).then(True),
    otherwise = False
)

## Any other antidiabetic drug use before pandemicstart_date (could also be a combo with metformin)
dataset.cov_date_sulfo_last = last_matching_med_dmd_before(sulfonylurea_dmd, pandemicstart_date).date 
dataset.cov_date_dpp4_last = last_matching_med_dmd_before(dpp4_dmd, pandemicstart_date).date 
#dataset.cov_date_dpp4_mono_last = last_matching_med_dmd_before(dpp4_mono_dmd, pandemicstart_date).date 
dataset.cov_date_tzd_last = last_matching_med_dmd_before(tzd_dmd, pandemicstart_date).date 
#dataset.cov_date_tzd_mono_last = last_matching_med_dmd_before(tzd_mono_dmd, pandemicstart_date).date 
dataset.cov_date_sglt2_last = last_matching_med_dmd_before(sglt2_dmd, pandemicstart_date).date 
#dataset.cov_date_sglt2_mono_last = last_matching_med_dmd_before(sglt2_mono_dmd, pandemicstart_date).date 
dataset.cov_date_glp1_last = last_matching_med_dmd_before(glp1_dmd, pandemicstart_date).date 
dataset.cov_date_megli_last = last_matching_med_dmd_before(meglitinides_dmd, pandemicstart_date).date 
dataset.cov_date_agi_last = last_matching_med_dmd_before(agi_dmd, pandemicstart_date).date 
dataset.cov_date_insulin_last = last_matching_med_dmd_before(insulin_dmd, pandemicstart_date).date

## Obesity, on or before pandemicstart_date
dataset.cov_bin_obesity = (
    last_matching_event_clinical_snomed_before(bmi_obesity_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(bmi_obesity_icd10, pandemicstart_date).exists_for_patient()
)

## Acute myocardial infarction, on or before pandemicstart_date
dataset.cov_bin_ami = (
    last_matching_event_clinical_snomed_before(ami_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(ami_prior_icd10 + ami_icd10, pandemicstart_date).exists_for_patient()
)

## All stroke, on or before pandemicstart_date
dataset.cov_bin_all_stroke = (
    last_matching_event_clinical_snomed_before(stroke_isch_snomed_clinical + stroke_sah_hs_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(stroke_isch_icd10 + stroke_sah_hs_icd10, pandemicstart_date).exists_for_patient()
)

## Other arterial embolism, on or before pandemicstart_date
dataset.cov_bin_other_arterial_embolism = (
    last_matching_event_clinical_snomed_before(other_arterial_embolism_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(other_arterial_embolism_icd10, pandemicstart_date).exists_for_patient()
)

## Venous thrombolism events, on or before pandemicstart_date
dataset.cov_bin_vte = (
  last_matching_event_clinical_snomed_before(portal_vein_thrombosis_snomed_clinical + dvt_dvt_snomed_clinical + dvt_icvt_snomed_clinical + dvt_pregnancy_snomed_clinical + other_dvt_snomed_clinical + pe_snomed_clinical, pandemicstart_date).exists_for_patient() |
  last_matching_event_apc_before(portal_vein_thrombosis_icd10 + dvt_dvt_icd10 + dvt_icvt_icd10 + dvt_pregnancy_icd10 + other_dvt_icd10 + icvt_pregnancy_icd10 + pe_icd10, pandemicstart_date).exists_for_patient()
)

## Heart failure, on or before pandemicstart_date
dataset.cov_bin_hf = (
    last_matching_event_clinical_snomed_before(hf_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(hf_icd10, pandemicstart_date).exists_for_patient()
)

## Angina, on or before pandemicstart_date
dataset.cov_bin_angina = (
    last_matching_event_clinical_snomed_before(angina_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(angina_icd10, pandemicstart_date).exists_for_patient()
)

## Dementia, on or before pandemicstart_date
dataset.cov_bin_dementia = (
    last_matching_event_clinical_snomed_before(dementia_snomed_clinical + dementia_vascular_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(dementia_icd10 + dementia_vascular_icd10, pandemicstart_date).exists_for_patient()
)

## Cancer, on or before pandemicstart_date
dataset.cov_bin_cancer = (
    last_matching_event_clinical_snomed_before(cancer_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(cancer_icd10, pandemicstart_date).exists_for_patient()
)

## Hypertension, on or before pandemicstart_date
dataset.cov_bin_hypertension = (
    last_matching_event_clinical_snomed_before(hypertension_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(hypertension_icd10, pandemicstart_date).exists_for_patient()
)

## Depression, on or before pandemicstart_date
dataset.cov_bin_depression = (
    last_matching_event_clinical_snomed_before(depression_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(depression_icd10, pandemicstart_date).exists_for_patient()
)

## Chronic obstructive pulmonary disease, on or before pandemicstart_date
dataset.cov_bin_copd = (
    last_matching_event_clinical_snomed_before(copd_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(copd_icd10, pandemicstart_date).exists_for_patient()
)

## Liver disease, on or before pandemicstart_date
dataset.cov_bin_liver_disease = (
    last_matching_event_clinical_snomed_before(liver_disease_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(liver_disease_icd10, pandemicstart_date).exists_for_patient()
)

## Chronic kidney disease, on or before pandemicstart_date
dataset.cov_bin_chronic_kidney_disease = (
    last_matching_event_clinical_snomed_before(ckd_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(ckd_icd10, pandemicstart_date).exists_for_patient()
)

## PCOS, on or before pandemicstart_date
dataset.cov_bin_pcos = (
    last_matching_event_clinical_snomed_before(pcos_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(pcos_icd10, pandemicstart_date).exists_for_patient()
)

## Prediabetes, on or before pandemicstart_date
# Any preDM diagnosis in primary care
tmp_cov_bin_prediabetes = last_matching_event_clinical_snomed_before(prediabetes_snomed, pandemicstart_date).exists_for_patient()
# Any HbA1c preDM in primary care
tmp_cov_bin_predm_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(pandemicstart_date))
    .where((clinical_events.numeric_value>=42) & (clinical_events.numeric_value<=47.9))
    .exists_for_patient()
)
# Any preDM diagnosis or Hb1Ac preDM range value (in period before pandemicstart_date)
dataset.cov_bin_prediabetes = tmp_cov_bin_prediabetes | tmp_cov_bin_predm_hba1c_mmol_mol

## Diabetes complications (foot, retino, neuro, nephro), on or before pandemicstart_date
dataset.cov_bin_diabetescomp = (
    last_matching_event_clinical_snomed_before(diabetescomp_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(diabetescomp_icd10, pandemicstart_date).exists_for_patient()
)

## Any HbA1c measurement, on or before pandemicstart_date ## does not make sense for landmark set up (stick to HbA1c level)
#dataset.cov_bin_hba1c_measurement = last_matching_event_clinical_snomed_before(hba1c_measurement_snomed, pandemicstart_date).exists_for_patient()

## Any OGTT done, on or before pandemicstart_date ## does not make sense for landmark set up (stick to HbA1c level)
#dataset.cov_bin_ogtt_measurement = last_matching_event_clinical_snomed_before(ogtt_measurement_snomed, pandemicstart_date).exists_for_patient()

## BMI, most recent value, within previous 2 years, on or before pandemicstart_date
bmi_measurement = most_recent_bmi(
    where=clinical_events.date.is_on_or_between(pandemicstart_date - days(2 * 366), pandemicstart_date),
    minimum_age_at_measurement=16,
)
dataset.cov_num_bmi = bmi_measurement.numeric_value
dataset.cov_cat_bmi_groups = case(
    when(dataset.cov_num_bmi < 18.5).then("Underweight"),
    when((dataset.cov_num_bmi >= 18.5) & (dataset.cov_num_bmi < 25.0)).then("Healthy weight (18.5-24.9)"),
    when((dataset.cov_num_bmi >= 25.0) & (dataset.cov_num_bmi < 30.0)).then("Overweight (25-29.9)"),
    when((dataset.cov_num_bmi >= 30.0) & (dataset.cov_num_bmi < 70.0)).then("Obese (>30)"), # Set maximum to avoid any impossibly extreme values being classified as obese
    otherwise = "missing", 
)

## HbA1c, most recent value, within previous 2 years, on or before pandemicstart_date
dataset.cov_num_hba1c_mmol_mol = last_matching_event_clinical_snomed_between(hba1c_snomed, pandemicstart_date - days(2*366), pandemicstart_date).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years. 

## Total Cholesterol, most recent value, within previous 2 years, on or before pandemicstart_date
dataset.tmp_cov_num_cholesterol = last_matching_event_clinical_snomed_between(cholesterol_snomed, pandemicstart_date - days(2*366), pandemicstart_date).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years. 

## HDL Cholesterol, most recent value, within previous 2 years, on or before pandemicstart_date
dataset.tmp_cov_num_hdl_cholesterol = last_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, pandemicstart_date - days(2*366), pandemicstart_date).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years.

## Healthcare worker at the time they received a COVID-19 vaccination
dataset.cov_bin_healthcare_worker = (
  occupation_on_covid_vaccine_record.where(
    occupation_on_covid_vaccine_record.is_healthcare_worker == True)
    .exists_for_patient()
)


#######################################################################################
# Table 6) Outcomes
#######################################################################################
# for descriptive purpose and for censoring reasons, add metformin initiation as an "outcome" AFTER landmark, as well as all other antidiabetics options
dataset.out_date_metfin_first = first_matching_med_dmd_between(metformin_dmd, pandemicstart_date, studyend_date).date
dataset.out_date_metfin_mono_first = first_matching_med_dmd_between(metformin_mono_dmd, pandemicstart_date, studyend_date).date

## Any other antidiabetic drug use after pandemicstart_date
dataset.out_date_sulfo_first = first_matching_med_dmd_between(sulfonylurea_dmd, pandemicstart_date, studyend_date).date 
dataset.out_date_dpp4_first = first_matching_med_dmd_between(dpp4_dmd, pandemicstart_date, studyend_date).date 
# dataset.out_date_dpp4_mono_first = first_matching_med_dmd_between(dpp4_mono_dmd, pandemicstart_date, studyend_date).date 
dataset.out_date_tzd_first = first_matching_med_dmd_between(tzd_dmd, pandemicstart_date, studyend_date).date 
# dataset.out_date_tzd_mono_first = first_matching_med_dmd_between(tzd_mono_dmd, pandemicstart_date, studyend_date).date 
dataset.out_date_sglt2_first = first_matching_med_dmd_between(sglt2_dmd, pandemicstart_date, studyend_date).date 
# dataset.out_date_sglt2_mono_first = first_matching_med_dmd_between(sglt2_mono_dmd, pandemicstart_date, studyend_date).date 
dataset.out_date_glp1_first = first_matching_med_dmd_between(glp1_dmd, pandemicstart_date, studyend_date).date 
dataset.out_date_megli_first = first_matching_med_dmd_between(meglitinides_dmd, pandemicstart_date, studyend_date).date 
dataset.out_date_agi_first = first_matching_med_dmd_between(agi_dmd, pandemicstart_date, studyend_date).date 
dataset.out_date_insulin_first = first_matching_med_dmd_between(insulin_dmd, pandemicstart_date, studyend_date).date

## Practice deregistration date 1: Based on registration at pandemicstart_date
# dataset.out_date_dereg_pandemicstart_date = registered.end_date

## Practice deregistration date 2: Based on registration in mid2018
dataset.out_date_dereg_mid2018 = registered_mid2018.end_date

## Practice deregistration date 3: Any dereg after mid2018 / but this does not take into account those who were registered at baseline - the above one is cleaner and for sure only takes into account the baseline registration
#dataset.out_date_dereg_any = (
#  practice_registrations.where(practice_registrations.end_date.is_on_or_between(mid2018_date, studyend_date))
#  .sort_by(practice_registrations.end_date)
#  .first_for_patient()
#  .end_date
#)

## First covid-19 related hospital admission, between pandemicstart_date and studyend_date (incl. those dates)
dataset.out_date_covid19_hes = first_matching_event_apc_between(covid_codes_incl_clin_diag, pandemicstart_date, studyend_date, only_prim_diagnoses=True).admission_date
## First covid-19 related emergency attendance, between pandemicstart_date and studyend_date (incl. those dates)
dataset.out_date_covid19_emergency = first_matching_event_ec_snomed_between(covid_emergency, pandemicstart_date, studyend_date).arrival_date
# combined: First covid-19 related hospitalization
dataset.out_date_covid19_hosp = minimum_of(dataset.out_date_covid19_hes, dataset.out_date_covid19_emergency)

## First COVID-19 diagnosis in primary care, between pandemicstart_date and studyend_date (incl. those dates)
tmp_covid19_primary_care_date = first_matching_event_clinical_ctv3_between(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae, pandemicstart_date, studyend_date).date
## First positive SARS-COV-2 PCR in primary care, between pandemicstart_date and studyend_date (incl. those dates)
tmp_covid19_sgss_date = (
  sgss_covid_all_tests.where(sgss_covid_all_tests.specimen_taken_date.is_on_or_between(pandemicstart_date, studyend_date))
  .where(sgss_covid_all_tests.is_positive)
  .sort_by(sgss_covid_all_tests.specimen_taken_date)
  .first_for_patient()
  .specimen_taken_date
)
## First COVID-19 diagnosis in primary care, or pos. test in primary care, or covid-19 hosp, between pandemicstart_date and studyend_date (incl. those dates)
dataset.out_date_covid19 = minimum_of(tmp_covid19_primary_care_date, tmp_covid19_sgss_date, dataset.out_date_covid19_hosp)

## COVID-related death, between pandemicstart_date and studyend_date (incl. those dates), stated anywhere on any of the 15 death certificate options
tmp_out_bin_death_covid = matching_death_between(covid_codes_incl_clin_diag, pandemicstart_date, studyend_date)
dataset.out_date_covid19_death = case(when(tmp_out_bin_death_covid).then(ons_deaths.date))
# combined: First covid-19 related hospitalization or death after pandemicstart_date, i.e. "severe COVID" outcome = secondary outcome
dataset.out_date_covid19_severe = minimum_of(dataset.out_date_covid19_death, dataset.out_date_covid19_hosp)

## First Long COVID code in primary care, between pandemicstart_date and studyend_date (incl. those dates), based on https://github.com/opensafely/long-covid/blob/main/analysis/codelists.py
#dataset.out_date_long_covid = first_matching_event_clinical_snomed_between(long_covid_diagnostic_snomed_clinical + long_covid_referral_snomed_clinical + long_covid_assessment_snomed_clinical, pandemicstart_date, studyend_date).date
## First Viral fatigue code in primary care, between pandemicstart_date and studyend_date (incl. those dates), based on https://github.com/opensafely/long-covid/blob/main/analysis/codelists.py
#dataset.out_date_viral_fatigue = first_matching_event_clinical_snomed_between(post_viral_fatigue_snomed_clinical, pandemicstart_date, studyend_date).date
# combined: First Long COVID code or Viral Fatigue code after pandemicstart_date = primary outcome
#dataset.out_date_long_fatigue = minimum_of(dataset.out_date_long_covid, dataset.out_date_viral_fatigue)

## Long COVID signs and symptoms, first code in primary care, between pandemicstart_date and studyend_date (incl. those dates)
#dataset.out_date_breathlessness = first_matching_event_clinical_snomed_between(breathlessness_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_cough = first_matching_event_clinical_snomed_between(cough_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_chest_pain = first_matching_event_clinical_snomed_between(chest_pain_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_chest_tightness = first_matching_event_clinical_snomed_between(chest_tightness_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_palpitations = first_matching_event_clinical_snomed_between(palpitations_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_fatigue = first_matching_event_clinical_snomed_between(fatigue_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_fever = first_matching_event_clinical_snomed_between(fever_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_pain = first_matching_event_clinical_snomed_between(pain_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_cog_impair = first_matching_event_clinical_snomed_between(cog_impair_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_headache = first_matching_event_clinical_snomed_between(headache_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_sleep = first_matching_event_clinical_snomed_between(sleep_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_pnp = first_matching_event_clinical_snomed_between(pnp_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_dizziness = first_matching_event_clinical_snomed_between(dizziness_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_delirium = first_matching_event_clinical_snomed_between(delirium_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_mob_impair = first_matching_event_clinical_snomed_between(mob_impair_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_visual = first_matching_event_clinical_snomed_between(visual_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_abdo_pain = first_matching_event_clinical_snomed_between(abdo_pain_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_nausea_vomiting = first_matching_event_clinical_snomed_between(nausea_vomiting_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_diarrhoea = first_matching_event_clinical_snomed_between(diarrhoea_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_weight_appetite = first_matching_event_clinical_snomed_between(weight_appetite_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_tinnitus = first_matching_event_clinical_snomed_between(tinnitus_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_earache = first_matching_event_clinical_snomed_between(earache_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_sore_throat = first_matching_event_clinical_snomed_between(sore_throat_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_smell_taste = first_matching_event_clinical_snomed_between(smell_taste_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_nasal = first_matching_event_clinical_snomed_between(nasal_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_hair_loss = first_matching_event_clinical_snomed_between(hair_loss_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_skin_rash = first_matching_event_clinical_snomed_between(skin_rash_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_anxiety = first_matching_event_clinical_snomed_between(anxiety_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_depression = first_matching_event_clinical_snomed_between(depression_snomed, pandemicstart_date, studyend_date).date
#dataset.out_date_ptsd = first_matching_event_clinical_snomed_between(ptsd_snomed, pandemicstart_date, studyend_date).date