#######################################################################################
# IMPORT
#######################################################################################
## ehrQL functions
from ehrql import (
    case,
    create_dataset,
    when,
    minimum_of,
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

# numpy for random seed - and set random seed
#import numpy as np 
#np.random.seed(1928374) # random seed


#######################################################################################
# DEFINE the feasibility study end date
#######################################################################################
with open("output/study_dates.json") as f:
  study_dates = json.load(f)
studyend_date = study_dates["studyend_date"]
landmark_date = study_dates["landmark_date"]


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
dataset.qa_bin_was_adult = (patients.age_on(landmark_date) >= 18) & (patients.age_on(landmark_date) <= 110) 
dataset.qa_bin_was_alive = (((patients.date_of_death.is_null()) | (patients.date_of_death.is_after(landmark_date))) & 
        ((ons_deaths.date.is_null()) | (ons_deaths.date.is_after(landmark_date))))
dataset.qa_bin_known_imd = addresses.for_patient_on(landmark_date).exists_for_patient() # known deprivation
dataset.qa_bin_was_registered = practice_registrations.spanning(landmark_date - days(914), landmark_date).exists_for_patient() # only include if registered on landmark_date spanning back 30 months (i.e. only deregistered later than landmark_date or has no deregistration date, see https://docs.opensafely.org/ehrql/reference/schemas/tpp/#practice_registrations.spanning). Calculated from 1 year = 365.25 days, taking into account leap years.

## Year of birth
dataset.qa_num_birth_year = patients.date_of_birth

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
prostate_cancer_hes = apcs.where(apcs.all_diagnoses.is_in(prostate_cancer_icd10)).exists_for_patient()

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

## DIABETES algo variables start ------------------------
## See https://github.com/opensafely/post-covid-diabetes/blob/main/analysis/common_variables.py 

## Type 1 Diabetes 
# First date from primary+secondary, but also primary care date separately for diabetes algo
dataset.tmp_cov_date_t1dm_ctv3 = first_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, landmark_date).date
dataset.cov_date_t1dm = minimum_of(
    (first_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, landmark_date).date),
    (first_matching_event_apc_before(diabetes_type1_icd10, landmark_date).admission_date)
)
# Count codes (individually and together, for diabetes algo)
dataset.tmp_cov_count_t1dm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, landmark_date)
dataset.tmp_cov_count_t1dm_hes = count_matching_event_apc_before(diabetes_type1_icd10, landmark_date)
dataset.tmp_cov_count_t1dm = dataset.tmp_cov_count_t1dm_ctv3 + dataset.tmp_cov_count_t1dm_hes

## Type 2 Diabetes
# First date from primary+secondary, but also primary care date separately for diabetes algo)
dataset.tmp_cov_date_t2dm_ctv3 = first_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, landmark_date).date
dataset.cov_date_t2dm = minimum_of(
    (first_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, landmark_date).date),
    (first_matching_event_apc_before(diabetes_type2_icd10, landmark_date).admission_date)
)
# Count codes (individually and together, for diabetes algo)
dataset.tmp_cov_count_t2dm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, landmark_date)
dataset.tmp_cov_count_t2dm_hes = count_matching_event_apc_before(diabetes_type2_icd10, landmark_date)
dataset.tmp_cov_count_t2dm = dataset.tmp_cov_count_t2dm_ctv3 + dataset.tmp_cov_count_t2dm_hes

## Diabetes unspecified/other
# First date
dataset.cov_date_otherdm = first_matching_event_clinical_ctv3_before(diabetes_other_ctv3_clinical, landmark_date).date
# Count codes
dataset.tmp_cov_count_otherdm = count_matching_event_clinical_ctv3_before(diabetes_other_ctv3_clinical, landmark_date)

## Gestational diabetes
# First date
dataset.cov_date_gestationaldm = first_matching_event_clinical_ctv3_before(diabetes_gestational_ctv3_clinical, landmark_date).date

## Diabetes diagnostic codes
# First date
dataset.tmp_cov_date_poccdm = first_matching_event_clinical_ctv3_before(diabetes_diagnostic_ctv3_clinical, landmark_date).date
# Count codes
dataset.tmp_cov_count_poccdm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_diagnostic_ctv3_clinical, landmark_date)

### Other variables needed to define diabetes
## HbA1c
# Maximum HbA1c measure (in the same period)
dataset.tmp_cov_num_max_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(landmark_date))
    .numeric_value.maximum_for_patient()
)
# Date of first maximum HbA1c measure
dataset.tmp_cov_date_max_hba1c = ( 
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(landmark_date)) # this line of code probably not needed again
    .where(clinical_events.numeric_value == dataset.tmp_cov_num_max_hba1c_mmol_mol)
    .sort_by(clinical_events.date)
    .first_for_patient() 
    .date
)

## Diabetes drugs
# First dates
dataset.tmp_cov_date_insulin_snomed = first_matching_med_dmd_before(insulin_dmd, landmark_date).date
dataset.tmp_cov_date_antidiabetic_drugs_snomed = first_matching_med_dmd_before(antidiabetic_drugs_snomed_clinical, landmark_date).date
dataset.tmp_cov_date_nonmetform_drugs_snomed = first_matching_med_dmd_before(non_metformin_dmd, landmark_date).date # this extra step makes sense for the diabetes algorithm (otherwise not)

# Identify first date (in same period) that any diabetes medication was prescribed
dataset.tmp_cov_date_diabetes_medication = minimum_of(dataset.tmp_cov_date_insulin_snomed, dataset.tmp_cov_date_antidiabetic_drugs_snomed) # why excluding tmp_cov_date_nonmetform_drugs_snomed? -> this extra step makes sense for the diabetes algorithm (otherwise not)

# Identify first date (in same period) that any diabetes diagnosis codes were recorded
dataset.tmp_cov_date_first_diabetes_diag = minimum_of(
  dataset.cov_date_t2dm, 
  dataset.cov_date_t1dm,
  dataset.cov_date_otherdm,
  dataset.cov_date_gestationaldm,
  dataset.tmp_cov_date_poccdm,
  dataset.tmp_cov_date_diabetes_medication,
  dataset.tmp_cov_date_nonmetform_drugs_snomed
)

## DIABETES algo variables end ------------------------

## Known hypersensitivity / intolerance to metformin, on or before baseline
dataset.tmp_elig_date_metfin_allergy = first_matching_event_clinical_snomed_before(metformin_allergy_snomed_clinical, landmark_date).date

## Moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; stage 4/5), on or before baseline
dataset.tmp_elig_date_ckd_45 = minimum_of(
    first_matching_event_clinical_snomed_before(ckd_snomed_clinical_45, landmark_date).date,
    first_matching_event_apc_before(ckd_stage4_icd10 + ckd_stage5_icd10, landmark_date).admission_date
)

## Advance decompensated liver cirrhosis, on or before baseline
dataset.tmp_elig_date_liver_cirrhosis = minimum_of(
    first_matching_event_clinical_snomed_before(advanced_decompensated_cirrhosis_snomed_clinical + ascitic_drainage_snomed_clinical, landmark_date).date,
    first_matching_event_apc_before(advanced_decompensated_cirrhosis_icd10, landmark_date).admission_date
)

## Use of the following medications in the last 14 days (drug-drug interaction with metformin)
dataset.tmp_elig_date_metfin_interaction = last_matching_med_dmd_before(metformin_interaction_dmd, landmark_date).date


#######################################################################################
# Table 4) INTERVENTION/EXPOSURE variables
#######################################################################################
## Any metformin use before baseline, i.e., extract first metformin use before landmark and use it a) for Intervention/Exposure assignment and b) to establish true eligibility variable (i.e. before T2DM diagnosis) in R => assign tmp_ to this variable
dataset.tmp_exp_date_metfin_first = first_matching_med_dmd_before(metformin_dmd, landmark_date).date # create exp_date_metfin_first and elig_bin_metfin_before_baseline from this tmp_ variable (-> in R)
dataset.exp_count_metfin = (
  medications.where(
    medications.dmd_code.is_in(metformin_dmd))
    .where(medications.date.is_on_or_before(landmark_date))
    .count_for_patient()
)


#######################################################################################
# Table 5) Potential confounders
#######################################################################################

## Sex
dataset.cov_cat_sex = patients.sex.is_in(["female", "male"])

## Age at landmark_date
dataset.cov_num_age = patients.age_on(landmark_date)

## Ethnicity in 6 categories (mainly for diabetes algo())
dataset.cov_cat_ethnicity = (
    clinical_events.where(clinical_events.ctv3_code.is_in(ethnicity_codes))
    .sort_by(clinical_events.date)
    .last_for_patient()
    .ctv3_code.to_category(ethnicity_codes)
)

## Index of Multiple Deprevation Rank (rounded down to nearest 100). 5 categories.
imd_rounded = addresses.for_patient_on(landmark_date).imd_rounded
dataset.cov_cat_deprivation_5 = case(
    when((imd_rounded >=0) & (imd_rounded < int(32844 * 1 / 5))).then("1 (most deprived)"),
    when(imd_rounded < int(32844 * 2 / 5)).then("2"),
    when(imd_rounded < int(32844 * 3 / 5)).then("3"),
    when(imd_rounded < int(32844 * 4 / 5)).then("4"),
    when(imd_rounded < int(32844 * 5 / 5)).then("5 (least deprived)"),
    otherwise="unknown"
)

## Practice registration info at landmark_date
registered = practice_registrations.for_patient_on(landmark_date)
dataset.cov_cat_region = registered.practice_nuts1_region_name ## Region
dataset.cov_cat_stp = registered.practice_stp ## Practice
dataset.cov_cat_rural_urban = addresses.for_patient_on(landmark_date).rural_urban_classification ## Rurality

## Smoking status at landmark_date
tmp_most_recent_smoking_cat = last_matching_event_clinical_ctv3_before(smoking_clear, landmark_date).ctv3_code.to_category(smoking_clear)
tmp_ever_smoked = last_matching_event_clinical_ctv3_before(ever_smoking, landmark_date).exists_for_patient() # uses a different codelist with ONLY smoking codes
dataset.cov_cat_smoking_status = case(
    when(tmp_most_recent_smoking_cat == "S").then("S"),
    when(tmp_most_recent_smoking_cat == "E").then("E"),
    when((tmp_most_recent_smoking_cat == "N") & (tmp_ever_smoked == True)).then("E"),
    when(tmp_most_recent_smoking_cat == "N").then("N"),
    when((tmp_most_recent_smoking_cat == "M") & (tmp_ever_smoked == True)).then("E"),
    when(tmp_most_recent_smoking_cat == "M").then("M"),
    otherwise = "M"
)

## Care home resident at landmark_date, see https://github.com/opensafely/opioids-covid-research/blob/main/analysis/define_dataset_table.py
# Flag care home based on primis (patients in long-stay nursing and residential care)
tmp_care_home_code = last_matching_event_clinical_snomed_before(carehome, landmark_date).exists_for_patient()
# Flag care home based on TPP
tmp_care_home_tpp1 = addresses.for_patient_on(landmark_date).care_home_is_potential_match
tmp_care_home_tpp2 = addresses.for_patient_on(landmark_date).care_home_requires_nursing
tmp_care_home_tpp3 = addresses.for_patient_on(landmark_date).care_home_does_not_require_nursing
# combine
dataset.cov_bin_carehome_status = case(
    when(tmp_care_home_code).then(True),
    when(tmp_care_home_tpp1).then(True),
    when(tmp_care_home_tpp2).then(True),
    when(tmp_care_home_tpp3).then(True),
    otherwise = False
)

## Any sulfonylurea use before landmark_date (could also be a combo with metformin)
dataset.cov_date_sulfo = last_matching_med_dmd_before(sulfonylurea_dmd, landmark_date).date 
dataset.cov_date_dpp4 = last_matching_med_dmd_before(dpp4_dmd, landmark_date).date 
dataset.cov_date_tzd = last_matching_med_dmd_before(tzd_dmd, landmark_date).date 
dataset.cov_date_sglt2 = last_matching_med_dmd_before(sglt2_dmd, landmark_date).date 
dataset.cov_date_glp1 = last_matching_med_dmd_before(glp1_dmd, landmark_date).date 
dataset.cov_date_megli = last_matching_med_dmd_before(meglitinides_dmd, landmark_date).date 
dataset.cov_date_agi = last_matching_med_dmd_before(agi_dmd, landmark_date).date 
dataset.cov_date_insulin = last_matching_med_dmd_before(insulin_dmd, landmark_date).date

## Obesity, on or before landmark_date
dataset.cov_bin_obesity = (
    last_matching_event_clinical_snomed_before(bmi_obesity_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(bmi_obesity_icd10, landmark_date).exists_for_patient()
)

## Acute myocardial infarction, on or before landmark_date
dataset.cov_bin_ami = (
    last_matching_event_clinical_snomed_before(ami_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(ami_prior_icd10 + ami_icd10, landmark_date).exists_for_patient()
)

## All stroke, on or before landmark_date
dataset.cov_bin_all_stroke = (
    last_matching_event_clinical_snomed_before(stroke_isch_snomed_clinical + stroke_sah_hs_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(stroke_isch_icd10 + stroke_sah_hs_icd10, landmark_date).exists_for_patient()
)

## Other arterial embolism, on or before landmark_date
dataset.cov_bin_other_arterial_embolism = (
    last_matching_event_clinical_snomed_before(other_arterial_embolism_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(other_arterial_embolism_icd10, landmark_date).exists_for_patient()
)

## Venous thrombolism events, on or before landmark_date
dataset.cov_bin_vte = (
  last_matching_event_clinical_snomed_before(portal_vein_thrombosis_snomed_clinical + dvt_dvt_snomed_clinical + dvt_icvt_snomed_clinical + dvt_pregnancy_snomed_clinical + other_dvt_snomed_clinical + pe_snomed_clinical, landmark_date).exists_for_patient() |
  last_matching_event_apc_before(portal_vein_thrombosis_icd10 + dvt_dvt_icd10 + dvt_icvt_icd10 + dvt_pregnancy_icd10 + other_dvt_icd10 + icvt_pregnancy_icd10 + pe_icd10, landmark_date).exists_for_patient()
)

## Heart failure, on or before landmark_date
dataset.cov_bin_hf = (
    last_matching_event_clinical_snomed_before(hf_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(hf_icd10, landmark_date).exists_for_patient()
)

## Angina, on or before landmark_date
dataset.cov_bin_angina = (
    last_matching_event_clinical_snomed_before(angina_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(angina_icd10, landmark_date).exists_for_patient()
)

## Dementia, on or before landmark_date
dataset.cov_bin_dementia = (
    last_matching_event_clinical_snomed_before(dementia_snomed_clinical + dementia_vascular_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(dementia_icd10 + dementia_vascular_icd10, landmark_date).exists_for_patient()
)

## Cancer, on or before landmark_date
dataset.cov_bin_cancer = (
    last_matching_event_clinical_snomed_before(cancer_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(cancer_icd10, landmark_date).exists_for_patient()
)

## Hypertension, on or before landmark_date
dataset.cov_bin_hypertension = (
    last_matching_event_clinical_snomed_before(hypertension_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(hypertension_icd10, landmark_date).exists_for_patient()
)

## Depression, on or before landmark_date
dataset.cov_bin_depression = (
    last_matching_event_clinical_snomed_before(depression_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(depression_icd10, landmark_date).exists_for_patient()
)

## Chronic obstructive pulmonary disease, on or before landmark_date
dataset.cov_bin_copd = (
    last_matching_event_clinical_snomed_before(copd_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(copd_icd10, landmark_date).exists_for_patient()
)

## Liver disease, on or before landmark_date
dataset.cov_bin_liver_disease = (
    last_matching_event_clinical_snomed_before(liver_disease_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(liver_disease_icd10, landmark_date).exists_for_patient()
)

## Chronic kidney disease, on or before landmark_date
dataset.cov_bin_chronic_kidney_disease = (
    last_matching_event_clinical_snomed_before(ckd_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(ckd_icd10, landmark_date).exists_for_patient()
)

## PCOS, on or before landmark_date
dataset.cov_bin_pcos = (
    last_matching_event_clinical_snomed_before(pcos_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(pcos_icd10, landmark_date).exists_for_patient()
)

## Prediabetes, on or before landmark_date
# Any preDM diagnosis in primary care
tmp_cov_bin_prediabetes = last_matching_event_clinical_snomed_before(prediabetes_snomed, landmark_date).exists_for_patient()
# Any HbA1c preDM in primary care
tmp_cov_bin_predm_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(landmark_date))
    .where((clinical_events.numeric_value>=42) & (clinical_events.numeric_value<=47.9))
    .exists_for_patient()
)
# Any preDM diagnosis or Hb1Ac preDM range value (in period before landmark_date)
dataset.cov_bin_prediabetes = tmp_cov_bin_prediabetes | tmp_cov_bin_predm_hba1c_mmol_mol

## Diabetes complications (foot, retino, neuro, nephro), on or before landmark_date
dataset.cov_bin_diabetescomp = (
    last_matching_event_clinical_snomed_before(diabetescomp_snomed_clinical, landmark_date).exists_for_patient() |
    last_matching_event_apc_before(diabetescomp_icd10, landmark_date).exists_for_patient()
)

## Any HbA1c measurement, on or before landmark_date
dataset.cov_bin_hba1c_measurement = last_matching_event_clinical_snomed_before(hba1c_measurement_snomed, landmark_date).exists_for_patient()

## Any OGTT done, on or before landmark_date
dataset.cov_bin_ogtt_measurement = last_matching_event_clinical_snomed_before(ogtt_measurement_snomed, landmark_date).exists_for_patient()

## BMI, most recent value, within previous 2 years, on or before landmark_date
bmi_measurement = most_recent_bmi(
    where=clinical_events.date.is_on_or_between(landmark_date - days(2 * 366), landmark_date),
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

## HbA1c, most recent value, within previous 2 years, on or before landmark_date
dataset.cov_num_hba1c_mmol_mol = last_matching_event_clinical_snomed_between(hba1c_snomed, landmark_date - days(2*366), landmark_date).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years. 

## Total Cholesterol, most recent value, within previous 2 years, on or before landmark_date
dataset.tmp_cov_num_cholesterol = last_matching_event_clinical_snomed_between(cholesterol_snomed, landmark_date - days(2*366), landmark_date).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years. 

## HDL Cholesterol, most recent value, within previous 2 years, on or before landmark_date
dataset.tmp_cov_num_hdl_cholesterol = last_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, landmark_date - days(2*366), landmark_date).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years.

## Healthcare worker at the time they received a COVID-19 vaccination
dataset.cov_bin_healthcare_worker = (
  occupation_on_covid_vaccine_record.where(
    occupation_on_covid_vaccine_record.is_healthcare_worker == True)
    .exists_for_patient()
)


#######################################################################################
# Table 6) Outcomes
#######################################################################################

## Practice deregistration (no support from any practices anymore), after landmark_date: not based on landmark_date
# the latest dereg date
#tmp_dereg_date = (
#  practice_registrations.where(
#    practice_registrations.end_date.is_not_null())
#    .sort_by(practice_registrations.end_date)
#    .last_for_patient()
#    .end_date
#)
#tmp_not_dereg = practice_registrations.where(practice_registrations.end_date.is_null()).exists_for_patient() # has not left the practice
#dataset.out_date_dereg = case(when(tmp_not_dereg == False).then(tmp_dereg_date))

## Practice deregistration date 2: Based on registration on landmark_date
dataset.out_date_dereg = registered.end_date

## First covid-19 related hospital admission, between landmark_date and studyend_date (incl. those dates)
dataset.out_date_covid19_hes = first_matching_event_apc_between(covid_codes_incl_clin_diag, landmark_date, studyend_date).admission_date
## First covid-19 related emergency attendance, between landmark_date and studyend_date (incl. those dates)
dataset.out_date_covid19_emergency = first_matching_event_ec_snomed_between(covid_emergency, landmark_date, studyend_date).arrival_date
# combined: First covid-19 related hospitalization
dataset.out_date_covid19_hosp = minimum_of(dataset.out_date_covid19_hes, dataset.out_date_covid19_emergency)

## First COVID-19 diagnosis in primary care, between landmark_date and studyend_date (incl. those dates)
tmp_covid19_primary_care_date = first_matching_event_clinical_ctv3_between(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae, landmark_date, studyend_date).date
## First positive SARS-COV-2 PCR in primary care, between landmark_date and studyend_date (incl. those dates)
tmp_covid19_sgss_date = (
  sgss_covid_all_tests.where(sgss_covid_all_tests.specimen_taken_date.is_on_or_between(landmark_date, studyend_date))
  .where(sgss_covid_all_tests.is_positive)
  .sort_by(sgss_covid_all_tests.specimen_taken_date)
  .first_for_patient()
  .specimen_taken_date
)
## First COVID-19 diagnosis in primary care, or pos. test in primary care, or covid-19 hosp, between landmark_date and studyend_date (incl. those dates)
dataset.out_date_covid19 = minimum_of(tmp_covid19_primary_care_date, tmp_covid19_sgss_date, dataset.out_date_covid19_hosp)

## COVID-related death, between landmark_date and studyend_date (incl. those dates), stated anywhere on any of the 15 death certificate options
tmp_out_bin_death_covid = matching_death_between(covid_codes_incl_clin_diag, landmark_date, studyend_date)
dataset.out_date_covid19_death = case(when(tmp_out_bin_death_covid).then(ons_deaths.date))
# combined: First covid-19 related hospitalization or death after landmark_date, i.e. "severe COVID" outcome = secondary outcome
dataset.out_date_covid19_severe = minimum_of(dataset.out_date_covid19_death, dataset.out_date_covid19_hosp)

## First Long COVID code in primary care, between landmark_date and studyend_date (incl. those dates), based on https://github.com/opensafely/long-covid/blob/main/analysis/codelists.py
dataset.out_date_long_covid = first_matching_event_clinical_snomed_between(long_covid_diagnostic_snomed_clinical + long_covid_referral_snomed_clinical + long_covid_assessment_snomed_clinical, landmark_date, studyend_date).date
## First Viral fatigue code in primary care, between landmark_date and studyend_date (incl. those dates), based on https://github.com/opensafely/long-covid/blob/main/analysis/codelists.py
dataset.out_date_fatigue = first_matching_event_clinical_snomed_between(post_viral_fatigue_snomed_clinical, landmark_date, studyend_date).date
# combined: First Long COVID code or Viral Fatigue code after landmark_date = primary outcome
dataset.out_date_long_fatigue = minimum_of(dataset.out_date_long_covid, dataset.out_date_fatigue)
