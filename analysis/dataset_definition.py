#######################################################################################
# IMPORT
#######################################################################################

## ehrQL functions
from ehrql import (
    case,
    create_dataset,
    days,
    when,
    minimum_of,
    maximum_of
)

## TPP tables
from ehrql.tables.tpp import (
    addresses,
    clinical_events,
    apcs,
    medications,
    patients,
    practice_registrations,
    occupation_on_covid_vaccine_record,
    ons_deaths,
    sgss_covid_all_tests,
    vaccinations,
)

## All codelists from codelists.py
from codelists import *

## variable helper functions 
from variable_helper_functions import *

## json (for the dates)
import json

# numpy for random seed - and set random seed
import numpy as np 
np.random.seed(1928374) # random seed

#######################################################################################
# DEFINE the dates: Import study dates defined in "study-dates.R" and exported to JSON
#######################################################################################
with open("output/study_dates.json") as f:
  study_dates = json.load(f)

studystart_date = study_dates["studystart_date"]
studyend_date = study_dates["studyend_date"]

#######################################################################################
# DEFINE the baseline date (= index date) based on SARS-CoV-2 infection
#######################################################################################
## First COVID-19 event in primary care in recruitment period
tmp_covid19_primary_care_date = first_matching_event_clinical_ctv3_between(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae, studystart_date, studyend_date).date

## First positive SARS-COV-2 PCR in recruitment period
tmp_covid19_sgss_date = (
  sgss_covid_all_tests.where(sgss_covid_all_tests.specimen_taken_date.is_on_or_between(studystart_date, studyend_date))
  .where(sgss_covid_all_tests.is_positive)
  .sort_by(sgss_covid_all_tests.specimen_taken_date)
  .first_for_patient()
  .specimen_taken_date
)

"""
## First covid-19 related hospital admission in recruitment period // exclude since we are only (?) interested in recruitment in primary care? / all_diagnoses versus primary/sec diagnosis?
tmp_covid19_hes_date = (
    hospital_admissions.where(hospital_admissions.all_diagnoses.is_in(covid_codes)) # double-check with https://github.com/opensafely/comparative-booster-spring2023/blob/main/analysis/codelists.py uses a different codelist: codelists/opensafely-covid-identification.csv
    .where(hospital_admissions.admission_date.is_on_or_between(studystart_date,studyend_date))
    .sort_by(hospital_admissions.admission_date)
    .first_for_patient()
    .admission_date
)
"""

### Define (first) baseline date within recruitment period
baseline_date = minimum_of(tmp_covid19_primary_care_date, tmp_covid19_sgss_date)


#######################################################################################
# INITIALISE the dataset, define the baseline date and event, and set the dummy dataset size
#######################################################################################
dataset = create_dataset()
dataset.configure_dummy_data(population_size=1000)
dataset.baseline_date = baseline_date
dataset.define_population(patients.exists_for_patient())


#######################################################################################
# QUALITY ASSURANCES and completeness criteria
#######################################################################################
# population variables for dataset definition 
dataset.qa_bin_is_female_or_male = patients.sex.is_in(["female", "male"]) 
dataset.qa_bin_was_adult = (patients.age_on(baseline_date) >= 18) & (patients.age_on(baseline_date) <= 110) 
dataset.qa_bin_was_alive = (((patients.date_of_death.is_null()) | (patients.date_of_death.is_after(baseline_date))) & 
        ((ons_deaths.date.is_null()) | (ons_deaths.date.is_after(baseline_date))))
dataset.qa_bin_known_imd = addresses.for_patient_on(baseline_date).exists_for_patient() # known deprivation
dataset.qa_bin_was_registered = practice_registrations.spanning(baseline_date - days(366), baseline_date).exists_for_patient() # only include if registered on baseline date spanning back 1 year. Calculated from 1 year = 365.25 days, taking into account leap years.
# double-check line above against code from Will, line 98: https://github.com/opensafely/comparative-booster-spring2023/blob/main/analysis/dataset_definition.py 

## Year of birth
dataset.qa_num_birth_year = patients.date_of_birth

## Date of death
dataset.qa_date_of_death = ons_deaths.date

## Pregnancy (over entire study period, not based on baseline date -> don't have helper function for this)
dataset.qa_bin_pregnancy = clinical_events.where(clinical_events.snomedct_code.is_in(pregnancy_snomed_clinical)).exists_for_patient()

## Combined oral contraceptive pill
dataset.qa_bin_cocp = last_matching_med_dmd_before(cocp_dmd, baseline_date).exists_for_patient()

## Hormone replacement therapy
dataset.qa_bin_hrt = last_matching_med_dmd_before(hrt_dmd, baseline_date).exists_for_patient()

## Prostate cancer (over entire study period, not based on baseline date -> don't have helper function for this)
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
# DEMOGRAPHIC variables
#######################################################################################
## Sex
dataset.cov_cat_sex = patients.sex

## Age
dataset.cov_num_age = patients.age_on(baseline_date)

## Ethnicity in 6 categories based on codelists/opensafely-ethnicity.csv only. https://github.com/opensafely/comparative-booster-spring2023/blob/main/analysis/codelists.py  
dataset.cov_cat_ethnicity = (
    clinical_events.where(clinical_events.ctv3_code.is_in(ethnicity_codes))
    .sort_by(clinical_events.date)
    .last_for_patient()
    .ctv3_code.to_category(ethnicity_codes)
)

## Index of Multiple Deprevation Rank (rounded down to nearest 100). 5 categories.
imd_rounded = addresses.for_patient_on(baseline_date).imd_rounded
dataset.cov_cat_deprivation_5 = case(
    when((imd_rounded >=0) & (imd_rounded < int(32844 * 1 / 5))).then("1 (most deprived)"),
    when(imd_rounded < int(32844 * 2 / 5)).then("2"),
    when(imd_rounded < int(32844 * 3 / 5)).then("3"),
    when(imd_rounded < int(32844 * 4 / 5)).then("4"),
    when(imd_rounded < int(32844 * 5 / 5)).then("5 (least deprived)"),
    otherwise="unknown"
)

## Registration info at baseline date
registered = practice_registrations.for_patient_on(baseline_date)

## Region
dataset.cov_cat_region = registered.practice_nuts1_region_name

## Rurality
dataset.cov_cat_rural_urban = addresses.for_patient_on(baseline_date).rural_urban_classification

## Practice registration
dataset.cov_cat_stp = registered.practice_stp


#######################################################################################
# ELIGIBILITY variables
#######################################################################################
### DIABETES start ---------
### BASED on https://github.com/opensafely/post-covid-diabetes/blob/main/analysis/common_variables.py 

## Type 1 Diabetes 
# Latest date from primary+secondary, but also primary care date separately for diabetes algo)
dataset.tmp_cov_date_t1dm_ctv3 = last_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, baseline_date).date
dataset.cov_date_t1dm = minimum_of(
    (last_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, baseline_date).date),
    (last_matching_event_apc_before(diabetes_type1_icd10, baseline_date).admission_date)
)
# Count codes (individually and together, for diabetes algo)
dataset.tmp_cov_count_t1dm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, baseline_date)
dataset.tmp_cov_count_t1dm_hes = count_matching_event_apc_before(diabetes_type1_icd10, baseline_date)
dataset.tmp_cov_count_t1dm = dataset.tmp_cov_count_t1dm_ctv3 + dataset.tmp_cov_count_t1dm_hes

## Type 2 Diabetes
# Latest date from primary+secondary, but also primary care date separately for diabetes algo)
dataset.tmp_cov_date_t2dm_ctv3 = last_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, baseline_date).date
dataset.cov_date_t2dm = minimum_of(
    (last_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, baseline_date).date),
    (last_matching_event_apc_before(diabetes_type2_icd10, baseline_date).admission_date)
)
# Count codes (individually and together, for diabetes algo)
dataset.tmp_cov_count_t2dm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, baseline_date)
dataset.tmp_cov_count_t2dm_hes = count_matching_event_apc_before(diabetes_type2_icd10, baseline_date)
dataset.tmp_cov_count_t2dm = dataset.tmp_cov_count_t2dm_ctv3 + dataset.tmp_cov_count_t2dm_hes

## Diabetes unspecified/other
# Latest date
dataset.cov_date_otherdm = last_matching_event_clinical_ctv3_before(diabetes_other_ctv3_clinical, baseline_date).date
# Count codes
dataset.tmp_cov_count_otherdm = count_matching_event_clinical_ctv3_before(diabetes_other_ctv3_clinical, baseline_date)

## Gestational diabetes
# Latest date
dataset.cov_date_gestationaldm = last_matching_event_clinical_ctv3_before(diabetes_gestational_ctv3_clinical, baseline_date).date

## Diabetes diagnostic codes
# Latest date
dataset.tmp_cov_date_poccdm = last_matching_event_clinical_ctv3_before(diabetes_diagnostic_ctv3_clinical, baseline_date).date
# Count codes
dataset.tmp_cov_count_poccdm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_diagnostic_ctv3_clinical, baseline_date)

### Other variables needed to define diabetes
## HbA1c
# Maximum HbA1c measure (in period before baseline_date) -> don't have helper function for this
dataset.tmp_cov_num_max_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(baseline_date))
    .numeric_value.maximum_for_patient()
)
# Date of latest maximum HbA1c measure
dataset.tmp_cov_date_max_hba1c = ( 
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(baseline_date)) # this line of code probably not needed again
    .where(clinical_events.numeric_value == dataset.tmp_cov_num_max_hba1c_mmol_mol)
    .sort_by(clinical_events.date)
    .last_for_patient() # translates in cohortextractor to "on_most_recent_day_of_measurement=True"
    .date
)

## Diabetes drugs
# Latest dates
dataset.tmp_cov_date_insulin_snomed = last_matching_med_dmd_before(insulin_snomed_clinical, baseline_date).date
dataset.tmp_cov_date_antidiabetic_drugs_snomed = last_matching_med_dmd_before(antidiabetic_drugs_snomed_clinical, baseline_date).date
dataset.tmp_cov_date_nonmetform_drugs_snomed = last_matching_med_dmd_before(non_metformin_dmd, baseline_date).date # this extra step makes sense for the diabetes algorithm (otherwise not)

# Identify latest date (in period before baseline_date) that any diabetes medication was prescribed
dataset.tmp_cov_date_diabetes_medication = maximum_of(dataset.tmp_cov_date_insulin_snomed, dataset.tmp_cov_date_antidiabetic_drugs_snomed) # why excluding tmp_cov_date_nonmetform_drugs_snomed? -> this extra step makes sense for the diabetes algorithm (otherwise not)

# Identify latest date (in period before baseline_date) that any diabetes diagnosis codes were recorded
dataset.tmp_cov_date_latest_diabetes_diag = maximum_of(
  dataset.cov_date_t2dm, 
  dataset.cov_date_t1dm,
  dataset.cov_date_otherdm,
  dataset.cov_date_gestationaldm,
  dataset.tmp_cov_date_poccdm,
  dataset.tmp_cov_date_diabetes_medication,
  dataset.tmp_cov_date_nonmetform_drugs_snomed
)

### DIABETES end ---------


## Prediabetes, on or before baseline
# Date of preDM code in primary care
tmp_cov_date_prediabetes = last_matching_event_clinical_snomed_before(prediabetes_snomed, baseline_date).date
# Date of preDM HbA1c measure in period before baseline_date in preDM range (mmol/mol): 42-47.9
tmp_cov_date_predm_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(baseline_date))
    .where((clinical_events.numeric_value>=42) & (clinical_events.numeric_value<=47.9))
    .sort_by(clinical_events.date)
    .last_for_patient()
    .date
)
# Latest date (in period before baseline_date) that any prediabetes was diagnosed or HbA1c in preDM range
dataset.cov_date_prediabetes = maximum_of(tmp_cov_date_prediabetes, tmp_cov_date_predm_hba1c_mmol_mol) 

# Any preDM diagnosis in primary care
tmp_cov_bin_prediabetes = last_matching_event_clinical_snomed_before(prediabetes_snomed, baseline_date).exists_for_patient()
# Any HbA1c preDM in primary care
tmp_cov_bin_predm_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(baseline_date))
    .where((clinical_events.numeric_value>=42) & (clinical_events.numeric_value<=47.9))
    .exists_for_patient()
)
# Any preDM diagnosis or Hb1Ac preDM range value (in period before baseline_date)
dataset.cov_bin_prediabetes = tmp_cov_bin_prediabetes | tmp_cov_bin_predm_hba1c_mmol_mol


## Hospitalization at baseline (incl. 1 day prior)
dataset.cov_bin_hosp_baseline = (
    apcs.where(apcs.admission_date.is_on_or_between(baseline_date - days(1), baseline_date))
    .exists_for_patient()
)

## Metformin use at baseline, defined as receiving a metformin prescription up until 6 months prior to baseline date (assuming half-yearly prescription for stable diabetes across GPs in the UK)
dataset.cov_bin_metfin_before_baseline = last_matching_med_dmd_between(metformin_codes_dmd, (baseline_date - days(183)), baseline_date).exists_for_patient() # https://www.opencodelists.org/codelist/user/john-tazare/metformin-dmd/48e43356/ / # Calculated from 1 year = 365.25 days, taking into account leap years. 
dataset.cov_date_metfin_before_baseline = last_matching_med_dmd_between(metformin_codes_dmd, (baseline_date - days(183)), baseline_date).date

## Known hypersensitivity / intolerance to metformin, on or before baseline
dataset.cov_bin_metfin_allergy = last_matching_event_clinical_snomed_before(metformin_allergy_snomed_clinical, baseline_date).exists_for_patient() 

## Moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; stage 4/5), on or before baseline
dataset.cov_bin_ckd_45 = (
    last_matching_event_clinical_snomed_before(ckd_snomed_clinical_45, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(ckd_stage4_icd10 + ckd_stage5_icd10, baseline_date).exists_for_patient()
)

## Advance decompensated liver cirrhosis, on or before baseline
dataset.cov_bin_liver_cirrhosis = (
    last_matching_event_clinical_snomed_before(advanced_decompensated_cirrhosis_snomed_clinical + ascitic_drainage_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(advanced_decompensated_cirrhosis_icd10, baseline_date).exists_for_patient()
)

## Use of the following medications in the last 14 days (drug-drug interaction with metformin)
dataset.cov_bin_metfin_interaction = last_matching_med_dmd_between(metformin_interaction_dmd, (baseline_date - days(14)), baseline_date).exists_for_patient()
dataset.cov_date_metfin_interaction = last_matching_med_dmd_between(metformin_interaction_dmd, (baseline_date - days(14)), baseline_date).date

## Prior Long COVID diagnosis, based on https://github.com/opensafely/long-covid/blob/main/analysis/codelists.py
dataset.cov_bin_long_covid = last_matching_event_clinical_snomed_before(long_covid_diagnostic_snomed_clinical + long_covid_referral_snomed_clinical + long_covid_assessment_snomed_clinical, baseline_date).exists_for_patient()
dataset.cov_date_long_covid = last_matching_event_clinical_snomed_before(long_covid_diagnostic_snomed_clinical + long_covid_referral_snomed_clinical + long_covid_assessment_snomed_clinical, baseline_date).date



#######################################################################################
# Potential CONFOUNDER variables / Covariates for IPTW and IPCW
#######################################################################################
## Smoking status at baseline
tmp_most_recent_smoking_cat = last_matching_event_clinical_ctv3_before(smoking_clear, baseline_date).ctv3_code.to_category(smoking_clear)
tmp_ever_smoked = last_matching_event_clinical_ctv3_before(ever_smoking, baseline_date).exists_for_patient() # uses a different codelist with ONLY smoking codes
dataset.cov_cat_smoking_status = case(
    when(tmp_most_recent_smoking_cat == "S").then("S"),
    when(tmp_most_recent_smoking_cat == "E").then("E"),
    when((tmp_most_recent_smoking_cat == "N") & (tmp_ever_smoked == True)).then("E"),
    when(tmp_most_recent_smoking_cat == "N").then("N"),
    when((tmp_most_recent_smoking_cat == "M") & (tmp_ever_smoked == True)).then("E"),
    when(tmp_most_recent_smoking_cat == "M").then("M"),
    otherwise = "M"
)

## Care home resident at baseline, like e.g. https://github.com/opensafely/opioids-covid-research/blob/main/analysis/define_dataset_table.py
# Flag care home based on primis (patients in long-stay nursing and residential care)
tmp_care_home_code = last_matching_event_clinical_snomed_before(carehome, baseline_date).exists_for_patient()
# Flag care home based on TPP
tmp_care_home_tpp1 = addresses.for_patient_on(baseline_date).care_home_is_potential_match
tmp_care_home_tpp2 = addresses.for_patient_on(baseline_date).care_home_requires_nursing
tmp_care_home_tpp3 = addresses.for_patient_on(baseline_date).care_home_does_not_require_nursing
# combine
dataset.cov_bin_carehome_status = case(
    when(tmp_care_home_code).then(True),
    when(tmp_care_home_tpp1).then(True),
    when(tmp_care_home_tpp2).then(True),
    when(tmp_care_home_tpp3).then(True),
    otherwise = False
)

## Obesity, on or before baseline
dataset.cov_bin_obesity = (
    last_matching_event_clinical_snomed_before(bmi_obesity_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(bmi_obesity_icd10, baseline_date).exists_for_patient()
)

## Acute myocardial infarction, on or before baseline
dataset.cov_bin_ami = (
    last_matching_event_clinical_snomed_before(ami_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(ami_prior_icd10 + ami_icd10, baseline_date).exists_for_patient()
)

## All stroke, on or before baseline
dataset.cov_bin_all_stroke = (
    last_matching_event_clinical_snomed_before(stroke_isch_snomed_clinical + stroke_sah_hs_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(stroke_isch_icd10 + stroke_sah_hs_icd10, baseline_date).exists_for_patient()
)

## Other arterial embolism, on or before baseline
dataset.cov_bin_other_arterial_embolism = (
    last_matching_event_clinical_snomed_before(other_arterial_embolism_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(other_arterial_embolism_icd10, baseline_date).exists_for_patient()
)

## Venous thrombolism events, on or before baseline
dataset.cov_bin_vte = (
  last_matching_event_clinical_snomed_before(portal_vein_thrombosis_snomed_clinical + dvt_dvt_snomed_clinical + dvt_icvt_snomed_clinical + dvt_pregnancy_snomed_clinical + other_dvt_snomed_clinical + pe_snomed_clinical, baseline_date).exists_for_patient() |
  last_matching_event_apc_before(portal_vein_thrombosis_icd10 + dvt_dvt_icd10 + dvt_icvt_icd10 + dvt_pregnancy_icd10 + other_dvt_icd10 + icvt_pregnancy_icd10 + pe_icd10, baseline_date).exists_for_patient()
)

## Heart failure, on or before baseline
dataset.cov_bin_hf = (
    last_matching_event_clinical_snomed_before(hf_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(hf_icd10, baseline_date).exists_for_patient()
)

## Angina, on or before baseline
dataset.cov_bin_angina = (
    last_matching_event_clinical_snomed_before(angina_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(angina_icd10, baseline_date).exists_for_patient()
)

## Dementia, on or before baseline
dataset.cov_bin_dementia = (
    last_matching_event_clinical_snomed_before(dementia_snomed_clinical + dementia_vascular_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(dementia_icd10 + dementia_vascular_icd10, baseline_date).exists_for_patient()
)

## Cancer, on or before baseline
dataset.cov_bin_cancer = (
    last_matching_event_clinical_snomed_before(cancer_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(cancer_icd10, baseline_date).exists_for_patient()
)

## Hypertension, on or before baseline
dataset.cov_bin_hypertension = (
    last_matching_event_clinical_snomed_before(hypertension_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(hypertension_icd10, baseline_date).exists_for_patient() | 
    last_matching_med_dmd_before(hypertension_drugs_dmd, baseline_date).exists_for_patient() # not sure this is a good idea..
)

## Depression, on or before baseline
dataset.cov_bin_depression = (
    last_matching_event_clinical_snomed_before(depression_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(depression_icd10, baseline_date).exists_for_patient()
)

## Chronic obstructive pulmonary disease, on or before baseline
dataset.cov_bin_copd = (
    last_matching_event_clinical_snomed_before(copd_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(copd_icd10, baseline_date).exists_for_patient()
)

## Liver disease, on or before baseline
dataset.cov_bin_liver_disease = (
    last_matching_event_clinical_snomed_before(liver_disease_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(liver_disease_icd10, baseline_date).exists_for_patient()
)

## Chronic kidney disease, on or before baseline
dataset.cov_bin_chronic_kidney_disease = (
    last_matching_event_clinical_snomed_before(ckd_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(ckd_icd10, baseline_date).exists_for_patient()
)

## Gestational diabetes
dataset.cov_bin_gestationaldm = (
    last_matching_event_clinical_ctv3_before(diabetes_gestational_ctv3_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(gestationaldm_icd10, baseline_date).exists_for_patient()
)

## PCOS
dataset.cov_bin_pcos = (
    last_matching_event_clinical_snomed_before(pcos_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(pcos_icd10, baseline_date).exists_for_patient()
)

## Type 1 Diabetes
dataset.cov_bin_t1dm = (
    last_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(diabetes_type1_icd10, baseline_date).exists_for_patient()
)

## Diabetes complications (foot, retino, neuro, nephro)
dataset.cov_bin_diabetescomp = (
    last_matching_event_clinical_snomed_before(diabetescomp_snomed_clinical, baseline_date).exists_for_patient() |
    last_matching_event_apc_before(diabetescomp_icd10, baseline_date).exists_for_patient()
)

### Any HbA1c measurement
dataset.cov_bin_hba1c_measurement = last_matching_event_clinical_snomed_before(hba1c_measurement_snomed, baseline_date).exists_for_patient()

### Any OGTT done
dataset.cov_bin_ogtt_measurement = last_matching_event_clinical_snomed_before(ogtt_measurement_snomed, baseline_date).exists_for_patient()

### Covid-19 vaccination history
dataset.cov_count_covid_vaccines = (
  vaccinations
  .where(vaccinations.target_disease.is_in(["SARS-2 CORONAVIRUS"]))
  .where(vaccinations.date.is_on_or_before(baseline_date))
  .count_for_patient()
)
dataset.cov_date_recent_covid_vaccines = (
  vaccinations
  .where(vaccinations.target_disease.is_in(["SARS-2 CORONAVIRUS"]))
  .where(vaccinations.date.is_on_or_before(baseline_date))
  .sort_by(vaccinations.date)
  .last_for_patient()
  .date
)

## BMI, most recent value, within previous 2 years
bmi_measurement = most_recent_bmi(
    where=clinical_events.date.is_after(baseline_date - days(2 * 366)),
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

## HbA1c, most recent value, within previous 2 years
dataset.cov_num_hba1c_mmol_mol = last_matching_event_clinical_snomed_between(hba1c_snomed, baseline_date - days(2*366), baseline_date).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years. 

## Total Cholesterol, most recent value, within previous 2 years
dataset.tmp_cov_num_cholesterol = last_matching_event_clinical_snomed_between(cholesterol_snomed, baseline_date - days(2*366), baseline_date).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years. 

## HDL Cholesterol, most recent value, within previous 2 years
dataset.tmp_cov_num_hdl_cholesterol = last_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, baseline_date - days(2*366), baseline_date).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years.

## Number of consultations in year prior to pandemic (2019)
#dataset.cov_num_consultation_rate = (
#  appointments.where(
#    appointments.status.is_in(["Arrived", "In Progress", "Finished", "Visit", "Waiting", "Patient Walked Out",]))
#    .where(appointments.seen_date.is_on_or_between(studystart_date - days(366),studystart_date)) # the year before the pandemic
#    .count_for_patient()
#)

## Healthcare worker at the time they received a COVID-19 vaccination
dataset.cov_bin_healthcare_worker = (
  occupation_on_covid_vaccine_record.where(
    occupation_on_covid_vaccine_record.is_healthcare_worker == True)
    .exists_for_patient()
)


#######################################################################################
# INTERVENTION/EXPOSURE variables
#######################################################################################
# METFORMIN, based on codelist https://www.opencodelists.org/codelist/user/john-tazare/metformin-dmd/48e43356/
dataset.exp_date_first_metfin = first_matching_med_dmd_between(metformin_codes_dmd, baseline_date, studyend_date).date 
dataset.exp_count_metfin = (
  medications.where(
    medications.dmd_code.is_in(metformin_codes_dmd))
    .where(medications.date.is_on_or_after(baseline_date))
    .count_for_patient()
)
# we will code the treatment window in R, but to compare
dataset.exp_bin_10d_metfin = first_matching_med_dmd_between(metformin_codes_dmd, baseline_date, baseline_date + days(10)).exists_for_patient()


#######################################################################################
# OUTCOME variables
#######################################################################################
# OUTCOMES irrespective of outcome window, this will be coded in R

## Practice deregistration
## the latest
tmp_dereg_date = (
  practice_registrations.where(
    practice_registrations.end_date.is_not_null())
    .sort_by(practice_registrations.end_date)
    .last_for_patient()
    .end_date
)
# has not left the practice
tmp_not_dereg = practice_registrations.where(practice_registrations.end_date.is_null()).exists_for_patient()
dataset.out_date_dereg = case(when(tmp_not_dereg == False).then(tmp_dereg_date))

# SARS-CoV-2, based on https://github.com/opensafely/long-covid/blob/main/analysis/codelists.py
# First covid-19 related hospital admission, after baseline date
dataset.out_date_covid19_hes = first_matching_event_apc_between(covid_codes_incl_clin_diag, baseline_date, studyend_date).admission_date
# First emergency attendance for covid, after baseline date
dataset.out_date_covid19_emergency = first_matching_event_ec_snomed_between(covid_emergency, baseline_date, studyend_date).arrival_date
## First hospitalisation or emergency attendance for covid, after baseline date
dataset.out_date_covid_hosp = minimum_of(dataset.out_date_covid19_hes, dataset.out_date_covid19_emergency)

# First COVID-19 code (diagnosis, positive test or sequelae) in primary care, after baseline date
dataset.out_date_covid19_primary_care = first_matching_event_clinical_ctv3_between(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae, baseline_date, studyend_date).date
# First positive SARS-COV-2 PCR, after baseline date
dataset.out_date_covid19_sgss = (
  sgss_covid_all_tests.where(
    sgss_covid_all_tests.specimen_taken_date.is_on_or_between(baseline_date, studyend_date))
    .where(sgss_covid_all_tests.is_positive)
    .sort_by(sgss_covid_all_tests.specimen_taken_date)
    .first_for_patient()
    .specimen_taken_date
)

# First covid-19 event (hospital admission, pos test, clinical diagnosis), after baseline date
dataset.out_date_covid19 = minimum_of(dataset.out_date_covid19_primary_care, dataset.out_date_covid19_sgss, dataset.out_date_covid19_hes, dataset.out_date_covid19_emergency)

## Covid-related Death
# dataset.out_death_date = ons_deaths.date # already defined as QA
# covid-related death (stated anywhere on any of the 15 death certificate options) # https://github.com/opensafely/comparative-booster-spring2023/blob/main/analysis/codelists.py uses a different codelist: codelists/opensafely-covid-identification.csv
tmp_out_bin_death_cause_covid = matching_death_between(covid_codes_incl_clin_diag, baseline_date, studyend_date)
# add default F
dataset.out_bin_death_cause_covid = case(
    when(tmp_out_bin_death_cause_covid).then(True),
    otherwise = False
)
dataset.out_date_death_cause_covid = case(when(tmp_out_bin_death_cause_covid).then(ons_deaths.date)) # double-check if used later or not
### First covid-19 related hosp, emergency hosp, or death, after baseline date
dataset.out_date_covid_hosp_death = minimum_of(dataset.out_date_death_cause_covid, dataset.out_date_covid19_hes, dataset.out_date_covid19_emergency)

## Long COVID, based on https://github.com/opensafely/long-covid/blob/main/analysis/codelists.py
dataset.out_bin_long_covid = first_matching_event_clinical_snomed_between(long_covid_diagnostic_snomed_clinical + long_covid_referral_snomed_clinical + long_covid_assessment_snomed_clinical, baseline_date, studyend_date).exists_for_patient()
dataset.out_date_long_covid_first = first_matching_event_clinical_snomed_between(long_covid_diagnostic_snomed_clinical + long_covid_referral_snomed_clinical + long_covid_assessment_snomed_clinical, baseline_date, studyend_date).date
# Any viral fatigue code in primary care after baseline date
dataset.out_bin_viral_fatigue = first_matching_event_clinical_snomed_between(post_viral_fatigue_snomed_clinical, baseline_date, studyend_date).exists_for_patient()
dataset.out_date_viral_fatigue_first = first_matching_event_clinical_snomed_between(post_viral_fatigue_snomed_clinical, baseline_date, studyend_date).date
# combined
dataset.out_date_long_fatigue = minimum_of(dataset.out_date_long_covid_first, dataset.out_date_viral_fatigue_first)