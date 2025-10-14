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
    occupation_on_covid_vaccine_record,
    appointments
)

## for import of diabetes algo created data
from ehrql.query_language import (
    table_from_file,
    PatientFrame,
    Series
)

## All codelists from codelists.py
from codelists import *

## variable helper functions 
from variable_helper_functions import *

import project_permissions

## json (for the dates)
import json

## for import of diabetes algo created data
from datetime import date

## import diabetes algo created data
@table_from_file("output/data_processed.csv.gz")
class data_processed(PatientFrame):
    #qa_num_birth_year = Series(int) # could import it, too, but creates friction with data formatting function
    ethnicity_cat = Series(str)
    t2dm_date = Series(date)
    tmp_first_diabetes_diag_date = Series(date)
    tmp_otherdm_date = Series(date)
    tmp_poccdm_date = Series(date)
    tmp_t1dm_date = Series(date)
    tmp_gestationaldm_date = Series(date)
    tmp_diabetes_medication_date = Series(date)
    tmp_nonmetform_drugs_dmd_date = Series(date)
    step_1 = Series(str)
    step_1a = Series(str)
    step_2 = Series(str)
    step_3 = Series(str)
    step_4 = Series(str)

# random seed (ideally use numpy, but currently not working on my local environment)
#import numpy as np 
#np.random.seed(19283) # random seed
import random
random.seed(19283) # random seed

#######################################################################################
# INITIALISE the dataset and set the dummy dataset size
#######################################################################################
dataset = create_dataset()
dataset.configure_dummy_data(population_size=7000)
dataset.define_population(patients.exists_for_patient())

#######################################################################################
# DEFINE the dates
#######################################################################################
dataset.elig_date_t2dm = data_processed.t2dm_date
with open("output/study_dates.json") as f:
  study_dates = json.load(f)
studyend_date = study_dates["studyend_date"]
pandemicstart_date = study_dates["pandemicstart_date"]
mid2018_date = study_dates["mid2018_date"]

#######################################################################################
# DM algorithm helper variables, to explore the DM algo
#######################################################################################
dataset.step_1 = data_processed.step_1
dataset.step_1a = data_processed.step_1a
dataset.step_2 = data_processed.step_2
dataset.step_3 = data_processed.step_3
dataset.step_4 = data_processed.step_4

dataset.tmp_first_diabetes_diag_date = data_processed.tmp_first_diabetes_diag_date
dataset.tmp_otherdm_date = data_processed.tmp_otherdm_date
dataset.tmp_poccdm_date = data_processed.tmp_poccdm_date
dataset.tmp_t1dm_date = data_processed.tmp_t1dm_date
dataset.tmp_gestationaldm_date = data_processed.tmp_gestationaldm_date
dataset.tmp_diabetes_medication_date = data_processed.tmp_diabetes_medication_date
dataset.tmp_nonmetform_drugs_dmd_date = data_processed.tmp_nonmetform_drugs_dmd_date

#######################################################################################
# Table 2) QUALITY ASSURANCES and completeness criteria
#######################################################################################
## Year of birth
dataset.qa_num_birth_year = patients.date_of_birth
# population variables for dataset definition 
dataset.qa_bin_is_female_or_male = patients.sex.is_in(["female", "male"]) 
dataset.qa_bin_was_adult = (patients.age_on(dataset.elig_date_t2dm) >= 18) & (patients.age_on(dataset.elig_date_t2dm) <= 85) 
dataset.qa_bin_was_alive = patients.is_alive_on(dataset.elig_date_t2dm)
dataset.qa_bin_known_imd = addresses.for_patient_on(dataset.elig_date_t2dm).exists_for_patient() # known deprivation
dataset.qa_bin_was_registered = practice_registrations.spanning_with_systmone(dataset.elig_date_t2dm - days(366), dataset.elig_date_t2dm).exists_for_patient() # see https://docs.opensafely.org/ehrql/reference/schemas/tpp/#practice_registrations.spanning. Calculated from 1 year = 365.25 days, taking into account leap year.

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

## Hospitalized at baseline (= T2DM diagnosis date)
dataset.elig_bin_hosp = (
   apcs.where(
      apcs.admission_date.is_on_or_before(dataset.elig_date_t2dm)
      & apcs.discharge_date.is_after(dataset.elig_date_t2dm)
      )
   .exists_for_patient()
)

## Care home resident at elig_date_t2dm, see https://github.com/opensafely/opioids-covid-research/blob/main/analysis/define_dataset_table.py
# Flag care home based on primis (patients in long-stay nursing and residential care)
tmp_care_home_code = last_matching_event_clinical_snomed_before(carehome, dataset.elig_date_t2dm).exists_for_patient()
# Flag care home based on TPP
tmp_care_home_tpp1 = addresses.for_patient_on(dataset.elig_date_t2dm).care_home_is_potential_match
tmp_care_home_tpp2 = addresses.for_patient_on(dataset.elig_date_t2dm).care_home_requires_nursing
tmp_care_home_tpp3 = addresses.for_patient_on(dataset.elig_date_t2dm).care_home_does_not_require_nursing
# combine
dataset.elig_bin_carehome_status = case(
    when(tmp_care_home_code).then(True),
    when(tmp_care_home_tpp1).then(True),
    when(tmp_care_home_tpp2).then(True),
    when(tmp_care_home_tpp3).then(True),
    otherwise = False
)

## Palliative care within 6 months prior to T2DM diagnosis
dataset.elig_date_palliative = last_matching_event_clinical_snomed_between(palliative_snomed, dataset.elig_date_t2dm - days(183), dataset.elig_date_t2dm).date

## Any metformin prescription BEFORE T2DM diagnosis
dataset.elig_date_metfin = first_matching_med_dmd_before(metformin_dmd, dataset.elig_date_t2dm).date
## Any other antidiabetic drug exposure BEFORE T2DM diagnosis
dataset.elig_date_sulfo = first_matching_med_dmd_before(sulfonylurea_dmd, dataset.elig_date_t2dm).date 
dataset.elig_date_dpp4 = first_matching_med_dmd_before(dpp4_dmd, dataset.elig_date_t2dm).date
dataset.elig_date_tzd = first_matching_med_dmd_before(tzd_dmd, dataset.elig_date_t2dm).date 
dataset.elig_date_sglt2 = first_matching_med_dmd_before(sglt2_dmd, dataset.elig_date_t2dm).date 
dataset.elig_date_glp1 = first_matching_med_dmd_before(glp1_dmd, dataset.elig_date_t2dm).date 
dataset.elig_date_megli = first_matching_med_dmd_before(meglitinides_dmd, dataset.elig_date_t2dm).date
dataset.elig_date_agi = first_matching_med_dmd_before(agi_dmd, dataset.elig_date_t2dm).date 
dataset.elig_date_insulin = first_matching_med_dmd_before(insulin_dmd, dataset.elig_date_t2dm).date

## Known hypersensitivity / intolerance to metformin, on or before elig_date_t2dm
dataset.elig_date_metfin_allergy_first = first_matching_event_clinical_snomed_before(metformin_allergy_snomed_clinical, dataset.elig_date_t2dm).date

## Moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; stage 4/5), on or before elig_date_t2dm
dataset.elig_date_ckd_45_first = minimum_of(
    first_matching_event_clinical_snomed_before(ckd_snomed_clinical_45, dataset.elig_date_t2dm).date,
    first_matching_event_apc_before(ckd_stage4_icd10 + ckd_stage5_icd10, dataset.elig_date_t2dm).admission_date
)

## Advance decompensated liver cirrhosis, on or before elig_date_t2dm
dataset.elig_date_liver_cirrhosis_first = minimum_of(
    first_matching_event_clinical_snomed_before(advanced_decompensated_cirrhosis_snomed_clinical + ascitic_drainage_snomed_clinical, dataset.elig_date_t2dm).date,
    first_matching_event_apc_before(advanced_decompensated_cirrhosis_icd10, dataset.elig_date_t2dm).admission_date
)

## Use of the following medications in the last 14 days before elig_date_t2dm (drug-drug interaction with metformin)
dataset.elig_date_metfin_interaction_last = last_matching_med_dmd_before(metformin_interaction_dmd, dataset.elig_date_t2dm).date
dataset.elig_date_metfin_interaction_landmark_last = last_matching_med_dmd_before(metformin_interaction_dmd, dataset.elig_date_t2dm + days(183)).date

## HbA1c, most recent value, within previous 2 years, on or before elig_date_t2dm + days(183), i.e., landmark date -> only used as an updated eligibility at landmark, not a baseline covariate! 
## cov_num_hba1c_mmol_mol (below) will be used for the baseline HbA1c exclusion
dataset.elig_num_hba1c_landmark_mmol_mol = last_matching_event_clinical_snomed_between(hba1c_snomed, dataset.elig_date_t2dm + days(183) - days(2*366), dataset.elig_date_t2dm + days(183)).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years. 

#######################################################################################
# Table 4) INTERVENTION/EXPOSURE variables
#######################################################################################
## First metformin prescription on/after T2DM diagnosis
dataset.exp_date_metfin_first = first_matching_med_dmd_between(metformin_dmd, dataset.elig_date_t2dm, studyend_date).date
dataset.exp_date_metfin_mono_first = first_matching_med_dmd_between(metformin_mono_dmd, dataset.elig_date_t2dm, studyend_date).date

## Other antidiabetic drug exposure on/after T2DM diagnosis (to define treatment strategy, or/and as intercurrent events for PP analysis)
dataset.exp_date_sulfo_first = first_matching_med_dmd_between(sulfonylurea_dmd, dataset.elig_date_t2dm, studyend_date).date 
dataset.exp_date_dpp4_first = first_matching_med_dmd_between(dpp4_dmd, dataset.elig_date_t2dm, studyend_date).date 
dataset.exp_date_tzd_first = first_matching_med_dmd_between(tzd_dmd, dataset.elig_date_t2dm, studyend_date).date 
dataset.exp_date_sglt2_first = first_matching_med_dmd_between(sglt2_dmd, dataset.elig_date_t2dm, studyend_date).date 
dataset.exp_date_glp1_first = first_matching_med_dmd_between(glp1_dmd, dataset.elig_date_t2dm, studyend_date).date 
dataset.exp_date_megli_first = first_matching_med_dmd_between(meglitinides_dmd, dataset.elig_date_t2dm, studyend_date).date 
dataset.exp_date_agi_first = first_matching_med_dmd_between(agi_dmd, dataset.elig_date_t2dm, studyend_date).date 
dataset.exp_date_insulin_first = first_matching_med_dmd_between(insulin_dmd, dataset.elig_date_t2dm, studyend_date).date

## Last metformin prescription before pandemic start: Use to establish who stopped before pandemic start (for PP analysis)
dataset.exp_date_metfin_last = last_matching_med_dmd_before(metformin_dmd, pandemicstart_date).date
dataset.exp_date_metfin_mono_last = last_matching_med_dmd_before(metformin_mono_dmd, pandemicstart_date).date 


#######################################################################################
# Table 5) Demographics, covariates and potential confounders
#######################################################################################
## Sex
dataset.cov_cat_sex = patients.sex

## Age at elig_date_t2dm
dataset.cov_num_age = patients.age_on(dataset.elig_date_t2dm)

## Ethnicity (import from the diabetes algo data)
dataset.cov_cat_ethnicity = data_processed.ethnicity_cat

## Index of Multiple Deprevation Rank (rounded down to nearest 100). 5 categories.
imd_rounded = addresses.for_patient_on(dataset.elig_date_t2dm).imd_rounded
dataset.cov_cat_deprivation_5 = case(
    when((imd_rounded >=0) & (imd_rounded < int(32844 * 1 / 5))).then("1 (most deprived)"),
    when(imd_rounded < int(32844 * 2 / 5)).then("2"),
    when(imd_rounded < int(32844 * 3 / 5)).then("3"),
    when(imd_rounded < int(32844 * 4 / 5)).then("4"),
    when(imd_rounded < int(32844 * 5 / 5)).then("5 (least deprived)"),
    otherwise="Unknown"
)
# dataset.cov_cat_deprivation_10 = case(
#     when((imd_rounded >= 0) & (imd_rounded <= int(32844 * 1 / 10))).then("1 (most deprived)"),
#     when(imd_rounded <= int(32844 * 2 / 10)).then("2"),
#     when(imd_rounded <= int(32844 * 3 / 10)).then("3"),
#     when(imd_rounded <= int(32844 * 4 / 10)).then("4"),
#     when(imd_rounded <= int(32844 * 5 / 10)).then("5"),
#     when(imd_rounded <= int(32844 * 6 / 10)).then("6"),
#     when(imd_rounded <= int(32844 * 7 / 10)).then("7"),
#     when(imd_rounded <= int(32844 * 8 / 10)).then("8"),
#     when(imd_rounded <= int(32844 * 9 / 10)).then("9"),
#     when(imd_rounded <= int(32844 * 10 / 10)).then("10 (least deprived)"),
#     otherwise="Unknown"
# )

## Practice registration info at elig_date_t2dm
# but use a mix between spanning (as per eligibility criteria) and for_patient_on() to sort the multiple rows: https://docs.opensafely.org/ehrql/reference/schemas/tpp/#practice_registrations.for_patient_on
spanning_regs = practice_registrations.spanning_with_systmone(dataset.elig_date_t2dm - days(366), dataset.elig_date_t2dm)
registered = spanning_regs.sort_by(
    practice_registrations.end_date,
    practice_registrations.practice_pseudo_id,
).last_for_patient()
dataset.strat_cat_region = registered.practice_nuts1_region_name ## Region
dataset.practice_date_systmone = registered.practice_systmone_go_live_date ## Date on which the practice started using the SystmOne EHR platform.
dataset.cov_cat_rural_urban = addresses.for_patient_on(dataset.elig_date_t2dm).rural_urban_classification ## Rurality

## Smoking status at elig_date_t2dm
tmp_most_recent_smoking_cat = last_matching_event_clinical_ctv3_before(smoking_clear, dataset.elig_date_t2dm).ctv3_code.to_category(smoking_clear)
tmp_ever_smoked = last_matching_event_clinical_ctv3_before(ever_smoking, dataset.elig_date_t2dm).exists_for_patient() # uses a different codelist with ONLY smoking codes
dataset.cov_cat_smoking_status = case(
    when(tmp_most_recent_smoking_cat == "S").then("S"),
    when(tmp_most_recent_smoking_cat == "E").then("E"),
    when((tmp_most_recent_smoking_cat == "N") & (tmp_ever_smoked == True)).then("E"),
    when(tmp_most_recent_smoking_cat == "N").then("N"),
    when((tmp_most_recent_smoking_cat == "M") & (tmp_ever_smoked == True)).then("E"),
    when(tmp_most_recent_smoking_cat == "M").then("M"),
    otherwise = "M"
)

## Healthcare worker at the time they received a COVID-19 vaccination
dataset.cov_bin_healthcare_worker = (
  occupation_on_covid_vaccine_record.where(
    occupation_on_covid_vaccine_record.is_healthcare_worker == True)
    .exists_for_patient()
)

## Consultation rate in previous year (mid2017 to mid2018) as a proxy for health seeking behaviour
### Consultation rate in 2019
tmp_cov_num_consrate = appointments.where(
    appointments.status.is_in([
        "Arrived",
        "In Progress",
        "Finished",
        "Visit",
        "Waiting",
        "Patient Walked Out",
        ]) & appointments.start_date.is_on_or_between(dataset.elig_date_t2dm - days(366), dataset.elig_date_t2dm)
        ).count_for_patient()    

dataset.cov_num_consrate = case(
    when(tmp_cov_num_consrate <= 365).then(tmp_cov_num_consrate),
    otherwise=365, # quality assurance
)

## Obesity, on or before elig_date_t2dm
dataset.cov_bin_obesity = (
    last_matching_event_clinical_snomed_before(bmi_obesity_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(bmi_obesity_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Acute myocardial infarction, on or before elig_date_t2dm
dataset.cov_bin_ami = (
    last_matching_event_clinical_snomed_before(ami_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(ami_prior_icd10 + ami_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## All stroke, on or before elig_date_t2dm
dataset.cov_bin_all_stroke = (
    last_matching_event_clinical_snomed_before(stroke_isch_snomed_clinical + stroke_sah_hs_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(stroke_isch_icd10 + stroke_sah_hs_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Other arterial embolism, on or before elig_date_t2dm
dataset.cov_bin_other_arterial_embolism = (
    last_matching_event_clinical_snomed_before(other_arterial_embolism_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(other_arterial_embolism_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Venous thrombolism events, on or before elig_date_t2dm
dataset.cov_bin_vte = (
  last_matching_event_clinical_snomed_before(portal_vein_thrombosis_snomed_clinical + dvt_dvt_snomed_clinical + dvt_icvt_snomed_clinical + dvt_pregnancy_snomed_clinical + other_dvt_snomed_clinical + pe_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
  last_matching_event_apc_before(portal_vein_thrombosis_icd10 + dvt_dvt_icd10 + dvt_icvt_icd10 + dvt_pregnancy_icd10 + other_dvt_icd10 + icvt_pregnancy_icd10 + pe_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Heart failure, on or before elig_date_t2dm
dataset.cov_bin_hf = (
    last_matching_event_clinical_snomed_before(hf_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(hf_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Angina, on or before elig_date_t2dm
dataset.cov_bin_angina = (
    last_matching_event_clinical_snomed_before(angina_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(angina_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Dementia, on or before elig_date_t2dm
dataset.cov_bin_dementia = (
    last_matching_event_clinical_snomed_before(dementia_snomed_clinical + dementia_vascular_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(dementia_icd10 + dementia_vascular_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Cancer, on or before elig_date_t2dm
dataset.cov_bin_cancer = (
    last_matching_event_clinical_snomed_before(cancer_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(cancer_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Hypertension, on or before elig_date_t2dm
dataset.cov_bin_hypertension = (
    last_matching_event_clinical_snomed_before(hypertension_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(hypertension_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Depression, on or before elig_date_t2dm
dataset.cov_bin_depression = (
    last_matching_event_clinical_snomed_before(depression_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(depression_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Chronic obstructive pulmonary disease, on or before elig_date_t2dm
dataset.cov_bin_copd = (
    last_matching_event_clinical_snomed_before(copd_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(copd_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Liver disease, on or before elig_date_t2dm
dataset.cov_bin_liver_disease = (
    last_matching_event_clinical_snomed_before(liver_disease_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(liver_disease_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Chronic kidney disease, on or before elig_date_t2dm
dataset.cov_bin_chronic_kidney_disease = (
    last_matching_event_clinical_snomed_before(ckd_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(ckd_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## PCOS, on or before elig_date_t2dm
dataset.cov_bin_pcos = (
    last_matching_event_clinical_snomed_before(pcos_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(pcos_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## Prediabetes, on or before elig_date_t2dm
# Any preDM diagnosis in primary care
tmp_cov_bin_prediabetes = last_matching_event_clinical_snomed_before(prediabetes_snomed, dataset.elig_date_t2dm).exists_for_patient()
# Any HbA1c preDM in primary care
tmp_cov_bin_predm_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(dataset.elig_date_t2dm))
    .where((clinical_events.numeric_value>=42) & (clinical_events.numeric_value<=47.9))
    .exists_for_patient()
)
# Any preDM diagnosis or Hb1Ac preDM range value (in period before elig_date_t2dm)
dataset.cov_bin_prediabetes = tmp_cov_bin_prediabetes | tmp_cov_bin_predm_hba1c_mmol_mol

## Diabetes complications (foot, retino, neuro, nephro), on or before elig_date_t2dm
dataset.cov_bin_diabetescomp = (
    last_matching_event_clinical_snomed_before(diabetescomp_snomed_clinical, dataset.elig_date_t2dm).exists_for_patient() |
    last_matching_event_apc_before(diabetescomp_icd10, dataset.elig_date_t2dm).exists_for_patient()
)

## BMI, most recent value, within previous 2 years, on or before elig_date_t2dm
bmi_measurement = most_recent_bmi(
    where=clinical_events.date.is_on_or_between(dataset.elig_date_t2dm - days(2 * 366), dataset.elig_date_t2dm),
    minimum_age_at_measurement=16,
)
dataset.cov_num_bmi_b = bmi_measurement.numeric_value
dataset.cov_cat_bmi_groups = case(
    when((dataset.cov_num_bmi_b < 18.5) & (dataset.cov_num_bmi_b >= 12.0)).then("Underweight"), # Set minimum to avoid any impossibly extreme values being classified as underweight
    when((dataset.cov_num_bmi_b >= 18.5) & (dataset.cov_num_bmi_b < 25.0)).then("Healthy weight (18.5-24.9)"),
    when((dataset.cov_num_bmi_b >= 25.0) & (dataset.cov_num_bmi_b < 30.0)).then("Overweight (25-29.9)"),
    when((dataset.cov_num_bmi_b >= 30.0) & (dataset.cov_num_bmi_b <= 70.0)).then("Obese (>30)"), # Set maximum to avoid any impossibly extreme values being classified as obese
    otherwise = "missing", 
)

## HbA1c, most recent value, within previous 2 years, on or before elig_date_t2dm
dataset.cov_num_hba1c_b = last_matching_event_clinical_snomed_between(hba1c_snomed, dataset.elig_date_t2dm - days(2*366), dataset.elig_date_t2dm).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years. 

## Total Cholesterol, most recent value, within previous 2 years, on or before elig_date_t2dm
dataset.cov_num_chol_b = last_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.elig_date_t2dm - days(2*366), dataset.elig_date_t2dm).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years. 

## HDL Cholesterol, most recent value, within previous 2 years, on or before elig_date_t2dm
dataset.cov_num_hdl_chol_b = last_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.elig_date_t2dm - days(2*366), dataset.elig_date_t2dm).numeric_value # Calculated from 1 year = 365.25 days, taking into account leap years.

## Count HbA1c measurement, within previous 2 years, on or before elig_date_t2dm
dataset.cov_num_counthba1c = count_matching_event_clinical_snomed_between(hba1c_measurement_snomed, dataset.elig_date_t2dm - days(2*366), dataset.elig_date_t2dm)

## Count lifestyle advice discussions, within previous 2 years, on or before elig_date_t2dm
dataset.cov_num_countlifestyle = count_matching_event_clinical_snomed_between(lifestyle_advice_snomed, dataset.elig_date_t2dm - days(2*366), dataset.elig_date_t2dm)


#######################################################################################
# Table 6) Outcomes and censoring events
#######################################################################################
### COVID-related HOSPITALIZAION
## First covid-19 related hospital admission, between elig_date_t2dm and studyend_date (incl. those dates)
dataset.out_date_covid_hes = first_matching_event_apc_between(covid_codes_incl_clin_diag, dataset.elig_date_t2dm, studyend_date, only_prim_diagnoses=True).admission_date
## First covid-19 related emergency attendance, between elig_date_t2dm and studyend_date (incl. those dates)
dataset.out_date_covid_emergency = first_matching_event_ec_snomed_between(covid_emergency, dataset.elig_date_t2dm, studyend_date).arrival_date
# combined: First covid-19 related hospitalization
dataset.out_date_covid_hosp = minimum_of(dataset.out_date_covid_hes, dataset.out_date_covid_emergency)

### COVID-related DEATH
## Between elig_date_t2dm and studyend_date (incl. those dates), stated anywhere on any of the 15 death certificate options
tmp_out_bin_death_covid = matching_death_between(covid_codes_incl_clin_diag, dataset.elig_date_t2dm, studyend_date)
dataset.out_date_covid_death = case(when(tmp_out_bin_death_covid).then(ons_deaths.date))

### COVID-related HOSPITALIZAION or DEATH
dataset.out_date_severecovid = minimum_of(dataset.out_date_covid_death, dataset.out_date_covid_hosp)

### COVID-related EVENT/DIAGNOSIS
## First COVID-19 diagnosis in primary care, between elig_date_t2dm and studyend_date (incl. those dates)
tmp_covid19_primary_care_date = first_matching_event_clinical_ctv3_between(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae, dataset.elig_date_t2dm, studyend_date).date
## First positive SARS-COV-2 PCR in primary care, between elig_date_t2dm and studyend_date (incl. those dates)
tmp_covid19_sgss_date = (
  sgss_covid_all_tests.where(sgss_covid_all_tests.specimen_taken_date.is_on_or_between(dataset.elig_date_t2dm, studyend_date))
  .where(sgss_covid_all_tests.is_positive)
  .sort_by(sgss_covid_all_tests.specimen_taken_date)
  .first_for_patient()
  .specimen_taken_date
)
## First COVID-19 diagnosis in primary care, or pos. test in primary care, or covid-19 hosp, or covid-19 death, between elig_date_t2dm and studyend_date (incl. those dates)
dataset.out_date_covid = minimum_of(tmp_covid19_primary_care_date, tmp_covid19_sgss_date, dataset.out_date_covid_hosp, dataset.out_date_covid_death)


### LONG COVID
## First Long COVID code in primary care, between elig_date_t2dm and studyend_date (incl. those dates), based on https://github.com/opensafely/long-covid/blob/main/analysis/codelists.py
dataset.out_date_longcovid = first_matching_event_clinical_snomed_between(long_covid_diagnostic_snomed_clinical + long_covid_referral_snomed_clinical + long_covid_assessment_snomed_clinical, dataset.elig_date_t2dm, studyend_date).date
## First Viral fatigue code in primary care, between elig_date_t2dm and studyend_date (incl. those dates), based on https://github.com/opensafely/long-covid/blob/main/analysis/codelists.py
dataset.out_date_virfat = first_matching_event_clinical_snomed_between(post_viral_fatigue_snomed_clinical, dataset.elig_date_t2dm, studyend_date).date
# combined: First Long COVID code or Viral Fatigue code after elig_date_t2dm = primary outcome
dataset.out_date_longcovid_virfat = minimum_of(dataset.out_date_longcovid, dataset.out_date_virfat)

## Long COVID signs and symptoms, first code in primary care, between elig_date_t2dm and studyend_date (incl. those dates)
# dataset.out_date_breathlessness = first_matching_event_clinical_snomed_between(breathlessness_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_cough = first_matching_event_clinical_snomed_between(cough_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_chest_pain = first_matching_event_clinical_snomed_between(chest_pain_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_chest_tightness = first_matching_event_clinical_snomed_between(chest_tightness_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_palpitations = first_matching_event_clinical_snomed_between(palpitations_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_fatigue = first_matching_event_clinical_snomed_between(fatigue_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_fever = first_matching_event_clinical_snomed_between(fever_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_pain = first_matching_event_clinical_snomed_between(pain_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_cog_impair = first_matching_event_clinical_snomed_between(cog_impair_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_headache = first_matching_event_clinical_snomed_between(headache_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_sleep = first_matching_event_clinical_snomed_between(sleep_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_pnp = first_matching_event_clinical_snomed_between(pnp_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_dizziness = first_matching_event_clinical_snomed_between(dizziness_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_delirium = first_matching_event_clinical_snomed_between(delirium_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_mob_impair = first_matching_event_clinical_snomed_between(mob_impair_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_visual = first_matching_event_clinical_snomed_between(visual_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_abdo_pain = first_matching_event_clinical_snomed_between(abdo_pain_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_nausea_vomiting = first_matching_event_clinical_snomed_between(nausea_vomiting_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_diarrhoea = first_matching_event_clinical_snomed_between(diarrhoea_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_weight_appetite = first_matching_event_clinical_snomed_between(weight_appetite_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_tinnitus = first_matching_event_clinical_snomed_between(tinnitus_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_earache = first_matching_event_clinical_snomed_between(earache_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_sore_throat = first_matching_event_clinical_snomed_between(sore_throat_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_smell_taste = first_matching_event_clinical_snomed_between(smell_taste_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_nasal = first_matching_event_clinical_snomed_between(nasal_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_hair_loss = first_matching_event_clinical_snomed_between(hair_loss_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_skin_rash = first_matching_event_clinical_snomed_between(skin_rash_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_anxiety = first_matching_event_clinical_snomed_between(anxiety_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_depression = first_matching_event_clinical_snomed_between(depression_snomed, dataset.elig_date_t2dm, studyend_date).date
# dataset.out_date_ptsd = first_matching_event_clinical_snomed_between(ptsd_snomed, dataset.elig_date_t2dm, studyend_date).date


### UPDATED eligibility and intercurrent events for potential censoring
## Practice deregistration date: Based on main registration at t2dm diagnosis date
# However, it does count those who only switch TPP practices
deregistered = spanning_regs.sort_by(
    practice_registrations.end_date,
    practice_registrations.practice_pseudo_id,
).first_for_patient()
dataset.cens_date_dereg = deregistered.end_date

## Known hypersensitivity / intolerance to metformin, after elig_date_t2dm
# dataset.cens_date_metfin_allergy_first = first_matching_event_clinical_snomed_between(metformin_allergy_snomed_clinical, dataset.elig_date_t2dm + days(1), studyend_date).date

## Moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; stage 4/5), after elig_date_t2dm
# dataset.cens_date_ckd_45_first = minimum_of(
#     first_matching_event_clinical_snomed_between(ckd_snomed_clinical_45, dataset.elig_date_t2dm + days(1), studyend_date).date,
#     first_matching_event_apc_between(ckd_stage4_icd10 + ckd_stage5_icd10, dataset.elig_date_t2dm + days(1), studyend_date).admission_date
# )

## Advance decompensated liver cirrhosis, after elig_date_t2dm
# dataset.cens_date_liver_cirrhosis_first = minimum_of(
#     first_matching_event_clinical_snomed_between(advanced_decompensated_cirrhosis_snomed_clinical + ascitic_drainage_snomed_clinical, dataset.elig_date_t2dm + days(1), studyend_date).date,
#     first_matching_event_apc_between(advanced_decompensated_cirrhosis_icd10, dataset.elig_date_t2dm + days(1), studyend_date).admission_date
# )

## Use of the following medications (drug-drug interaction with metformin)
# dataset.cens_date_metfin_interaction_first = first_matching_med_dmd_between(metformin_interaction_dmd, dataset.elig_date_t2dm + days(1), studyend_date).date


### Sensitivity analyses, neg & pos control
## Pos control: Diabetes complications (foot, retino, neuro, nephro), after elig_date_t2dm
# dataset.out_date_diabetescomp = minimum_of(
#     first_matching_event_clinical_snomed_between(diabetescomp_snomed_clinical, dataset.elig_date_t2dm + days(1), studyend_date).date,
#     first_matching_event_apc_between(diabetescomp_icd10, dataset.elig_date_t2dm + days(1), studyend_date).admission_date
# )
## Pos control: Diabetes-related death, after elig_date_t2dm, stated anywhere on any of the 15 death certificate options
dataset.tmp_out_bin_dm_death = matching_death_between(diabetes_type2_icd10, dataset.elig_date_t2dm, studyend_date)
dataset.out_date_dm_death = case(when(dataset.tmp_out_bin_dm_death).then(ons_deaths.date))
## Neg control: Fracture, after elig_date_t2dm
dataset.out_date_fracture = first_matching_event_apc_between(fracture_icd10, dataset.elig_date_t2dm, studyend_date).admission_date

## Flag patients on diet intervention only, no drug intervention, between elig_date_t2dm and studyend_date (incl. those dates)
dataset.tmp_date_diet_only = first_matching_event_clinical_snomed_between(diet_only_snomed, dataset.elig_date_t2dm, studyend_date).date