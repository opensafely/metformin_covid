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
    appointments,
    vaccinations
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

## json (for the dates)
import json

## for import of diabetes algo created data
from datetime import date

## import diabetes algo created data
@table_from_file("output/data_processed.csv.gz")
class data_processed(PatientFrame):
    ethnicity_cat = Series(str)
    t2dm_date = Series(date)

# random seed
import random
random.seed(19283) # random seed

#######################################################################################
# INITIALISE the dataset and set the dummy dataset size
#######################################################################################
dataset = create_dataset()
dataset.configure_dummy_data(population_size=8000)
dataset.define_population(patients.exists_for_patient())

#######################################################################################
# DEFINE the dates
#######################################################################################
dataset.elig_date_t2dm = data_processed.t2dm_date
with open("output/study_dates.json") as f:
  study_dates = json.load(f)
studyend_date = study_dates["studyend_date"]
pandemicstart_date = study_dates["pandemicstart_date"]

#######################################################################################
# Table 2) QUALITY ASSURANCES and completeness criteria
#######################################################################################
dataset.qa_num_birth_year = patients.date_of_birth
dataset.qa_bin_is_female_or_male = patients.sex.is_in(["female", "male"]) 
dataset.qa_bin_was_adult = (patients.age_on(pandemicstart_date) >= 18) & (patients.age_on(pandemicstart_date) <= 110) 
dataset.qa_bin_was_alive = patients.is_alive_on(pandemicstart_date)
dataset.qa_bin_known_imd = addresses.for_patient_on(pandemicstart_date).exists_for_patient() # known deprivation
dataset.qa_bin_was_registered = practice_registrations.spanning(pandemicstart_date - days(366), pandemicstart_date).exists_for_patient() # see https://docs.opensafely.org/ehrql/reference/schemas/tpp/#practice_registrations.spanning. Calculated from 1 year = 365.25 days, taking into account leap year.

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
# type 2 diabetes defined through diabetes-algo
###

### CATEGORY (1): "first ever & time-fixed" eigibility

## first known hypersensitivity / intolerance to metformin, on or before studyend_date
dataset.elig_date_metfin_allergy_first = first_matching_event_clinical_snomed_before(metformin_allergy_snomed_clinical, studyend_date).date

## first moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; stage 4/5), on or before studyend_date
dataset.elig_date_ckd_45_first = minimum_of(
    first_matching_event_clinical_snomed_before(ckd_snomed_clinical_45, studyend_date).date,
    first_matching_event_apc_before(ckd_stage4_icd10 + ckd_stage5_icd10, studyend_date).admission_date
)

## first advanced decompensated liver cirrhosis, on or before studyend_date
dataset.elig_date_liver_cirrhosis_first = minimum_of(
    first_matching_event_clinical_snomed_before(advanced_decompensated_cirrhosis_snomed_clinical + ascitic_drainage_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(advanced_decompensated_cirrhosis_icd10, studyend_date).admission_date
)

## first use of the following medications (drug-drug interaction with metformin), on or before studyend_date
## we assume someone starting any of these medications remains on the medication, i.e. continues to be ineligible after first ever event
## if we only care about e.g. use in 2 weeks window prior to trial start, then extract ELD, to apply windows at the start of each trial (see below)
dataset.elig_date_metfin_interaction_first = first_matching_med_dmd_before(metformin_interaction_dmd, studyend_date).date

### CATEGORY (2): "time-updated" eligibility -> use ELD
# medication table
medication_events = (
   medications
   .where(medications.date.is_after(pandemicstart_date))
   .where(medications.date.is_on_or_before(studyend_date))
   .sort_by(medications.date)
)
metfin_interactions = medication_events.where(medication_events.dmd_code.is_in(metformin_interaction_dmd))
dataset.add_event_table("metfin_interactions", date = metfin_interactions.date, eld_elig_bin_metfin_interaction = metfin_interactions.dmd_code.is_in(metformin_interaction_dmd))


#######################################################################################
# Table 4) INTERVENTION/EXPOSURE variables
#######################################################################################
## First metformin prescription ever
dataset.exp_date_metfin_first = first_matching_med_dmd_before(metformin_dmd, studyend_date).date
dataset.exp_date_metfin_mono_first = first_matching_med_dmd_before(metformin_dmd, studyend_date).date
## First other antidiabetic prescription ever
dataset.exp_date_sulfo_first = first_matching_med_dmd_before(sulfonylurea_dmd, studyend_date).date 
dataset.exp_date_dpp4_first = first_matching_med_dmd_before(dpp4_dmd, studyend_date).date
dataset.exp_date_tzd_first = first_matching_med_dmd_before(tzd_dmd, studyend_date).date 
dataset.exp_date_sglt2_first = first_matching_med_dmd_before(sglt2_dmd, studyend_date).date 
dataset.exp_date_glp1_first = first_matching_med_dmd_before(glp1_dmd, studyend_date).date 
dataset.exp_date_megli_first = first_matching_med_dmd_before(meglitinides_dmd, studyend_date).date
dataset.exp_date_agi_first = first_matching_med_dmd_before(agi_dmd, studyend_date).date 
dataset.exp_date_insulin_first = first_matching_med_dmd_before(insulin_dmd, studyend_date).date
## we assume someone starting any of these medications remains on the medication
## if we need more detailed on/off prescription patterns (though hard to define), then extract medications as ELD. Not for now.


#######################################################################################
# Table 5) Demographics, covariates and potential confounders
#######################################################################################

### CATEGORY (1): "first ever & time-fixed" covariates (or/and anchor them at pandemic start date)
## Age at pandemicstart_date
dataset.cov_num_age = patients.age_on(pandemicstart_date)

## Sex
dataset.cov_cat_sex = patients.sex

## Ethnicity (import from the diabetes algo data)
dataset.cov_cat_ethnicity = data_processed.ethnicity_cat

## Index of Multiple Deprevation Rank (rounded down to nearest 100). 5 categories.
imd_rounded = addresses.for_patient_on(pandemicstart_date).imd_rounded
dataset.cov_cat_deprivation_5 = case(
    when((imd_rounded >=0) & (imd_rounded < int(32844 * 1 / 5))).then("1 (most deprived)"),
    when(imd_rounded < int(32844 * 2 / 5)).then("2"),
    when(imd_rounded < int(32844 * 3 / 5)).then("3"),
    when(imd_rounded < int(32844 * 4 / 5)).then("4"),
    when(imd_rounded < int(32844 * 5 / 5)).then("5 (least deprived)"),
    otherwise="Unknown"
)

## Practice registration info at pandemicstart_date
# use a mix between spanning (as per eligibility criteria) and for_patient_on() to sort the multiple rows: https://docs.opensafely.org/ehrql/reference/schemas/tpp/#practice_registrations.for_patient_on
spanning_regs = practice_registrations.spanning(pandemicstart_date - days(366), pandemicstart_date)
registered = spanning_regs.sort_by(
    practice_registrations.end_date,
    practice_registrations.practice_pseudo_id,
).last_for_patient()
dataset.strat_cat_region = registered.practice_nuts1_region_name 
dataset.cov_cat_rural_urban = addresses.for_patient_on(pandemicstart_date).rural_urban_classification 

## Healthcare worker at the time they received a COVID-19 vaccination (which is anyway not really an updated variable -> we assume stable)
dataset.cov_bin_healthcare_worker = (
  occupation_on_covid_vaccine_record.where(
    occupation_on_covid_vaccine_record.is_healthcare_worker == True)
    .exists_for_patient()
)

## First (if any) acute myocardial infarction, on or before study end date
dataset.cov_bin_ami = (
    first_matching_event_clinical_snomed_before(ami_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(ami_prior_icd10 + ami_icd10, studyend_date).exists_for_patient()
)
## First (if any) stroke, on or before study end date
dataset.cov_bin_all_stroke = (
    first_matching_event_clinical_snomed_before(stroke_isch_snomed_clinical + stroke_sah_hs_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(stroke_isch_icd10 + stroke_sah_hs_icd10, studyend_date).exists_for_patient()
)

## Other arterial embolism, on or before study end date
dataset.cov_bin_other_arterial_embolism = (
   first_matching_event_clinical_snomed_before(other_arterial_embolism_snomed_clinical, studyend_date).exists_for_patient() |
   first_matching_event_apc_before(other_arterial_embolism_icd10, studyend_date).exists_for_patient()
)

## Venous thrombolism events, on or before study end date
dataset.cov_bin_vte = (
  first_matching_event_clinical_snomed_before(portal_vein_thrombosis_snomed_clinical + dvt_dvt_snomed_clinical + dvt_icvt_snomed_clinical + dvt_pregnancy_snomed_clinical + other_dvt_snomed_clinical + pe_snomed_clinical, studyend_date).exists_for_patient() |
  first_matching_event_apc_before(portal_vein_thrombosis_icd10 + dvt_dvt_icd10 + dvt_icvt_icd10 + dvt_pregnancy_icd10 + other_dvt_icd10 + icvt_pregnancy_icd10 + pe_icd10, studyend_date).exists_for_patient()
)

## Heart failure, on or before study end date
dataset.cov_bin_hf = (
    first_matching_event_clinical_snomed_before(hf_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(hf_icd10, studyend_date).exists_for_patient()
)

## Angina, on or before study end date
dataset.cov_bin_angina = (
    first_matching_event_clinical_snomed_before(angina_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(angina_icd10, studyend_date).exists_for_patient()
)

## Dementia, on or before study end date
dataset.cov_bin_dementia = (
    first_matching_event_clinical_snomed_before(dementia_snomed_clinical + dementia_vascular_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(dementia_icd10 + dementia_vascular_icd10, studyend_date).exists_for_patient()
)

## Cancer, on or before study end date
dataset.cov_bin_cancer = (
    first_matching_event_clinical_snomed_before(cancer_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(cancer_icd10, studyend_date).exists_for_patient()
)

## Hypertension, on or before study end date
dataset.cov_bin_hypertension = (
    first_matching_event_clinical_snomed_before(hypertension_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(hypertension_icd10, studyend_date).exists_for_patient()
)

## Depression, on or before study end date
dataset.cov_bin_depression = (
    first_matching_event_clinical_snomed_before(depression_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(depression_icd10, studyend_date).exists_for_patient()
)

## Chronic obstructive pulmonary disease, on or before study end date
dataset.cov_bin_copd = (
    first_matching_event_clinical_snomed_before(copd_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(copd_icd10, studyend_date).exists_for_patient()
)

## Liver disease,on or before study end date
dataset.cov_bin_liver_disease = (
    first_matching_event_clinical_snomed_before(liver_disease_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(liver_disease_icd10, studyend_date).exists_for_patient()
)

## Chronic kidney disease, on or before study end date
dataset.cov_bin_chronic_kidney_disease = (
    first_matching_event_clinical_snomed_before(ckd_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(ckd_icd10, studyend_date).exists_for_patient()
)

## PCOS, on or before study end date
dataset.cov_bin_pcos = (
    first_matching_event_clinical_snomed_before(pcos_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(pcos_icd10, studyend_date).exists_for_patient()
)

## Prediabetes, on or before study end date 
# Any preDM diagnosis in primary care
tmp_cov_bin_prediabetes = first_matching_event_clinical_snomed_before(prediabetes_snomed, studyend_date).exists_for_patient()
# Any HbA1c preDM in primary care || reconsider! Does not make much sense, would need to combine with ELD
#tmp_cov_bin_predm_hba1c_mmol_mol = (
#  clinical_events.where(
#    clinical_events.snomedct_code.is_in(hba1c_snomed))
#    .where(clinical_events.date.is_on_or_before(studyend_date))
#    .where((clinical_events.numeric_value>=42) & (clinical_events.numeric_value<=47.9))
#    .exists_for_patient()
#)
# Any preDM diagnosis or Hb1Ac preDM range value (in period before elig_date_t2dm)
dataset.cov_bin_prediabetes = tmp_cov_bin_prediabetes

## Diabetes complications (foot, retino, neuro, nephro), on or before study end date
dataset.cov_bin_diabetescomp = (
    first_matching_event_clinical_snomed_before(diabetescomp_snomed_clinical, studyend_date).exists_for_patient() |
    first_matching_event_apc_before(diabetescomp_icd10, studyend_date).exists_for_patient()
)


### CATEGORY (2): "time-fixed at pandemic start" covariates, but we might update them later, either with or without ELD
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

## Consultation rate in previous year as a proxy for health seeking behaviour
# Consultation rate in 2019, pre-pandemic
tmp_cov_num_consrate = appointments.where(
    appointments.status.is_in([
        "Arrived",
        "In Progress",
        "Finished",
        "Visit",
        "Waiting",
        "Patient Walked Out",
        ]) & appointments.start_date.is_on_or_between("2019-02-01", "2020-02-01")
        ).count_for_patient()    

dataset.cov_num_consrate = case(
    when(tmp_cov_num_consrate <= 365).then(tmp_cov_num_consrate),
    otherwise=365, # quality assurance
)


### CATEGORY (3): TIME-VARYING covariates (they may reduce/increase) between pandemic start (= study start) and study end
## FIRST, define the baseline value at pandemic start to use for first trial
# Obesity, on or before pandemic start
dataset.cov_bin_obesity_b = (
    last_matching_event_clinical_snomed_before(bmi_obesity_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(bmi_obesity_icd10, pandemicstart_date).exists_for_patient()
)

# BMI, most recent value, within previous 2 years, on or before pandemic start
bmi_measurement = most_recent_bmi(
    where=clinical_events.date.is_on_or_between(pandemicstart_date - days(2 * 366), pandemicstart_date),
    minimum_age_at_measurement=16,
)
dataset.cov_num_bmi_b = bmi_measurement.numeric_value
dataset.cov_cat_bmi_groups_b = case(
    when((dataset.cov_num_bmi_b < 18.5) & (dataset.cov_num_bmi_b >= 12.0)).then("Underweight"), # Set minimum to avoid any impossibly extreme values being classified as underweight
    when((dataset.cov_num_bmi_b >= 18.5) & (dataset.cov_num_bmi_b < 25.0)).then("Healthy weight (18.5-24.9)"),
    when((dataset.cov_num_bmi_b >= 25.0) & (dataset.cov_num_bmi_b < 30.0)).then("Overweight (25-29.9)"),
    when((dataset.cov_num_bmi_b >= 30.0) & (dataset.cov_num_bmi_b <= 70.0)).then("Obese (>30)"), # Set maximum to avoid any impossibly extreme values being classified as obese
    otherwise = "missing", 
)

# HbA1c, most recent value, within previous 2 years, on or before pandemic start
dataset.cov_num_hba1c_mmol_mol_b = last_matching_event_clinical_snomed_between(hba1c_snomed, pandemicstart_date - days(2*366), pandemicstart_date).numeric_value 

# Total Cholesterol, most recent value, within previous 2 years, on or before pandemic start
dataset.tmp_cov_num_cholesterol_b = last_matching_event_clinical_snomed_between(cholesterol_snomed, pandemicstart_date - days(2*366), pandemicstart_date).numeric_value 

# HDL Cholesterol, most recent value, within previous 2 years, on or before pandemic start
dataset.tmp_cov_num_hdl_cholesterol_b = last_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, pandemicstart_date - days(2*366), pandemicstart_date).numeric_value 

## SECOND, use event-level data to capture follow-up data for certain covariates, after pandemic start. And add vaccination data.
# All SARS-CoV-2 vaccinations someone received during study period; to dynamically adjust the treatment weights for no vs any vs multiple vaccinations at the start of each seq trial
covid_vaccinations = (
  vaccinations
  .where(vaccinations.target_disease.is_in(["SARS-2 CORONAVIRUS"]))
  .where(vaccinations.date.is_on_or_between(pandemicstart_date, studyend_date))
  .sort_by(vaccinations.date)
)
dataset.add_event_table(
  "covid_vaccinations",
  date=covid_vaccinations.date,
  eld_cov_cat_vacc = covid_vaccinations.product_name # not sure I need it, currently for testing purpose
)

## Obesity, HbA1c, lipids recordings during study period
## to dynamically adjust the treatment weights at the start of each sequential trial and inform the censoring weights for the per protocol analysis (start of metformin in control)
# 1) clinical_events table, i.e. primary care
pc_events = (
   clinical_events
   .where(clinical_events.date.is_after(pandemicstart_date))
   .where(clinical_events.date.is_on_or_before(studyend_date))
   .sort_by(clinical_events.date)
)
# 2) apcs table, i.e. secondary care
sc_events = (
   apcs
   .where(apcs.admission_date.is_after(pandemicstart_date))
   .where(apcs.admission_date.is_on_or_before(studyend_date))
   .sort_by(apcs.admission_date)
)

# Obesity diagnosis, from primary and secondary care
# Consider combining obesity with BMI measurement as above?
obesity_events_pc = pc_events.where(pc_events.snomedct_code.is_in(bmi_obesity_snomed_clinical))
dataset.add_event_table("obesity_pc", date = obesity_events_pc.date, eld_cov_bin_obesity_pc = obesity_events_pc.snomedct_code.is_in(bmi_obesity_snomed_clinical))
obesity_events_sc = sc_events.where(sc_events.all_diagnoses.contains_any_of(bmi_obesity_icd10))
dataset.add_event_table("obesity_sc", date = obesity_events_sc.admission_date, eld_cov_bin_obesity_sc = obesity_events_sc.all_diagnoses.contains_any_of(bmi_obesity_icd10))

# BMI, from primary care
#bmi_events = bmi_eld(
#   where=clinical_events.date.is_on_or_between(pandemicstart_date, studyend_date), 
#   minimum_age_at_measurement=16)
bmi_events = bmi_eld(
   where=(clinical_events.date.is_after(pandemicstart_date) & clinical_events.date.is_on_or_before(studyend_date)),
   minimum_age_at_measurement=16
   )
dataset.add_event_table("bmi", date = bmi_events.date, eld_cov_num_bmi = bmi_events.numeric_value)

# HbA1c, from primary care
hba1c_events = pc_events.where(pc_events.snomedct_code.is_in(hba1c_snomed))
dataset.add_event_table("hba1c", date = hba1c_events.date, eld_cov_num_hba1c = hba1c_events.numeric_value)

# Total Cholesterol, from primary care
chol_events = pc_events.where(pc_events.snomedct_code.is_in(cholesterol_snomed))
dataset.add_event_table("chol", date = chol_events.date, eld_cov_num_chol = chol_events.numeric_value)

# HDL Cholesterol, from primary care
hdl_events = pc_events.where(pc_events.snomedct_code.is_in(hdl_cholesterol_snomed))
dataset.add_event_table("hdl", date = hdl_events.date, eld_cov_num_hdl = hdl_events.numeric_value)


#######################################################################################
# Table 5) Outcomes (also important for time-updated eligibility), and censoring events
#######################################################################################

### CATEGORY (1): "first ever & time-fixed" outcomes
## COVID-related DEATH (all-cause deaths already defined above in QA)
# Between pandemicstart_date and studyend_date (incl. those dates), stated anywhere on any of the 15 death certificate options
tmp_out_bin_death_covid = matching_death_between(covid_codes_incl_clin_diag, pandemicstart_date, studyend_date)
dataset.out_date_covid_death = case(when(tmp_out_bin_death_covid).then(ons_deaths.date))

## Composite: first COVID-related DEATH or COVID-related hospitalization
tmp_out_date_covid_hes = first_matching_event_apc_between(covid_codes_incl_clin_diag, pandemicstart_date, studyend_date, only_prim_diagnoses=True).admission_date
tmp_out_date_covid_emergency = first_matching_event_ec_snomed_between(covid_emergency, pandemicstart_date, studyend_date).arrival_date
dataset.out_date_severecovid = minimum_of(dataset.out_date_covid_death, tmp_out_date_covid_hes, tmp_out_date_covid_emergency)

## Deregistration, for censoring for LTFU
dataset.cens_date_dereg = registered.end_date


### CATEGORY (2): TIME-VARYING outcomes (they may reduce/increase) between pandemic start (= study start) and study end
### they are also relevant for time-updated eligibility to be included in SeqTrials => use ELD

## COVID-related hospitalizations, diagnosis recorded only in primary diagnosis field
covid_hosp = sc_events.where(sc_events.primary_diagnosis.is_in(covid_codes_incl_clin_diag))
dataset.add_event_table("covid_hosp", date = covid_hosp.admission_date, eld_out_bin_covid_hosp = covid_hosp.primary_diagnosis.is_in(covid_codes_incl_clin_diag))

## COVID-related emergency care admissions
covid_ec = ec_eld(covid_emergency, pandemicstart_date, studyend_date)
dataset.add_event_table(
  "covid_ec",
  date = covid_ec.arrival_date
)

## COVID pos. tests
covid_sgss = (
   sgss_covid_all_tests
   .where(sgss_covid_all_tests.specimen_taken_date.is_after(pandemicstart_date))
   .where(sgss_covid_all_tests.specimen_taken_date.is_on_or_before(studyend_date))
   .where(sgss_covid_all_tests.is_positive)
   .sort_by(sgss_covid_all_tests.specimen_taken_date)
)
dataset.add_event_table(
  "covid_sgss",
  date = covid_sgss.specimen_taken_date,
  eld_out_bin_covid_test = covid_sgss.is_positive
)

## COVID diagnosis/events in primary care
covid_pc = pc_events.where(pc_events.ctv3_code.is_in(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae))
dataset.add_event_table("covid_pc", date = covid_pc.date, eld_out_bin_covid_pc = covid_pc.ctv3_code.is_in(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae))