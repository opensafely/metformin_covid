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
mid2018_date = study_dates["mid2018_date"]

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

### FIRST, define the "first ever" (stable) events
## first known hypersensitivity / intolerance to metformin, on or before studyend_date
dataset.elig_date_metfin_allergy_first = first_matching_event_clinical_snomed_before(metformin_allergy_snomed_clinical, studyend_date).date

## first moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; stage 4/5), on or before studyend_date
dataset.elig_date_ckd_45_first = minimum_of(
    first_matching_event_clinical_snomed_before(ckd_snomed_clinical_45, studyend_date).date,
    first_matching_event_apc_before(ckd_stage4_icd10 + ckd_stage5_icd10, studyend_date).admission_date
)

## first advance decompensated liver cirrhosis, on or before studyend_date
dataset.elig_date_liver_cirrhosis_first = minimum_of(
    first_matching_event_clinical_snomed_before(advanced_decompensated_cirrhosis_snomed_clinical + ascitic_drainage_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(advanced_decompensated_cirrhosis_icd10, studyend_date).admission_date
)

## First use of the following medications (drug-drug interaction with metformin), on or before studyend_date
## we assume someone starting any of these medications remains on the medication, i.e. continues to be ineligible after first ever event
## if we only care about e.g. use in 2 weeks window prior to trial start, then extract ELD, to apply windows at the start of each trial
dataset.elig_date_metfin_interaction_first = first_matching_med_dmd_before(metformin_interaction_dmd, studyend_date).date


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
## if we need more detailed on/off prescription patterns (though hard to define), then extract medications as ELD


#######################################################################################
# Table 5) Demographics, covariates and potential confounders
#######################################################################################
### CATEGORY 1) Demographics measured at pandemic start and not updated thereafter
# Age at pandemicstart_date
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

## Healthcare worker at the time they received a COVID-19 vaccination
dataset.cov_bin_healthcare_worker = (
  occupation_on_covid_vaccine_record.where(
    occupation_on_covid_vaccine_record.is_healthcare_worker == True)
    .exists_for_patient()
)

### CATEGORY 2) Measured at pandemic start but we might update them later (to have time-updated baseline info for sequential trial baseline adjustment), not sure ELD is needed, though
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
# Consultation rate in 2019
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

### CATEGORY 3) ONLY ONE-TIME (i.e., first ever ) recording is important, since these conditions are absorbing conditions (they might increase, however, but that's not important for now).
### PS: reduced to 2 clinical diagnoses only, there would be many more, but to reduce overloading testing with ELD
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

### CATEGORY 4) TIME-VARYING conditions (they may reduce/increase) between pandemic start (= study start) and study end
## FIRST, just measured value at pandemic start to use for first trial
# Obesity, on or before pandemic start
dataset.cov_bin_obesity = (
    last_matching_event_clinical_snomed_before(bmi_obesity_snomed_clinical, pandemicstart_date).exists_for_patient() |
    last_matching_event_apc_before(bmi_obesity_icd10, pandemicstart_date).exists_for_patient()
)

# BMI, most recent value, within previous 2 years, on or before pandemic start
bmi_measurement = most_recent_bmi(
    where=clinical_events.date.is_on_or_between(pandemicstart_date - days(2 * 366), pandemicstart_date),
    minimum_age_at_measurement=16,
)
dataset.cov_num_bmi = bmi_measurement.numeric_value
dataset.cov_cat_bmi_groups = case(
    when((dataset.cov_num_bmi < 18.5) & (dataset.cov_num_bmi >= 12.0)).then("Underweight"), # Set minimum to avoid any impossibly extreme values being classified as underweight
    when((dataset.cov_num_bmi >= 18.5) & (dataset.cov_num_bmi < 25.0)).then("Healthy weight (18.5-24.9)"),
    when((dataset.cov_num_bmi >= 25.0) & (dataset.cov_num_bmi < 30.0)).then("Overweight (25-29.9)"),
    when((dataset.cov_num_bmi >= 30.0) & (dataset.cov_num_bmi <= 70.0)).then("Obese (>30)"), # Set maximum to avoid any impossibly extreme values being classified as obese
    otherwise = "missing", 
)

# HbA1c, most recent value, within previous 2 years, on or before pandemic start
dataset.cov_num_hba1c_mmol_mol = last_matching_event_clinical_snomed_between(hba1c_snomed, pandemicstart_date - days(2*366), pandemicstart_date).numeric_value 

# Total Cholesterol, most recent value, within previous 2 years, on or before pandemic start
dataset.tmp_cov_num_cholesterol = last_matching_event_clinical_snomed_between(cholesterol_snomed, pandemicstart_date - days(2*366), pandemicstart_date).numeric_value 

# HDL Cholesterol, most recent value, within previous 2 years, on or before pandemic start
dataset.tmp_cov_num_hdl_cholesterol = last_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, pandemicstart_date - days(2*366), pandemicstart_date).numeric_value 

#######
## SECOND, use event-level data for the above CATEGORY 4) TIME-VARYING conditions - and now also include vaccination data. And integrated outcomes in ELD tables.
#######

## All SARS-CoV-2 vaccinations someone received during study period; to dynamically adjust the treatment weights for no vs any vs multiple vaccinations at the start of each seq trial
covid_vaccinations = (
  vaccinations
  .where(vaccinations.target_disease.is_in(["SARS-2 CORONAVIRUS"]))
  .where(vaccinations.date.is_on_or_between(pandemicstart_date, studyend_date))
  .sort_by(vaccinations.date)
)
dataset.add_event_table(
  "covid_vaccinations",
  date=covid_vaccinations.date,
  vaccine_product = covid_vaccinations.product_name # not sure I need it, just for testing purpose for now
)

## Obesity, HbA1c, lipids recordings during study period; to dynamically adjust the treatment weights at the start of each sequential trial and inform the censoring weights for the PPA (start of metformin in control)
## Include the time-varying outcomes here, too, if I search already in corresponding table
# Consider combining obesity with BMI measurement as above?

# 1) clinical_events table
pc_events = (
   clinical_events
   .where(clinical_events.date.is_after(pandemicstart_date))
   .where(clinical_events.date.is_on_or_before(studyend_date))
      .sort_by(clinical_events.date)
)

obesity_events = pc_events.where(pc_events.snomedct_code.is_in(bmi_obesity_snomed_clinical))
dataset.add_event_table("pc_obesity", date=obesity_events.date, pc_obesity = obesity_events.snomedct_code.is_in(bmi_obesity_snomed_clinical))

covid_events = pc_events.where(pc_events.ctv3_code.is_in(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae))
dataset.add_event_table("pc_covid", date=covid_events.date, pc_covid = covid_events.ctv3_code.is_in(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae))

hba1c_events = pc_events.where(pc_events.snomedct_code.is_in(hba1c_snomed))
dataset.add_event_table("hba1c", date=hba1c_events.date, hba1c_num=hba1c_events.numeric_value)

chol_events = pc_events.where(pc_events.snomedct_code.is_in(cholesterol_snomed))
dataset.add_event_table("chol", date=chol_events.date, chol_num=chol_events.numeric_value)

hdl_events = pc_events.where(pc_events.snomedct_code.is_in(hdl_cholesterol_snomed))
dataset.add_event_table("hdl", date=hdl_events.date, hdl_num=hdl_events.numeric_value)

#clinical_events_pc = (
#   clinical_events
#   .where(clinical_events.date.is_after(pandemicstart_date))
#   .where(clinical_events.date.is_on_or_before(studyend_date))
#   .where((clinical_events.snomedct_code.is_in(
#      bmi_obesity_snomed_clinical
#      + hba1c_snomed
#      + cholesterol_snomed
#      + hdl_cholesterol_snomed
#      )) | (clinical_events.ctv3_code.is_in(
#      covid_primary_care_code
#      + covid_primary_care_positive_test
#      + covid_primary_care_sequelae))
#      )
#      .sort_by(clinical_events.date)
#)
#dataset.add_event_table(
#  "clinical_events_pc",
#  date = clinical_events_pc.date,
#  obesity = clinical_events_pc.snomedct_code.is_in(bmi_obesity_snomed_clinical),
#  covid = clinical_events_pc.ctv3_code.is_in(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae),
#  hba1c = clinical_events_pc.snomedct_code.is_in(hba1c_snomed),
#  chol = clinical_events_pc.snomedct_code.is_in(cholesterol_snomed),
#  hdl = clinical_events_pc.snomedct_code.is_in(hdl_cholesterol_snomed),
#  num = clinical_events_pc.numeric_value
#)

# 2) apcs table, diagnosis recorded in any diagnosis field
apcs_obesity = (
   apcs
   .where(apcs.admission_date.is_after(pandemicstart_date))
   .where(apcs.admission_date.is_on_or_before(studyend_date))
   .where(apcs.all_diagnoses.contains_any_of(bmi_obesity_icd10))
   .sort_by(apcs.admission_date)
)
dataset.add_event_table(
  "sc_obesity",
  date=apcs_obesity.admission_date,
  sc_obesity=apcs_obesity.all_diagnoses.contains_any_of(bmi_obesity_icd10)
)

# 3) apcs table, diagnosis recorded only in primary diagnosis field
apcs_covid = (
   apcs
   .where(apcs.admission_date.is_after(pandemicstart_date))
   .where(apcs.admission_date.is_on_or_before(studyend_date))
   .where(apcs.primary_diagnosis.is_in(covid_codes_incl_clin_diag))
   .sort_by(apcs.admission_date)
)
dataset.add_event_table(
  "sc_covid",
  date=apcs_covid.admission_date,
  sc_covid=apcs_covid.primary_diagnosis.is_in(covid_codes_incl_clin_diag)
)

# 4) emergency_care_attendances table, diagnosis recorded in any diagnosis field
def ec_eld(codelist, start_date, end_date, where=True):
    conditions = [
        getattr(emergency_care_attendances, column_name).is_in(codelist)
        for column_name in ([f"diagnosis_{i:02d}" for i in range(1, 25)])
    ]
    return(
        emergency_care_attendances.where(where)
        .where(any_of(conditions))
        .where(emergency_care_attendances.arrival_date.is_after(start_date))
        .where(emergency_care_attendances.arrival_date.is_on_or_before(end_date))
        .sort_by(emergency_care_attendances.arrival_date)
    )

ec_covid = ec_eld(covid_emergency, pandemicstart_date, studyend_date)
dataset.add_event_table(
  "ec_covid",
  date=ec_covid.arrival_date
)

# 5) sgss_covid_all_tests table
sgss_covid = (
   sgss_covid_all_tests
   .where(sgss_covid_all_tests.specimen_taken_date.is_after(pandemicstart_date))
   .where(sgss_covid_all_tests.specimen_taken_date.is_on_or_before(studyend_date))
   .where(sgss_covid_all_tests.is_positive)
   .sort_by(sgss_covid_all_tests.specimen_taken_date)
)
dataset.add_event_table(
  "sgss_covid",
  date=sgss_covid.specimen_taken_date,
  sgss_covid=sgss_covid.is_positive
)

#######################################################################################
# Table 6) Outcomes (also important for time-updated eligibility), and censoring events
#######################################################################################
### CATEGORY 1) ONE-TIME conditions (absorbing conditions, no change after first measurement)
## COVID-related DEATH
# Between pandemicstart_date and studyend_date (incl. those dates), stated anywhere on any of the 15 death certificate options
tmp_out_bin_death_covid = matching_death_between(covid_codes_incl_clin_diag, pandemicstart_date, studyend_date)
dataset.out_date_covid_death = case(when(tmp_out_bin_death_covid).then(ons_deaths.date))

## Deregistration, for censoring
dataset.cens_date_dereg = registered.end_date

### CATEGORY 2) TIME-VARYING conditions; they may be present/absent - and be relevant for time-updated eligibility to be included in SeqTrials => use event-level data
#### All are above integrated in ELD tables ####


### COVID-related HOSPITALIZAION
## First covid-19 related hospital admission, between pandemicstart_date and studyend_date (incl. those dates)
#dataset.out_date_covid_hes = first_matching_event_apc_between(covid_codes_incl_clin_diag, pandemicstart_date, studyend_date, only_prim_diagnoses=True).admission_date
## First covid-19 related emergency attendance, between pandemicstart_date and studyend_date (incl. those dates)
#dataset.out_date_covid_emergency = first_matching_event_ec_snomed_between(covid_emergency, pandemicstart_date, studyend_date).arrival_date
# combined: First covid-19 related hospitalization
#dataset.out_date_covid_hosp = minimum_of(dataset.out_date_covid_hes, dataset.out_date_covid_emergency)

### COVID-related HOSPITALIZAION or DEATH
#dataset.out_date_severecovid = minimum_of(dataset.out_date_covid_death, dataset.out_date_covid_hosp)

### COVID-related EVENT/DIAGNOSIS
## First COVID-19 diagnosis in primary care, between elig_date_t2dm and studyend_date (incl. those dates)
#tmp_covid19_primary_care_date = first_matching_event_clinical_ctv3_between(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae,  pandemicstart_date, studyend_date).date
## First positive SARS-COV-2 PCR in primary care, between elig_date_t2dm and studyend_date (incl. those dates)
#tmp_covid19_sgss_date = (
#  sgss_covid_all_tests.where(sgss_covid_all_tests.specimen_taken_date.is_on_or_between(pandemicstart_date, studyend_date))
#  .where(sgss_covid_all_tests.is_positive)
#  .sort_by(sgss_covid_all_tests.specimen_taken_date)
#  .first_for_patient()
#  .specimen_taken_date
#)
## First COVID-19 diagnosis in primary care, or pos. test in primary care, or covid-19 hosp, or covid-19 death, between elig_date_t2dm and studyend_date (incl. those dates)
#dataset.out_date_covid = minimum_of(tmp_covid19_primary_care_date, tmp_covid19_sgss_date, dataset.out_date_covid_hosp, dataset.out_date_covid_death)