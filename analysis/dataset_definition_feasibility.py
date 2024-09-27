#######################################################################################
# IMPORT
#######################################################################################

## ehrQL functions
from ehrql import (
    case,
    create_dataset,
    when,
    minimum_of,
    maximum_of
)

## TPP tables
from ehrql.tables.tpp import (
    clinical_events,
    medications,
    patients,
    practice_registrations,
    ons_deaths,
    sgss_covid_all_tests,
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
# DEFINE the feasibility study end date
#######################################################################################
with open("output/study_dates.json") as f:
  study_dates = json.load(f)
feasibilityend_date = study_dates["feasibilityend_date"]

#######################################################################################
# INITIALISE the dataset and set the dummy dataset size
#######################################################################################
dataset = create_dataset()
dataset.configure_dummy_data(population_size=10000)
dataset.define_population(patients.exists_for_patient())

#######################################################################################
# DEMOGRAPHIC variables
#######################################################################################
## Sex
dataset.cov_cat_sex = patients.sex.is_in(["female", "male"])
## Age
dataset.cov_num_age = patients.age_on(feasibilityend_date)
## Registration info at feasibilityend_date
registered = practice_registrations.for_patient_on(feasibilityend_date)
## Region
dataset.cov_cat_region = registered.practice_nuts1_region_name
## Practice
dataset.cov_cat_stp = registered.practice_stp

#######################################################################################
# COVID-19 variables: Only COVID diagnoses (clinical + tests) and 2 outcomes (COVID deaths and Long COVID)
#######################################################################################
## First COVID-19 diagnosis in primary care before feasibility end date 
tmp_covid19_primary_care_date = first_matching_event_clinical_ctv3_before(covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae, feasibilityend_date).date
## First positive SARS-COV-2 PCR before feasibility end date
tmp_covid19_sgss_date = (
  sgss_covid_all_tests.where(sgss_covid_all_tests.specimen_taken_date.is_on_or_before(feasibilityend_date))
  .where(sgss_covid_all_tests.is_positive)
  .sort_by(sgss_covid_all_tests.specimen_taken_date)
  .first_for_patient()
  .specimen_taken_date
)
dataset.cov_date_covid19 = minimum_of(tmp_covid19_primary_care_date, tmp_covid19_sgss_date)

## COVID-related deaths (stated anywhere on any of the 15 death certificate options), based on https://github.com/opensafely/comparative-booster-spring2023/blob/main/analysis/codelists.py uses a different codelist: codelists/opensafely-covid-identification.csv
tmp_out_bin_death_covid = matching_death_before(covid_codes_incl_clin_diag, feasibilityend_date)
dataset.out_date_death_covid = case(when(tmp_out_bin_death_covid).then(ons_deaths.date))

## Long COVID code in primary care, based on https://github.com/opensafely/long-covid/blob/main/analysis/codelists.py
dataset.out_date_long_covid_first = first_matching_event_clinical_snomed_before(long_covid_diagnostic_snomed_clinical + long_covid_referral_snomed_clinical + long_covid_assessment_snomed_clinical, feasibilityend_date).date
## Viral fatigue code in primary care
dataset.out_date_fatigue_first = first_matching_event_clinical_snomed_before(post_viral_fatigue_snomed_clinical, feasibilityend_date).date
# combined
dataset.out_date_long_fatigue_first = minimum_of(dataset.out_date_long_covid_first, dataset.out_date_fatigue_first)

#######################################################################################
# INTERVENTION/EXPOSURE variables
#######################################################################################
# METFORMIN, based on codelist https://www.opencodelists.org/codelist/user/john-tazare/metformin-dmd/48e43356/
dataset.exp_date_first_metfin = first_matching_med_dmd_before(metformin_codes_dmd, feasibilityend_date).date 
dataset.exp_count_metfin = (
  medications.where(
    medications.dmd_code.is_in(metformin_codes_dmd))
    .where(medications.date.is_on_or_before(feasibilityend_date))
    .count_for_patient()
)

#######################################################################################
# DIABETES variables: Simple code-based, and for diabetes_algo()
#######################################################################################

## Type 2 Diabetes (cov_date_t2dm) & Type 1 Diabetes (cov_date_t1dm) & Gestational Diabetes (cov_date_gestationaldm) defined in algo below
## PCOS
dataset.cov_date_pcos = maximum_of(
    last_matching_event_clinical_snomed_before(pcos_snomed_clinical, feasibilityend_date).date,
    last_matching_event_apc_before(pcos_icd10, feasibilityend_date).admission_date
)
## Prediabetes, on or before feasibilityend_date
# Date of preDM code in primary care
tmp_cov_date_prediabetes = last_matching_event_clinical_snomed_before(prediabetes_snomed, feasibilityend_date).date
# Date of preDM HbA1c measure in period before feasibilityend_date in preDM range (mmol/mol): 42-47.9
tmp_cov_date_predm_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(feasibilityend_date))
    .where((clinical_events.numeric_value>=42) & (clinical_events.numeric_value<=47.9))
    .sort_by(clinical_events.date)
    .last_for_patient()
    .date
)
# Latest date (in period before feasibilityend_date) that any prediabetes was diagnosed or HbA1c in preDM range
dataset.cov_date_prediabetes = maximum_of(tmp_cov_date_prediabetes, tmp_cov_date_predm_hba1c_mmol_mol) 

## DIABETES algo variables start ------------------------
## BASED on https://github.com/opensafely/post-covid-diabetes/blob/main/analysis/common_variables.py 

## Type 1 Diabetes 
# Latest date from primary+secondary, but also primary care date separately for diabetes algo
dataset.tmp_cov_date_t1dm_ctv3 = last_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, feasibilityend_date).date
dataset.cov_date_t1dm = minimum_of(
    (last_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, feasibilityend_date).date),
    (last_matching_event_apc_before(diabetes_type1_icd10, feasibilityend_date).admission_date)
)
# Count codes (individually and together, for diabetes algo)
dataset.tmp_cov_count_t1dm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, feasibilityend_date)
dataset.tmp_cov_count_t1dm_hes = count_matching_event_apc_before(diabetes_type1_icd10, feasibilityend_date)
dataset.tmp_cov_count_t1dm = dataset.tmp_cov_count_t1dm_ctv3 + dataset.tmp_cov_count_t1dm_hes

## Type 2 Diabetes
# Latest date from primary+secondary, but also primary care date separately for diabetes algo)
dataset.tmp_cov_date_t2dm_ctv3 = last_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, feasibilityend_date).date
dataset.cov_date_t2dm = minimum_of(
    (last_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, feasibilityend_date).date),
    (last_matching_event_apc_before(diabetes_type2_icd10, feasibilityend_date).admission_date)
)
# Count codes (individually and together, for diabetes algo)
dataset.tmp_cov_count_t2dm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, feasibilityend_date)
dataset.tmp_cov_count_t2dm_hes = count_matching_event_apc_before(diabetes_type2_icd10, feasibilityend_date)
dataset.tmp_cov_count_t2dm = dataset.tmp_cov_count_t2dm_ctv3 + dataset.tmp_cov_count_t2dm_hes

## Diabetes unspecified/other
# Latest date
dataset.cov_date_otherdm = last_matching_event_clinical_ctv3_before(diabetes_other_ctv3_clinical, feasibilityend_date).date
# Count codes
dataset.tmp_cov_count_otherdm = count_matching_event_clinical_ctv3_before(diabetes_other_ctv3_clinical, feasibilityend_date)

## Gestational diabetes
# Latest date
dataset.cov_date_gestationaldm = last_matching_event_clinical_ctv3_before(diabetes_gestational_ctv3_clinical, feasibilityend_date).date

## Diabetes diagnostic codes
# Latest date
dataset.tmp_cov_date_poccdm = last_matching_event_clinical_ctv3_before(diabetes_diagnostic_ctv3_clinical, feasibilityend_date).date
# Count codes
dataset.tmp_cov_count_poccdm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_diagnostic_ctv3_clinical, feasibilityend_date)

### Other variables needed to define diabetes
## HbA1c
# Maximum HbA1c measure (in period before feasibilityend_date) -> don't have helper function for this
dataset.tmp_cov_num_max_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(feasibilityend_date))
    .numeric_value.maximum_for_patient()
)
# Date of latest maximum HbA1c measure
dataset.tmp_cov_date_max_hba1c = ( 
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(feasibilityend_date)) # this line of code probably not needed again
    .where(clinical_events.numeric_value == dataset.tmp_cov_num_max_hba1c_mmol_mol)
    .sort_by(clinical_events.date)
    .last_for_patient() # translates in cohortextractor to "on_most_recent_day_of_measurement=True"
    .date
)

## Diabetes drugs
# Latest dates
dataset.tmp_cov_date_insulin_snomed = last_matching_med_dmd_before(insulin_snomed_clinical, feasibilityend_date).date
dataset.tmp_cov_date_antidiabetic_drugs_snomed = last_matching_med_dmd_before(antidiabetic_drugs_snomed_clinical, feasibilityend_date).date
dataset.tmp_cov_date_nonmetform_drugs_snomed = last_matching_med_dmd_before(non_metformin_dmd, feasibilityend_date).date # this extra step makes sense for the diabetes algorithm (otherwise not)

# Identify latest date (in period before feasibilityend_date) that any diabetes medication was prescribed
dataset.tmp_cov_date_diabetes_medication = maximum_of(dataset.tmp_cov_date_insulin_snomed, dataset.tmp_cov_date_antidiabetic_drugs_snomed) # why excluding tmp_cov_date_nonmetform_drugs_snomed? -> this extra step makes sense for the diabetes algorithm (otherwise not)

# Identify latest date (in period before feasibilityend_date) that any diabetes diagnosis codes were recorded
dataset.tmp_cov_date_latest_diabetes_diag = maximum_of(
  dataset.cov_date_t2dm, 
  dataset.cov_date_t1dm,
  dataset.cov_date_otherdm,
  dataset.cov_date_gestationaldm,
  dataset.tmp_cov_date_poccdm,
  dataset.tmp_cov_date_diabetes_medication,
  dataset.tmp_cov_date_nonmetform_drugs_snomed
)
## DIABETES algo variables end ------------------------