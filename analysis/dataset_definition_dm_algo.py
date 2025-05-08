#######################################################################################
# IMPORT
#######################################################################################
## ehrQL functions
from ehrql import (
    create_dataset,
    minimum_of,
    case,
    when
)

## TPP tables
from ehrql.tables.tpp import (
    clinical_events,
    patients,
    ethnicity_from_sus
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
# DEFINE the exposure end date (pandemic start)
#######################################################################################
with open("output/study_dates.json") as f:
  study_dates = json.load(f)
pandemicstart_date = study_dates["pandemicstart_date"]

#######################################################################################
# INITIALISE the dataset and set the dummy dataset size
#######################################################################################
dataset = create_dataset()
dataset.configure_dummy_data(population_size=5000)
dataset.define_population(patients.exists_for_patient())

#######################################################################################
# Table 3) ELIGIBILITY criteria
#######################################################################################
## Year of birth
dataset.qa_num_birth_year = patients.date_of_birth

## Ethnicity in 6 categories (mainly for diabetes algo())
ethnicity_snomed = (
    clinical_events.where(clinical_events.snomedct_code.is_in(ethnicity_codes))
    .sort_by(clinical_events.date)
    .last_for_patient()
    .snomedct_code.to_category(ethnicity_codes, default="Missing")
)

ethnicity_sus = ethnicity_from_sus.code

dataset.cov_cat_ethnicity = case(
  when((ethnicity_snomed == "1") | ((ethnicity_snomed == "Missing") & (ethnicity_sus.is_in(["A", "B", "C"])))).then("White"),
  when((ethnicity_snomed == "2") | ((ethnicity_snomed == "Missing") & (ethnicity_sus.is_in(["D", "E", "F", "G"])))).then("Mixed"),
  when((ethnicity_snomed == "3") | ((ethnicity_snomed == "Missing") & (ethnicity_sus.is_in(["H", "J", "K", "L"])))).then("Asian"),
  when((ethnicity_snomed == "4") | ((ethnicity_snomed == "Missing") & (ethnicity_sus.is_in(["M", "N", "P"])))).then("Black"),
  when((ethnicity_snomed == "5") | ((ethnicity_snomed == "Missing") & (ethnicity_sus.is_in(["R", "S"])))).then("Other"),
  otherwise="Unknown", 
)

## DIABETES algo variables start ------------------------

## Type 1 Diabetes 
# First date from primary+secondary, but also primary care date separately for diabetes algo
dataset.tmp_elig_date_t1dm_ctv3 = first_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, pandemicstart_date).date
dataset.elig_date_t1dm = minimum_of(
    (first_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, pandemicstart_date).date),
    (first_matching_event_apc_before(diabetes_type1_icd10, pandemicstart_date).admission_date)
)
# Count codes (individually and together, for diabetes algo)
tmp_elig_count_t1dm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_type1_ctv3_clinical, pandemicstart_date)
tmp_elig_count_t1dm_hes = count_matching_event_apc_before(diabetes_type1_icd10, pandemicstart_date)
dataset.tmp_elig_count_t1dm = tmp_elig_count_t1dm_ctv3 + tmp_elig_count_t1dm_hes

## Type 2 Diabetes
# First date from primary+secondary, but also primary care date separately for diabetes algo)
dataset.tmp_elig_date_t2dm_ctv3 = first_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, pandemicstart_date).date
dataset.elig_date_t2dm = minimum_of(
    (first_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, pandemicstart_date).date),
    (first_matching_event_apc_before(diabetes_type2_icd10, pandemicstart_date).admission_date)
)

# Count codes (individually and together, for diabetes algo)
tmp_elig_count_t2dm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_type2_ctv3_clinical, pandemicstart_date)
tmp_elig_count_t2dm_hes = count_matching_event_apc_before(diabetes_type2_icd10, pandemicstart_date)
dataset.tmp_elig_count_t2dm = tmp_elig_count_t2dm_ctv3 + tmp_elig_count_t2dm_hes

## Diabetes unspecified/other
# First date
dataset.elig_date_otherdm = first_matching_event_clinical_ctv3_before(diabetes_other_ctv3_clinical, pandemicstart_date).date
# Count codes
dataset.tmp_elig_count_otherdm = count_matching_event_clinical_ctv3_before(diabetes_other_ctv3_clinical, pandemicstart_date)

## Gestational diabetes ## Comment 10/12/2024: Search in both primary and secondary
# First date from primary+secondary
dataset.elig_date_gestationaldm = minimum_of(
    (first_matching_event_clinical_ctv3_before(diabetes_gestational_ctv3_clinical, pandemicstart_date).date),
    (first_matching_event_apc_before(diabetes_gestational_icd10, pandemicstart_date).admission_date)
)

## Diabetes diagnostic codes
# First date
dataset.tmp_elig_date_poccdm = first_matching_event_clinical_ctv3_before(diabetes_diagnostic_ctv3_clinical, pandemicstart_date).date
# Count codes
dataset.tmp_elig_count_poccdm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_diagnostic_ctv3_clinical, pandemicstart_date)

### Other variables needed to define diabetes
## HbA1c
# Maximum HbA1c measure (in the same period)
dataset.tmp_elig_num_max_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(pandemicstart_date))
    .numeric_value.maximum_for_patient()
)
# Date of first maximum HbA1c measure
dataset.tmp_elig_date_max_hba1c = ( 
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(pandemicstart_date)) # this line of code probably not needed again
    .where(clinical_events.numeric_value == dataset.tmp_elig_num_max_hba1c_mmol_mol)
    .sort_by(clinical_events.date)
    .first_for_patient() 
    .date
)

## Diabetes drugs
# First dates
dataset.tmp_elig_date_insulin_snomed = first_matching_med_dmd_before(insulin_dmd, pandemicstart_date).date
dataset.tmp_elig_date_antidiabetic_drugs_snomed = first_matching_med_dmd_before(antidiabetic_drugs_snomed_clinical, pandemicstart_date).date
dataset.tmp_elig_date_nonmetform_drugs_snomed = first_matching_med_dmd_before(non_metformin_dmd, pandemicstart_date).date # this extra step makes sense for the diabetes algorithm (otherwise not)

# Identify first date (in same period) that any diabetes medication was prescribed
dataset.tmp_elig_date_diabetes_medication = minimum_of(dataset.tmp_elig_date_insulin_snomed, dataset.tmp_elig_date_antidiabetic_drugs_snomed) # why excluding tmp_elig_date_nonmetform_drugs_snomed? -> this extra step makes sense for the diabetes algorithm (otherwise not)

# Identify first date (in same period) that any diabetes diagnosis codes were recorded
dataset.tmp_elig_date_first_diabetes_diag = minimum_of(
  dataset.elig_date_t2dm, 
  dataset.elig_date_t1dm,
  dataset.elig_date_otherdm,
  dataset.elig_date_gestationaldm,
  dataset.tmp_elig_date_poccdm,
  dataset.tmp_elig_date_diabetes_medication,
  dataset.tmp_elig_date_nonmetform_drugs_snomed
)

## DIABETES algo variables end ------------------------