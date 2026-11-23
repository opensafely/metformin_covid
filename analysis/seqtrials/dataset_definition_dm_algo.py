#######################################################################################
# SEQTRIALS PIPELINE - dataset for the diabetes algorithm
# Differences to the landmark pipeline (analysis/dataset_definition_dm_algo.py):
# (a) all diabetes-related events are searched looking back from the study end date
#     (not the pandemic start), since the sequential trials enrol people with a T2DM
#     diagnosis throughout the pandemic
# (b) uses the updated codelists from diabetes-algo PR #31 (dmalgo_* in codelists.py):
#     SNOMED instead of CTV3 for T1DM/T2DM/other/gestational DM, plus NOS DM (SNOMED + ICD-10);
#     the non-diagnostic diabetes codes remain CTV3
#######################################################################################

#######################################################################################
# IMPORT
#######################################################################################
## Make the shared modules in analysis/ (codelists.py, variable_helper_functions.py)
## importable: ehrQL only adds this file's own folder (analysis/seqtrials/) to the path
## (without this, `from codelists import *` silently imports the top-level codelists/ folder)
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

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
import random
random.seed(19283) # random seed

#######################################################################################
# DEFINE the date up to which diabetes events are searched (study end date)
#######################################################################################
with open("output/seqtrials/study_dates.json") as f:
  study_dates = json.load(f)
studyend_date = study_dates["studyend_date"]

#######################################################################################
# INITIALISE the dataset and set the dummy dataset size
#######################################################################################
dataset = create_dataset()
dataset.configure_dummy_data(population_size=7000)
dataset.define_population(patients.exists_for_patient())

#######################################################################################
# Variables needed for the diabetes algorithm
#######################################################################################
## Year of birth
dataset.qa_num_birth_year = patients.date_of_birth

## Ethnicity in 6 categories
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

## Type 1 Diabetes
# First date from primary+secondary, but also primary care date separately for diabetes algo
dataset.tmp_elig_date_t1dm_primarycare = first_matching_event_clinical_snomed_before(dmalgo_type1_snomed, studyend_date).date
dataset.tmp_elig_date_t1dm = minimum_of(
    (first_matching_event_clinical_snomed_before(dmalgo_type1_snomed, studyend_date).date),
    (first_matching_event_apc_before(dmalgo_type1_icd10, studyend_date).admission_date)
)
# Count codes (individually and together, for diabetes algo)
tmp_elig_count_t1dm_primarycare = count_matching_event_clinical_snomed_before(dmalgo_type1_snomed, studyend_date)
tmp_elig_count_t1dm_hes = count_matching_event_apc_before(dmalgo_type1_icd10, studyend_date)
dataset.tmp_elig_count_t1dm = tmp_elig_count_t1dm_primarycare + tmp_elig_count_t1dm_hes

## Type 2 Diabetes
# First date from primary+secondary, but also primary care date separately for diabetes algo
dataset.tmp_elig_date_t2dm_primarycare = first_matching_event_clinical_snomed_before(dmalgo_type2_snomed, studyend_date).date
dataset.elig_date_t2dm = minimum_of(
    (first_matching_event_clinical_snomed_before(dmalgo_type2_snomed, studyend_date).date),
    (first_matching_event_apc_before(dmalgo_type2_icd10, studyend_date).admission_date)
)
# Count codes (individually and together, for diabetes algo)
tmp_elig_count_t2dm_primarycare = count_matching_event_clinical_snomed_before(dmalgo_type2_snomed, studyend_date)
tmp_elig_count_t2dm_hes = count_matching_event_apc_before(dmalgo_type2_icd10, studyend_date)
dataset.tmp_elig_count_t2dm = tmp_elig_count_t2dm_primarycare + tmp_elig_count_t2dm_hes

## Diabetes of other specified aetiology or not otherwise specified (other DM)
# First date from primary+secondary
dataset.tmp_elig_date_otherdm = minimum_of(
    (first_matching_event_clinical_snomed_before(dmalgo_other_snomed + dmalgo_nos_snomed, studyend_date).date),
    (first_matching_event_apc_before(dmalgo_nos_icd10, studyend_date).admission_date)
)
# Count codes (primary+secondary)
tmp_elig_count_otherdm_primarycare = count_matching_event_clinical_snomed_before(dmalgo_other_snomed + dmalgo_nos_snomed, studyend_date)
tmp_elig_count_otherdm_hes = count_matching_event_apc_before(dmalgo_nos_icd10, studyend_date)
dataset.tmp_elig_count_otherdm = tmp_elig_count_otherdm_primarycare + tmp_elig_count_otherdm_hes

## Gestational diabetes
# First date from primary+secondary
dataset.tmp_elig_date_gestationaldm = minimum_of(
    (first_matching_event_clinical_snomed_before(dmalgo_gestational_snomed, studyend_date).date),
    (first_matching_event_apc_before(dmalgo_gestational_icd10, studyend_date).admission_date)
)

## Diabetes non-diagnostic codes (still CTV3)
# First date
dataset.tmp_elig_date_poccdm = first_matching_event_clinical_ctv3_before(diabetes_diagnostic_ctv3_clinical, studyend_date).date
# Count codes
dataset.tmp_elig_count_poccdm_ctv3 = count_matching_event_clinical_ctv3_before(diabetes_diagnostic_ctv3_clinical, studyend_date)

### Other variables needed to define diabetes
## HbA1c
# Maximum HbA1c measure (in the same period)
dataset.tmp_elig_num_max_hba1c_mmol_mol = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(studyend_date))
    .numeric_value.maximum_for_patient()
)
# Date of first maximum HbA1c measure
dataset.tmp_elig_date_max_hba1c = (
  clinical_events.where(
    clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.date.is_on_or_before(studyend_date))
    .where(clinical_events.numeric_value == dataset.tmp_elig_num_max_hba1c_mmol_mol)
    .sort_by(clinical_events.date)
    .first_for_patient()
    .date
)

## Diabetes drugs
# First dates
dataset.tmp_elig_date_insulin = first_matching_med_dmd_before(insulin_dmd, studyend_date).date
dataset.tmp_elig_date_antidiabetic_drugs = first_matching_med_dmd_before(antidiabetic_drugs_dmd, studyend_date).date
dataset.tmp_elig_date_nonmetform_drugs = first_matching_med_dmd_before(dmalgo_non_metformin_dmd, studyend_date).date

# Identify first date (in same period) that any diabetes medication was prescribed
dataset.tmp_elig_date_diabetes_medication = minimum_of(dataset.tmp_elig_date_insulin, dataset.tmp_elig_date_antidiabetic_drugs)

# Identify first date (in same period) that any diabetes diagnosis codes were recorded
dataset.tmp_elig_date_first_diabetes_diag = minimum_of(
  dataset.elig_date_t2dm,
  dataset.tmp_elig_date_t1dm,
  dataset.tmp_elig_date_otherdm,
  dataset.tmp_elig_date_gestationaldm,
  dataset.tmp_elig_date_poccdm,
  dataset.tmp_elig_date_diabetes_medication,
  dataset.tmp_elig_date_nonmetform_drugs
)
