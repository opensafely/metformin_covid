#######################################################################################
# SEQTRIALS PIPELINE - dataset 2: people with a T2DM diagnosis in the T2DM window +
# data quality assurance (QA) and completeness criteria variables ONLY (no other
# eligibility criteria, covariates, treatment or outcome variables; these follow in dataset 3)
#
# Population: T2DM diagnosis between 6 months before study start (the first trial, T0,
# enrols people diagnosed within the previous 6 months) and study end.
#
# Anchor date: T2DM diagnosis date (elig_date_t2dm), used only where a variable depends
# on time (deprivation and STP). Since all diagnoses are recent, it is at most 6 months
# before the baseline of any trial a person can enter. The QA variables are not anchored:
# they check the plausibility of the entire record. Sex is fixed in the patients table.
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
    case,
    create_dataset,
    when
)

## TPP tables
from ehrql.tables.tpp import (
    addresses,
    apcs,
    clinical_events,
    medications,
    ons_deaths,
    patients,
    practice_registrations
)

## for import of diabetes algo created data
from ehrql.query_language import (
    table_from_file,
    PatientFrame,
    Series
)
from datetime import date

## All codelists from codelists.py
from codelists import *

## variable helper functions
from variable_helper_functions import *

## json (for the dates)
import json

# random seed (ideally use numpy, but currently not working on my local environment)
import random
random.seed(19283) # random seed

## import diabetes algo created data: only the T2DM diagnosis date
## (diabetes-algo sets t2dm_date only for people classified as T2DM, otherwise missing)
@table_from_file("output/seqtrials/data_processed_dm_algo.csv.gz")
class data_processed(PatientFrame):
    t2dm_date = Series(date)

#######################################################################################
# DEFINE the dates
#######################################################################################
with open("output/seqtrials/study_dates.json") as f:
  study_dates = json.load(f)
studystart_date = study_dates["studystart_date"]
studyend_date = study_dates["studyend_date"]

## Start of the T2DM window: 6 calendar months before study start (2020-02-01 -> 2019-08-01)
_start = date.fromisoformat(studystart_date)
_month0 = _start.year * 12 + (_start.month - 1) - 6
t2dm_window_start = date(_month0 // 12, _month0 % 12 + 1, _start.day).isoformat()

#######################################################################################
# INITIALISE the dataset and set the dummy dataset size
#######################################################################################
dataset = create_dataset()
dataset.configure_dummy_data(population_size=7000)

## T2DM diagnosis date (= anchor date)
dataset.elig_date_t2dm = data_processed.t2dm_date
## Population: T2DM diagnosis within the T2DM window (both dates inclusive)
dataset.define_population(
    dataset.elig_date_t2dm.is_on_or_between(t2dm_window_start, studyend_date)
)

#######################################################################################
# DATA QUALITY ASSURANCE variables (entire record, no anchor date)
#######################################################################################
## Rules 1-3: Year of birth (missing; after year of death; before 1838 or after study end)
# Kept as date of birth: fn_extract_data.R converts *_birth_year columns to the year
dataset.qa_num_birth_year = patients.date_of_birth

## Rules 2 and 4: Date of death (ONS death registry)
dataset.qa_date_of_death = ons_deaths.date

## Sex: needed for rules 5-7 (and completeness criterion 1: known sex)
dataset.cov_cat_sex = patients.sex

## Rule 5: Pregnancy codes (primary care)
dataset.qa_bin_pregnancy = (
    clinical_events.where(clinical_events.snomedct_code.is_in(pregnancy_snomed_clinical))
    .exists_for_patient()
)

## Rule 6: Hormone replacement therapy or combined oral contraceptive pill (prescriptions)
dataset.qa_bin_hrt = medications.where(medications.dmd_code.is_in(hrt_dmd)).exists_for_patient()
dataset.qa_bin_cocp = medications.where(medications.dmd_code.is_in(cocp_dmd)).exists_for_patient()

## Rule 7: Prostate cancer codes, in primary care, hospital admissions (any diagnosis
## position) or anywhere on the death certificate.
# ICD-10 3-character codes are additionally matched with "X" padding (e.g. C61X), since
# cause of death is matched exactly (as in the diabetes-algo PR #31 codelists)
prostate_cancer_icd10_padded = prostate_cancer_icd10 + [code + "X" for code in prostate_cancer_icd10 if len(code) == 3]
dataset.qa_bin_prostate_cancer = case(
    when(clinical_events.where(clinical_events.snomedct_code.is_in(prostate_cancer_snomed_clinical)).exists_for_patient()).then(True),
    when(apcs.where(apcs.all_diagnoses.contains_any_of(prostate_cancer_icd10_padded)).exists_for_patient()).then(True),
    when(cause_of_death_matches(prostate_cancer_icd10_padded)).then(True),
    otherwise=False
)

#######################################################################################
# DATA COMPLETENESS variables (at anchor date = T2DM diagnosis)
#######################################################################################
## Criterion 1: Known sex (female/male) -> derived from cov_cat_sex above

## Criterion 2: Known deprivation (IMD) at T2DM diagnosis
# imd_rounded can be missing even if an address exists (e.g. no postcode), hence not
# only checking for the existence of an address
dataset.qa_bin_known_imd = addresses.for_patient_on(dataset.elig_date_t2dm).imd_rounded.is_not_null()

## Criterion 3: Known sustainability and transformation partnership (STP) at T2DM
## diagnosis, from the practice registration on that date
# (registration only needs to be active on that date; the eligibility criterion of
# >= 12 months registration is applied at each trial baseline in dataset 3)
dataset.qa_bin_known_stp = practice_registrations.for_patient_on(dataset.elig_date_t2dm).practice_stp.is_not_null()
