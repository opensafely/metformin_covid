#######################################################################################
# SEQTRIALS PIPELINE - main dataset: eligibility, covariates, treatment and outcomes
# for the people with T2DM in the study period who passed the QA and completeness criteria (data_t2dm_qa.arrow)
#
# Design: trial k starts at the T2DM diagnosis date + 30*k days (k = 0, ..., 6)
# Eligibility is re-assessed at each trial start. 
# No event-level tables possible in this repo, therefore four kinds of variables:
#   A. Fixed variables (no date needed)
#   B. First-ever dates up to study end (compared with each trial start in R:
#      e.g. "before the trial start" = exclusion, "within 30 days after" = treatment)
#   C. Status at each trial start, suffix _t0 ... _t6
#   D. Lab values on the 30-day grid, suffix _g0 ... _g31: the most recent value at T2DM diagnosis + 30*i days, from diagnosis to study end
#######################################################################################

#######################################################################################
# IMPORT
#######################################################################################
## Make the shared modules in analysis/ (codelists.py, variable_helper_functions.py) importable: ehrQL only adds this file's own folder (analysis/seqtrials/) to the path
## (without this, `from codelists import *` silently imports the top-level codelists/ folder)
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

## ehrQL functions
from ehrql import (
    case,
    create_dataset,
    days,
    minimum_of,
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

## for import of the eligible people from dataset 2
from ehrql.query_language import (
    table_from_file,
    PatientFrame,
    Series
)
from datetime import date, timedelta

## All codelists from codelists.py
from codelists import *

## variable helper functions
from variable_helper_functions import *

## json (for the dates)
import json

# random seed
import random
random.seed(19283) # random seed

## import the people who passed the QA and completeness criteria (output of data_process_qa.R)
@table_from_file("output/seqtrials/data/data_t2dm_qa.arrow")
class data_t2dm_qa(PatientFrame):
    elig_date_t2dm = Series(date)

#######################################################################################
# DEFINE the dates and the trial grid
#######################################################################################
with open("output/seqtrials/study_dates.json") as f:
  study_dates = json.load(f)
studystart_date = study_dates["studystart_date"]
studyend_date = study_dates["studyend_date"]

## Trials: start at T2DM diagnosis + TRIAL_STEP*k days, k = 0, ..., N_TRIALS - 1
# Corresponds to the T2DM window in dataset_definition_t2dm_qa.py and data_process_qa.R:
# window start = study start - TRIAL_STEP * (N_TRIALS - 1) = study start - 180 days)
TRIAL_STEP = 30
N_TRIALS = 7

## Grid for the lab values: T2DM diagnosis + TRIAL_STEP*i days, i = 0, ..., N_GRID - 1,
# up to study end for the earliest possible diagnosis (start of the T2DM window)
t2dm_window_start = date.fromisoformat(studystart_date) - timedelta(days=TRIAL_STEP * (N_TRIALS - 1))
N_GRID = (date.fromisoformat(studyend_date) - t2dm_window_start).days // TRIAL_STEP + 1

## Lookback window for the most recent lab value: 2 years (1 year = 366 days, as in landmark)
LOOKBACK_LAB = 2 * 366

#######################################################################################
# INITIALISE the dataset and set the dummy dataset size
#######################################################################################
dataset = create_dataset()
dataset.configure_dummy_data(population_size=7000)
dataset.define_population(data_t2dm_qa.exists_for_patient())

## T2DM diagnosis date (trial k starts at elig_date_t2dm + 30*k days)
dataset.elig_date_t2dm = data_t2dm_qa.elig_date_t2dm

#######################################################################################
# B. FIRST-EVER DATES (up to study end) used for eligibility
#######################################################################################
## Date of death (ONS): alive at trial start; later also the competing event
dataset.out_date_death = ons_deaths.date

## Metformin and other antidiabetic drugs: first-ever prescription (no lower date bound)
# Eligibility: no prescription BEFORE the trial start (strictly before, i.e. starting on the trial start day counts as treatment, as in landmark). The same dates define the treatment (start within 30 days after the trial start) and deviations from it, hence exp_*
dataset.exp_date_metfin_first = first_matching_med_dmd_before(metformin_dmd, studyend_date).date
dataset.exp_date_sulfo_first = first_matching_med_dmd_before(sulfonylurea_dmd, studyend_date).date
dataset.exp_date_dpp4_first = first_matching_med_dmd_before(dpp4_dmd, studyend_date).date
dataset.exp_date_tzd_first = first_matching_med_dmd_before(tzd_dmd, studyend_date).date
dataset.exp_date_sglt2_first = first_matching_med_dmd_before(sglt2_dmd, studyend_date).date
dataset.exp_date_glp1_first = first_matching_med_dmd_before(glp1_dmd, studyend_date).date
dataset.exp_date_megli_first = first_matching_med_dmd_before(meglitinides_dmd, studyend_date).date
dataset.exp_date_agi_first = first_matching_med_dmd_before(agi_dmd, studyend_date).date
dataset.exp_date_insulin_first = first_matching_med_dmd_before(insulin_dmd, studyend_date).date

## Contraindications to metformin: first-ever record (eligibility: on or before the trial start)
# Known hypersensitivity / intolerance to metformin
dataset.cov_date_metfin_allergy_first = first_matching_event_clinical_snomed_before(metformin_allergy_snomed_clinical, studyend_date).date
# Moderate to severe renal impairment (CKD stage 4/5), primary care or hospital (any diagnosis position)
dataset.cov_date_ckd45_first = minimum_of(
    first_matching_event_clinical_snomed_before(ckd_snomed_clinical_45, studyend_date).date,
    first_matching_event_apc_before(ckd_stage4_icd10 + ckd_stage5_icd10, studyend_date).admission_date
)
# Advanced decompensated liver cirrhosis, primary care or hospital (any diagnosis position)
dataset.cov_date_liver_cirrhosis_first = minimum_of(
    first_matching_event_clinical_snomed_before(advanced_decompensated_cirrhosis_snomed_clinical + ascitic_drainage_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(advanced_decompensated_cirrhosis_icd10, studyend_date).admission_date
)

#######################################################################################
# C. STATUS AT EACH TRIAL START (suffix _t0 ... _t6)
#######################################################################################
for k in range(N_TRIALS):
    trial_start = dataset.elig_date_t2dm + days(TRIAL_STEP * k)

    ## Age at trial start (eligibility: 18-85 years)
    dataset.add_column(f"cov_num_age_t{k}", patients.age_on(trial_start))

    ## Registered with a TPP practice for at least 1 year (366 days) up to the trial start
    dataset.add_column(
        f"cov_bin_registered_t{k}",
        practice_registrations.spanning_with_systmone(trial_start - days(366), trial_start).exists_for_patient()
    )

    ## Hospitalised on the trial start: admitted on or before and discharged after it (a missing discharge date is not counted, as in landmark)
    dataset.add_column(
        f"cov_bin_hosp_t{k}",
        apcs.where(
            apcs.admission_date.is_on_or_before(trial_start)
            & apcs.discharge_date.is_after(trial_start)
        ).exists_for_patient()
    )

    ## Care home resident on the trial start: PRIMIS code (any time before) or TPP address flags on that date (as in landmark, from opioids-covid-research)
    care_home_address = addresses.for_patient_on(trial_start)
    dataset.add_column(
        f"cov_bin_carehome_t{k}",
        case(
            when(last_matching_event_clinical_snomed_before(carehome, trial_start).exists_for_patient()).then(True),
            when(care_home_address.care_home_is_potential_match).then(True),
            when(care_home_address.care_home_requires_nursing).then(True),
            when(care_home_address.care_home_does_not_require_nursing).then(True),
            otherwise=False
        )
    )

    ## Palliative care code within 6 months (183 days) before or on the trial start
    dataset.add_column(
        f"cov_bin_palliative_t{k}",
        last_matching_event_clinical_snomed_between(palliative_snomed, trial_start - days(183), trial_start).exists_for_patient()
    )

    ## Prescription of a drug interacting with metformin within 14 days before or on the trial start
    dataset.add_column(
        f"cov_bin_metfin_interaction_t{k}",
        last_matching_med_dmd_between(metformin_interaction_dmd, trial_start - days(14), trial_start).exists_for_patient()
    )

#######################################################################################
# D. LAB VALUES ON THE 30-DAY GRID (suffix _g0 ... _g31)
#######################################################################################
## HbA1c (IFCC, mmol/mol): most recent value within 2 years before or on the grid point
# Eligibility at trial k: value at grid point k exists and is <= 75 mmol/mol (as in landmark).
# Missing and non-positive values are skipped (then the previous valid value is used, LVCF).
# The remaining plausibility checks are done in R
hba1c_valid = (
    clinical_events.where(clinical_events.snomedct_code.is_in(hba1c_snomed))
    .where(clinical_events.numeric_value > 0)
)
for i in range(N_GRID):
    grid_date = dataset.elig_date_t2dm + days(TRIAL_STEP * i)
    hba1c_last = (
        hba1c_valid.where(clinical_events.date.is_on_or_between(grid_date - days(LOOKBACK_LAB), grid_date))
        .sort_by(clinical_events.date)
        .last_for_patient()
    )
    dataset.add_column(f"cov_num_hba1c_g{i}", hba1c_last.numeric_value)
    dataset.add_column(f"cov_date_hba1c_g{i}", hba1c_last.date)
