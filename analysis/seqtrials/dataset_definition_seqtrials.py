#######################################################################################
# SEQTRIALS PIPELINE - main dataset: eligibility, covariates, treatment and outcomes
# for the people with T2DM in the study period who passed the QA and completeness criteria (data_t2dm_qa.arrow)
#
# Design: trial k starts at the T2DM diagnosis date + 30*k days (k = 0, ..., 6)
# Eligibility is re-assessed at each trial start. 
# No event-level tables possible in this repo, therefore the following kinds of variables:
#   A. Fixed variables (no date needed)
#   B. First-ever dates up to study end (compared with each trial start in R:
#      e.g. "before the trial start" = exclusion, "within 30 days after" = treatment)
#   C. Status at each trial start, suffix _t0 ... _t6
#   D. Lab values on the 30-day grid, suffix _g0 ... _g31: the most recent value at T2DM diagnosis + 30*i days, from diagnosis to study end
#   E. Outcomes: first date between study start and study end (any prior COVID-19 = ineligible)
#   F. Loss to follow-up per trial (suffix _t0 ... _t6)
#   G. Negative control outcome (fracture) per 30-day period (suffix _g0 ... _g31)
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
    appointments,
    clinical_events,
    medications,
    occupation_on_covid_vaccine_record,
    ons_deaths,
    patients,
    practice_registrations,
    sgss_covid_all_tests,
    vaccinations
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

## permissions for the restricted tables (appointments, occupation_on_covid_vaccine_record), as in landmark
import project_permissions

## json (for the dates)
import json

# random seed
import random
random.seed(19283) # random seed

## import the people who passed the QA and completeness criteria (output of data_process_qa.R)
@table_from_file("output/seqtrials/data/data_t2dm_qa.arrow")
class data_t2dm_qa(PatientFrame):
    elig_date_t2dm = Series(date)

## import ethnicity from the diabetes algorithm output (the same variable the algorithm used for the T2DM classification)
@table_from_file("output/seqtrials/data_processed_dm_algo.csv.gz")
class data_dm_algo(PatientFrame):
    ethnicity_cat = Series(str)

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
# A. FIXED VARIABLES
#######################################################################################
## Sex (only female/male left after the completeness criteria)
dataset.cov_cat_sex = patients.sex

## Ethnicity (5 groups + Unknown, from the diabetes algorithm output)
dataset.cov_cat_ethnicity = data_dm_algo.ethnicity_cat

## Healthcare worker, as recorded at any COVID-19 vaccination (recorded from Dec 2020 onwards, i.e. sometimes after the trial start; used as a proxy for the occupation)
dataset.cov_bin_healthcare_worker = (
    occupation_on_covid_vaccine_record.where(occupation_on_covid_vaccine_record.is_healthcare_worker == True)
    .exists_for_patient()
)

#######################################################################################
# B. FIRST-EVER DATES (up to study end)
#######################################################################################
#### B1. Used for eligibility
## Date of death (ONS): alive at trial start; later also the competing event
dataset.out_date_death = ons_deaths.date

## Metformin and other antidiabetic drugs: first-ever prescription (no lower date bound)
# Eligibility: no prescription BEFORE the trial start (strictly before, i.e. starting on the trial start day counts as treatment, as in landmark). The same dates define the treatment (start within 30 days after the trial start) and deviations from it, hence exp_*
dataset.exp_date_metfin_first = first_matching_med_dmd_before(metformin_dmd, studyend_date).date
## Metformin monotherapy (without other antidiabetic drugs in the same product): only this counts as treatment start (as in landmark).
# In R, per trial: metformin mono started within 30 days after the trial start = treated; no antidiabetic drug = control; metformin combo or another drug = excluded from that trial.
# The class codelists below include the combination products with metformin, i.e. a combo prescription sets both exp_date_metfin_first and the date of the other class
dataset.exp_date_metfin_mono_first = first_matching_med_dmd_before(metformin_mono_dmd, studyend_date).date
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

#### B2. Comorbidities (confounders): first-ever record in primary care or hospital (any diagnosis position); in R: present at trial k = date on or before the trial start
## Obesity (codes; the BMI value itself is on the grid, see D)
dataset.cov_date_obesity = minimum_of(
    first_matching_event_clinical_snomed_before(bmi_obesity_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(bmi_obesity_icd10, studyend_date).admission_date
)
## Acute myocardial infarction
dataset.cov_date_ami = minimum_of(
    first_matching_event_clinical_snomed_before(ami_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(ami_prior_icd10 + ami_icd10, studyend_date).admission_date
)
## All stroke
dataset.cov_date_all_stroke = minimum_of(
    first_matching_event_clinical_snomed_before(stroke_isch_snomed_clinical + stroke_sah_hs_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(stroke_isch_icd10 + stroke_sah_hs_icd10, studyend_date).admission_date
)
## Other arterial embolism
dataset.cov_date_other_arterial_embolism = minimum_of(
    first_matching_event_clinical_snomed_before(other_arterial_embolism_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(other_arterial_embolism_icd10, studyend_date).admission_date
)
## Venous thromboembolism
dataset.cov_date_vte = minimum_of(
    first_matching_event_clinical_snomed_before(portal_vein_thrombosis_snomed_clinical + dvt_dvt_snomed_clinical + dvt_icvt_snomed_clinical + dvt_pregnancy_snomed_clinical + other_dvt_snomed_clinical + pe_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(portal_vein_thrombosis_icd10 + dvt_dvt_icd10 + dvt_icvt_icd10 + dvt_pregnancy_icd10 + other_dvt_icd10 + icvt_pregnancy_icd10 + pe_icd10, studyend_date).admission_date
)
## Heart failure
dataset.cov_date_hf = minimum_of(
    first_matching_event_clinical_snomed_before(hf_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(hf_icd10, studyend_date).admission_date
)
## Angina
dataset.cov_date_angina = minimum_of(
    first_matching_event_clinical_snomed_before(angina_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(angina_icd10, studyend_date).admission_date
)
## Dementia
dataset.cov_date_dementia = minimum_of(
    first_matching_event_clinical_snomed_before(dementia_snomed_clinical + dementia_vascular_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(dementia_icd10 + dementia_vascular_icd10, studyend_date).admission_date
)
## Cancer
dataset.cov_date_cancer = minimum_of(
    first_matching_event_clinical_snomed_before(cancer_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(cancer_icd10, studyend_date).admission_date
)
## Hypertension
dataset.cov_date_hypertension = minimum_of(
    first_matching_event_clinical_snomed_before(hypertension_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(hypertension_icd10, studyend_date).admission_date
)
## Depression
dataset.cov_date_depression = minimum_of(
    first_matching_event_clinical_snomed_before(depression_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(depression_icd10, studyend_date).admission_date
)
## Chronic obstructive pulmonary disease
dataset.cov_date_copd = minimum_of(
    first_matching_event_clinical_snomed_before(copd_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(copd_icd10, studyend_date).admission_date
)
## Liver disease
dataset.cov_date_liver_disease = minimum_of(
    first_matching_event_clinical_snomed_before(liver_disease_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(liver_disease_icd10, studyend_date).admission_date
)
## Chronic kidney disease (any stage)
dataset.cov_date_chronic_kidney_disease = minimum_of(
    first_matching_event_clinical_snomed_before(ckd_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(ckd_icd10, studyend_date).admission_date
)
## Polycystic ovary syndrome
dataset.cov_date_pcos = minimum_of(
    first_matching_event_clinical_snomed_before(pcos_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(pcos_icd10, studyend_date).admission_date
)
## Diabetes complications (foot, retino, neuro, nephro)
dataset.cov_date_diabetescomp = minimum_of(
    first_matching_event_clinical_snomed_before(diabetescomp_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(diabetescomp_icd10, studyend_date).admission_date
)
## Prediabetes: first prediabetes code, or first HbA1c in the prediabetes range (42-47.9 mmol/mol, as in landmark), in primary care
# In R: compare with the T2DM DIAGNOSIS date (as in landmark), not the trial start: an HbA1c of 42-47.9 after the diagnosis means well-controlled diabetes, not prediabetes
dataset.cov_date_prediabetes = minimum_of(
    first_matching_event_clinical_snomed_before(prediabetes_snomed, studyend_date).date,
    first_matching_event_clinical_snomed_before(
        hba1c_snomed, studyend_date,
        where=(clinical_events.numeric_value >= 42) & (clinical_events.numeric_value <= 47.9)
    ).date
)

#### B3. COVID-19 vaccination (confounder; time-updated in R: number of doses = dose dates on or before the trial start / grid date)
# Based on https://github.com/opensafely/covid-vaccine-history: vaccinations table, target disease SARS-CoV-2, doses from the start of the vaccination programme (2020-12-07; earlier records = vaccine trials, excluded) up to study end.
# Doses less than 14 days after the previous (kept) dose are not counted (duplicates/data quality), applied directly in the chain below: dose n = first record at least 14 days after dose n-1.
# Up to study end (2022-03-01) at most ~4 doses (primary course, 3rd primary dose for the immunosuppressed, autumn 2021 booster); dose 6 should be (almost) empty (will be checked in R)
VAX_PROGRAMME_START = "2020-12-07"
VAX_MIN_GAP = 14
N_VAX = 6
covid_vaccinations = vaccinations.where(
    (vaccinations.target_disease == "SARS-2 CORONAVIRUS")
    & vaccinations.date.is_on_or_between(VAX_PROGRAMME_START, studyend_date)
)
previous_vax_date = None
for n in range(1, N_VAX + 1):
    if previous_vax_date is None:
        current_vax = covid_vaccinations.sort_by(vaccinations.date).first_for_patient()
    else:
        current_vax = (
            covid_vaccinations.where(vaccinations.date.is_on_or_after(previous_vax_date + days(VAX_MIN_GAP)))
            .sort_by(vaccinations.date)
            .first_for_patient()
        )
    dataset.add_column(f"cov_date_covid_vax_{n}", current_vax.date)
    previous_vax_date = current_vax.date

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

    ## Index of Multiple Deprivation (IMD) on the trial start, in fifths (rank rounded to the nearest 100; 32844 LSOAs in England)
    imd_rounded = addresses.for_patient_on(trial_start).imd_rounded
    dataset.add_column(
        f"cov_cat_deprivation_5_t{k}",
        case(
            when((imd_rounded >= 0) & (imd_rounded < int(32844 * 1 / 5))).then("1 (most deprived)"),
            when(imd_rounded < int(32844 * 2 / 5)).then("2"),
            when(imd_rounded < int(32844 * 3 / 5)).then("3"),
            when(imd_rounded < int(32844 * 4 / 5)).then("4"),
            when(imd_rounded <= 32844).then("5 (least deprived)"),
            otherwise="Unknown"
        )
    )

    ## Rural/urban classification of the address on the trial start
    dataset.add_column(f"cov_cat_rural_urban_t{k}", addresses.for_patient_on(trial_start).rural_urban_classification)

    ## STP and region of the practice registration on the trial start
    registration = practice_registrations.for_patient_on(trial_start)
    dataset.add_column(f"cov_cat_stp_t{k}", registration.practice_stp)
    dataset.add_column(f"cov_cat_region_t{k}", registration.practice_nuts1_region_name)

    ## Smoking status on the trial start (as in landmark): most recent smoking code (S, E, N, M); never/missing are reclassified as ex-smoker if any smoking code was ever recorded
    most_recent_smoking_cat = last_matching_event_clinical_ctv3_before(smoking_clear, trial_start).ctv3_code.to_category(smoking_clear)
    ever_smoked = last_matching_event_clinical_ctv3_before(ever_smoking, trial_start).exists_for_patient()
    dataset.add_column(
        f"cov_cat_smoking_status_t{k}",
        case(
            when(most_recent_smoking_cat == "S").then("S"),
            when(most_recent_smoking_cat == "E").then("E"),
            when((most_recent_smoking_cat == "N") & ever_smoked).then("E"),
            when(most_recent_smoking_cat == "N").then("N"),
            when((most_recent_smoking_cat == "M") & ever_smoked).then("E"),
            otherwise="M"
        )
    )

    ## Consultation rate: number of attended appointments in the year (366 days) before or on the trial start, capped at 365 (as in landmark); proxy for health-seeking behaviour
    consultations = appointments.where(
        appointments.status.is_in([
            "Arrived",
            "In Progress",
            "Finished",
            "Visit",
            "Waiting",
            "Patient Walked Out",
        ])
        & appointments.start_date.is_on_or_between(trial_start - days(366), trial_start)
    ).count_for_patient()
    dataset.add_column(
        f"cov_num_consrate_t{k}",
        case(when(consultations <= 365).then(consultations), otherwise=365)
    )

    ## Number of HbA1c tests (any HbA1c code, with or without a value) in the 2 years before or on the trial start
    dataset.add_column(
        f"cov_num_counthba1c_t{k}",
        count_matching_event_clinical_snomed_between(hba1c_measurement_snomed, trial_start - days(2 * 366), trial_start)
    )

    ## Number of lifestyle advice codes in the 2 years before or on the trial start
    dataset.add_column(
        f"cov_num_countlifestyle_t{k}",
        count_matching_event_clinical_snomed_between(lifestyle_advice_snomed, trial_start - days(2 * 366), trial_start)
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

## BMI, total cholesterol and HDL cholesterol: most recent value (and its date) within 2 years before or on the grid point
# Missing and non-positive values are skipped; the plausibility checks are done in R (e.g. BMI 12-70, as in landmark)
for i in range(N_GRID):
    grid_date = dataset.elig_date_t2dm + days(TRIAL_STEP * i)

    ## BMI: recorded BMI (CTV3 22K..), measured at age 16 or older (helper from landmark)
    bmi_last = most_recent_bmi(
        where=clinical_events.date.is_on_or_between(grid_date - days(LOOKBACK_LAB), grid_date) & (clinical_events.numeric_value > 0),
        minimum_age_at_measurement=16
    )
    dataset.add_column(f"cov_num_bmi_g{i}", bmi_last.numeric_value)
    dataset.add_column(f"cov_date_bmi_g{i}", bmi_last.date)

    ## Total cholesterol
    chol_last = last_matching_event_clinical_snomed_between(cholesterol_snomed, grid_date - days(LOOKBACK_LAB), grid_date, where=clinical_events.numeric_value > 0)
    dataset.add_column(f"cov_num_chol_g{i}", chol_last.numeric_value)
    dataset.add_column(f"cov_date_chol_g{i}", chol_last.date)

    ## HDL cholesterol
    hdl_chol_last = last_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, grid_date - days(LOOKBACK_LAB), grid_date, where=clinical_events.numeric_value > 0)
    dataset.add_column(f"cov_num_hdl_chol_g{i}", hdl_chol_last.numeric_value)
    dataset.add_column(f"cov_date_hdl_chol_g{i}", hdl_chol_last.date)

#######################################################################################
# E. OUTCOMES (first date between study start and study end, both inclusive)
#######################################################################################
# Any prior COVID-19 (out_date_covid, and Long COVID/post-viral fatigue codes) on or before the trial start = ineligible (applied in R).
# Hence, for the eligible people, the first date in the study period is the first event after the trial start, and first dates are sufficient for all COVID-19 outcomes.
# COVID-19 codes and tests do not exist before 2020, hence no events are lost by starting at study start.

### COVID-19 HOSPITALISATION
## Hospital admission with COVID-19 as the primary diagnosis
dataset.out_date_covid_hes = first_matching_event_apc_between(covid_codes_incl_clin_diag, studystart_date, studyend_date, only_prim_diagnoses=True).admission_date
## Emergency care attendance with a COVID-19 diagnosis
dataset.out_date_covid_emergency = first_matching_event_ec_snomed_between(covid_emergency, studystart_date, studyend_date).arrival_date
## Combined: first COVID-19 hospitalisation
dataset.out_date_covid_hosp = minimum_of(dataset.out_date_covid_hes, dataset.out_date_covid_emergency)

### COVID-19 DEATH: COVID-19 anywhere on the death certificate (underlying or contributing cause)
dataset.out_date_covid_death = case(
    when(matching_death_between(covid_codes_incl_clin_diag, studystart_date, studyend_date)).then(ons_deaths.date)
)

### SEVERE COVID-19 (primary outcome): first COVID-19 hospitalisation or death
dataset.out_date_severecovid = minimum_of(dataset.out_date_covid_hosp, dataset.out_date_covid_death)

### ANY COVID-19 (outcome, and eligibility: no prior COVID-19): first COVID-19 code in primary care, positive SARS-CoV-2 test, COVID-19 hospitalisation or death
## COVID-19 code in primary care (clinical code, positive test or sequelae; CTV3)
tmp_covid_primary_care_date = first_matching_event_clinical_ctv3_between(
    covid_primary_care_code + covid_primary_care_positive_test + covid_primary_care_sequelae, studystart_date, studyend_date
).date
## Positive SARS-CoV-2 test (SGSS, date of specimen)
tmp_covid_sgss_date = (
    sgss_covid_all_tests.where(sgss_covid_all_tests.specimen_taken_date.is_on_or_between(studystart_date, studyend_date))
    .where(sgss_covid_all_tests.is_positive)
    .sort_by(sgss_covid_all_tests.specimen_taken_date)
    .first_for_patient()
    .specimen_taken_date
)
## Combined
dataset.out_date_covid = minimum_of(tmp_covid_primary_care_date, tmp_covid_sgss_date, dataset.out_date_covid_hosp, dataset.out_date_covid_death)

### LONG COVID (outcome, and eligibility: no prior Long COVID/post-viral fatigue code, as it implies a prior infection)
## Long COVID code in primary care (diagnosis, referral or assessment)
dataset.out_date_longcovid = first_matching_event_clinical_snomed_between(
    long_covid_diagnostic_snomed_clinical + long_covid_referral_snomed_clinical + long_covid_assessment_snomed_clinical, studystart_date, studyend_date
).date
## Long COVID diagnosis in hospital (any diagnosis position)
dataset.out_date_longcovid_sec_care = first_matching_event_apc_between(long_covid_diagnostic_icd10, studystart_date, studyend_date).admission_date
## Post-viral fatigue code in primary care
dataset.out_date_virfat = first_matching_event_clinical_snomed_between(post_viral_fatigue_snomed_clinical, studystart_date, studyend_date).date
## Combined: first Long COVID or post-viral fatigue code
dataset.out_date_longcovid_virfat = minimum_of(dataset.out_date_longcovid, dataset.out_date_longcovid_sec_care, dataset.out_date_virfat)

### POSITIVE CONTROL: diabetes-related death (T2DM anywhere on the death certificate; ICD-10 with "X" padding, as cause of death is matched exactly)
# CAVEAT (kept for now): T2DM is often a contributing cause on deaths from any cause, incl. COVID-19 deaths (overlap with the primary outcome); its recording may depend on treatment; and the expected effect of metformin within 1-2 years is unknown.
# A second positive control, the HbA1c change after the trial start, is derived in R from the HbA1c grid (D)
dataset.out_date_dm_death = case(
    when(matching_death_between(dmalgo_type2_icd10, studystart_date, studyend_date)).then(ons_deaths.date)
)

### All-cause death: see out_date_death above (B1; competing event)

#######################################################################################
# F. LOSS TO FOLLOW-UP, per trial (suffix _t0 ... _t6)
#######################################################################################
for k in range(N_TRIALS):
    trial_start = dataset.elig_date_t2dm + days(TRIAL_STEP * k)

    ## Deregistration: end date of the practice registration active on the trial start (as in landmark; missing = still registered)
    # for_patient_on() picks the most recent registration if several overlap on that date
    registration_end = practice_registrations.for_patient_on(trial_start).end_date
    dataset.add_column(f"cens_date_dereg_t{k}", registration_end)

    ## Practice move: registered again (with any TPP practice) on the day after the deregistration date, i.e. the person moved between TPP practices and is not actually lost to follow-up
    # Used to report how many deregistrations are practice moves (if many, we may change this approach and also extract on the grid instead)
    dataset.add_column(
        f"cens_bin_dereg_move_t{k}",
        practice_registrations.for_patient_on(registration_end + days(1)).exists_for_patient()
    )

#######################################################################################
# G. NEGATIVE CONTROL OUTCOME ON THE 30-DAY GRID (suffix _g0 ... _g31)
#######################################################################################
# Fracture: first hospital admission for a fracture (any diagnosis position) within period i = [elig_date_t2dm + 30*i, elig_date_t2dm + 30*(i+1)).
# Per period because fractures recur and a prior fracture does not make a person ineligible, so each trial counts the first fracture after its start
for i in range(N_GRID):
    period_start = dataset.elig_date_t2dm + days(TRIAL_STEP * i)
    period_end = dataset.elig_date_t2dm + days(TRIAL_STEP * (i + 1) - 1)
    dataset.add_column(
        f"out_date_fracture_g{i}",
        first_matching_event_apc_between(fracture_icd10, period_start, period_end).admission_date
    )
