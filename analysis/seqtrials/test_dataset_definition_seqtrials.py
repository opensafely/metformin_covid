#######################################################################################
# ASSURANCE TESTS for dataset_definition_seqtrials.py
#
# Run locally:
#   opensafely exec ehrql:v1 assure analysis/seqtrials/test_dataset_definition_seqtrials.py
#
# table_from_file always reads its file (not the test data), so this script writes small
# input files with the test patients (output/seqtrials/assure_test_input/) and points the dataset definition to them via the environment variable SEQTRIALS_ASSURE_INPUT_DIR
#######################################################################################
import csv
import os
from datetime import date, timedelta

D = date(2020, 6, 1)  # T2DM diagnosis date of all test patients


def d(n):
    """Date n days after (n < 0: before) the T2DM diagnosis date"""
    return D + timedelta(days=n)


## Codes (one example per codelist)
HBA1C = "999791000000106"  # hba1c_snomed
CHOL = "1005671000000105"  # cholesterol_snomed
HDL = "1005681000000107"  # hdl_cholesterol_snomed
BMI_CTV3 = "22K.."  # most_recent_bmi()
INTERACTION_DRUG = "42351911000001100"  # metformin_interaction_dmd (cimetidine)
PALLIATIVE = "1025611000000105"  # palliative_snomed
FRACTURE = "M484"  # fracture_icd10
COVID_ICD10 = "U071"  # covid_codes_incl_clin_diag
COVID_VACCINE = "SARS-2 CORONAVIRUS"

EVENT_TABLES = [
    "addresses", "apcs", "appointments", "clinical_events", "emergency_care_attendances",
    "medications", "occupation_on_covid_vaccine_record", "practice_registrations",
    "sgss_covid_all_tests", "vaccinations",
]


def make_patient(**tables):
    """Test patient with default (empty) tables; keyword arguments replace the defaults"""
    patient = {name: [] for name in EVENT_TABLES}
    patient["patients"] = {"date_of_birth": date(1960, 1, 1), "sex": "female"}
    patient["ons_deaths"] = {}
    patient["data_t2dm_qa"] = {"elig_date_t2dm": D}
    patient["data_dm_algo"] = {"ethnicity_cat": "White"}
    patient.update(tables)
    return patient


def registration(start, end):
    return {
        "start_date": start, "end_date": end, "practice_pseudo_id": 1,
        "practice_stp": "E54000001", "practice_nuts1_region_name": "London",
        "practice_systmone_go_live_date": date(2000, 1, 1),
    }


test_data = {
    ###################################################################################
    # D. LAB GRID
    ###################################################################################
    ## 1: HbA1c, BMI, cholesterol, HDL: window boundaries, invalid values skipped, no leaking from later dates
    1: make_patient(
        clinical_events=[
            # HbA1c
            {"date": d(-732), "snomedct_code": HBA1C, "numeric_value": 55},    # exactly 2 years (732 days) before g0 -> included
            {"date": d(-10), "snomedct_code": HBA1C, "numeric_value": 0},      # invalid -> skipped
            {"date": d(30), "snomedct_code": HBA1C, "numeric_value": 60},      # on g1 -> included in g1, not in g0
            {"date": d(31), "snomedct_code": HBA1C, "numeric_value": 70},      # day after g1 -> first in g2
            {"date": d(59), "snomedct_code": HBA1C, "numeric_value": None},    # missing value -> skipped
            # BMI
            {"date": d(-5), "ctv3_code": BMI_CTV3, "numeric_value": 30.5},
            {"date": d(-1), "ctv3_code": BMI_CTV3, "numeric_value": 0},        # invalid -> skipped
            # Cholesterol and HDL
            {"date": d(-100), "snomedct_code": CHOL, "numeric_value": 5.5},
            {"date": d(40), "snomedct_code": CHOL, "numeric_value": 0},        # invalid -> skipped
            {"date": d(-100), "snomedct_code": HDL, "numeric_value": 1.5},
        ],
        expected_in_population=True,
        expected_columns={
            "cov_num_hba1c_g0": 55, "cov_date_hba1c_g0": d(-732),
            "cov_num_hba1c_g1": 60, "cov_date_hba1c_g1": d(30),
            "cov_num_hba1c_g2": 70, "cov_date_hba1c_g2": d(31),
            "cov_num_hba1c_g25": 70, "cov_date_hba1c_g25": d(31),  # g25 = D+750: D+31 still within 732 days
            "cov_num_hba1c_g26": None, "cov_date_hba1c_g26": None,  # g26 = D+780: D+31 too old, D+59 has no value
            "cov_num_bmi_g0": 30.5, "cov_date_bmi_g0": d(-5),
            "cov_num_chol_g0": 5.5, "cov_date_chol_g0": d(-100),
            "cov_num_chol_g2": 5.5, "cov_date_chol_g2": d(-100),
            "cov_num_hdl_chol_g0": 1.5, "cov_date_hdl_chol_g0": d(-100),
        },
    ),
    ## 2: HbA1c 733 days before g0 -> outside the window
    2: make_patient(
        clinical_events=[{"date": d(-733), "snomedct_code": HBA1C, "numeric_value": 50}],
        expected_in_population=True,
        expected_columns={"cov_num_hba1c_g0": None, "cov_date_hba1c_g0": None},
    ),
    ## 3: BMI measured before the age of 16 -> ignored
    3: make_patient(
        patients={"date_of_birth": date(2005, 6, 1), "sex": "male"},
        clinical_events=[{"date": d(-1), "ctv3_code": BMI_CTV3, "numeric_value": 25}],
        expected_in_population=True,
        expected_columns={"cov_num_bmi_g0": None, "cov_num_age_t0": 15},
    ),

    ###################################################################################
    # B3. COVID-19 VACCINATION
    ###################################################################################
    ## 4: before programme start ignored; same-day duplicate and dose < 14 days later dropped; other vaccines and doses after study end ignored
    4: make_patient(
        vaccinations=[
            {"date": date(2020, 11, 1), "target_disease": COVID_VACCINE},   # before 2020-12-07 -> ignored
            {"date": date(2021, 1, 10), "target_disease": COVID_VACCINE},   # dose 1
            {"date": date(2021, 1, 10), "target_disease": COVID_VACCINE},   # same-day duplicate -> dropped
            {"date": date(2021, 1, 20), "target_disease": COVID_VACCINE},   # 10 days after dose 1 -> dropped
            {"date": date(2021, 1, 24), "target_disease": COVID_VACCINE},   # 14 days after dose 1 -> dose 2
            {"date": date(2021, 4, 1), "target_disease": COVID_VACCINE},    # dose 3
            {"date": date(2021, 6, 1), "target_disease": "INFLUENZA"},      # other vaccine -> ignored
            {"date": date(2022, 3, 5), "target_disease": COVID_VACCINE},    # after study end -> ignored
        ],
        expected_in_population=True,
        expected_columns={
            "cov_date_covid_vax_1": date(2021, 1, 10),
            "cov_date_covid_vax_2": date(2021, 1, 24),
            "cov_date_covid_vax_3": date(2021, 4, 1),
            "cov_date_covid_vax_4": None,
        },
    ),

    ###################################################################################
    # G. FRACTURE PER 30-DAY PERIOD
    ###################################################################################
    ## 5: period i = [D+30i, D+30i+29]; first fracture per period
    5: make_patient(
        apcs=[
            {"admission_date": d(-1), "all_diagnoses": FRACTURE},   # before g0 -> not in any period
            {"admission_date": d(29), "all_diagnoses": FRACTURE},   # last day of period 0
            {"admission_date": d(30), "all_diagnoses": FRACTURE},   # first day of period 1
            {"admission_date": d(35), "all_diagnoses": FRACTURE},   # second fracture in period 1 -> not the first
        ],
        expected_in_population=True,
        expected_columns={
            "out_date_fracture_g0": d(29),
            "out_date_fracture_g1": d(30),
            "out_date_fracture_g2": None,
        },
    ),

    ###################################################################################
    # C/F. REGISTRATION, DEREGISTRATION AND PRACTICE MOVE
    ###################################################################################
    ## 6: ongoing registration -> no deregistration
    6: make_patient(
        practice_registrations=[registration(date(2015, 1, 1), None)],
        expected_in_population=True,
        expected_columns={
            "cov_bin_registered_t0": True,
            "cens_date_dereg_t0": None, "cens_bin_dereg_move_t0": False,
            "cov_cat_region_t0": "London",
        },
    ),
    ## 7: new registration the day after the deregistration -> practice move; at trial 4 the new registration is active
    7: make_patient(
        practice_registrations=[registration(date(2015, 1, 1), d(100)), registration(d(101), None)],
        expected_in_population=True,
        expected_columns={
            "cens_date_dereg_t0": d(100), "cens_bin_dereg_move_t0": True,
            "cens_date_dereg_t4": None, "cens_bin_dereg_move_t4": False,
            "cov_bin_registered_t4": False,  # new registration < 1 year at trial 4 (t4 = D+120)
        },
    ),
    ## 8: new registration 2 days after the deregistration -> not a practice move (gap)
    8: make_patient(
        practice_registrations=[registration(date(2015, 1, 1), d(100)), registration(d(102), None)],
        expected_in_population=True,
        expected_columns={"cens_date_dereg_t0": d(100), "cens_bin_dereg_move_t0": False},
    ),
    ## 9: overlapping registrations -> the one started most recently is used; the older one is still active the day after
    9: make_patient(
        practice_registrations=[registration(date(2015, 1, 1), d(200)), registration(date(2019, 1, 1), d(50))],
        expected_in_population=True,
        expected_columns={"cens_date_dereg_t0": d(50), "cens_bin_dereg_move_t0": True},
    ),

    ###################################################################################
    # C. WINDOWS AT EACH TRIAL START
    ###################################################################################
    ## 10: interaction drug (14 days), palliative care (183 days), hospitalisation on the trial start
    10: make_patient(
        medications=[
            {"date": d(-14), "dmd_code": INTERACTION_DRUG},   # t0 - 14 days -> included in t0
            {"date": d(15), "dmd_code": INTERACTION_DRUG},    # t1 - 15 days -> not in t1
        ],
        clinical_events=[
            {"date": d(-183), "snomedct_code": PALLIATIVE},   # t0 - 183 days -> included in t0
            {"date": d(-154), "snomedct_code": PALLIATIVE},   # t1 - 184 days -> not in t1
        ],
        apcs=[
            {"admission_date": d(55), "discharge_date": d(65)},    # covers t2 (D+60) only
            {"admission_date": d(120), "discharge_date": d(121)},  # admitted on t4 -> hospitalised on t4
            {"admission_date": d(140), "discharge_date": d(150)},  # discharged on t5 -> not hospitalised on t5
        ],
        expected_in_population=True,
        expected_columns={
            "cov_bin_metfin_interaction_t0": True, "cov_bin_metfin_interaction_t1": False,
            "cov_bin_palliative_t0": True, "cov_bin_palliative_t1": False,
            "cov_bin_hosp_t1": False, "cov_bin_hosp_t2": True, "cov_bin_hosp_t3": False,
            "cov_bin_hosp_t4": True, "cov_bin_hosp_t5": False,
        },
    ),

    ###################################################################################
    # E. OUTCOMES AND B2. PREDIABETES
    ###################################################################################
    ## 11: COVID-19 test, hospitalisation (primary diagnosis) and death (contributing cause); prediabetes HbA1c 41.9 vs 42
    11: make_patient(
        sgss_covid_all_tests=[{"specimen_taken_date": d(80), "is_positive": True}],
        apcs=[{"admission_date": d(90), "primary_diagnosis": COVID_ICD10, "all_diagnoses": COVID_ICD10}],
        ons_deaths={"date": d(100), "underlying_cause_of_death": "I219", "cause_of_death_01": COVID_ICD10},
        clinical_events=[
            {"date": d(-400), "snomedct_code": HBA1C, "numeric_value": 41.9},  # below the prediabetes range
            {"date": d(-300), "snomedct_code": HBA1C, "numeric_value": 42},    # lower limit -> prediabetes
        ],
        expected_in_population=True,
        expected_columns={
            "out_date_covid": d(80),
            "out_date_covid_hes": d(90),
            "out_date_covid_hosp": d(90),
            "out_date_covid_death": d(100),
            "out_date_severecovid": d(90),
            "out_date_death": d(100),
            "out_date_dm_death": None,
            "cov_date_prediabetes": d(-300),
        },
    ),
    ## 12: diabetes-related death (padded code E11X as contributing cause); prediabetes HbA1c 48 vs 47.9
    12: make_patient(
        ons_deaths={"date": d(200), "underlying_cause_of_death": "I219", "cause_of_death_02": "E11X"},
        clinical_events=[
            {"date": d(-500), "snomedct_code": HBA1C, "numeric_value": 48},    # above the prediabetes range
            {"date": d(-200), "snomedct_code": HBA1C, "numeric_value": 47.9},  # upper limit -> prediabetes
        ],
        expected_in_population=True,
        expected_columns={
            "out_date_dm_death": d(200),
            "out_date_covid_death": None,
            "cov_date_prediabetes": d(-200),
        },
    ),
}

#######################################################################################
# Write the input files for table_from_file and import the dataset definition
#######################################################################################
# (must be within the workspace for ehrQL; output/ is not committed)
input_dir = "output/seqtrials/assure_test_input"
os.makedirs(input_dir, exist_ok=True)
with open(f"{input_dir}/data_t2dm_qa.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["patient_id", "elig_date_t2dm"])
    for patient_id, patient in test_data.items():
        writer.writerow([patient_id, patient["data_t2dm_qa"]["elig_date_t2dm"].isoformat()])
with open(f"{input_dir}/data_dm_algo.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["patient_id", "ethnicity_cat"])
    for patient_id, patient in test_data.items():
        writer.writerow([patient_id, patient["data_dm_algo"]["ethnicity_cat"]])
os.environ["SEQTRIALS_ASSURE_INPUT_DIR"] = input_dir

from dataset_definition_seqtrials import dataset  # noqa: E402
