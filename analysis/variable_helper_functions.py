#######################################################################################
# FUNCTIONS for the dataset_definition
# Adapted from https://github.com/grpEHR/post-covid-events-ehrQL/blob/dataset-definition/analysis/variable_helper_functions.py
# Credits to https://github.com/ZoeMZou 
#######################################################################################
from ehrql.tables.tpp import (
    apcs, 
    opa_diag, 
    clinical_events, 
    medications, 
    ons_deaths,
    emergency_care_attendances,
    patients
)
from ehrql import days # for BMI function
from ehrql.codes import CTV3Code # for BMI function

### HELPER functions, based on https://github.com/opensafely/comparative-booster-spring2023/blob/main/analysis/dataset_definition.py
import operator
from functools import reduce
def any_of(conditions):
    return reduce(operator.or_, conditions)

#######################################################################################
### ANY HISTORY of ... and give latest ... (including baseline_date) 
#######################################################################################
## In PRIMARY CARE
# CTV3/Read
def last_matching_event_clinical_ctv3_before(codelist, baseline_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.ctv3_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_before(baseline_date))
        .sort_by(clinical_events.date)
        .last_for_patient()
    )
# Snomed
def last_matching_event_clinical_snomed_before(codelist, baseline_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.snomedct_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_before(baseline_date))
        .sort_by(clinical_events.date)
        .last_for_patient()
    )
# Medication
def last_matching_med_dmd_before(codelist, baseline_date, where=True):
    return(
        medications.where(where)
        .where(medications.dmd_code.is_in(codelist))
        .where(medications.date.is_on_or_before(baseline_date))
        .sort_by(medications.date)
        .last_for_patient()
    )

## In SECONDARY CARE (Hospital Episodes)
def last_matching_event_apc_before(codelist, baseline_date, only_prim_diagnoses=False, where=True):
    query = apcs.where(where).where(apcs.admission_date.is_on_or_before(baseline_date))
    if only_prim_diagnoses:
        # If set to True, then check only primary diagnosis field
        query = query.where(
            apcs.primary_diagnosis.is_in(codelist)
        )
    else:
        # Else, check all diagnoses (default, i.e. when only_prim_diagnoses argument not defined)
        query = query.where(apcs.all_diagnoses.contains_any_of(codelist))
    return query.sort_by(apcs.admission_date).last_for_patient()

## In OUTPATIENT CARE
def last_matching_event_opa_before(codelist, baseline_date, where=True):
    return(
        opa_diag.where(where)
        .where(opa_diag.primary_diagnosis_code.is_in(codelist) | opa_diag.secondary_diagnosis_code_1.is_in(codelist))
        .where(opa_diag.appointment_date.is_on_or_before(baseline_date))
        .sort_by(opa_diag.appointment_date)
        .last_for_patient()
    )

## In EMERGENCY CARE
def last_matching_event_ec_snomed_before(codelist, baseline_date, where=True):
    conditions = [
        getattr(emergency_care_attendances, column_name).is_in(codelist)
        for column_name in ([f"diagnosis_{i:02d}" for i in range(1, 25)])
    ]
    return(
        emergency_care_attendances.where(where)
        .where(any_of(conditions))
        .where(emergency_care_attendances.arrival_date.is_before(baseline_date))
        .sort_by(emergency_care_attendances.arrival_date)
        .last_for_patient()
    )

## DEATH
def matching_death_before(codelist, baseline_date, where=True):
    conditions = [
        getattr(ons_deaths, column_name).is_in(codelist)
        for column_name in (["underlying_cause_of_death"]+[f"cause_of_death_{i:02d}" for i in range(1, 16)])
    ]
    return any_of(conditions) & ons_deaths.date.is_before(baseline_date)

#######################################################################################
### HISTORY of ... in past ... days/months/years ... and give latest (including baseline_date)
#######################################################################################
## In PRIMARY CARE
# CTV3/Read
def last_matching_event_clinical_ctv3_between(codelist, start_date, baseline_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.ctv3_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_between(start_date, baseline_date))
        .sort_by(clinical_events.date)
        .last_for_patient()
    )
# Snomed
def last_matching_event_clinical_snomed_between(codelist, start_date, baseline_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.snomedct_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_between(start_date, baseline_date))
        .sort_by(clinical_events.date)
        .last_for_patient()
    )
# Medication
def last_matching_med_dmd_between(codelist, start_date, baseline_date, where=True):
    return(
        medications.where(where)
        .where(medications.dmd_code.is_in(codelist))
        .where(medications.date.is_on_or_between(start_date, baseline_date))
        .sort_by(medications.date)
        .last_for_patient()
    )

## In SECONDARY CARE (Hospital Episodes)
def last_matching_event_apc_between(codelist, start_date, baseline_date, only_prim_diagnoses=False, where=True):
    query = apcs.where(where).where(apcs.admission_date.is_on_or_between(start_date, baseline_date))
    if only_prim_diagnoses:
        # If set to True, then check only primary diagnosis field
        query = query.where(
            apcs.primary_diagnosis.is_in(codelist)
        )
    else:
        # Else, check all diagnoses (default, i.e. when only_prim_diagnoses argument not defined)
        query = query.where(apcs.all_diagnoses.contains_any_of(codelist))
    return query.sort_by(apcs.admission_date).last_for_patient()

## In OUTPATIENT CARE
def last_matching_event_opa_between(codelist, start_date, baseline_date, where=True):
    return(
        opa_diag.where(where)
        .where(opa_diag.primary_diagnosis_code.is_in(codelist) | opa_diag.secondary_diagnosis_code_1.is_in(codelist))
        .where(opa_diag.appointment_date.is_on_or_between(start_date, baseline_date))
        .sort_by(opa_diag.appointment_date)
        .last_for_patient()
    )

## In EMERGENCY CARE
def last_matching_event_ec_snomed_between(codelist, start_date, baseline_date, where=True):
    conditions = [
        getattr(emergency_care_attendances, column_name).is_in(codelist)
        for column_name in ([f"diagnosis_{i:02d}" for i in range(1, 25)])
    ]
    return(
        emergency_care_attendances.where(where)
        .where(any_of(conditions))
        .where(emergency_care_attendances.arrival_date.is_on_or_between(start_date, baseline_date))
        .sort_by(emergency_care_attendances.arrival_date)
        .last_for_patient()
    )

#######################################################################################
### COUNT all prior events (including baseline_date)
#######################################################################################
## In PRIMARY CARE
# CTV3/Read
def count_matching_event_clinical_ctv3_before(codelist, baseline_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.ctv3_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_before(baseline_date))
        .count_for_patient()
    )
# Snomed
def count_matching_event_clinical_snomed_before(codelist, baseline_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.snomedct_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_before(baseline_date))
        .count_for_patient()
    )

## In SECONDARY CARE (Hospital Episodes)
def count_matching_event_apc_before(codelist, baseline_date, only_prim_diagnoses=False, where=True):
    query = apcs.where(where).where(apcs.admission_date.is_on_or_before(baseline_date))
    if only_prim_diagnoses:
        # If set to True, then check only primary diagnosis field
        query = query.where(
            apcs.primary_diagnosis.is_in(codelist)
        )
    else:
        # Else, check all diagnoses (default, i.e. when only_prim_diagnoses argument not defined)
        query = query.where(apcs.all_diagnoses.contains_any_of(codelist))
    return query.count_for_patient()

## In OUTPATIENT CARE
def count_matching_event_opa_before(codelist, baseline_date, where=True):
    return(
        opa_diag.where(where)
        .where(opa_diag.primary_diagnosis_code.is_in(codelist) | opa_diag.secondary_diagnosis_code_1.is_in(codelist))
        .where(opa_diag.appointment_date.is_on_or_before(baseline_date))
        .count_for_patient()
    )

#######################################################################################
### COUNT all events in between ... and ... (including baseline_date)
#######################################################################################
## In PRIMARY CARE
# CTV3/Read
def count_matching_event_clinical_ctv3_between(codelist, start_date, baseline_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.ctv3_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_between(start_date, baseline_date))
        .count_for_patient()
    )
# Snomed
def count_matching_event_clinical_snomed_between(codelist, start_date, baseline_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.snomedct_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_between(start_date, baseline_date))
        .count_for_patient()
    )

## In SECONDARY CARE (Hospital Episodes)
def count_matching_event_apc_between(codelist, start_date, baseline_date, only_prim_diagnoses=False, where=True):
    query = apcs.where(where).where(apcs.admission_date.is_on_or_between(start_date, baseline_date))
    if only_prim_diagnoses:
        # If set to True, then check only primary diagnosis field
        query = query.where(
            apcs.primary_diagnosis.is_in(codelist)
        )
    else:
        # Else, check all diagnoses (default, i.e. when only_prim_diagnoses argument not defined)
        query = query.where(apcs.all_diagnoses.contains_any_of(codelist))
    return query.count_for_patient()

## In OUTPATIENT CARE
def count_matching_event_opa_between(codelist, start_date, baseline_date, where=True):
    return(
        opa_diag.where(where)
        .where(opa_diag.primary_diagnosis_code.is_in(codelist) | opa_diag.secondary_diagnosis_code_1.is_in(codelist))
        .where(opa_diag.appointment_date.is_on_or_between(start_date, baseline_date))
        .count_for_patient()
    )

#######################################################################################
### ANY HISTORY of ... and give first ... (including baseline_date) 
#######################################################################################
## In PRIMARY CARE
# CTV3/Read
def first_matching_event_clinical_ctv3_before(codelist, baseline_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.ctv3_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_before(baseline_date))
        .sort_by(clinical_events.date)
        .first_for_patient()
    )
# Snomed
def first_matching_event_clinical_snomed_before(codelist, baseline_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.snomedct_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_before(baseline_date))
        .sort_by(clinical_events.date)
        .first_for_patient()
    )
# Medication
def first_matching_med_dmd_before(codelist, baseline_date, where=True):
    return(
        medications.where(where)
        .where(medications.dmd_code.is_in(codelist))
        .where(medications.date.is_on_or_before(baseline_date))
        .sort_by(medications.date)
        .first_for_patient()
    )

## In SECONDARY CARE (Hospital Episodes)
def first_matching_event_apc_before(codelist, baseline_date, only_prim_diagnoses=False, where=True):
    query = apcs.where(where).where(apcs.admission_date.is_on_or_before(baseline_date))
    if only_prim_diagnoses:
         # If set to True, then check only primary diagnosis field
        query = query.where(
            apcs.primary_diagnosis.is_in(codelist)
        )
    else:
        # Else, check all diagnoses (default, i.e. when only_prim_diagnoses argument not defined)
        query = query.where(apcs.all_diagnoses.contains_any_of(codelist))
    return query.sort_by(apcs.admission_date).first_for_patient()

## In EMERGENCY CARE
def first_matching_event_ec_snomed_before(codelist, baseline_date, where=True):
    conditions = [
        getattr(emergency_care_attendances, column_name).is_in(codelist)
        for column_name in ([f"diagnosis_{i:02d}" for i in range(1, 25)])
    ]
    return(
        emergency_care_attendances.where(where)
        .where(any_of(conditions))
        .where(emergency_care_attendances.arrival_date.is_before(baseline_date))
        .sort_by(emergency_care_attendances.arrival_date)
        .first_for_patient()
    )

#######################################################################################
### Any future events (including baseline_date and study end_date)
#######################################################################################
## In PRIMARY CARE
# CTV3/Read
def first_matching_event_clinical_ctv3_between(codelist, baseline_date, end_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.ctv3_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_between(baseline_date, end_date))
        .sort_by(clinical_events.date)
        .first_for_patient()
    )
# Snomed
def first_matching_event_clinical_snomed_between(codelist, baseline_date, end_date, where=True):
    return(
        clinical_events.where(where)
        .where(clinical_events.snomedct_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_between(baseline_date, end_date))
        .sort_by(clinical_events.date)
        .first_for_patient()
    )
# Medication
def first_matching_med_dmd_between(codelist, baseline_date, end_date, where=True):
    return(
        medications.where(where)
        .where(medications.dmd_code.is_in(codelist))
        .where(medications.date.is_on_or_between(baseline_date, end_date))
        .sort_by(medications.date)
        .first_for_patient()
    )

## In SECONDARY CARE
def first_matching_event_apc_between(codelist, baseline_date, end_date, only_prim_diagnoses=False, where=True):
    query = apcs.where(where).where(apcs.admission_date.is_on_or_between(baseline_date, end_date))
    if only_prim_diagnoses:
        # If set to True, then check only primary diagnosis field
        query = query.where(
            apcs.primary_diagnosis.is_in(codelist)
        )
    else:
        # Else, check all diagnoses (default, i.e. when only_prim_diagnoses argument not defined)
        query = query.where(apcs.all_diagnoses.contains_any_of(codelist))
    return query.sort_by(apcs.admission_date).first_for_patient()

## In OUTPATIENT CARE
def first_matching_event_opa_between(codelist, baseline_date, end_date, where=True):
    return(
        opa_diag.where(where)
        .where(opa_diag.primary_diagnosis_code.is_in(codelist) | opa_diag.secondary_diagnosis_code_1.is_in(codelist))
        .where(opa_diag.appointment_date.is_on_or_between(baseline_date, end_date))
        .sort_by(opa_diag.appointment_date)
        .first_for_patient()
    )

## In EMERGENCY CARE
def first_matching_event_ec_snomed_between(codelist, baseline_date, end_date, where=True):
    conditions = [
        getattr(emergency_care_attendances, column_name).is_in(codelist)
        for column_name in ([f"diagnosis_{i:02d}" for i in range(1, 25)])
    ]
    return(
        emergency_care_attendances.where(where)
        .where(any_of(conditions))
        .where(emergency_care_attendances.arrival_date.is_on_or_between(baseline_date, end_date))
        .sort_by(emergency_care_attendances.arrival_date)
        .first_for_patient()
    )

## DEATH
def matching_death_between(codelist, baseline_date, end_date, where=True):
    conditions = [
        getattr(ons_deaths, column_name).is_in(codelist)
        for column_name in (["underlying_cause_of_death"]+[f"cause_of_death_{i:02d}" for i in range(1, 16)])
    ]
    return any_of(conditions) & ons_deaths.date.is_on_or_between(baseline_date, end_date)

### CAUSES of DEATH without any date restrictions
def cause_of_death_matches(codelist):
    conditions = [
        getattr(ons_deaths, column_name).is_in(codelist)
        for column_name in (["underlying_cause_of_death"]+[f"cause_of_death_{i:02d}" for i in range(1, 16)])
    ]
    return any_of(conditions)


### BMI calculation
def most_recent_bmi(*, minimum_age_at_measurement, where=True):
    age_threshold = patients.date_of_birth + days(
        # This is obviously inexact but, given that the dates of birth are rounded to
        # the first of the month anyway, there's no point trying to be more accurate
        int(365.25 * minimum_age_at_measurement)
    )
    return (
        # This captures just explicitly recorded BMI observations rather than attempting
        # to calculate it from height and weight measurements. Investigation has shown
        # this to have no real benefit it terms of coverage or accuracy.
        clinical_events.where(where)
        .where(clinical_events.ctv3_code == CTV3Code("22K.."))
        .where(clinical_events.date >= age_threshold)
        .sort_by(clinical_events.date)
        .last_for_patient()
    )