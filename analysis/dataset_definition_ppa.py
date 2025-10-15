#######################################################################################
# IMPORT
#######################################################################################
## ehrQL functions
from ehrql import (
    create_dataset,
    minimum_of
)

## TPP tables
from ehrql.tables.tpp import (
    clinical_events
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

import project_permissions

## json (for the dates)
import json

## for import of diabetes algo created data
from datetime import date

## import eligible participants and all the data which we have defined already
elig_participants = "output/data/data_processed.arrow"
@table_from_file(elig_participants)
class elig_patients(PatientFrame):
    landmark_date = Series(date)
    elig_date_t2dm = Series(date)

# random seed (ideally use numpy, but currently not working on my local environment)
#import numpy as np 
#np.random.seed(19283) # random seed
import random
random.seed(19283) # random seed

#######################################################################################
# DEFINE the dates
#######################################################################################
with open("output/study_dates.json") as f:
  study_dates = json.load(f)
studyend_date = study_dates["studyend_date"]

#######################################################################################
# INITIALISE the dataset of eligible participants
#######################################################################################
dataset = create_dataset()
dataset.define_population(
    elig_patients.exists_for_patient()
)
dataset.landmark_date = elig_patients.landmark_date
dataset.elig_date_t2dm = elig_patients.elig_date_t2dm

#######################################################################################
# ELIGIBILITY variables
#######################################################################################
## Defined at baseline and not updated thereafter since we run 1 trial only => take from main dataset

#######################################################################################
# INTERVENTION/EXPOSURE variables
#######################################################################################
## I have the first metformin and other antidiabetic prescription in the main dataset (all with such prescriptions before baseline are already excluded)
## We assume someone starting any of these medications remains on the medication
## if we need more detailed on/off prescription patterns (though hard to define), then extract medications as ELD. Not for now.

#######################################################################################
# Covariates and potential confounders
#######################################################################################
### CATEGORY (1): Covariates which are only defined at baseline date (t2dm diagnosis) => take from main dataset:
## Age (increase age in R script)
## Sex
## Ethnicity
## Index of Multiple Deprevation
## Practice registration info
## Healthcare worker
## Smoking status (reconsider to update it)
## Care home
## Consultation rate

### CATEGORY (2): stable time-updated covariates 

## First (if any) acute myocardial infarction, on or before study end date
dataset.cov_date_ami = minimum_of(
    first_matching_event_clinical_snomed_before(ami_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(ami_prior_icd10 + ami_icd10, studyend_date).admission_date
)
## First (if any) stroke, on or before study end date
dataset.cov_date_all_stroke = minimum_of(
    first_matching_event_clinical_snomed_before(stroke_isch_snomed_clinical + stroke_sah_hs_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(stroke_isch_icd10 + stroke_sah_hs_icd10, studyend_date).admission_date
)

## Other arterial embolism, on or before study end date
dataset.cov_date_other_arterial_embolism = minimum_of(
   first_matching_event_clinical_snomed_before(other_arterial_embolism_snomed_clinical, studyend_date).date,
   first_matching_event_apc_before(other_arterial_embolism_icd10, studyend_date).admission_date
)

## Venous thrombolism events, on or before study end date
dataset.cov_date_vte = minimum_of(
  first_matching_event_clinical_snomed_before(portal_vein_thrombosis_snomed_clinical + dvt_dvt_snomed_clinical + dvt_icvt_snomed_clinical + dvt_pregnancy_snomed_clinical + other_dvt_snomed_clinical + pe_snomed_clinical, studyend_date).date,
  first_matching_event_apc_before(portal_vein_thrombosis_icd10 + dvt_dvt_icd10 + dvt_icvt_icd10 + dvt_pregnancy_icd10 + other_dvt_icd10 + icvt_pregnancy_icd10 + pe_icd10, studyend_date).admission_date
)

## Heart failure, on or before study end date
dataset.cov_date_hf = minimum_of(
    first_matching_event_clinical_snomed_before(hf_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(hf_icd10, studyend_date).admission_date
)

## Angina, on or before study end date
dataset.cov_date_angina = minimum_of(
    first_matching_event_clinical_snomed_before(angina_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(angina_icd10, studyend_date).admission_date
)

## Dementia, on or before study end date
dataset.cov_date_dementia = minimum_of(
    first_matching_event_clinical_snomed_before(dementia_snomed_clinical + dementia_vascular_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(dementia_icd10 + dementia_vascular_icd10, studyend_date).admission_date
)

## Cancer, on or before study end date
dataset.cov_date_cancer = minimum_of(
    first_matching_event_clinical_snomed_before(cancer_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(cancer_icd10, studyend_date).admission_date
)

## Hypertension, on or before study end date
dataset.cov_date_hypertension = minimum_of(
    first_matching_event_clinical_snomed_before(hypertension_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(hypertension_icd10, studyend_date).admission_date
)

## Depression, on or before study end date
dataset.cov_date_depression = minimum_of(
    first_matching_event_clinical_snomed_before(depression_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(depression_icd10, studyend_date).admission_date
)

## Chronic obstructive pulmonary disease, on or before study end date
dataset.cov_date_copd = minimum_of(
    first_matching_event_clinical_snomed_before(copd_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(copd_icd10, studyend_date).admission_date
)

## Liver disease,on or before study end date
dataset.cov_date_liver_disease = minimum_of(
    first_matching_event_clinical_snomed_before(liver_disease_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(liver_disease_icd10, studyend_date).admission_date
)

## Chronic kidney disease, on or before study end date
dataset.cov_date_chronic_kidney_disease = minimum_of(
    first_matching_event_clinical_snomed_before(ckd_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(ckd_icd10, studyend_date).admission_date
)

## PCOS, on or before study end date
dataset.cov_date_pcos = minimum_of(
    first_matching_event_clinical_snomed_before(pcos_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(pcos_icd10, studyend_date).admission_date
)

## Prediabetes, on or before study end date 
# Any preDM diagnosis in primary care
dataset.cov_date_prediabetes = first_matching_event_clinical_snomed_before(prediabetes_snomed, studyend_date).date
# Any HbA1c preDM in primary care || reconsider! Does not make much sense, would need to combine with ELD and will have HbA1c's

## Diabetes complications (foot, retino, neuro, nephro), on or before study end date
dataset.cov_date_diabetescomp = minimum_of(
    first_matching_event_clinical_snomed_before(diabetescomp_snomed_clinical, studyend_date).date,
    first_matching_event_apc_before(diabetescomp_icd10, studyend_date).admission_date
)

### CATEGORY (3): dynamic time-updated covariates, between baseline and study end
## Obesity/BMI, HbA1c, lipids (TotChol & HDL)

## BMI: Up to 8 measurements per person should cover 2 years
dataset.cov_date_bmi_1 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.landmark_date, studyend_date),
   minimum_age_at_measurement=16,
).date
dataset.cov_num_bmi_1 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.landmark_date, studyend_date),
   minimum_age_at_measurement=16,
).numeric_value
dataset.cov_date_bmi_2 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_1, studyend_date),
   minimum_age_at_measurement=16,
).date
dataset.cov_num_bmi_2 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_1, studyend_date),
   minimum_age_at_measurement=16,
).numeric_value
dataset.cov_date_bmi_3 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_2, studyend_date),
   minimum_age_at_measurement=16,
).date
dataset.cov_num_bmi_3 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_2, studyend_date),
   minimum_age_at_measurement=16,
).numeric_value
dataset.cov_date_bmi_4 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_3, studyend_date),
   minimum_age_at_measurement=16,
).date
dataset.cov_num_bmi_4 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_3, studyend_date),
   minimum_age_at_measurement=16,
).numeric_value
dataset.cov_date_bmi_5 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_4, studyend_date),
   minimum_age_at_measurement=16,
).date
dataset.cov_num_bmi_5 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_4, studyend_date),
   minimum_age_at_measurement=16,
).numeric_value
dataset.cov_date_bmi_6 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_5, studyend_date),
   minimum_age_at_measurement=16,
).date
dataset.cov_num_bmi_6 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_5, studyend_date),
   minimum_age_at_measurement=16,
).numeric_value
dataset.cov_date_bmi_7 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_6, studyend_date),
   minimum_age_at_measurement=16,
).date
dataset.cov_num_bmi_7 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_6, studyend_date),
   minimum_age_at_measurement=16,
).numeric_value
dataset.cov_date_bmi_8 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_7, studyend_date),
   minimum_age_at_measurement=16,
).date
dataset.cov_num_bmi_8 = first_bmi(
   where=clinical_events.date.is_on_or_between(dataset.cov_date_bmi_7, studyend_date),
   minimum_age_at_measurement=16,
).numeric_value

## HbA1c: Up to 8 measurements per person should cover 2 years
dataset.cov_date_hba1c_1 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.landmark_date, studyend_date).date
dataset.cov_num_hba1c_1 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.landmark_date, studyend_date).numeric_value
dataset.cov_date_hba1c_2 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_1, studyend_date).date
dataset.cov_num_hba1c_2 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_1, studyend_date).numeric_value
dataset.cov_date_hba1c_3 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_2, studyend_date).date
dataset.cov_num_hba1c_3 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_2, studyend_date).numeric_value
dataset.cov_date_hba1c_4 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_3, studyend_date).date
dataset.cov_num_hba1c_4 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_3, studyend_date).numeric_value
dataset.cov_date_hba1c_5 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_4, studyend_date).date
dataset.cov_num_hba1c_5 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_4, studyend_date).numeric_value
dataset.cov_date_hba1c_6 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_5, studyend_date).date
dataset.cov_num_hba1c_6 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_5, studyend_date).numeric_value
dataset.cov_date_hba1c_7 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_6, studyend_date).date
dataset.cov_num_hba1c_7 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_6, studyend_date).numeric_value
dataset.cov_date_hba1c_8 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_7, studyend_date).date
dataset.cov_num_hba1c_8 = first_matching_event_clinical_snomed_between(hba1c_snomed, dataset.cov_date_hba1c_7, studyend_date).numeric_value

## Total Cholesterol: Up to 8 measurements per person should cover 2 years
dataset.cov_date_chol_1 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.landmark_date, studyend_date).date
dataset.cov_num_chol_1 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.landmark_date, studyend_date).numeric_value
dataset.cov_date_chol_2 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_1, studyend_date).date
dataset.cov_num_chol_2 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_1, studyend_date).numeric_value
dataset.cov_date_chol_3 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_2, studyend_date).date
dataset.cov_num_chol_3 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_2, studyend_date).numeric_value
dataset.cov_date_chol_4 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_3, studyend_date).date
dataset.cov_num_chol_4 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_3, studyend_date).numeric_value
dataset.cov_date_chol_5 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_4, studyend_date).date
dataset.cov_num_chol_5 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_4, studyend_date).numeric_value
dataset.cov_date_chol_6 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_5, studyend_date).date
dataset.cov_num_chol_6 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_5, studyend_date).numeric_value
dataset.cov_date_chol_7 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_6, studyend_date).date
dataset.cov_num_chol_7 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_6, studyend_date).numeric_value
dataset.cov_date_chol_8 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_7, studyend_date).date
dataset.cov_num_chol_8 = first_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.cov_date_chol_7, studyend_date).numeric_value

## HDL Cholesterol: Up to 8 measurements per person should cover 2 years
dataset.cov_date_hdl_chol_1 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.landmark_date, studyend_date).date
dataset.cov_num_hdl_chol_1 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.landmark_date, studyend_date).numeric_value
dataset.cov_date_hdl_chol_2 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_1, studyend_date).date
dataset.cov_num_hdl_chol_2 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_1, studyend_date).numeric_value
dataset.cov_date_hdl_chol_3 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_2, studyend_date).date
dataset.cov_num_hdl_chol_3 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_2, studyend_date).numeric_value
dataset.cov_date_hdl_chol_4 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_3, studyend_date).date
dataset.cov_num_hdl_chol_4 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_3, studyend_date).numeric_value
dataset.cov_date_hdl_chol_5 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_4, studyend_date).date
dataset.cov_num_hdl_chol_5 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_4, studyend_date).numeric_value
dataset.cov_date_hdl_chol_6 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_5, studyend_date).date
dataset.cov_num_hdl_chol_6 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_5, studyend_date).numeric_value
dataset.cov_date_hdl_chol_7 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_6, studyend_date).date
dataset.cov_num_hdl_chol_7 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_6, studyend_date).numeric_value
dataset.cov_date_hdl_chol_8 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_7, studyend_date).date
dataset.cov_num_hdl_chol_8 = first_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.cov_date_hdl_chol_7, studyend_date).numeric_value


# ### OTHER (just for our specific case): baseline values of the dynamic time-updated covariates
# ## Total Cholesterol, most recent value, within previous 2 years, on or before elig_date_t2dm
# dataset.cov_num_cholesterol_b = last_matching_event_clinical_snomed_between(cholesterol_snomed, dataset.elig_date_t2dm - days(2*366), dataset.elig_date_t2dm).numeric_value 
# ## HDL Cholesterol, most recent value, within previous 2 years, on or before elig_date_t2dm
# dataset.cov_num_hdl_cholesterol_b = last_matching_event_clinical_snomed_between(hdl_cholesterol_snomed, dataset.elig_date_t2dm - days(2*366), dataset.elig_date_t2dm).numeric_value 


#######################################################################################
# Outcomes, and censoring events
#######################################################################################
# We don't need them for updated eligibility, we may need covid-hosp if we run a PPA on the Long COVID outcome (unless first ever hosp is enough and extracted in main dataset)