from ehrql import (
    case,
    codelist_from_csv,
    create_dataset,
    days,
    when,
)


#######################################################################################
# DEFINE the baseline date based on SARS-CoV-2 infection
#######################################################################################
## COVID-19
covid_primary_care_positive_test = codelist_from_csv("codelists/opensafely-covid-identification-in-primary-care-probable-covid-positive-test.csv", column="CTV3ID")
covid_primary_care_code = codelist_from_csv("codelists/opensafely-covid-identification-in-primary-care-probable-covid-clinical-code.csv", column="CTV3ID")
covid_primary_care_sequelae = codelist_from_csv("codelists/opensafely-covid-identification-in-primary-care-probable-covid-sequelae.csv", column="CTV3ID")
covid_codes = codelist_from_csv("codelists/user-RochelleKnight-confirmed-hospitalised-covid-19.csv", column="code") # only PCR-confirmed! => U071 (covid19 virus identified)


#######################################################################################
# QUALITY ASSURANCES variables
#######################################################################################
# prostate
prostate_cancer_icd10 = codelist_from_csv("codelists/user-RochelleKnight-prostate_cancer_icd10.csv",column="code")
prostate_cancer_snomed_clinical = codelist_from_csv("codelists/user-RochelleKnight-prostate_cancer_snomed.csv",column="code")

# pregnancy
pregnancy_snomed_clinical = codelist_from_csv("codelists/user-RochelleKnight-pregnancy_and_birth_snomed.csv",column="code")

# combined oral contraceptive pill
cocp_dmd = codelist_from_csv("codelists/user-elsie_horne-cocp_dmd.csv",column="dmd_id")

# hormone replacement therapy
hrt_dmd = codelist_from_csv("codelists/user-elsie_horne-hrt_dmd.csv",column="dmd_id")


#######################################################################################
# DEMOGRAPHIC variables
#######################################################################################
# ethnicity
ethnicity_codes = codelist_from_csv(
    "codelists/opensafely-ethnicity.csv",
    column="Code",
    category_column="Grouping_6",
)
primis_covid19_vacc_update_ethnicity = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-eth2001.csv",
    column="code",
    category_column="grouping_6_id",
)


#######################################################################################
# ELIGIBILITY variables
#######################################################################################
## DIABETES
# T1DM
diabetes_type1_ctv3_clinical = codelist_from_csv("codelists/user-hjforbes-type-1-diabetes.csv",column="code")
# T2DM
diabetes_type2_ctv3_clinical = codelist_from_csv("codelists/user-hjforbes-type-2-diabetes.csv",column="code")
# Other or non-specific diabetes
diabetes_other_ctv3_clinical = codelist_from_csv("codelists/user-hjforbes-other-or-nonspecific-diabetes.csv",column="code")
# Gestational diabetes
diabetes_gestational_ctv3_clinical = codelist_from_csv("codelists/user-hjforbes-gestational-diabetes.csv",column="code")
# Type 1 diabetes secondary care
diabetes_type1_icd10 = codelist_from_csv("codelists/opensafely-type-1-diabetes-secondary-care.csv",column="icd10_code")
# Type 2 diabetes secondary care
diabetes_type2_icd10 = codelist_from_csv("codelists/user-r_denholm-type-2-diabetes-secondary-care-bristol.csv",column="code")
# Non-diagnostic diabetes codes
diabetes_diagnostic_ctv3_clinical = codelist_from_csv("codelists/user-hjforbes-nondiagnostic-diabetes-codes.csv",column="code")
# HbA1c
hba1c_new_codes = codelist_from_csv("codelists/user-alainamstutz-hba1c-bristol.csv",column="code")
# Antidiabetic drugs
insulin_snomed_clinical = codelist_from_csv("codelists/opensafely-insulin-medication.csv",column="id")
antidiabetic_drugs_snomed_clinical = codelist_from_csv("codelists/opensafely-antidiabetic-drugs.csv",column="id")
non_metformin_dmd = codelist_from_csv("codelists/user-r_denholm-non-metformin-antidiabetic-drugs_bristol.csv",column="id")

## Prediabetes
prediabetes_snomed = codelist_from_csv("codelists/opensafely-prediabetes-snomed.csv",column="code")

## metformin
metformin_codes = codelist_from_csv("codelists/user-john-tazare-metformin-dmd.csv",column="code")

## metformin allergy
metformin_allergy = codelist_from_csv("codelists/user-alainamstutz-metformin-intolerance-bristol.csv",column="code")

## moderate to severe renal impairment (eGFR of <30ml/min/1.73 m2; stage 4/5)
ckd_snomed_clinical_45 = codelist_from_csv("codelists/nhsd-primary-care-domain-refsets-ckdatrisk1_cod.csv",column="code")
ckd_stage4_icd10 = ["N184"]
ckd_stage5_icd10 = ["N185"]

## advanced decompensated liver cirrhosis
advanced_decompensated_cirrhosis_snomed_codes = codelist_from_csv("codelists/opensafely-condition-advanced-decompensated-cirrhosis-of-the-liver.csv",column="code")
advanced_decompensated_cirrhosis_icd10_codes = codelist_from_csv("codelists/opensafely-condition-advanced-decompensated-cirrhosis-of-the-liver-and-associated-conditions-icd-10.csv",column="code")
# ascitic drainage
ascitic_drainage_snomed_codes = codelist_from_csv("codelists/opensafely-procedure-ascitic-drainage.csv",column="code")

## drug-drug interaction with metformin
metformin_interaction_codes = codelist_from_csv("codelists/user-alainamstutz-metformin-drug-drug-interaction-bristol-dmd.csv",column="code")

## Prior Long COVID diagnosis
long_covid_diagnostic_codes = codelist_from_csv("codelists/opensafely-nice-managing-the-long-term-effects-of-covid-19.csv",column="code")
long_covid_referral_codes = codelist_from_csv("codelists/opensafely-referral-and-signposting-for-long-covid.csv",column="code")
long_covid_assessment_codes = codelist_from_csv("codelists/opensafely-assessment-instruments-and-outcome-measures-for-long-covid.csv",column="code")

post_viral_fatigue_codes = codelist_from_csv("codelists/user-alex-walker-post-viral-syndrome.csv",column="code")


#######################################################################################
# Potential CONFOUNDER variables
#######################################################################################
# smoking
smoking_clear = codelist_from_csv("codelists/opensafely-smoking-clear.csv",
    column="CTV3Code",
    category_column="Category",
)
ever_smoking = codelist_from_csv("codelists/user-alainamstutz-ever-smoking-bristol.csv",column="code")

# Patients in long-stay nursing and residential care
carehome = codelist_from_csv("codelists/primis-covid19-vacc-uptake-longres.csv",column="code")

# obesity
bmi_obesity_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-bmi_obesity_snomed.csv",column="code")
bmi_obesity_icd10 = codelist_from_csv("codelists/user-elsie_horne-bmi_obesity_icd10.csv",column="code")

# acute myocardial infarction
ami_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-ami_snomed.csv",column="code")
ami_icd10 = codelist_from_csv("codelists/user-RochelleKnight-ami_icd10.csv",column="code")
ami_prior_icd10 = codelist_from_csv("codelists/user-elsie_horne-ami_prior_icd10.csv",column="code")

# all strokes
stroke_isch_icd10 = codelist_from_csv("codelists/user-RochelleKnight-stroke_isch_icd10.csv",column="code")
stroke_isch_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-stroke_isch_snomed.csv",column="code")
stroke_sah_hs_icd10 = codelist_from_csv("codelists/user-RochelleKnight-stroke_sah_hs_icd10.csv",column="code")
stroke_sah_hs_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-stroke_sah_hs_snomed.csv",column="code")

# other arterial embolism
other_arterial_embolism_snomed_clinical = codelist_from_csv("codelists/user-tomsrenin-other_art_embol.csv",column="code")
other_arterial_embolism_icd10 = codelist_from_csv("codelists/user-elsie_horne-other_arterial_embolism_icd10.csv",column="code")

# All VTEs in SNOMED
# Portal vein thrombosis
portal_vein_thrombosis_snomed_clinical = codelist_from_csv("codelists/user-tomsrenin-pvt.csv",column="code")
# DVT
dvt_dvt_snomed_clinical = codelist_from_csv("codelists/user-tomsrenin-dvt_main.csv",column="code")
# ICVT
dvt_icvt_snomed_clinical = codelist_from_csv("codelists/user-tomsrenin-dvt_icvt.csv",column="code")
# DVT in pregnancy
dvt_pregnancy_snomed_clinical = codelist_from_csv("codelists/user-tomsrenin-dvt-preg.csv",column="code")
# Other DVT
other_dvt_snomed_clinical = codelist_from_csv("codelists/user-tomsrenin-dvt-other.csv",column="code")
# pulmonary embolism
pe_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-pe_snomed.csv",column="code")
# All VTEs ICD10
portal_vein_thrombosis_icd10 = codelist_from_csv("codelists/user-elsie_horne-portal_vein_thrombosis_icd10.csv",column="code")
dvt_dvt_icd10 = codelist_from_csv("codelists/user-RochelleKnight-dvt_dvt_icd10.csv",column="code")
dvt_icvt_icd10 = codelist_from_csv("codelists/user-elsie_horne-dvt_icvt_icd10.csv",column="code")
dvt_pregnancy_icd10 = codelist_from_csv("codelists/user-elsie_horne-dvt_pregnancy_icd10.csv",column="code")
other_dvt_icd10 = codelist_from_csv("codelists/user-elsie_horne-other_dvt_icd10.csv",column="code")
icvt_pregnancy_icd10 = codelist_from_csv("codelists/user-elsie_horne-icvt_pregnancy_icd10.csv",column="code")
pe_icd10 = codelist_from_csv("codelists/user-RochelleKnight-pe_icd10.csv",column="code")

# heart failure
hf_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-hf_snomed.csv",column="code")
hf_icd10 = codelist_from_csv("codelists/user-RochelleKnight-hf_icd10.csv",column="code")

# angina
angina_snomed_clinical = codelist_from_csv("codelists/user-hjforbes-angina_snomed.csv",column="code")
angina_icd10 = codelist_from_csv("codelists/user-RochelleKnight-angina_icd10.csv",column="code")

# dementia
dementia_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-dementia_snomed.csv",column="code")
dementia_icd10 = codelist_from_csv("codelists/user-elsie_horne-dementia_icd10.csv",column="code")
dementia_vascular_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-dementia_vascular_snomed.csv",column="code")
dementia_vascular_icd10 = codelist_from_csv("codelists/user-elsie_horne-dementia_vascular_icd10.csv",column="code")

# cancer
cancer_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-cancer_snomed.csv",column="code")
cancer_icd10 = codelist_from_csv("codelists/user-elsie_horne-cancer_icd10.csv",column="code")

# hypertension
hypertension_snomed_clinical = codelist_from_csv("codelists/nhsd-primary-care-domain-refsets-hyp_cod.csv",column="code")
hypertension_icd10 = codelist_from_csv("codelists/user-elsie_horne-hypertension_icd10.csv",column="code")
hypertension_drugs_dmd = codelist_from_csv("codelists/user-elsie_horne-hypertension_drugs_dmd.csv",column="dmd_id")

# depression
depression_snomed_clinical = codelist_from_csv("codelists/user-hjforbes-depression-symptoms-and-diagnoses.csv",column="code")
depression_icd10 = codelist_from_csv("codelists/user-kurttaylor-depression_icd10.csv",column="code")

# COPD
copd_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-copd_snomed.csv",column="code")
copd_icd10 = codelist_from_csv("codelists/user-elsie_horne-copd_icd10.csv",column="code")

# liver disease
liver_disease_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-liver_disease_snomed.csv",column="code")
liver_disease_icd10 = codelist_from_csv("codelists/user-elsie_horne-liver_disease_icd10.csv",column="code")

# chronic kidney disease
ckd_snomed_clinical = codelist_from_csv("codelists/user-elsie_horne-ckd_snomed.csv",column="code")
ckd_icd10 = codelist_from_csv("codelists/user-elsie_horne-ckd_icd10.csv",column="code")

# gestational diabetes ICD10
gestationaldm_icd10 = codelist_from_csv("codelists/user-alainamstutz-gestational-diabetes-icd10-bristol.csv",column="code")

# PCOS
pcos_snomed_clinical = codelist_from_csv("codelists/user-alainamstutz-pcos-bristol.csv",column="code")
pcos_icd10 = codelist_from_csv("codelists/user-alainamstutz-pcos-icd10-bristol.csv",column="code")

# key diabetes complications (foot, retino, neuro, nephro)
diabetescomp_snomed_clinical = codelist_from_csv("codelists/user-alainamstutz-diabetes-complications-bristol.csv",column="code")
diabetescomp_icd10 = codelist_from_csv("codelists/user-alainamstutz-diabetes-complications-icd10-bristol.csv",column="code")

# Any HbA1c measurement
hba1c_measurement_snomed = codelist_from_csv("codelists/opensafely-glycated-haemoglobin-hba1c-tests.csv",column="code")

# Any OGTT done
ogtt_measurement_snomed = codelist_from_csv("codelists/user-alainamstutz-ogtt-bristol.csv",column="code")

# Total Cholesterol
cholesterol_snomed = codelist_from_csv("codelists/opensafely-cholesterol-tests-numerical-value.csv",column="code")

# HDL Cholesterol
hdl_cholesterol_snomed = codelist_from_csv("codelists/bristol-hdl-cholesterol.csv",column="code")


#######################################################################################
# OUTCOME variables
#######################################################################################

# covid infection at hosp incl. clin diagnosis without PCR
# covid_codes_incl_clin_diag = codelist_from_csv("codelists/opensafely-covid-identification.csv",column="icd10_code")

# overwrite imported codelist to add 2 additional codes, see blog post here: https://github.com/opensafely/documentation/discussions/1480 
covid_codes_incl_clin_diag = ["U071", "U072", "U109", "U099"]

# covid_emergency = codelist_from_csv(
#     "codelists-opensafely-covid-19-ae-diagnosis-codes.csv",
#     column="Code",
# )
# option without "post-covid syndrome" (> 3 months after infection) based on https://github.com/opensafely/comparative-booster-spring2023/blob/main/analysis/codelists.py 
covid_emergency = ["1240751000000100", "1325171000000109", "1325181000000106"]

# long covid (see in eligibility criteria)
