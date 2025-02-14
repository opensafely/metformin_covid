fn_diabetes_algorithm <- function(data_extracted){
  data_extracted <- data_extracted %>%
    
    # Preparation Step: Define helper variables needed for step 5
    mutate(tmp_elig_year_first_diabetes_diag = as.integer(format(tmp_elig_date_first_diabetes_diag,"%Y")),
           tmp_age_1st_diag = tmp_elig_year_first_diabetes_diag - qa_num_birth_year,
           tmp_age_1st_diag = replace(tmp_age_1st_diag, which(tmp_age_1st_diag < 0), NA),
           tmp_age_under_35_30_1st_diag = ifelse(!is.na(tmp_age_1st_diag) & 
                                            (tmp_age_1st_diag < 35 & (cov_cat_ethnicity == "White" | cov_cat_ethnicity == "Mixed" | cov_cat_ethnicity == "Other")) | 
                                            (tmp_age_1st_diag < 30), "Yes", "No"),
           # earliest HbA1c date for only those with >=47.5
           tmp_hba1c_date_step7 = as_date(case_when(tmp_elig_num_max_hba1c_mmol_mol >= 47.5 ~ pmin(tmp_elig_date_max_hba1c, na.rm = TRUE))),
           # take the first process code date in those individuals that have 5 or more process codes
           tmp_over5_pocc_step7 = as_date(case_when(tmp_elig_count_poccdm_ctv3 >= 5 ~ pmin(tmp_elig_date_poccdm, na.rm = TRUE)))) %>%    
    
    # Step 1. Any gestational diabetes code?
    mutate(step_1 = ifelse(!is.na(elig_date_gestationaldm), "Yes", "No")) %>%

    # Step 1a. Any T1/ T2 diagnostic codes present? Denominator for step 1a is those with yes to step 1
    mutate(step_1a = ifelse(step_1 == "Yes" &
                              (!is.na(elig_date_t1dm) | !is.na(elig_date_t2dm)), "Yes",
                            ifelse(step_1 == "Yes" &
                                     is.na(elig_date_t1dm) &
                                     is.na(elig_date_t2dm), "No", NA))) %>% # NA will never be fulfilled as long as step_1=="Yes"

    # Step 2. Non-metformin antidiabetic denominator for step 2: no to step 1 or yes to step 1a
    mutate(step_2 = ifelse((step_1 == "No" | step_1a == "Yes" ) &
                             !is.na(tmp_elig_date_nonmetform_drugs_snomed), "Yes",
                           ifelse((step_1 == "No" | step_1a == "Yes") &
                                    is.na(tmp_elig_date_nonmetform_drugs_snomed), "No", NA))) %>% # NA will never be fulfilled

    # Step 3. Type 1 code in the absence of type 2 code? denominator for step 3: no to step 2
    mutate(step_3 = ifelse(step_2=="No" &
                             !is.na(elig_date_t1dm) &
                             is.na(elig_date_t2dm), "Yes",
                           ifelse(step_2 == "No", "No", NA))) %>% # NA will never be fulfilled as long as step_2=="No"

    # Step 4. Type 2 code in the absence of type 1 code denominator for step 3: no to step 3
    mutate(step_4 = ifelse(step_3 == "No" &
                             is.na(elig_date_t1dm) &
                             !is.na(elig_date_t2dm), "Yes",
                           ifelse(step_3 == "No", "No", NA))) %>% # NA will never be fulfilled

    # Step 5. Aged <35yrs (or <30 yrs for SAs and AFCS) at first diagnostic code? denominator for step 5: no to step 4
    mutate(step_5 = ifelse(step_4 == "No" &
                             tmp_age_under_35_30_1st_diag == "Yes", "Yes", ### includes many codes (incl. gestational DM, but excluded in Step 1): elig_date_t2dm, elig_date_t1dm, elig_date_otherdm, elig_date_gestationaldm, tmp_elig_date_poccdm, tmp_elig_date_diabetes_medication, tmp_elig_date_nonmetform_drugs_snomed
                           ifelse(step_4 == "No" &
                                    tmp_age_under_35_30_1st_diag == "No", "No", NA))) %>%
    mutate(step_5 = ifelse(step_5 == "No" |
                             is.na(step_5) & step_4 == "No", "No", "Yes")) %>% # => step_5 will never be NA, as an extra security step. However, these two lines are not really needed; "tmp_age_under_35_30_1st_diag" cannot be NA.

    # Step 6. Type 1 and type 2 codes present? denominator for step 6: no to step 5
    mutate(step_6 = ifelse(step_5 == "No" &
                             !is.na(elig_date_t1dm) &
                             !is.na(elig_date_t2dm), "Yes", # step_6 == Yes does not contain any NA
                           ifelse(step_5 == "No" &
                                    (is.na(elig_date_t1dm) |
                                       is.na(elig_date_t2dm)), "No", NA))) %>% # NA will never be fulfilled

    # Step 6a. Type 1 only reported in primary care? denominator for step 6a: yes to step 6
    mutate(step_6a = ifelse(step_6 == "Yes" &
                              !is.na(tmp_elig_date_t1dm_ctv3) &
                              is.na(tmp_elig_date_t2dm_ctv3), "Yes",
                            ifelse(step_6 == "Yes", "No", NA))) %>% # NA will never be fulfilled

    # Step 6b. Type 2 only reported in primary care? denominator for step 6b: no to step 6a
    mutate(step_6b = ifelse(step_6a == "No" &
                              is.na(tmp_elig_date_t1dm_ctv3) &
                              !is.na(tmp_elig_date_t2dm_ctv3), "Yes",
                            ifelse(step_6a == "No", "No", NA))) %>% # NA will never be fulfilled

    # Step 6c. Number of type 1 codes>number of type 2 codes? denominator for step 6c: no to step 6b
    mutate(step_6c = ifelse(step_6b == "No" &
                              tmp_elig_count_t1dm > tmp_elig_count_t2dm, "Yes", # count variable cannot be NA, according to ehrQL: count_for_patient() "Note this will be 0 rather than NULL if the patient has no rows at all in the frame." https://docs.opensafely.org/ehrql/reference/language/#PatientFrame.count_for_patient
                            ifelse(step_6b == "No" &
                                     tmp_elig_count_t1dm <= tmp_elig_count_t2dm, "No", NA))) %>% # NA will never be fulfilled

    # Step 6d. Number of type 2 codes>number of type 1 codes denominator for step 6d: no to step 6c
    mutate(step_6d = ifelse(step_6c == "No" &
                              tmp_elig_count_t2dm > tmp_elig_count_t1dm, "Yes",
                            ifelse(step_6c == "No" &
                                     tmp_elig_count_t2dm <= tmp_elig_count_t1dm, "No", NA))) %>% # NA will never be fulfilled

    # Step 6e. Type 2 code most recent? denominator for step 6e: no to step 6d
    mutate(step_6e = ifelse(step_6d == "No" &
                              elig_date_t2dm > elig_date_t1dm, "Yes",
                            ifelse(step_6d == "No" &
                                     elig_date_t2dm <= elig_date_t1dm, "No", NA))) %>% ### NA will not be fulfilled. Type 1 overrules Type 2 in case of same date.

    # Step 7. Diabetes medication or >5 process of care codes or HbA1c >= 6.5 (47.5 mmol/mol)? denominator for step 7: no to step 6
    mutate(step_7 = ifelse(step_6 == "No" & # includes missing in elig_date_t1dm or elig_date_t2dm in step_6
                             ((!is.na(tmp_elig_date_diabetes_medication)) |
                                (tmp_elig_num_max_hba1c_mmol_mol >= 47.5) |
                                (tmp_elig_count_poccdm_ctv3 >= 5)), "Yes",
                           ifelse(step_6=="No" , "No", NA))) %>% # NA will never be fulfilled. Those with missing elig_date_t1dm/t2dm and missing any other evidence for DM (medication, hba1c, procedure codes) will be classified as no diabetes (step_7 == "No")

    # Create Diabetes Variable
    mutate(elig_cat_diabetes = ifelse(step_1 == "No" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
                                       step_5 == "No" & step_6 == "No" & step_7 == "No" |
                                       step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
                                       step_5 == "No" & step_6 == "No" & step_7 == "No" ,
                                     "DM unlikely",
                                     ifelse(step_1 == "No" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
                                              step_5 == "No" & step_6 == "No" & step_7 == "Yes" |
                                              step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
                                              step_5 == "No" & step_6 == "No" & step_7 == "Yes",
                                            "DM_other",
                                            ifelse(step_1 == "No" & step_2 == "Yes" |
                                                     step_1 == "Yes" & step_1a == "Yes" & step_2 == "Yes" |
                                                     step_1 == "No" & step_2 == "No" & step_3 == "No" & step_4 == "Yes" |
                                                     step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3 == "No" & step_4 == "Yes" |
                                                     step_1 == "No" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
                                                     step_5 == "No" & step_6 == "Yes" & step_6a == "No" & step_6b=="Yes" |
                                                     step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
                                                     step_5 == "No" & step_6 == "Yes" & step_6a == "No" & step_6b=="Yes" |
                                                     step_1 == "No" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
                                                     step_5 == "No" & step_6 == "Yes" & step_6a == "No" & step_6b=="No" &
                                                     step_6c == "No" & step_6d == "Yes" |
                                                     step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
                                                     step_5 == "No" & step_6 == "Yes" & step_6a == "No" & step_6b=="No" &
                                                     step_6c == "No" & step_6d == "Yes" |
                                                     step_1 == "No" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
                                                     step_5 == "No" & step_6 == "Yes" & step_6a == "No" & step_6b=="No" &
                                                     step_6c == "No" & step_6d == "No" & step_6e == "Yes" |
                                                     step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
                                                     step_5 == "No" & step_6 == "Yes" & step_6a == "No" & step_6b=="No" &
                                                     step_6c == "No" & step_6d == "No" & step_6e == "Yes",
                                                   "T2DM",
                                                   ifelse(step_1 == "No" & step_2 == "No" & step_3=="Yes" |
                                                            step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3=="Yes" |
                                                            step_1 == "No" & step_2 == "No" & step_3 =="No" & step_4 == "No" & step_5 == "Yes" |
                                                            step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3 =="No" & step_4 == "No" &
                                                            step_5 == "Yes" |
                                                            step_1 == "No" & step_2 == "No" & step_3 =="No" & step_4 == "No" & step_5 == "No" &
                                                            step_6 == "Yes" & step_6a == "Yes" |
                                                            step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3 =="No" & step_4 == "No" &
                                                            step_5 == "No" &
                                                            step_6 == "Yes" & step_6a == "Yes" |
                                                            step_1 == "No" & step_2 == "No" & step_3 =="No" & step_4 == "No" & step_5 == "No" &
                                                            step_6 == "Yes" & step_6a == "No" & step_6b == "No" & step_6c == "Yes" |
                                                            step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3 =="No" & step_4 == "No" &
                                                            step_5 == "No" &
                                                            step_6 == "Yes" & step_6a == "No" & step_6b == "No" & step_6c == "Yes" |
                                                            step_1 == "No" & step_2 == "No" & step_3 =="No" & step_4 == "No" & step_5 == "No" &
                                                            step_6 == "Yes" & step_6a == "No" & step_6b == "No" & step_6c == "No" &
                                                            step_6d == "No" & step_6e == "No" |
                                                            step_1 == "Yes" & step_1a == "Yes" & step_2 == "No" & step_3 =="No" & step_4 == "No" & step_5 == "No" &
                                                            step_6 == "Yes" & step_6a == "No" & step_6b == "No" & step_6c == "No" &
                                                            step_6d == "No" & step_6e == "No",
                                                          "T1DM",
                                                          ifelse(step_1 == "Yes" & step_1a == "No", "GDM", NA)))))) %>%
    # replace NAs with None (no diabetes)
    mutate_at(vars(elig_cat_diabetes), ~replace_na(., "None")) %>%

    # Define incident diabetes date variables
    # Uses diabetes category from algorithm above and date of first diabetes related code

    # remove original diabetes variables to avoid duplication
    dplyr::select(- elig_date_t1dm, - elig_date_t2dm, - elig_date_otherdm, - elig_date_gestationaldm) %>%
           # GESTATIONAL
    mutate(elig_date_gestationaldm = as_date(case_when(elig_cat_diabetes == "GDM" ~ tmp_elig_date_first_diabetes_diag,
                                                       TRUE ~ as.Date(NA))),
           # T2DM
           elig_date_t2dm = as_date(case_when(elig_cat_diabetes == "T2DM" ~ tmp_elig_date_first_diabetes_diag,
                                              TRUE ~ as.Date(NA))),
           # T1DM
           elig_date_t1dm = as_date(case_when(elig_cat_diabetes == "T1DM" ~ tmp_elig_date_first_diabetes_diag,
                                              TRUE ~ as.Date(NA))),
           # OTHER
           elig_date_otherdm = as_date(case_when(elig_cat_diabetes == "DM_other" ~ pmin(tmp_hba1c_date_step7, tmp_over5_pocc_step7, na.rm = TRUE),
                                                 TRUE ~ as.Date(NA)))) %>%
    # drop unnecessary helper variables
    dplyr::select(-contains("tmp"), -contains("step"))
}
