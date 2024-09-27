diabetes_algo <- function(data_extracted){
  data_extracted <- data_extracted %>%

    # Step 1. Any gestational diabetes code?
    mutate(step_1 = ifelse(!is.na(cov_date_gestationaldm), "Yes", "No")) %>%

    # Step 1a. Any T1/ T2 diagnostic codes present? Denominator for step 1a is those with yes to step 1
    mutate(step_1a = ifelse(step_1 == "Yes" &
                              (!is.na(cov_date_t1dm) | !is.na(cov_date_t2dm)), "Yes",
                            ifelse(step_1 == "Yes" &
                                     is.na(cov_date_t1dm) &
                                     is.na(cov_date_t2dm), "No", NA))) %>% # NA will never be fulfilled as long as step_1=="Yes"

    # Step 2. Non-metformin antidiabetic denominator for step 2: no to step 1 or yes to step 1a
    mutate(step_2 = ifelse((step_1 == "No" | step_1a == "Yes" ) &
                             !is.na(tmp_cov_date_nonmetform_drugs_snomed), "Yes",
                           ifelse((step_1 == "No" | step_1a == "Yes") &
                                    is.na(tmp_cov_date_nonmetform_drugs_snomed), "No", NA))) %>% # NA will never be fulfilled

    # Step 3. Type 1 code in the absence of type 2 code? denominator for step 3: no to step 2
    mutate(step_3 = ifelse(step_2=="No" &
                             !is.na(cov_date_t1dm) &
                             is.na(cov_date_t2dm), "Yes",
                           ifelse(step_2 == "No", "No", NA))) %>% # NA will never be fulfilled as long as step_2=="No"

    # Step 4. Type 2 code in the absence of type 1 code denominator for step 3: no to step 3
    mutate(step_4 = ifelse(step_3 == "No" &
                             is.na(cov_date_t1dm) &
                             !is.na(cov_date_t2dm), "Yes",
                           ifelse(step_3 == "No", "No", NA))) %>% # NA will never be fulfilled

    # Step 5. Aged <35yrs (or <30 yrs for SAs and AFCS) at first diagnostic code? denominator for step 5: no to step 4
    mutate(step_5 = ifelse(step_4 == "No" &
                             tmp_age_under_35_30_1st_diag == "Yes", "Yes", ### includes NA and many codes (incl. gestational DM, but excluded in Step 1): cov_date_t2dm, cov_date_t1dm, cov_date_otherdm, cov_date_gestationaldm, tmp_cov_date_poccdm, tmp_cov_date_diabetes_medication, tmp_cov_date_nonmetform_drugs_snomed
                           ifelse(step_4 == "No" &
                                    tmp_age_under_35_30_1st_diag == "No", "No", NA))) %>%
    mutate(step_5 = ifelse(step_5 == "No" |
                             is.na(step_5) & step_4 == "No", "No", "Yes")) %>% # => step_5 will never be NA

    # Step 6. Type 1 and type 2 codes present? denominator for step 6: no to step 5
    mutate(step_6 = ifelse(step_5 == "No" &
                             !is.na(cov_date_t1dm) &
                             !is.na(cov_date_t2dm), "Yes", # step_6 == Yes does not contain any missing in cov_date_t1dm & _t2dm
                           ifelse(step_5 == "No" &
                                    (is.na(cov_date_t1dm) |
                                       is.na(cov_date_t2dm)), "No", NA))) %>% # NA will never be fulfilled

    # Step 6a. Type 1 only reported in primary care. denominator for step 6a: "yes" to step 6 ("adapted") ###
    mutate(step_6a = ifelse(step_6 == "Yes" &
                              !is.na(tmp_cov_date_t1dm_ctv3) &
                              is.na(tmp_cov_date_t2dm_ctv3), "Yes",
                            ifelse(step_6 == "Yes", "No", NA))) %>% # NA will never be fulfilled

    # Step 6b. Type 2 only reported in primary care. denominator for step 6b: no to step 6"a" ("adapted") ###
    mutate(step_6b = ifelse(step_6a == "No" &
                              is.na(tmp_cov_date_t1dm_ctv3) &
                              !is.na(tmp_cov_date_t2dm_ctv3), "Yes",
                            ifelse(step_6a == "No", "No", NA))) %>% # NA will never be fulfilled

    # Step 6c. Number of type 1 codes>number of type 2 codes? denominator for step 6c: no to step 6b
    mutate(step_6c = ifelse(step_6b == "No" &
                              tmp_cov_count_t1dm > tmp_cov_count_t2dm, "Yes",
                            ifelse(step_6b == "No" &
                                     tmp_cov_count_t1dm <= tmp_cov_count_t2dm, "No", NA))) %>% # NA will never be fulfilled

    # Step 6d. Number of type 2 codes>number of type 1 codes denominator for step 6d: no to step 6c
    mutate(step_6d = ifelse(step_6c == "No" &
                              tmp_cov_count_t2dm > tmp_cov_count_t1dm, "Yes",
                            ifelse(step_6c == "No" &
                                     tmp_cov_count_t2dm <= tmp_cov_count_t1dm, "No", NA))) %>% # NA will never be fulfilled

    # Step 6e. Type 2 code most recent? denominator for step 6e: no to step 6d
    mutate(step_6e = ifelse(step_6d == "No" &
                              cov_date_t2dm > cov_date_t1dm, "Yes",
                            ifelse(step_6d == "No" &
                                     cov_date_t2dm < cov_date_t1dm, "No", NA))) %>% ### NA will not be fulfilled except if cov_date_t2dm = cov_date_t1dm

    # Step 7. Diabetes medication or >5 process of care codes or HbA1c>=6.5? denominator for step 7: no to step 6
    mutate(step_7 = ifelse(step_6 == "No" & # includes missing in cov_date_t1dm or cov_date_t2dm in step_6
                             ((!is.na(tmp_cov_date_diabetes_medication)) |
                                (tmp_cov_num_max_hba1c_mmol_mol >= 47.5) |
                                (tmp_cov_count_poccdm_ctv3 >= 5)), "Yes",
                           ifelse(step_6=="No" , "No", NA))) %>% # NA will never be fulfilled. Those with missing cov_date_t1dm/t2dm and missing any other evidence for DM (medication, hba1c, procedure codes) will be classified as no diabetes (step_7 == "No")

    # Create Diabetes Variable
    mutate(cov_cat_diabetes = ifelse(step_1 == "No" & step_2 == "No" & step_3 == "No" & step_4 == "No" &
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
    mutate_at(vars(cov_cat_diabetes), ~replace_na(., "None")) %>%

    # Define incident diabetes date variables needed for cox analysis
    # Uses diabetes category from algorithm above and date of "latest" diabetes related code ("adapted") ###

    # remove old diabetes variables to avoid duplication / confusion
    dplyr::select(- cov_date_t1dm, - cov_date_t2dm, - cov_date_otherdm, - cov_date_gestationaldm) %>%
    # GESTATIONAL
    mutate(cov_date_gestationaldm = as_date(case_when(cov_cat_diabetes == "GDM" ~ tmp_cov_date_latest_diabetes_diag)),
           # T2DM
           cov_date_t2dm = as_date(case_when(cov_cat_diabetes == "T2DM" ~ tmp_cov_date_latest_diabetes_diag)),
           # T1DM
           cov_date_t1dm = as_date(case_when(cov_cat_diabetes == "T1DM" ~ tmp_cov_date_latest_diabetes_diag)),
           # OTHER
           cov_date_otherdm = as_date(case_when(cov_cat_diabetes == "DM_other" ~ pmin(tmp_hba1c_date_step7, tmp_over5_pocc_step7, na.rm = TRUE))))
}
