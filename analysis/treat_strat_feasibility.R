####
## This script does the following:
# 1. Import the full dataset that was processed through data_process.R
# 2. Add various treatment regimen patterns, cum inc plot of treatments, and ICEs
# 3. Save all output
####


# Import libraries and user functions -------------------------------------
library('arrow')
library('readr')
library('here')
library('lubridate')
library('dplyr')
library('tidyr')
source(here::here("analysis", "functions", "utility.R"))


# Create directories for output -------------------------------------------
fs::dir_create(here::here("output", "data_description"))


# Import dates ------------------------------------------------------------
source(here::here("analysis", "metadates.R"))
study_dates <- lapply(study_dates, function(x) as.Date(x))


# Define redaction threshold ----------------------------------------------
threshold <- 6


# Import the data ---------------------------------------------------------
df <- read_feather(here("output", "data", "data_processed_full.arrow"))


# Investigate treatment patterns -------------------------------------------
df <- df %>% 
  mutate(
    ## (1) Starting any metformin, from T2DM diagnosis until study end date
    # metformin combo
    exp_bin_metfin_anytime = !is.na(exp_date_metfin_first),
    exp_tb_T2DMdiag_metfin_anytime = case_when(exp_bin_metfin_anytime == TRUE ~ as.numeric(difftime(exp_date_metfin_first, elig_date_t2dm, units = "days")),
                                               TRUE ~ NA_real_),
    exp_bin_metfin_anytime_3m = !is.na(exp_tb_T2DMdiag_metfin_anytime) & exp_tb_T2DMdiag_metfin_anytime <= 90,
    # exp_bin_metfin_anytime_6m = !is.na(exp_tb_T2DMdiag_metfin_anytime) & exp_tb_T2DMdiag_metfin_anytime <= 183, # is exp_bin_metfin
    # metformin mono
    exp_bin_metfin_mono_anytime = !is.na(exp_date_metfin_mono_first),
    exp_tb_T2DMdiag_metfin_mono_anytime = case_when(exp_bin_metfin_mono_anytime == TRUE ~ as.numeric(difftime(exp_date_metfin_mono_first, elig_date_t2dm, units = "days")),
                                                    TRUE ~ NA_real_),
    exp_bin_metfin_mono_anytime_3m = !is.na(exp_tb_T2DMdiag_metfin_mono_anytime) & exp_tb_T2DMdiag_metfin_mono_anytime <= 90,
    # exp_bin_metfin_mono_anytime_6m = !is.na(exp_tb_T2DMdiag_metfin_mono_anytime) & exp_tb_T2DMdiag_metfin_mono_anytime <= 183, # is exp_bin_metfin_mono
    
    ## (2) Let's investigate those who did not start any metfin, OVER ENTIRE STUDY PERIOD, i.e. exp_bin_metfin_anytime == FALSE
    # DPP4 mono (or combo with SGLT2)
    # metfin combo codelist entails also combo with dpp4 (e.g. Janumet), but they are already set to exp_bin_metfin_anytime == TRUE; same for all below (where combo is relevant)
    exp_bin_dpp4_mono_anytime = exp_bin_metfin_anytime == FALSE & !is.na(exp_date_dpp4_first),
    exp_date_dpp4_mono_anytime = case_when(exp_bin_dpp4_mono_anytime == TRUE ~ exp_date_dpp4_first, 
                                           TRUE ~ as.Date(NA)),
    # TZD mono
    exp_bin_tzd_mono_anytime = exp_bin_metfin_anytime == FALSE & !is.na(exp_date_tzd_first),
    exp_date_tzd_mono_anytime = case_when(exp_bin_tzd_mono_anytime == TRUE ~ exp_date_tzd_first, 
                                          TRUE ~ as.Date(NA)),
    # SGLT2 mono (or combo with DPP4)
    exp_bin_sglt2_mono_anytime = exp_bin_metfin_anytime == FALSE & !is.na(exp_date_sglt2_first),
    exp_date_sglt2_mono_anytime = case_when(exp_bin_sglt2_mono_anytime == TRUE ~ exp_date_sglt2_first, 
                                            TRUE ~ as.Date(NA)),
    # sulfo mono (glucovance/glibenclamid + metfin not in use anymore)
    exp_bin_sulfo_mono_anytime = exp_bin_metfin_anytime == FALSE & !is.na(exp_date_sulfo_first),
    exp_date_sulfo_mono_anytime = case_when(exp_bin_sulfo_mono_anytime == TRUE ~ exp_date_sulfo_first, 
                                            TRUE ~ as.Date(NA)),
    # glp1 mono
    exp_bin_glp1_mono_anytime = exp_bin_metfin_anytime == FALSE & !is.na(exp_date_glp1_first),
    exp_date_glp1_mono_anytime = case_when(exp_bin_glp1_mono_anytime == TRUE ~ exp_date_glp1_first, 
                                           TRUE ~ as.Date(NA)),
    # megli mono
    exp_bin_megli_mono_anytime = exp_bin_metfin_anytime == FALSE & !is.na(exp_date_megli_first),
    exp_date_megli_mono_anytime = case_when(exp_bin_megli_mono_anytime == TRUE ~ exp_date_megli_first, 
                                            TRUE ~ as.Date(NA)),
    # agi mono
    exp_bin_agi_mono_anytime = exp_bin_metfin_anytime == FALSE & !is.na(exp_date_agi_first),
    exp_date_agi_mono_anytime = case_when(exp_bin_agi_mono_anytime == TRUE ~ exp_date_agi_first, 
                                          TRUE ~ as.Date(NA)),
    # insulin mono
    exp_bin_insulin_mono_anytime = exp_bin_metfin_anytime == FALSE & !is.na(exp_date_insulin_first),
    exp_date_insulin_mono_anytime = case_when(exp_bin_insulin_mono_anytime == TRUE ~ exp_date_insulin_first, 
                                              TRUE ~ as.Date(NA)),
    ## (3) No prescription at all OVER ENTIRE STUDY PERIOD after T2DM diagnosis
    exp_bin_treat_nothing_anytime = case_when(exp_bin_metfin_anytime == FALSE
                                              & exp_bin_dpp4_mono_anytime == FALSE
                                              & exp_bin_tzd_mono_anytime == FALSE
                                              & exp_bin_sglt2_mono_anytime == FALSE 
                                              & exp_bin_sulfo_mono_anytime == FALSE
                                              & exp_bin_glp1_mono_anytime == FALSE 
                                              & exp_bin_megli_mono_anytime == FALSE
                                              & exp_bin_agi_mono_anytime == FALSE
                                              & exp_bin_insulin_mono_anytime == FALSE ~ TRUE,
                                              TRUE ~ FALSE),
    
    ## (4) Let's investigate those who did initiate anything else than metfin UNTIL 6M after T2DM diagnosis
    # the *_mono variables were created in data_process
    exp_bin_oad = case_when(exp_bin_metfin == FALSE
                            & (exp_bin_dpp4_mono == TRUE
                               | exp_bin_tzd_mono == TRUE 
                               | exp_bin_sglt2_mono == TRUE 
                               | exp_bin_sulfo_mono == TRUE
                               | exp_bin_glp1_mono == TRUE 
                               | exp_bin_megli_mono == TRUE
                               | exp_bin_agi_mono == TRUE
                               | exp_bin_insulin_mono == TRUE) ~ TRUE,
                            TRUE ~ FALSE),
    
    ## (5) Let's investigate those who did initiate anything else than metfin UNTIL 6M after T2DM diagnosis OR nothing
    # should be more than exp_bin_treat_nothing
    exp_bin_oad_nothing = case_when(exp_bin_metfin == FALSE
                                    | exp_bin_dpp4_mono == TRUE
                                    | exp_bin_tzd_mono == TRUE
                                    | exp_bin_sglt2_mono == TRUE
                                    | exp_bin_sulfo_mono == TRUE
                                    | exp_bin_glp1_mono == TRUE
                                    | exp_bin_megli_mono == TRUE
                                    | exp_bin_agi_mono == TRUE
                                    | exp_bin_insulin_mono == TRUE ~ TRUE,
                                    TRUE ~ FALSE),
    
    ## (6) Let's investigate those who did not start any metfin MONO, UNTIL 6M LANDMARK, i.e. exp_bin_metfin_mono == FALSE
    # DPP4 (or combo with SGLT2) +/- metformin
    # Now we may have people who initiated DPP4 + metformin in exp_bin_dpp4
    exp_bin_dpp4 = exp_bin_metfin_mono == FALSE & exp_date_dpp4_first <= landmark_date,
    # TZD +/- metformin
    exp_bin_tzd = exp_bin_metfin_mono == FALSE & exp_date_tzd_first <= landmark_date,
    # SGLT2 (or combo with DPP4) +/- metformin
    exp_bin_sglt2 = exp_bin_metfin_mono == FALSE & exp_date_sglt2_first <= landmark_date,
    # sulfo +/- metformin
    exp_bin_sulfo = exp_bin_metfin_mono == FALSE & exp_date_sulfo_first <= landmark_date,
    # glp1 +/- metformin
    exp_bin_glp1 = exp_bin_metfin_mono == FALSE & exp_date_glp1_first <= landmark_date,
    # megli +/- metformin
    exp_bin_megli = exp_bin_metfin_mono == FALSE & exp_date_megli_first <= landmark_date,
    # agi +/- metformin
    exp_bin_agi = exp_bin_metfin_mono == FALSE & exp_date_agi_first <= landmark_date,
    # insulin +/- metformin
    exp_bin_insulin = exp_bin_metfin_mono == FALSE & exp_date_insulin_first <= landmark_date,
    
    ## (7) No prescription at all UNTIL 6M after T2DM diagnosis | Should give the same result as exp_bin_treat_nothing, just a different distribution across the arms
    exp_bin_treat_nothing2 = case_when(exp_bin_metfin_mono == FALSE
                                       & exp_bin_dpp4 == FALSE
                                       & exp_bin_tzd == FALSE 
                                       & exp_bin_sglt2 == FALSE 
                                       & exp_bin_sulfo == FALSE
                                       & exp_bin_glp1 == FALSE 
                                       & exp_bin_megli == FALSE
                                       & exp_bin_agi == FALSE
                                       & exp_bin_insulin == FALSE ~ TRUE,
                                       TRUE ~ FALSE),
    
    ## (8) OAD prescription (except metformin mono) UNTIL 6M after T2DM diagnosis (i.e. might have some metfin combo in control arm)
    # should be more than exp_bin_oad
    exp_bin_oad_metfincombo = case_when(exp_bin_metfin_mono == FALSE
                                        & (exp_bin_dpp4 == TRUE
                                           | exp_bin_tzd == TRUE 
                                           | exp_bin_sglt2 == TRUE 
                                           | exp_bin_sulfo == TRUE
                                           | exp_bin_glp1 == TRUE 
                                           | exp_bin_megli == TRUE
                                           | exp_bin_agi == TRUE
                                           | exp_bin_insulin == TRUE) ~ TRUE,
                                        TRUE ~ FALSE),
    
    ## (9) OAD prescription (except metformin mono) UNTIL 6M after T2DM diagnosis (i.e. might have some metfin combo in control arm) OR nothing
    # should be more than exp_bin_oad_nothing
    exp_bin_oad_metfincombo_nothing = case_when(exp_bin_metfin_mono == FALSE
                                                | exp_bin_dpp4 == TRUE
                                                | exp_bin_tzd == TRUE 
                                                | exp_bin_sglt2 == TRUE 
                                                | exp_bin_sulfo == TRUE
                                                | exp_bin_glp1 == TRUE 
                                                | exp_bin_megli == TRUE
                                                | exp_bin_agi == TRUE
                                                | exp_bin_insulin == TRUE ~ TRUE,
                                                TRUE ~ FALSE),
    
    ## (10) Additional treatment strategies
    # exp_bin_metfin_mono and exp_bin_treat_nothing defined in data_process.R
    exp_cat_treat_all = case_when(exp_bin_metfin_mono == TRUE ~ 1, 
                                  exp_bin_treat_nothing == TRUE ~ 2,
                                  exp_bin_oad == TRUE ~ 3,
                                  exp_bin_oad_metfincombo == TRUE ~ 4,
                                  TRUE ~ NA_real_),
    exp_cat_treat_3groups = case_when(exp_bin_metfin_mono == TRUE ~ 1, 
                                      exp_bin_treat_nothing == TRUE ~ 2,
                                      exp_bin_oad == TRUE | exp_bin_metfin == TRUE ~ 3,
                                      TRUE ~ NA_real_)
  )


# Summarise as aggregate table --------------------------------------------
n_exp_out <- df %>% 
  summarise(
    # these are from data_process.R
    n_exp_bin_metfin = sum(exp_bin_metfin), 
    n_exp_bin_metfin_mono = sum(exp_bin_metfin_mono), 
    # (1)
    n_exp_bin_metfin_anytime = sum(exp_bin_metfin_anytime),
    n_exp_bin_metfin_anytime_3m = sum(exp_bin_metfin_anytime_3m), 
    n_exp_bin_metfin_mono_anytime = sum(exp_bin_metfin_mono_anytime),
    n_exp_bin_metfin_mono_anytime_3m = sum(exp_bin_metfin_mono_anytime_3m),
    # (2)
    n_exp_bin_dpp4_mono_anytime = sum(exp_bin_dpp4_mono_anytime),
    n_exp_bin_tzd_mono_anytime = sum(exp_bin_tzd_mono_anytime),
    n_exp_bin_sglt2_mono_anytime = sum(exp_bin_sglt2_mono_anytime),
    n_exp_bin_sulfo_mono_anytime = sum(exp_bin_sulfo_mono_anytime),
    n_exp_bin_glp1_mono_anytime = sum(exp_bin_glp1_mono_anytime),
    n_exp_bin_megli_mono_anytime = sum(exp_bin_megli_mono_anytime),
    n_exp_bin_agi_mono_anytime = sum(exp_bin_agi_mono_anytime),
    n_exp_bin_insulin_mono_anytime = sum(exp_bin_insulin_mono_anytime),
    # (3)
    n_exp_bin_treat_nothing_anytime = sum(exp_bin_treat_nothing_anytime),
    # these are from data_process.R
    n_exp_bin_dpp4_mono = sum(exp_bin_dpp4_mono),
    n_exp_bin_tzd_mono = sum(exp_bin_tzd_mono),
    n_exp_bin_sglt2_mono = sum(exp_bin_sglt2_mono),
    n_exp_bin_sulfo_mono = sum(exp_bin_sulfo_mono),
    n_exp_bin_glp1_mono = sum(exp_bin_glp1_mono),
    n_exp_bin_megli_mono = sum(exp_bin_megli_mono),
    n_exp_bin_agi_mono = sum(exp_bin_agi_mono),
    n_exp_bin_insulin_mono = sum(exp_bin_insulin_mono),
    n_exp_bin_treat_nothing = sum(exp_bin_treat_nothing),
    # (4)
    n_exp_bin_oad = sum(exp_bin_oad),
    # (5)
    n_exp_bin_oad_nothing = sum(exp_bin_oad_nothing),
    # (6)
    n_exp_bin_dpp4 = sum(exp_bin_dpp4),
    n_exp_bin_tzd = sum(exp_bin_tzd),
    n_exp_bin_sglt2 = sum(exp_bin_sglt2),
    n_exp_bin_sulfo = sum(exp_bin_sulfo),
    n_exp_bin_glp1 = sum(exp_bin_glp1),
    n_exp_bin_megli = sum(exp_bin_megli),
    n_exp_bin_agi = sum(exp_bin_agi),
    n_exp_bin_insulin = sum(exp_bin_insulin),
    # (7)
    n_exp_bin_treat_nothing2 = sum(exp_bin_treat_nothing2),
    # (8)
    n_exp_bin_oad_metfincombo = sum(exp_bin_oad_metfincombo),
    # (9)
    n_exp_bin_oad_metfincombo_nothing = sum(exp_bin_oad_metfincombo_nothing),
    
    median_tb_T2DMdiag_metfin_anytime = median(exp_tb_T2DMdiag_metfin_anytime, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_anytime = quantile(exp_tb_T2DMdiag_metfin_anytime, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_anytime = quantile(exp_tb_T2DMdiag_metfin_anytime, 0.75, na.rm = TRUE),
    
    median_tb_T2DMdiag_metfin_mono_anytime = median(exp_tb_T2DMdiag_metfin_mono_anytime, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_mono_anytime = quantile(exp_tb_T2DMdiag_metfin_mono_anytime, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_mono_anytime = quantile(exp_tb_T2DMdiag_metfin_mono_anytime, 0.75, na.rm = TRUE)
    
  ) %>% 
  
  # pivot (for easier data review in L4)
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Value"
  )


# midpoint6 rounded
n_exp_out_midpoint6 <- df %>% 
  summarise(
    # these are from data_process.R
    n_exp_bin_metfin_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono, na.rm = TRUE), threshold),
    # (1)
    n_exp_bin_metfin_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_anytime_3m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_anytime_3m, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_metfin_mono_anytime_3m_midpoint6 = fn_roundmid_any(sum(exp_bin_metfin_mono_anytime_3m, na.rm = TRUE), threshold), 
    # (2)
    n_exp_bin_dpp4_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_sulfo_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_glp1_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_megli_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_megli_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_agi_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_agi_mono_anytime, na.rm = TRUE), threshold), 
    n_exp_bin_insulin_mono_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin_mono_anytime, na.rm = TRUE), threshold), 
    # (3)
    n_exp_bin_treat_nothing_anytime_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_nothing_anytime, na.rm = TRUE), threshold), 
    # these are from data_process.R
    n_exp_bin_dpp4_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4_mono, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd_mono, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2_mono, na.rm = TRUE), threshold), 
    n_exp_bin_sulfo_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo_mono, na.rm = TRUE), threshold), 
    n_exp_bin_glp1_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1_mono, na.rm = TRUE), threshold), 
    n_exp_bin_megli_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_megli_mono, na.rm = TRUE), threshold), 
    n_exp_bin_agi_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_agi_mono, na.rm = TRUE), threshold), 
    n_exp_bin_insulin_mono_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin_mono, na.rm = TRUE), threshold), 
    n_exp_bin_treat_nothing_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_nothing, na.rm = TRUE), threshold), 
    # (4)
    n_exp_bin_oad_midpoint6 = fn_roundmid_any(sum(exp_bin_oad, na.rm = TRUE), threshold),     
    # (5)
    n_exp_bin_oad_nothing_midpoint6 = fn_roundmid_any(sum(exp_bin_oad_nothing, na.rm = TRUE), threshold), 
    # (6)
    n_exp_bin_dpp4_midpoint6 = fn_roundmid_any(sum(exp_bin_dpp4, na.rm = TRUE), threshold), 
    n_exp_bin_tzd_midpoint6 = fn_roundmid_any(sum(exp_bin_tzd, na.rm = TRUE), threshold), 
    n_exp_bin_sglt2_midpoint6 = fn_roundmid_any(sum(exp_bin_sglt2, na.rm = TRUE), threshold), 
    n_exp_bin_sulfo_midpoint6 = fn_roundmid_any(sum(exp_bin_sulfo, na.rm = TRUE), threshold), 
    n_exp_bin_glp1_midpoint6 = fn_roundmid_any(sum(exp_bin_glp1, na.rm = TRUE), threshold), 
    n_exp_bin_megli_midpoint6 = fn_roundmid_any(sum(exp_bin_megli, na.rm = TRUE), threshold), 
    n_exp_bin_agi_midpoint6 = fn_roundmid_any(sum(exp_bin_agi, na.rm = TRUE), threshold), 
    n_exp_bin_insulin_midpoint6 = fn_roundmid_any(sum(exp_bin_insulin, na.rm = TRUE), threshold), 
    # (7)
    n_exp_bin_treat_nothing2_midpoint6 = fn_roundmid_any(sum(exp_bin_treat_nothing2, na.rm = TRUE), threshold), 
    # (8)
    n_exp_bin_oad_metfincombo_midpoint6 = fn_roundmid_any(sum(exp_bin_oad_metfincombo, na.rm = TRUE), threshold),
    # (9)
    n_exp_bin_oad_metfincombo_nothing_midpoint6 = fn_roundmid_any(sum(exp_bin_oad_metfincombo_nothing, na.rm = TRUE), threshold),
    
    median_tb_T2DMdiag_metfin_anytime = median(exp_tb_T2DMdiag_metfin_anytime, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_anytime = quantile(exp_tb_T2DMdiag_metfin_anytime, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_anytime = quantile(exp_tb_T2DMdiag_metfin_anytime, 0.75, na.rm = TRUE),
    
    median_tb_T2DMdiag_metfin_mono_anytime = median(exp_tb_T2DMdiag_metfin_mono_anytime, na.rm = TRUE),
    IQR_lower_tb_T2DMdiag_metfin_mono_anytime = quantile(exp_tb_T2DMdiag_metfin_mono_anytime, 0.25, na.rm = TRUE),
    IQR_upper_tb_T2DMdiag_metfin_mono_anytime = quantile(exp_tb_T2DMdiag_metfin_mono_anytime, 0.75, na.rm = TRUE)
    
  ) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# Define labels
labels <- c(
  # these are from data_process.R
  n_exp_bin_metfin_midpoint6 = "Metformin (combo) within 6m",
  n_exp_bin_metfin_mono_midpoint6 = "Metformin mono within 6m",
  # (1)
  n_exp_bin_metfin_anytime_midpoint6 = "Metformin (combo) anytime",
  n_exp_bin_metfin_anytime_3m_midpoint6 = "Metformin (combo) within 3m",
  n_exp_bin_metfin_mono_anytime_midpoint6 = "Metformin mono anytime",
  n_exp_bin_metfin_mono_anytime_3m_midpoint6 = "Metformin mono within 3m",
  # (2)
  n_exp_bin_dpp4_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, DPP4",
  n_exp_bin_tzd_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, TZD",
  n_exp_bin_sglt2_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, SGLT2",
  n_exp_bin_sulfo_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, sulfonylurea",
  n_exp_bin_glp1_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, GLP1",
  n_exp_bin_megli_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, meglitinide",
  n_exp_bin_agi_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, alpha-glucosidase",
  n_exp_bin_insulin_mono_anytime_midpoint6 = "Among those without metformin (combo) anytime, insulin",
  # (3)
  n_exp_bin_treat_nothing_anytime_midpoint6 = "No metformin (combo) or any other antidiabetic anytime",
  # these are from data_process.R
  n_exp_bin_dpp4_mono_midpoint6 = "Among those without metformin (combo) 6m, DPP4",
  n_exp_bin_tzd_mono_midpoint6 = "Among those without metformin (combo) 6m, TZD",
  n_exp_bin_sglt2_mono_midpoint6 = "Among those without metformin (combo) 6m, SGLT2",
  n_exp_bin_sulfo_mono_midpoint6 = "Among those without metformin (combo) 6m, sulfonylurea",
  n_exp_bin_glp1_mono_midpoint6 = "Among those without metformin (combo) 6m, GLP1",
  n_exp_bin_megli_mono_midpoint6 = "Among those without metformin (combo) 6m, meglitinide",
  n_exp_bin_agi_mono_midpoint6 = "Among those without metformin (combo) 6m, alpha-glucosidase",
  n_exp_bin_insulin_mono_midpoint6 = "Among those without metformin (combo) 6m, insulin",
  n_exp_bin_treat_nothing_midpoint6 = "No metformin (combo) or any other antidiabetic within 6m",
  # (4)
  n_exp_bin_oad_midpoint6 = "Among those without metformin (combo) 6m, any antidiabetic (mono)",
  # (5)
  n_exp_bin_oad_nothing_midpoint6 = "Among those without metformin (combo) 6m, any antidiabetic or nothing",
  # (6)
  n_exp_bin_dpp4_midpoint6 = "Among those without metformin mono 6m, DPP4 (+/- metformin)",
  n_exp_bin_tzd_midpoint6 = "Among those without metformin mono 6m, TZD (+/- metformin)",
  n_exp_bin_sglt2_midpoint6 = "Among those without metformin mono 6m, SGLT2 (+/- metformin)",
  n_exp_bin_sulfo_midpoint6 = "Among those without metformin mono 6m, sulfonylurea (+/- metformin)",
  n_exp_bin_glp1_midpoint6 = "Among those without metformin mono 6m, GLP1 (+/- metformin)",
  n_exp_bin_megli_midpoint6 = "Among those without metformin mono 6m, meglitinide (+/- metformin)",
  n_exp_bin_agi_midpoint6 = "Among those without metformin mono 6m, alpha-glucosidase (+/- metformin)",
  n_exp_bin_insulin_midpoint6 = "Among those without metformin mono 6m, insulin (+/- metformin)",
  # (7)
  n_exp_bin_treat_nothing2_midpoint6 = "No metformin mono or any other antidiabetic (+/- metformin) within 6m",
  # (8)
  n_exp_bin_oad_metfincombo_midpoint6 = "Among those without metformin mono 6m, any antidiabetic (+/- metformin)",
  # (9)
  n_exp_bin_oad_metfincombo_nothing_midpoint6 = "Among those without metformin mono 6m, any antidiabetic (+/- metformin) or nothing",
  
  median_tb_T2DMdiag_metfin_anytime = "Median time from T2DM diagnosis to metformin (combo) start",
  IQR_lower_tb_T2DMdiag_metfin_anytime = "IQR lower bound: T2DM diagnosis to metformin (combo)",
  IQR_upper_tb_T2DMdiag_metfin_anytime = "IQR upper bound: T2DM diagnosis to metformin (combo)",
  
  median_tb_T2DMdiag_metfin_mono_anytime = "Median time from T2DM diagnosis to metformin mono start",
  IQR_lower_tb_T2DMdiag_metfin_mono_anytime = "IQR lower bound: T2DM diagnosis to metformin mono",
  IQR_upper_tb_T2DMdiag_metfin_mono_anytime = "IQR upper bound: T2DM diagnosis to metformin mono"
)

# Apply them
n_exp_out_midpoint6 <- n_exp_out_midpoint6 %>%
  mutate(Variable = labels[Variable])


# Output for cum inc treatment plot ---------------------------------------
data_plots <- df %>%
  dplyr::select(patient_id, elig_date_t2dm, exp_date_metfin_first, exp_bin_metfin_anytime, exp_date_metfin_mono_first, exp_bin_metfin_mono_anytime,
                exp_date_dpp4_mono_anytime, exp_bin_dpp4_mono_anytime, exp_date_tzd_mono_anytime, exp_bin_tzd_mono_anytime,
                exp_date_sglt2_mono_anytime, exp_bin_sglt2_mono_anytime, exp_date_sulfo_mono_anytime, exp_bin_sulfo_mono_anytime,
                exp_date_glp1_mono_anytime, exp_bin_glp1_mono_anytime, exp_date_megli_mono_anytime, exp_bin_megli_mono_anytime,
                exp_date_agi_mono_anytime, exp_bin_agi_mono_anytime, exp_date_insulin_mono_anytime, exp_bin_insulin_mono_anytime,
                out_date_severecovid, qa_date_of_death)
data_plots <- data_plots %>% # double-check for plot to avoid event_time is == 0
  dplyr::filter(elig_date_t2dm < exp_date_metfin_first | is.na(exp_date_metfin_first)) %>%
  dplyr::filter(elig_date_t2dm < exp_date_metfin_mono_first | is.na(exp_date_metfin_mono_first)) %>%
  dplyr::filter(elig_date_t2dm < exp_date_dpp4_mono_anytime | is.na(exp_date_dpp4_mono_anytime)) %>%
  dplyr::filter(elig_date_t2dm < exp_date_tzd_mono_anytime | is.na(exp_date_tzd_mono_anytime)) %>%
  dplyr::filter(elig_date_t2dm < exp_date_sglt2_mono_anytime | is.na(exp_date_sglt2_mono_anytime)) %>%
  dplyr::filter(elig_date_t2dm < exp_date_sulfo_mono_anytime | is.na(exp_date_sulfo_mono_anytime)) %>%
  dplyr::filter(elig_date_t2dm < exp_date_glp1_mono_anytime | is.na(exp_date_glp1_mono_anytime)) %>%
  dplyr::filter(elig_date_t2dm < exp_date_megli_mono_anytime | is.na(exp_date_megli_mono_anytime)) %>%
  dplyr::filter(elig_date_t2dm < exp_date_agi_mono_anytime | is.na(exp_date_agi_mono_anytime)) %>%
  dplyr::filter(elig_date_t2dm < exp_date_insulin_mono_anytime | is.na(exp_date_insulin_mono_anytime)) %>%
  dplyr::filter(elig_date_t2dm < qa_date_of_death | is.na(qa_date_of_death))


# Save output -------------------------------------------------------------
# data for cumulative incidence plots re treatment regimen pattern
write_feather(data_plots, here::here("output", "data", "data_plots.feather"))
# descriptive data re treatment patterns, events between index_date and pandemic start, and main outcome
write.csv(n_exp_out_midpoint6, file = here::here("output", "data_description", "n_exp_out_midpoint6.csv"))
write.csv(n_exp_out, file = here::here("output", "data_description", "n_exp_out.csv"))
