#### 
## This script modifies dummy data of df_months_months to: 
## - randomly select 80% of baseline lab values
## - add some random values (all of them are currently empty in the dummy data)
## - this will test the fill-up below
####

set.seed(123)  # for reproducibility

df_months <- df_months %>%
  mutate(
    cov_num_hba1c_mmol_mol = if_else(
      month == 0,
      # generate values only for baseline, not the rest of the monthly person-intervals
      round(pmax(pmin(rnorm(n(), mean = 50, sd = 15), 120), 0)),
      cov_num_hba1c_mmol_mol
    ),
    cov_num_cholesterol_b = if_else(
      month == 0,
      round(pmax(pmin(rnorm(n(), mean = 5, sd = 1.2), 20), 1.75)),
      cov_num_cholesterol_b
    ),
    cov_num_hdl_cholesterol_b = if_else(
      month == 0,
      round(pmax(pmin(rnorm(n(), mean = 1.3, sd = 0.3), 5), 0.4)),
      cov_num_hdl_cholesterol_b
    ),
    cov_num_bmi_b = if_else(
      month == 0,
      round(pmax(pmin(rnorm(n(), mean = 27, sd = 4), 50), 16)),
      cov_num_bmi_b
    )
  ) %>%
  # randomly set ~20% of baseline values to NA
  mutate(
    cov_num_hba1c_mmol_mol = if_else(month == 0 & runif(n()) < 0.2, NA_real_, cov_num_hba1c_mmol_mol),
    cov_num_cholesterol_b   = if_else(month == 0 & runif(n()) < 0.2, NA_real_, cov_num_cholesterol_b),
    cov_num_hdl_cholesterol_b = if_else(month == 0 & runif(n()) < 0.2, NA_real_, cov_num_hdl_cholesterol_b),
    cov_num_bmi_b           = if_else(month == 0 & runif(n()) < 0.2, NA_real_, cov_num_bmi_b)
  )