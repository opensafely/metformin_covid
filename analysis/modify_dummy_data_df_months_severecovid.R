#### 
# Modify dummy data to tackle the issue of rank-deficiency:
## - most binary variables are NA, but in real data they will either be 0 or 1

set.seed(123)
n <- nrow(df_months_severecovid)

df_months_severecovid <- df_months_severecovid %>%
  mutate(
    cov_bin_ami = rbinom(n, 1, 0.05),
    cov_bin_all_stroke = rbinom(n, 1, 0.10),
    cov_bin_other_arterial_embolism = rbinom(n, 1, 0.03),
    cov_bin_vte = rbinom(n, 1, 0.15),
    cov_bin_hf = rbinom(n, 1, 0.20),
    cov_bin_angina = rbinom(n, 1, 0.10),
    cov_bin_dementia = rbinom(n, 1, 0.09),
    cov_bin_hypertension = rbinom(n, 1, 0.15),
    cov_bin_copd = rbinom(n, 1, 0.11),
    cov_bin_liver_disease = rbinom(n, 1, 0.02),
    cov_bin_chronic_kidney_disease = rbinom(n, 1, 0.15),
    cov_bin_pcos = rbinom(n, 1, 0.02),
    cov_bin_prediabetes = rbinom(n, 1, 0.08),
    cov_bin_diabetescomp = rbinom(n, 1, 0.07)
  )