#### 
## This script modifies dummy data of df_months to: 
## - randomly select 50% of cov_date_ami 
## - change them to be very close to cens_date_metfin_start_cont (within 5 days before/after cens_date_metfin_start_cont)
## - this will test edge cases for fn_assign_time_fixed_cov
####

# Identify eligible persons (both dates non-NA)
eligible_persons <- df_months %>%
  distinct(patient_id, cov_date_ami, cens_date_metfin_start_cont) %>%
  filter(!is.na(cov_date_ami) & !is.na(cens_date_metfin_start_cont)) %>%
  pull(patient_id)

# Randomly select ~50% of them
n_change <- ceiling(0.50 * length(eligible_persons))
persons_to_change <- sample(eligible_persons, n_change)

# Compute new cov_date_ami for selected persons
new_dates <- df_months %>%
  distinct(patient_id, cov_date_ami, cens_date_metfin_start_cont) %>%
  filter(patient_id %in% persons_to_change) %>%
  mutate(
    new_cov_date_ami = cens_date_metfin_start_cont + sample(c(-5:-1, 1:5), n(), replace=TRUE)
  ) %>%
  select(patient_id, new_cov_date_ami)

# Merge back into original data
df_months <- df_months %>%
  left_join(new_dates, by="patient_id") %>%
  mutate(
    cov_date_ami = ifelse(patient_id %in% persons_to_change, new_cov_date_ami, cov_date_ami)
  ) %>%
  select(-new_cov_date_ami)

df_months$cov_date_ami <- as.Date(df_months$cov_date_ami, origin = "1970-01-01")