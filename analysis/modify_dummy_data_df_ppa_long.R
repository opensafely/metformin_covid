#### 
## This script modifies dummy data of df_ppa_long to: 
## - randomly select 5% of "date" (which is the date variable for all measurements in df_ppa_long)
## - change them to be very close to cens_date_metfin_start_cont (within 3 days before/after cens_date_metfin_start_cont)
## - this will test edge cases
####

# Identify eligible persons (both dates non-NA)
eligible_persons <- df_ppa_long %>%
  distinct(patient_id, date, variable, instance, cens_date_metfin_start_cont) %>%
  filter(!is.na(date) & !is.na(cens_date_metfin_start_cont)) %>%
  pull(patient_id)

# Randomly select ~5% of them
n_change <- ceiling(0.05 * length(eligible_persons))
persons_to_change <- sample(eligible_persons, n_change)

# Compute new date for selected persons
new_dates <- df_ppa_long %>%
  distinct(patient_id, date, variable, instance, cens_date_metfin_start_cont) %>%
  filter(patient_id %in% persons_to_change) %>%
  mutate(
    new_date = cens_date_metfin_start_cont + sample(c(-3:-1, 1:3), n(), replace=TRUE)
  ) %>%
  select(patient_id, new_date, variable, instance)
new_dates <- new_dates %>% arrange(patient_id, variable, new_date)

# Merge back into original data
df_ppa_long <- df_ppa_long %>% arrange(patient_id, variable, date)
df_ppa_long <- df_ppa_long %>%
  left_join(new_dates, by=c("patient_id", "variable", "instance")) %>%
  mutate(
    date = ifelse(patient_id %in% persons_to_change, new_date, date)
  ) %>%
  select(-new_date)

df_ppa_long$date <- as.Date(df_ppa_long$date, origin = "1970-01-01")
