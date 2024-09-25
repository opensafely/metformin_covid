######################################

# This script contains one function used in data_process.R:
# - add_status_and_fu_primary: adds colums status_primary and fu_primary to data
# linda.nab@thedatalab.org 20220627
######################################

library("lubridate")
library("tidyverse")

# Function 'add_status_and_fu_primary' add columns status_primary and fu_primary
# Input:
# - data: data.frame with the data extracted using study_definition.py
# Output:
# - data.frame with two columns added (status_primary and fu_primary, see below)
add_status_and_fu_primary <- function(data){
  data %>%
    mutate(
      # STATUS 'PRIMARY' ----
      # "none": followed for 28 days
      # "dereg": lost to follow up (deregistered)
      # "covid_hosp": covid hospitalisation
      # "covid_death": covid death
      # "noncovid_death": non-covid death
      min_date_primary = pmin(out_date_dereg_28,
                              out_date_death_28,
                              out_date_covid_death_28,
                              out_date_noncovid_death_28,
                              study_window,
                              na.rm = TRUE),
      status_primary = case_when(
        # pt should not have both noncovid and covid death, coded here to
        # circumvent mistakes if database errors exist
        min_date_primary == out_date_covid_death_28 ~ "covid_death",
        min_date_primary == out_date_noncovid_death_28 ~ "noncovid_death",
        min_date_primary == out_date_covid_hosp_28 ~ "covid_hosp",
        min_date_primary == out_date_dereg_28 ~ "dereg",
        TRUE ~ "none"
      ),
      # FOLLOW UP STATUS 'PRIMARY' ----
      fu_primary = difftime(min_date_primary,
                            baseline_date,
                            units = "days") %>% as.numeric(),
      # combine covid death and hospitalisation
      status_primary = if_else(status_primary == "covid_hosp" | status_primary == "covid_death", "covid_hosp_death", status_primary) %>%
        factor(levels = c("covid_hosp_death", "noncovid_death", "dereg", "none"))
    )
}
