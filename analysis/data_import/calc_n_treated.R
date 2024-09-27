################################################################################
# A custom made function to calculate the number of eligible participants treated within grace period window (and beyond)
################################################################################
calc_n_treated <- function(data_processed){
  n_before_processing <-
    data_processed %>%
    nrow()
  # treated within window
  n_treated_window  <-
    data_processed %>%
    filter(exp_treatment == "Treated") %>%
    nrow()
  # treated beyond window
  n_treated <-
    data_processed %>%
    filter(exp_any_treatment_cat == "Treated") %>%
    nrow()
  # treated beyond window, when?
  n_treated_time <-
    data_processed %>%
    filter(exp_any_treatment_cat == "Treated") %>%
    summarise(median = median(exp_num_tb_postest_treat),
              IQR = IQR(exp_num_tb_postest_treat),
              Q1 = quantile(exp_num_tb_postest_treat, probs = 0.25),
              Q3 = quantile(exp_num_tb_postest_treat, probs = 0.75))

  out <- tibble(n_before_processing,
                n_treated_window,
                n_treated,
                n_treated_time
                )
}
