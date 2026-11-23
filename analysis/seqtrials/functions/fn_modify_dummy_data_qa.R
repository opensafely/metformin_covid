################################################################################
# Modify the DUMMY data of the QA/completeness dataset (only run locally, never on
# real data): the ehrQL dummy data are random (e.g. ~60% with a prostate cancer code,
# ~15% with known deprivation), so almost everybody would be excluded.
# Make the variables plausible, with a few people failing each criterion, so that
# every step of the flowchart is tested.
################################################################################
fn_modify_dummy_data_qa <- function(data, studyend_date) {
  set.seed(36)
  n <- nrow(data)
  data %>%
    mutate(
      # Year of birth: adults, plus a few missing or implausible (rules 1 and 3)
      qa_num_birth_year = sample(1935:2003, n, replace = TRUE),
      qa_num_birth_year = case_when(
        runif(n) < 0.02 ~ NA_real_,
        runif(n) < 0.02 ~ 1830,
        TRUE ~ qa_num_birth_year
      ),
      # Date of death: mostly alive; deaths after the T2DM diagnosis, plus a few
      # implausible ones (rules 2 and 4)
      qa_date_of_death = if_else(
        runif(n) < 0.1,
        elig_date_t2dm + sample(1:1000, n, replace = TRUE),
        as.Date(NA)
      ),
      # (implausible = in the future; a date before 1900 would already be caught by rule 2)
      qa_date_of_death = case_when(
        runif(n) < 0.02 ~ Sys.Date() + 100,
        TRUE ~ qa_date_of_death
      ),
      # Sex: mostly female/male, a few unknown (completeness 1)
      cov_cat_sex = factor(
        sample(c("female", "male", "unknown"), n, replace = TRUE, prob = c(0.48, 0.5, 0.02)),
        levels = levels(cov_cat_sex)
      ),
      # QA codes: rare (rules 5-7)
      qa_bin_pregnancy = runif(n) < 0.03,
      qa_bin_hrt = runif(n) < 0.03,
      qa_bin_cocp = runif(n) < 0.03,
      qa_bin_prostate_cancer = runif(n) < 0.03,
      # Deprivation and STP: mostly known (completeness 2 and 3)
      qa_bin_known_imd = runif(n) < 0.95,
      qa_bin_known_stp = runif(n) < 0.95
    )
}
