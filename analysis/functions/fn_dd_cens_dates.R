####
# Custom made function for dummy data modification:
# Ensure random censoring dates after baseline date until studyend date
# Impute 5% of date values
####
fn_dd_cens_dates <- function(baseline_date, studyend_date) {
  if (runif(1) <= 0.05) {  # 5% chance
    if (as.Date(baseline_date) + 1 <= as.Date(studyend_date)) {
      return(sample(seq(as.Date(baseline_date) + 1, as.Date(studyend_date), by = "day"), 1))
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}