####
# Custom made function for dummy data modification:
# Ensure all exposure and outcome dates are between baseline date and studyend date 
# If the date is NA, leave it as NA
# If the date is before the baseline date, replace it with a random valid date
# If the date is valid, keep it unchanged
####
fn_dd_exp_out_dates <- function(date_var, baseline_date, studyend_date) {
  ifelse(is.na(date_var), NA,
         ifelse(
           date_var <= baseline_date, 
           if (as.Date(baseline_date) + 1 <= as.Date(studyend_date)) {
             sample(seq(as.Date(baseline_date) + 1, as.Date(studyend_date), by = "day"), 1)
           } else {
             NA # Avoid errors if range is invalid
           },
           date_var # Keep valid existing dates
         ))
}