###
# Custom-made function to perform disclosure-safe standardization/marginalization
###

# Args:
#   df: The data frame containing the patient data (e.g., df_months_severecovid).
#   K: The number of time points (e.g., total months of follow-up).
#   model: The logistic regression model for predicting discrete-time hazards.
#   group_col: A string specifying the name of the group name
#   time_col: A string specifying the name of the time column | CAVE: We model time as time and time_sqr. Needs to be double-checked if in line with analysis approach.
#   patient_id_col: The name of the patient ID column to group by. -> esp. important for use in bootstrap!
#   min_count: disclosure threshold, current OpenSAFELY default is 6
#   apply_rounding: Logical. If TRUE (default), applies fn_roundmid_any disclosure control to cumulative numerators and denominator before computing cumulative incidence.
#                   Set to FALSE (e.g. inside bootstrap to avoid adding artificial rounding noise to each replicate, which would inflate CI variance). In that case, rounding should only be applied to the point estimate from the original sample.

# Returns:
#   A list containing the final standardized risk data frame and the final
#   risk difference (rd) and risk ratio (rr) estimates.
# After binding:
# Row 1: month = 0 (zero row, time_0 = 1; from zero_row, see fn_standardize_risks)
# Row 2: month = 0 (time_0 = 1; from risk_summary, see fn_standardize_risks)
# Row 3: month = 1 (time_0 = 2)

# get fn_roundmid_any
source(here::here("analysis", "functions", "utility.R"))

fn_standardize_risks <- function(df, K, model, group_col = "exp_bin_treat", time_col = "month", patient_id_col = "patient_id", min_count = 6, apply_rounding = TRUE) {
  # Create a dataset with all time points for each individual under each treatment level
  # This sets up a "prediction" dataset where we can estimate outcomes for every patient at every time point under each treatment scenario.
  df_pred <- df %>%
    dplyr::filter(!!rlang::sym(time_col) == 0) %>%
    dplyr::select(-!!rlang::sym(time_col)) %>%
    tidyr::crossing(!!rlang::sym(time_col) := 0:(K - 1))
  
  # Add a squared month variable, as it's used in the model
  df_pred$monthsqr <- df_pred[[time_col]]^2
  
  # Predict hazards for treatment = 0 and treatment = 1
  # Instead of creating separate datasets, temporarily overwrite the group_col, trick predict() into using modified values without duplicating the whole data frame.
  df_pred <- df_pred %>%
    dplyr::mutate(
      # predicted hazard under control (0)
      p.event0 = predict(
        model,
        {
          tmp <- .
          tmp[[group_col]] <- 0
          tmp
        },
        type = "response"
      ),
      # predicted hazard under treatment (1)
      p.event1 = predict(
        model,
        {
          tmp <- .
          tmp[[group_col]] <- 1
          tmp
        },
        type = "response"
      )
    )
  
  # Compute survival curves (first survival probabilities, then the inverse -> risks)
  # The survival probability at time t is the cumulative product of (1 - hazard) up to that time point.
  df_pred <- df_pred %>%
    dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
    dplyr::mutate(
      surv0 = cumprod(1 - p.event0),
      surv1 = cumprod(1 - p.event1),
      risk0 = 1 - surv0,
      risk1 = 1 - surv1
    ) %>%
    dplyr::ungroup()
  
  # DISCLOSURE CONTROL
  # Get the mean risk in each treatment group at each time point
  # This aggregates the individual risks to get the standardized population-level risk.
  
  # N is fixed (same patients at every time point in g-computation)
  N <- dplyr::n_distinct(df_pred[[patient_id_col]])
  
  risk_summary <- df_pred %>%
    dplyr::group_by(!!rlang::sym(time_col)) %>%
    dplyr::summarise(
      # Sum of individual predicted risks = "expected event" numerator
      cml.event0 = sum(risk0),
      cml.event1 = sum(risk1),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Apply midpoint rounding when producing the point estimate for release.
      # When apply_rounding = FALSE (i.e. inside bootstrap), pass through raw values, to avoid adding artificial rounding noise to each replicate.
      cml.event0.rounded = if (apply_rounding) fn_roundmid_any(cml.event0, to = min_count) else cml.event0,
      cml.event1.rounded = if (apply_rounding) fn_roundmid_any(cml.event1, to = min_count) else cml.event1,
      # Denominator: total N, rounded once (only meaningful for release output)
      denom = if (apply_rounding) fn_roundmid_any(N, to = min_count) else N,
      # Back-calculate cumulative incidence from (rounded) counts
      risk0 = cml.event0.rounded / denom,
      risk1 = cml.event1.rounded / denom,
      # Edit data frame to reflect that risks are estimated at the END of each interval
      time_0 = .data[[time_col]] + 1
    )
  
  # add the initial zero row
  zero_row <- tibble::tibble(
    !!time_col := 0,
    cml.event0.rounded = 0,
    cml.event1.rounded = 0,
    denom = if (apply_rounding) fn_roundmid_any(N, to = min_count) else N,
    risk0 = 0,
    risk1 = 0,
    time_0 = 1
  )
  
  # build the graph and calculate risk difference and risk ratio
  # graph_data serves as both the release output (disclosure-safe when apply_rounding = TRUE) and the source for plotting. It contains rounded numerators and denominator so that OpenSAFELY output checkers can verify no percentage masks a small underlying count.
  graph <- dplyr::bind_rows(zero_row, risk_summary) %>%
    dplyr::mutate(
      rd = risk1 - risk0,
      rr = risk1 / risk0
    )
  
  # extract final values at end of follow-up
  final_row <- graph %>% dplyr::filter(.data[[time_col]] == K - 1)
  
  return(list(
    graph_data = graph, # disclosure-safe output for graph building and release
    final_risk0 = final_row$risk0,
    final_risk1 = final_row$risk1,
    final_rd = final_row$rd,
    final_rr = final_row$rr
  ))

}