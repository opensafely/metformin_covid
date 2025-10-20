###
# Custom-made function to perform standardization/marginalization
###

# Args:
#   df: The data frame containing the patient data (e.g., df_months_severecovid).
#   K: The number of time points (e.g., total months of follow-up).
#   model: The logistic regression model for predicting discrete-time hazards.
#   group_col: A string specifying the name of the group name
#   time_col: A string specifying the name of the time column
#   patient_id_col: The name of the patient ID column to group by. -> esp. important for use in bootstrap!
# Returns:
#   A list containing the final standardized risk data frame and the final
#   risk difference (rd) and risk ratio (rr) estimates.

fn_standardize_risks <- function(df, K, model, group_col = "exp_bin_treat", time_col = "month", patient_id_col = "patient_id") {
  # Create a dataset with all time points for each individual under each treatment level
  # This sets up a "prediction" dataset where we can estimate outcomes for every patient at every time point under each treatment scenario.
  df_pred <- df %>%
    dplyr::filter(!!rlang::sym(time_col) == 0) %>%
    dplyr::select(-!!rlang::sym(time_col)) %>%
    tidyr::crossing(!!rlang::sym(time_col) := 0:(K - 1))
  
  # Add a squared month variable, as it's likely used in the model
  df_pred$monthsqr <- df_pred[[time_col]]^2
  
  # --- Predict outcomes for the control group (everyone untreated) ---
  df_pred0 <- df_pred
  df_pred0[[group_col]] <- 0
  df_pred0$p.event0 <- predict(model, df_pred0, type = "response")
  
  # --- Predict outcomes for the treatment group (everyone treated) ---
  df_pred1 <- df_pred
  df_pred1[[group_col]] <- 1
  df_pred1$p.event1 <- predict(model, df_pred1, type = "response")
  
  # --- Obtain predicted survival probabilities from discrete-time hazards ---
  # The survival probability at time t is the cumulative product of (1 - hazard) up to that time point.
  df_pred0 <- df_pred0 %>%
    dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
    dplyr::mutate(surv0 = cumprod(1 - p.event0)) %>%
    dplyr::ungroup()
  
  df_pred1 <- df_pred1 %>%
    dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
    dplyr::mutate(surv1 = cumprod(1 - p.event1)) %>%
    dplyr::ungroup()
  
  # --- Estimate risks from survival probabilities ---
  # Risk is simply 1 - survival probability.
  df_pred0$risk0 <- 1 - df_pred0$surv0
  df_pred1$risk1 <- 1 - df_pred1$surv1
  
  # --- Get the mean risk in each treatment group at each time point ---
  # This aggregates the individual risks to get the standardized population-level risk.
  risk0 <- aggregate(df_pred0[c(group_col, time_col, "risk0")],
                     by = list(df_pred0[[time_col]]), FUN = mean)[c(group_col, time_col, "risk0")]
  
  risk1 <- aggregate(df_pred1[c(group_col, time_col, "risk1")],
                     by = list(df_pred1[[time_col]]), FUN = mean)[c(group_col, time_col, "risk1")]
  
  # --- Combine the risk estimates into a single data frame ---
  graph.pred <- merge(risk0, risk1, by = c(time_col))
  
  # Edit data frame to reflect that risks are estimated at the END of each interval
  graph.pred$time_0 <- graph.pred[[time_col]] + 1
  zero <- data.frame(val1 = 0, val2 = 0, val3 = 0, val4 = 1, val5 = 0, val6 = 0)
  zero <- setNames(zero, c(time_col, paste0(group_col, ".x"), "risk0", paste0(group_col, ".y"), "risk1", "time_0"))
  graph <- rbind(zero, graph.pred)
  
  # Add Risk Difference (RD) and Risk Ratio (RR) to the graph data frame
  graph$rd <- graph$risk1 - graph$risk0
  graph$rr <- graph$risk1 / graph$risk0
  
  # Extract final risk estimates at the end of follow-up (K-1)
  final_risk0 <- graph$risk0[which(graph[[time_col]] == K - 1)]
  final_risk1 <- graph$risk1[which(graph[[time_col]] == K - 1)]
  final_rd <- graph$rd[which(graph[[time_col]] == K - 1)]
  final_rr <- graph$rr[which(graph[[time_col]] == K - 1)]
  
  # Return a list of all relevant outputs
  return(list(
    graph_data = graph,
    final_risk0 = final_risk0,
    final_risk1 = final_risk1,
    final_rd = final_rd,
    final_rr = final_rr
  ))
}