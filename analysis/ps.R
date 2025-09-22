####
## This script does the following:
# 1. Import processed data
# 2. Create a propensity score for treatment
# 3. Extract SMDs before and after IPW based on PS
# 4. Density plot; untrimmed and trimmed
# 5. Save all output
####

# Import libraries and functions ------------------------------------------
print('Import libraries and functions')
library(arrow)
library(here)
library(tidyverse)
library(rms) # splines and strat
library(cobalt) # SMD density
library(ggplot2)
library(gtsummary)

# Create directories for output -------------------------------------------
print('Create directories for output')
fs::dir_create(here::here("output", "ps"))

# Import the data ---------------------------------------------------------
print('Import the data')
df <- read_feather(here("output", "data", "data_processed.arrow"))

# Add splines -------------------------------------------------------------
print('Add/compute splines')
# Compute knot locations based on percentiles, according to study protocol
age_knots <- quantile(df$cov_num_age, probs = c(0.10, 0.50, 0.90))
# print(age_knots)

# Define treatment variable, and covariates -------------------------------
print('Define treatment variable, and covariates')
# Set format and reference to ensure proper logistic regression modeling
df$exp_bin_treat <- ordered(df$exp_bin_treat, 
                                       levels = c(0, 1), # Reference first
                                       labels = c("nothing", "metformin"))

# Define the covariates included in the PS model
covariate_names <- names(df) %>%
  grep("^(cov_|strat_)", ., value = TRUE) %>% 
  # exclude those not needed in the model:
  ## cov_cat_bmi_groups is covering cov_bin_obesity & cov_num_bmi_b,
  ## cov_cat_hba1c_b is covering for cov_num_hba1c_b
  ## cov_cat_tc_hdl_ratio_b is covering for cov_num_tc_hdl_ratio_b
  ## spline(cov_num_age) is covering for cov_cat_age
  setdiff(c("cov_num_bmi_b", "cov_bin_obesity", "cov_num_hba1c_b", "cov_cat_age", "cov_num_tc_hdl_ratio_b", "cov_num_hdl_chol_b", "cov_num_chol_b"
            )) 
print(covariate_names)

# PS model ----------------------------------------------------------------
print('PS model and predict')
ps_formula <- as.formula(paste("exp_bin_treat ~ rcs(cov_num_age, age_knots) +", paste(covariate_names, collapse = " + ")))
# Fit the PS model to estimate the PS for being in the metformin-mono group
ps_model <- glm(ps_formula, family = binomial(link = "logit"), data = df)
# summary(ps_model)

# Extract propensity scores
df$ps <- predict(ps_model, type = "response")

# Export full PS model summary --------------------------------------------
print('Export PS model summary')

# Convert model summary coefficients to a data frame
ps_model_summary <- summary(ps_model)$coefficients %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable")

# Inverse probability weighting -------------------------------------------
print('Inverse probability weighting')
# Stabilized
p_treat <- mean(df$exp_bin_treat == "metformin") # prop treated, as per formular for stab weights
df$sw <- ifelse(df$exp_bin_treat == "metformin", 
                                 p_treat / df$ps, 
                                 (1 - p_treat) / (1 - df$ps))

# Assess covariate balance BEFORE IPW -------------------------------------
print('Assess covariate balance BEFORE IPW')
tbl1_unweighted <- bal.tab(df[covariate_names], 
                          treat = df$exp_bin_treat, 
                          binary = "std",
                          s.d.denom = "pooled") # standardization for binary variables
smd_unweighted <- as.data.frame(tbl1_unweighted$Balance)
smd_unweighted <- smd_unweighted %>% 
  select(Diff.Un) %>% 
  rename("SMD_Unweighted" = Diff.Un)

# Assess covariate balance AFTER IPW --------------------------------------
print('Assess covariate balance AFTER IPW')
tbl1_sw <- bal.tab(df[covariate_names], 
                    treat = df$exp_bin_treat, 
                    weights = df$sw,
                    binary = "std", 
                    s.d.denom = "pooled")
smd_sw <- as.data.frame(tbl1_sw$Balance)
smd_sw <- smd_sw %>% 
  select(Diff.Adj) %>% 
  rename("SMD_weighted_stabilized" = Diff.Adj)

# Trimming ----------------------------------------------------------------
print('Min./Max. PS per group and resulting common min/max after trimming, then trimming the dataset, and double-checking weights before/after')
# Determine the common min overlap range in PS
# Reduce to the max. of all group-wise min. PS (i.e. everyone has at least this PS) & the min. of all group-wise max. PS (i.e. no-one has a PS above) 
ps_trim <- df %>%
  group_by(exp_bin_treat) %>%
  summarise(min_ps = min(ps), max_ps = max(ps)) %>%
  ungroup() %>%
  summarise(min_common = max(min_ps), max_common = min(max_ps))
df_trimmed <- df %>%
  filter(ps >= ps_trim$min_common[1] & ps <= ps_trim$max_common[1])

# Document min/max PS before trimming and resulting common min/max after trimming
ps_minmax <- df %>%
  group_by(exp_bin_treat) %>%
  summarise(min_ps = min(ps), max_ps = max(ps))
ps_summary <- ps_minmax %>%
  mutate(
    min_common = ps_trim$min_common,
    max_common = ps_trim$max_common
  )

## Min, 25th percentile, median, mean, SD, 75th percentile, and max
summary(df$sw)
sd(df$sw)
summary(df_trimmed$sw)
sd(df_trimmed$sw)

# Density plot ------------------------------------------------------------
print('Density plot')
## Untrimmed
dens_treated <- density(df$ps[df$exp_bin_treat == "metformin"])
df_treated <- data.frame(dens_x = dens_treated$x, dens_y = dens_treated$y, group = "metformin (untrimmed)")
dens_untreated <- density(df$ps[df$exp_bin_treat == "nothing"])
df_untreated <- data.frame(dens_x = dens_untreated$x, dens_y = dens_untreated$y, group = "nothing (untrimmed)")
ps_density_data_untrimmed <- bind_rows(df_treated, df_untreated)

density_plot_untrimmed <- ggplot(ps_density_data_untrimmed, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1), 
    minor_breaks = seq(0, 1, by = 0.05), 
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity Score Density Plot; Untrimmed",
    x = "Propensity Score",
    y = "Density",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 10)  # optional: make labels clearer
  )

## Trimmed
dens_treated_trim <- density(df_trimmed$ps[df_trimmed$exp_bin_treat == "metformin"])
df_treated_trim <- data.frame(dens_x = dens_treated_trim$x, dens_y = dens_treated_trim$y, group = "metformin (trimmed)")
dens_untreated_trim <- density(df_trimmed$ps[df_trimmed$exp_bin_treat == "nothing"])
df_untreated_trim <- data.frame(dens_x = dens_untreated_trim$x, dens_y = dens_untreated_trim$y, group = "nothing (trimmed)")
ps_density_data_trimmed <- bind_rows(df_treated_trim, df_untreated_trim)

density_plot_trimmed <- ggplot(ps_density_data_trimmed, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1), 
    minor_breaks = seq(0, 1, by = 0.05), 
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity Score Density Plot; Untrimmed",
    x = "Propensity Score",
    y = "Density",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 10)  # optional: make labels clearer
  )

# Histograms  -------------------------------------------------------------
print('Density histograms')
## Untrimmed
df_untrimmed_hist <- df %>%
  select(ps, exp_bin_treat) %>%
  mutate(group = as.character(exp_bin_treat))
# Histogram
histogram_untrimmed <- ggplot(df_untrimmed_hist, aes(x = ps, fill = group)) +
  geom_histogram(binwidth = 0.025, position = "identity", alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("nothing" = "lightgrey", "metformin" = "black"))  +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    minor_breaks = seq(0, 1, by = 0.05),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity Score Histogram; Untrimmed",
    x = "Propensity Score",
    y = "Count",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "gray80"),
    panel.grid.minor.x = element_line(color = "gray90"),
    axis.text.x = element_text(size = 10)
  )

## Trimmed
df_trimmed_hist <- df_trimmed %>%
  select(ps, exp_bin_treat) %>%
  mutate(group = as.character(exp_bin_treat))
# Histogram
histogram_trimmed <- ggplot(df_trimmed_hist, aes(x = ps, fill = group)) +
  geom_histogram(binwidth = 0.025, position = "identity", alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("nothing" = "lightgrey", "metformin" = "black"))  +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    minor_breaks = seq(0, 1, by = 0.05),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity Score Histogram; Untrimmed",
    x = "Propensity Score",
    y = "Count",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "gray80"),
    panel.grid.minor.x = element_line(color = "gray90"),
    axis.text.x = element_text(size = 10)
  )

# Density plot for HbA1c subgroups -----------------------------------------
print('Density plot for HbA1c subgroups')

dens_treated_HbA1c59orabove <- density(df$ps[df$exp_bin_treat == "metformin" & df$cov_cat_hba1c_b == "59-75"])
df_treated_HbA1c59orabove <- data.frame(dens_x = dens_treated_HbA1c59orabove$x, dens_y = dens_treated_HbA1c59orabove$y, group = "metformin (untrimmed)")
dens_untreated_HbA1c59orabove <- density(df$ps[df$exp_bin_treat == "nothing" & df$cov_cat_hba1c_b == "59-75"])
df_untreated_HbA1c59orabove <- data.frame(dens_x = dens_untreated_HbA1c59orabove$x, dens_y = dens_untreated_HbA1c59orabove$y, group = "nothing (untrimmed)")
ps_density_data_untrimmed_HbA1c59orabove <- bind_rows(df_treated_HbA1c59orabove, df_untreated_HbA1c59orabove)
density_plot_untrimmed_HbA1c59orabove <- ggplot(ps_density_data_untrimmed_HbA1c59orabove, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1), 
    minor_breaks = seq(0, 1, by = 0.05), 
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity Score Density Plot; Untrimmed; _HbA1c59orabove",
    x = "Propensity Score",
    y = "Density",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 10)
  )

dens_treated_HbA1c42to58 <- density(df$ps[df$exp_bin_treat == "metformin" & df$cov_cat_hba1c_b == "42-58"])
df_treated_HbA1c42to58 <- data.frame(dens_x = dens_treated_HbA1c42to58$x, dens_y = dens_treated_HbA1c42to58$y, group = "metformin (untrimmed)")
dens_untreated_HbA1c42to58 <- density(df$ps[df$exp_bin_treat == "nothing" & df$cov_cat_hba1c_b == "42-58"])
df_untreated_HbA1c42to58 <- data.frame(dens_x = dens_untreated_HbA1c42to58$x, dens_y = dens_untreated_HbA1c42to58$y, group = "nothing (untrimmed)")
ps_density_data_untrimmed_HbA1c42to58 <- bind_rows(df_treated_HbA1c42to58, df_untreated_HbA1c42to58)
density_plot_untrimmed_HbA1c42to58 <- ggplot(ps_density_data_untrimmed_HbA1c42to58, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1), 
    minor_breaks = seq(0, 1, by = 0.05), 
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity Score Density Plot; Untrimmed; _HbA1c42to58",
    x = "Propensity Score",
    y = "Density",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 10)
  )

dens_treated_HbA1cbelow42 <- density(df$ps[df$exp_bin_treat == "metformin" & df$cov_cat_hba1c_b == "below 42"])
df_treated_HbA1cbelow42 <- data.frame(dens_x = dens_treated_HbA1cbelow42$x, dens_y = dens_treated_HbA1cbelow42$y, group = "metformin (untrimmed)")
dens_untreated_HbA1cbelow42 <- density(df$ps[df$exp_bin_treat == "nothing" & df$cov_cat_hba1c_b == "below 42"])
df_untreated_HbA1cbelow42 <- data.frame(dens_x = dens_untreated_HbA1cbelow42$x, dens_y = dens_untreated_HbA1cbelow42$y, group = "nothing (untrimmed)")
ps_density_data_untrimmed_HbA1cbelow42 <- bind_rows(df_treated_HbA1cbelow42, df_untreated_HbA1cbelow42)
density_plot_untrimmed_HbA1cbelow42 <- ggplot(ps_density_data_untrimmed_HbA1cbelow42, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1), 
    minor_breaks = seq(0, 1, by = 0.05), 
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity Score Density Plot; Untrimmed; _HbA1cbelow42",
    x = "Propensity Score",
    y = "Density",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 10)
  )


# Save output -------------------------------------------------------------
print('Save output')
# PS model summary
write.csv(ps_model_summary, file = here::here("output", "ps", "ps_model_summary.csv"))
# min/max PS per group
write.csv(ps_summary, file = here::here("output", "ps", "ps_summary.csv"))
# Standardized mean differences, unweighted and weighted (unstablized and stabilized)
write.csv(smd_unweighted, file = here::here("output", "ps", "smd_unweighted.csv"))
write.csv(smd_sw, file = here::here("output", "ps", "smd_weighted_stabilized.csv"))
# Density plot and underlying data
write.csv(ps_density_data_untrimmed, file = here::here("output", "ps", "density_plot_untrimmed.csv"))
write.csv(ps_density_data_trimmed, file = here::here("output", "ps", "density_plot_trimmed.csv"))
write.csv(ps_density_data_untrimmed_HbA1c59orabove, file = here::here("output", "ps", "ps_density_data_untrimmed_HbA1c59orabove.csv"))
write.csv(ps_density_data_untrimmed_HbA1c42to58, file = here::here("output", "ps", "ps_density_data_untrimmed_HbA1c42to58.csv"))
write.csv(ps_density_data_untrimmed_HbA1cbelow42, file = here::here("output", "ps", "ps_density_data_untrimmed_HbA1cbelow42.csv"))
ggsave(filename = here::here("output", "ps", "density_plot_untrimmed.png"), density_plot_untrimmed, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "ps", "density_plot_trimmed.png"), density_plot_trimmed, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "ps", "density_plot_untrimmed_HbA1c59orabove.png"), density_plot_untrimmed_HbA1c59orabove, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "ps", "density_plot_untrimmed_HbA1c42to58.png"), density_plot_untrimmed_HbA1c42to58, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "ps", "density_plot_untrimmed_HbA1cbelow42.png"), density_plot_untrimmed_HbA1cbelow42, width = 20, height = 20, units = "cm")
# Histograms
ggsave(filename = here::here("output", "ps", "histogram_untrimmed.png"), histogram_untrimmed, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "ps", "histogram_trimmed.png"), histogram_trimmed, width = 20, height = 20, units = "cm")
