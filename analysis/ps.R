####
## This script does the following:
# 1. Import processed data
# 2. Create a propensity score for treatment
# 3. Extract SMDs before and after IPW based on PS
# 4. Density plot; untrimmed and trimmed
# 5. Save all output
####

# Import libraries and functions ------------------------------------------
library(arrow)
library(here)
library(dplyr)
library(rms) # splines
library(cobalt) # SMD density
library(ggplot2)

# Create directories for output -------------------------------------------
fs::dir_create(here::here("output", "ps"))

# Import the data ---------------------------------------------------------
df <- read_feather(here("output", "data", "data_processed.arrow"))

# Add splines -------------------------------------------------------------
# Compute knot locations based on percentiles, according to study protocol
age_knots <- quantile(df$cov_num_age, probs = c(0.10, 0.50, 0.90))
# print(age_knots)

# PS model ----------------------------------------------------------------
# Set format and reference to ensure proper logistic regression modeling
df$exp_bin_treat <- factor(df$exp_bin_treat, 
                                       levels = c(0, 1), # Reference first
                                       labels = c("nothing", "metformin"))

# Define the covariates included in the PS model
covariate_names <- names(df) %>%
  grep("^cov_", ., value = TRUE) %>% 
  # exclude those not needed in the model:
  ## cov_bin_obesity is covering for cov_num_bmi & cov_cat_bmi_groups,
  ## cov_cat_hba1c_mmol_mol is covering for cov_num_hba1c_mmol_mol
  ## cov_cat_tc_hdl_ratio is covering for cov_num_tc_hdl_ratio
  ## spline(cov_num_age) is covering for cov_cat_age
  setdiff(c("cov_num_bmi", "cov_cat_bmi_groups", "cov_num_hba1c_mmol_mol", "cov_cat_age", "cov_num_tc_hdl_ratio"
            )) 
# print(covariate_names)

# Construct model formula dynamically
ps_formula <- as.formula(paste("exp_bin_treat ~ rcs(cov_num_age, age_knots) +", paste(covariate_names, collapse = " + "), "+ strata(strat_cat_region)"))

# Fit the PS model to estimate the PS for being in the metformin-mono group
ps_model <- glm(ps_formula, family = binomial(link = "logit"), data = df)
# summary(ps_model)

# Extract propensity scores
df$ps <- predict(ps_model, type = "response")

# Inverse probability weighting -------------------------------------------
# Stabilized
p_treat <- mean(df$exp_bin_treat == "metformin") # prop treated, as per formular for stab weights
df$sw <- ifelse(df$exp_bin_treat == "metformin", 
                                 p_treat / df$ps, 
                                 (1 - p_treat) / (1 - df$ps))

# Assess covariate balance BEFORE IPW -------------------------------------
tbl1_unweighted <- bal.tab(df[covariate_names], 
                          treat = df$exp_bin_treat, 
                          binary = "std",
                          s.d.denom = "pooled") # standardization for binary variables
smd_unweighted <- as.data.frame(tbl1_unweighted$Balance)
smd_unweighted <- smd_unweighted %>% 
  select(Diff.Un) %>% 
  rename("SMD_Unweighted" = Diff.Un)

# Assess covariate balance AFTER IPW --------------------------------------
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
# Determine the common min overlap range in PS
# Reduce to the max. of all group-wise min. PS (i.e. everyone has at least this PS) & the min. of all group-wise max. PS (i.e. no-one has a PS above) 
ps_trim <- df %>%
  group_by(exp_bin_treat) %>%
  summarise(min_ps = min(ps), max_ps = max(ps)) %>%
  ungroup() %>%
  summarise(min_common = max(min_ps), max_common = min(max_ps))
df_trimmed <- df %>%
  filter(ps >= ps_trim$min_common[1] & ps <= ps_trim$max_common[1])

## Min, 25th percentile, median, mean, SD, 75th percentile, and max
summary(df$sw)
sd(df$sw)
summary(df_trimmed$sw)
sd(df_trimmed$sw)

# Density plot ------------------------------------------------------------
## Untrimmed
dens_treated <- density(df$ps[df$exp_bin_treat == "metformin"])
df_treated <- data.frame(dens_x = dens_treated$x, dens_y = dens_treated$y, group = "metformin (untrimmed)")
dens_untreated <- density(df$ps[df$exp_bin_treat == "nothing"])
df_untreated <- data.frame(dens_x = dens_untreated$x, dens_y = dens_untreated$y, group = "nothing (untrimmed)")
ps_density_data_untrimmed <- bind_rows(df_treated, df_untreated)

density_plot_untrimmed <- ggplot(ps_density_data_untrimmed, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  labs(title = "Propensity Score Density Plot; Untrimmed",
       x = "Propensity Score",
       y = "Density",
       color = "Group") +
  theme_minimal()

## Trimmed
dens_treated_trim <- density(df_trimmed$ps[df_trimmed$exp_bin_treat == "metformin"])
df_treated_trim <- data.frame(dens_x = dens_treated_trim$x, dens_y = dens_treated_trim$y, group = "metformin (trimmed)")
dens_untreated_trim <- density(df_trimmed$ps[df_trimmed$exp_bin_treat == "nothing"])
df_untreated_trim <- data.frame(dens_x = dens_untreated_trim$x, dens_y = dens_untreated_trim$y, group = "nothing (trimmed)")
ps_density_data_trimmed <- bind_rows(df_treated_trim, df_untreated_trim)

density_plot_trimmed <- ggplot(ps_density_data_trimmed, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  labs(title = "Propensity Score Density Plot; Trimmed",
       x = "Propensity Score",
       y = "Density",
       color = "Group") +
  theme_minimal()

# Save output -------------------------------------------------------------
# Standardized mean differences, unweighted and weighted (unstablized and stabilized)
write.csv(smd_unweighted, file = here::here("output", "ps", "smd_unweighted.csv"))
write.csv(smd_sw, file = here::here("output", "ps", "smd_weighted_stabilized.csv"))
# Density plot and underlying data
write.csv(ps_density_data_untrimmed, file = here::here("output", "ps", "density_plot_untrimmed.csv"))
write.csv(ps_density_data_trimmed, file = here::here("output", "ps", "density_plot_trimmed.csv"))
ggsave(filename = here::here("output", "ps", "density_plot_untrimmed.png"), density_plot_untrimmed, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "ps", "density_plot_trimmed.png"), density_plot_trimmed, width = 20, height = 20, units = "cm")
