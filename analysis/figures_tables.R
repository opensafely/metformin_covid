####
## This script does the following:
# 1. Import released midpoint rounded disclosure-safe output data
# 2. Create various tables and figures
# Local/dummy-data run only, see the "figures_tables" action in project.yaml -
# inputs are manually-curated, already-released data_release_* snapshots
# (gitignored, local only), so this must never be requested as a real-data job.
####

# Import libraries and user functions -----------------------------
library('dplyr')
library('tidyr')
library('readr')
# library('stringr')
library('here')
# library('DiagrammeR')
library('ggplot2')
# library('RColorBrewer')
library('gt')
library(forestplot)
library(purrr)

# Create output directory ----
fs::dir_create(here::here("output", "figures_tables"))

# Import the data -------------------------------------------------
df_comp <- read_csv(here("output", "data_release_20260331", "data_description", "n_completeness_excluded_midpoint6.csv"))
df_qa <- read_csv(here("output", "data_release_20260331", "data_description", "n_qa_excluded_midpoint6.csv"))
df_elig <- read_csv(here("output", "data_release_20260331", "data_description", "n_elig_excluded_midpoint6.csv"))
df_elig_landmark <- read_csv(here("output", "data_release_20260331", "data_description", "n_restricted_midpoint6.csv"))

df_main <- read_csv(here("output", "data_release_20260331", "data_description", "table1_main_midpoint6.csv"))
df_death_ltfu1 <- read_csv(here("output", "data_release_20260331", "data_description", "table1_death_ltfu1_midpoint6.csv"))
df_death_ltfu2 <- read_csv(here("output", "data_release_20260610", "data_description", "table1_death_ltfu2_midpoint6.csv"))

df_ps <- read_csv(here("output", "data_release_20260331", "ps", "density_plot_untrimmed.csv"))
df_ps_HbA1c59orabove <- read_csv(here("output", "data_release_20260331", "ps", "ps_density_data_untrimmed_HbA1c59orabove.csv"))
df_ps_HbA1c42to58 <- read_csv(here("output", "data_release_20260331", "ps", "ps_density_data_untrimmed_HbA1c42to58.csv"))
df_ps_HbA1cbelow42 <- read_csv(here("output", "data_release_20260331", "ps", "ps_density_data_untrimmed_HbA1cbelow42.csv"))

df_smd <- read_csv(here("output", "data_release_20260610", "ps", "smd_combined.csv"))
df_smd_ps_overlap <- read_csv(here("output", "data_release_20260610", "ps", "smd_combined_ps_overlap.csv"))
df_smd_combined_HbA1cbelow42 <- read_csv(here("output", "data_release_20260610", "ps", "smd_combined_HbA1cbelow42.csv"))
df_smd_combined_HbA1c42to58 <- read_csv(here("output", "data_release_20260610", "ps", "smd_combined_HbA1c42to58.csv"))
df_smd_combined_HbA1c59orabove <- read_csv(here("output", "data_release_20260610", "ps", "smd_combined_HbA1c59orabove.csv"))

df_results_cox <- read_csv(here("output", "data_release_20260610", "make_output", "results_cox_midpoint6.csv"))

df_results_km_primary <- read_csv(here("output", "data_release_20260331", "te", "km_primary", "estimates.csv"))
df_results_km_covid_event <- read_csv(here("output", "data_release_20260331", "te", "km_covid_event", "estimates.csv"))
df_results_km_longvirfat <- read_csv(here("output", "data_release_20260331", "te", "km_longvirfat", "estimates.csv"))
df_results_km_primary_weighted <- read_csv(here("output", "data_release_20260610", "te", "km_primary_weighted", "estimates.csv"))

df_results_ppa_interaction <- read_csv(here("output", "data_release_20260610", "te", "pooled_log_reg", "cum_risk_treat_ltfu_ci_severecovid.csv"))
df_results_ppa_no_interaction <- read_csv(here("output", "data_release_20260610", "te", "pooled_log_reg", "cum_risk_treat_ltfu_ci_severecovid_no_interaction.csv"))
df_results_ppa_interaction_untrunc <- read_csv(here("output", "data_release_20260610", "te", "pooled_log_reg", "cum_risk_treat_ltfu_ci_severecovid_untrunc.csv"))

df_results_ppa_cens_ltfu <- read_csv(here("output", "data_release_20260610", "te", "pooled_log_reg", "censoring_ltfu_rates.csv"))
df_results_ppa_comp <- read_csv(here("output", "data_release_20260610", "te", "pooled_log_reg", "censoring_comp_rates.csv"))

# Baseline table, based on output from table1.R ----------------------

# Hide NAs where not needed
hide_NA <- function(df) {
  df %>%
    mutate(metformin = replace_na(as.character(metformin), "")) %>% 
    mutate(control = replace_na(as.character(control), ""))
}
df_main <- df_main %>% 
  mutate(metformin = replace_na(as.character(metformin), "")) %>% 
  mutate(control = replace_na(as.character(control), ""))
df_death_ltfu1 <- df_death_ltfu1 %>% 
  mutate(`Alive and in care at landmark` = replace_na(as.character(`Alive and in care at landmark`), "")) %>% 
  mutate(`Died or LTFU until landmark` = replace_na(as.character(`Died or LTFU until landmark`), ""))
df_death_ltfu2 <- df_death_ltfu2 %>% 
  mutate(`Alive and in care at pandemic start` = replace_na(as.character(`Alive and in care at pandemic start`), "")) %>% 
  mutate(`Died or LTFU between landmark and pandemic start` = replace_na(as.character(`Died or LTFU between landmark and pandemic start`), ""))

# Reduce to main variables
df_main <- df_main %>%
  filter(!grepl("INT: Any", var_label)) %>%
  filter(!grepl("CONT: Any", var_label)) %>%
  filter(!grepl("Diet", var_label)) %>%
  filter(!grepl("Calendar", var_label)) %>%
  filter(!grepl("COVID hosp or death", var_label)) %>%
  filter(!grepl("COVID hosp", var_label)) %>%
  filter(!grepl("COVID death", var_label)) %>%
  filter(!grepl("Any covid diagnosis", var_label)) %>%
  filter(!grepl("Any death", var_label)) %>%
  filter(!grepl("Any Viral Fatigue", var_label)) %>%
  filter(!grepl("Any LTFU", var_label)) %>%
  filter(!grepl("Any Long COVID", var_label))
df_death_ltfu1 <- df_death_ltfu1 %>%
  filter(!grepl("INT: Any", var_label)) %>%
  filter(!grepl("CONT: Any", var_label)) %>%
  filter(!grepl("Diet", var_label)) %>%
  filter(!grepl("Calendar", var_label)) %>%
  filter(!grepl("COVID hosp or death", var_label)) %>%
  filter(!grepl("COVID hosp", var_label)) %>%
  filter(!grepl("COVID death", var_label)) %>%
  filter(!grepl("Any covid diagnosis", var_label)) %>%
  filter(!grepl("Any death", var_label)) %>%
  filter(!grepl("Any Viral Fatigue", var_label)) %>%
  filter(!grepl("Any LTFU", var_label)) %>%
  filter(!grepl("Any Long COVID", var_label))
df_death_ltfu2 <- df_death_ltfu2 %>%
  filter(!grepl("INT: Any", var_label)) %>%
  filter(!grepl("CONT: Any", var_label)) %>%
  filter(!grepl("Diet", var_label)) %>%
  filter(!grepl("Calendar", var_label)) %>%
  filter(!grepl("COVID hosp or death", var_label)) %>%
  filter(!grepl("COVID hosp", var_label)) %>%
  filter(!grepl("COVID death", var_label)) %>%
  filter(!grepl("Any covid diagnosis", var_label)) %>%
  filter(!grepl("Any death", var_label)) %>%
  filter(!grepl("Any Viral Fatigue", var_label)) %>%
  filter(!grepl("Any LTFU", var_label)) %>%
  filter(!grepl("Any Long COVID", var_label))

# Create the gt tables
tbl_gt_main <- df_main %>%
  select(label, `metformin`, control) %>%
  gt() %>%
  tab_header(title = "Baseline Characteristics") %>%
  cols_label(
    label = "Characteristic",
    `metformin` = "Metformin",
    control = "Nothing"
  ) %>%
  tab_options(
    table.font.size = px(12)
  )
print(tbl_gt_main)
# Save as RTF (Word-readable)
gtsave(tbl_gt_main, here::here("output", "figures_tables", "baseline_characteristics_main.rtf"))

tbl_gt_death_ltfu1 <- df_death_ltfu1 %>%
  select(label, `Alive and in care at landmark`, `Died or LTFU until landmark`) %>%
  gt() %>%
  tab_header(title = "Baseline Characteristics") %>%
  cols_label(
    label = "Characteristic",
    `Alive and in care at landmark` = "Alive and in care at landmark",
    `Died or LTFU until landmark` = "Died or LTFU until landmark"
  ) %>%
  tab_options(
    table.font.size = px(12)
  )
print(tbl_gt_death_ltfu1)
# Save as RTF (Word-readable)
gtsave(tbl_gt_death_ltfu1, here::here("output", "figures_tables", "baseline_characteristics_death_ltfu1.rtf"))

tbl_gt_death_ltfu2 <- df_death_ltfu2 %>%
  select(label, `Alive and in care at pandemic start`, `Died or LTFU between landmark and pandemic start`) %>%
  gt() %>%
  tab_header(title = "Baseline Characteristics") %>%
  cols_label(
    label = "Characteristic",
    `Alive and in care at pandemic start` = "Alive and in care at pandemic start",
    `Died or LTFU between landmark and pandemic start` = "Died or LTFU between landmark and pandemic start"
  ) %>%
  tab_options(
    table.font.size = px(12)
  )
print(tbl_gt_death_ltfu2)
# Save as RTF (Word-readable)
gtsave(tbl_gt_death_ltfu2, here::here("output", "figures_tables", "baseline_characteristics_death_ltfu2.rtf"))


# PS figures ----------------------
density_plot_totpop <- ggplot(df_ps, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1), 
    minor_breaks = seq(0, 1, by = 0.05), 
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity score - density plot",
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
density_plot_HbA1c59orabove <- ggplot(df_ps_HbA1c59orabove, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1), 
    minor_breaks = seq(0, 1, by = 0.05), 
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity score - density plot - HbA1c 59 or above",
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
density_plot_HbA1c42to58 <- ggplot(df_ps_HbA1c42to58, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1), 
    minor_breaks = seq(0, 1, by = 0.05), 
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity score - density plot - HbA1c 42 to 58",
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
density_plot_HbA1cbelow42 <- ggplot(df_ps_HbA1cbelow42, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1), 
    minor_breaks = seq(0, 1, by = 0.05), 
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity score - density plot - HbA1c below 42",
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


# Love plots by HbA1c subgroup and PS overlap region ----------------------
print('Love plots by HbA1c subgroup and PS overlap region')

ggplot(df_smd, aes(y = Covariate)) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "red", linewidth = 0.5) +
    geom_point(aes(x = SMD_Unweighted, colour = "Unweighted"), shape = 1, size = 2.5) +
    geom_point(aes(x = SMD_Weighted,   colour = "Weighted (stabilized)"), shape = 16, size = 2.5) +
    scale_colour_manual(values = c("Unweighted" = "grey50", "Weighted (stabilized)" = "#2E9FDF")) +
    labs(x = "Standardized Mean Difference", y = NULL, colour = NULL,
         title = paste("Covariate Balance Before and After IPW")) +
    theme_minimal() +
    theme(
      axis.line.x = element_line(colour = "black"),
      legend.position = "bottom",
      panel.grid.major.y = element_line(colour = "gray90"),
      panel.grid.minor = element_blank()
    )
ggplot(df_smd_ps_overlap, aes(y = Covariate)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "red", linewidth = 0.5) +
  geom_point(aes(x = SMD_Unweighted, colour = "Unweighted"), shape = 1, size = 2.5) +
  geom_point(aes(x = SMD_Weighted,   colour = "Weighted (stabilized)"), shape = 16, size = 2.5) +
  scale_colour_manual(values = c("Unweighted" = "grey50", "Weighted (stabilized)" = "#2E9FDF")) +
  labs(x = "Standardized Mean Difference", y = NULL, colour = NULL,
       title = paste("Covariate Balance Before and After IPW; PS-overlap population")) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(colour = "black"),
    legend.position = "bottom",
    panel.grid.major.y = element_line(colour = "gray90"),
    panel.grid.minor = element_blank()
  )
ggplot(df_smd_combined_HbA1cbelow42, aes(y = Covariate)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "red", linewidth = 0.5) +
  geom_point(aes(x = SMD_Unweighted, colour = "Unweighted"), shape = 1, size = 2.5) +
  geom_point(aes(x = SMD_Weighted,   colour = "Weighted (stabilized)"), shape = 16, size = 2.5) +
  scale_colour_manual(values = c("Unweighted" = "grey50", "Weighted (stabilized)" = "#2E9FDF")) +
  labs(x = "Standardized Mean Difference", y = NULL, colour = NULL,
       title = paste("Covariate Balance Before and After IPW; HbA1c below 42 mmol/mol")) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(colour = "black"),
    legend.position = "bottom",
    panel.grid.major.y = element_line(colour = "gray90"),
    panel.grid.minor = element_blank()
  )
ggplot(df_smd_combined_HbA1c42to58, aes(y = Covariate)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "red", linewidth = 0.5) +
  geom_point(aes(x = SMD_Unweighted, colour = "Unweighted"), shape = 1, size = 2.5) +
  geom_point(aes(x = SMD_Weighted,   colour = "Weighted (stabilized)"), shape = 16, size = 2.5) +
  scale_colour_manual(values = c("Unweighted" = "grey50", "Weighted (stabilized)" = "#2E9FDF")) +
  labs(x = "Standardized Mean Difference", y = NULL, colour = NULL,
       title = paste("Covariate Balance Before and After IPW; HbA1c 42 to 58 mmol/mol")) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(colour = "black"),
    legend.position = "bottom",
    panel.grid.major.y = element_line(colour = "gray90"),
    panel.grid.minor = element_blank()
  )
ggplot(df_smd_combined_HbA1c59orabove, aes(y = Covariate)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "red", linewidth = 0.5) +
  geom_point(aes(x = SMD_Unweighted, colour = "Unweighted"), shape = 1, size = 2.5) +
  geom_point(aes(x = SMD_Weighted,   colour = "Weighted (stabilized)"), shape = 16, size = 2.5) +
  scale_colour_manual(values = c("Unweighted" = "grey50", "Weighted (stabilized)" = "#2E9FDF")) +
  labs(x = "Standardized Mean Difference", y = NULL, colour = NULL,
       title = paste("Covariate Balance Before and After IPW; HbA1c 59 mmol/mol or above")) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(colour = "black"),
    legend.position = "bottom",
    panel.grid.major.y = element_line(colour = "gray90"),
    panel.grid.minor = element_blank()
  )

# Results figures ----------------------
## Main results ----
df_results <- df_results_cox %>% 
  filter(term == "days0_730" | term == "days_pre")

df_results_main <- df_results %>% 
  filter(!grepl("sg_", name)) %>% 
  filter(!grepl("_landmark", name)) %>% 
  filter(grepl("mdl_max_adj", model))

df_results_main <- df_results_main %>%
  mutate(
    rn = row_number(),
    events_c = if_else(rn %% 2 == 1,
                           lead(N_events_midpoint6),
                           NA_real_)
  ) %>%
  filter(rn %% 2 == 1) %>% # keep only odd rows
  select(-rn)
df_results_main <- df_results_main %>%
  rename(events_i = N_events_midpoint6)

df_results_main <- df_results_main %>%
  rename(tot_i = N_exposed_midpoint6) %>% 
  mutate(tot_c = N_total_midpoint6-tot_i) %>% 
  mutate(pct_i = round(100*(events_i/tot_i),2)) %>% 
  mutate(pct_c = round(100*(events_c/tot_c),2))

df_results_main <- df_results_main %>%
  mutate(
    events_i = paste0(events_i, " (", sprintf("%.2f", pct_i), "%)"),
    events_c = paste0(events_c, " (", sprintf("%.2f", pct_c), "%)")
  )

df_results_main$name <- recode(df_results_main$name,
                               "primary" = "Primary: COVID-19-related hospitalization or death",
                               "covid_event" = "Secondary: Any COVID-19 diagnosis",
                               "longvirfat" = "Secondary: Long COVID",
                               "neg_control_pandemic" = "Negative control outcome: Bone fracture",
                               "pos_control_pandemic" = "Positive control outcome: Diabetes-related death",
                               "primary_ps_overlap" = "Sensitivity: Propensity score overlap population"
)
desired_order <- c(
  "Primary: COVID-19-related hospitalization or death",
  "Secondary: Any COVID-19 diagnosis",
  "Secondary: Long COVID",
  "Negative control outcome: Bone fracture",
  "Positive control outcome: Diabetes-related death",
  "Sensitivity: Propensity score overlap population"
)
df_results_main <- df_results_main %>%
  mutate(
    name = factor(name, levels = desired_order)
  ) %>%
  arrange(name)

# build forestplot
base_data <- tibble(
  mean = df_results_main$hr,
  lower = df_results_main$conf_low,
  upper = df_results_main$conf_high,
  name = as.character(df_results_main$name),
  events_i = as.character(df_results_main$events_i),
  events_c = as.character(df_results_main$events_c),
  estimates = paste0(formatC(df_results_main$hr, format = "f", digits = 2), 
                     " (", formatC(df_results_main$conf_low, format = "f", digits = 2), 
                     " - ", formatC(df_results_main$conf_high, format = "f", digits = 2), ")"))
header <- tibble(
  name = "Outcome",
  events_i = "Events Metformin\n n (%)",
  events_c = "Events No Antidiabetic\n n (%)",
  estimates = "aHR (95% CI)",
  mean = NA, lower = NA, upper = NA)

fp <- bind_rows(header, base_data)
font <- "sans"

# Create the forest plot
pdf(file = here::here("output", "figures_tables", "Fp_results_main.pdf"),
    width = 14,
    height = 7)

fp %>%
  forestplot(
    labeltext = c(name, 
                  events_i, events_c, 
                  estimates),
    mean = mean, lower = lower, upper = upper, # Numeric columns for the plot
    txt_gp = fpTxtGp(
      label = gpar(fontfamily = font, cex = 1),
      ticks = gpar(cex = 0.88),
      xlab = gpar(cex = 0.88)
    ),
    graph.pos = 4, 
    hrzl_lines = list("2" = gpar(lty = 2),
                      "3" = gpar(lty = 2),
                      "4" = gpar(lty = 2),
                      "5" = gpar(lty = 2),
                      "6" = gpar(lty = 2),
                      "7" = gpar(lty = 2),
                      "8" = gpar(lty = 2)),
    xlog = TRUE,
    # xticks = log(c(0.80, 1, 1.25)),
    lty.ci = 1,
    col = fpColors(
      box = "maroon4",
      line = "maroon1",
      summary = "magenta4",
      hrz_lines = "gray63"
    ),
    vertices = TRUE,
    xlab = NULL,
    zero = 1
  )

grid.text(
  "Favours Metformin group < > Favours No Antidiabetic group",
  x = unit(0.72, "npc"),
  y = unit(0.015, "npc"),
  gp = gpar(cex = 0.9)
)

dev.off()


## Minimally adjusted results ----
df_results_min <- df_results %>% 
  filter(!grepl("sg_", name)) %>% 
  filter(!grepl("_landmark", name)) %>% 
  filter(!grepl("mdl_max_adj", model))

df_results_min <- df_results_min %>%
  mutate(
    rn = row_number(),
    events_c = if_else(rn %% 2 == 1,
                       lead(N_events_midpoint6),
                       NA_real_)
  ) %>%
  filter(rn %% 2 == 1) %>% # keep only odd rows
  select(-rn)
df_results_min <- df_results_min %>%
  rename(events_i = N_events_midpoint6)

df_results_min <- df_results_min %>%
  rename(tot_i = N_exposed_midpoint6) %>% 
  mutate(tot_c = N_total_midpoint6-tot_i) %>% 
  mutate(pct_i = round(100*(events_i/tot_i),2)) %>% 
  mutate(pct_c = round(100*(events_c/tot_c),2))

df_results_min <- df_results_min %>%
  mutate(
    events_i = paste0(events_i, " (", sprintf("%.2f", pct_i), "%)"),
    events_c = paste0(events_c, " (", sprintf("%.2f", pct_c), "%)")
  )

df_results_min$name <- recode(df_results_min$name,
                              "primary" = "Primary: COVID-19-related hospitalization or death",
                              "covid_event" = "Secondary: Any COVID-19 diagnosis",
                              "longvirfat" = "Secondary: Long COVID",
                              "neg_control_pandemic" = "Negative control outcome: Bone fracture",
                              "pos_control_pandemic" = "Positive control outcome: Diabetes-related death",
                              "primary_ps_overlap" = "Sensitivity: Propensity score overlap population"
)
desired_order <- c(
  "Primary: COVID-19-related hospitalization or death",
  "Secondary: Any COVID-19 diagnosis",
  "Secondary: Long COVID",
  "Negative control outcome: Bone fracture",
  "Positive control outcome: Diabetes-related death",
  "Sensitivity: Propensity score overlap population"
)
df_results_min <- df_results_min %>%
  mutate(
    name = factor(name, levels = desired_order)
  ) %>%
  arrange(name)

# build forestplot
base_data <- tibble(
  mean = df_results_min$hr,
  lower = df_results_min$conf_low,
  upper = df_results_min$conf_high,
  name = as.character(df_results_min$name),
  events_i = as.character(df_results_min$events_i),
  events_c = as.character(df_results_min$events_c),
  estimates = paste0(formatC(df_results_min$hr, format = "f", digits = 2), 
                     " (", formatC(df_results_min$conf_low, format = "f", digits = 2), 
                     " - ", formatC(df_results_min$conf_high, format = "f", digits = 2), ")"))
header <- tibble(
  name = "Outcome",
  events_i = "Events Metformin\n n (%)",
  events_c = "Events No Antidiabetic\n n (%)",
  estimates = "aHR (95% CI)",
  mean = NA, lower = NA, upper = NA)

fp <- bind_rows(header, base_data)
font <- "sans"

# Create the forest plot
pdf(file = here::here("output", "figures_tables", "Fp_results_min.pdf"),
    width = 14,
    height = 7)

fp %>%
  forestplot(
    labeltext = c(name, 
                  events_i, events_c, 
                  estimates),
    mean = mean, lower = lower, upper = upper, # Numeric columns for the plot
    txt_gp = fpTxtGp(
      label = gpar(fontfamily = font, cex = 1),
      ticks = gpar(cex = 0.88),
      xlab = gpar(cex = 0.88)
    ),
    graph.pos = 4, 
    hrzl_lines = list("2" = gpar(lty = 2),
                      "3" = gpar(lty = 2),
                      "4" = gpar(lty = 2),
                      "5" = gpar(lty = 2),
                      "6" = gpar(lty = 2),
                      "7" = gpar(lty = 2),
                      "8" = gpar(lty = 2)),
    xlog = TRUE,
    # xticks = log(c(0.80, 1, 1.25)),
    lty.ci = 1,
    col = fpColors(
      box = "maroon4",
      line = "maroon1",
      summary = "magenta4",
      hrz_lines = "gray63"
    ),
    vertices = TRUE,
    xlab = NULL,
    zero = 1
  )

grid.text(
  "Favours Metformin group < > Favours No Antidiabetic group",
  x = unit(0.72, "npc"),
  y = unit(0.015, "npc"),
  gp = gpar(cex = 0.9)
)

dev.off()


## Subgroup analyses ----
### Run Qtest subgroups on main output
cochran_Q <- function(lnhr, se_lnhr) {
  w  <- 1 / se_lnhr^2
  b  <- sum(w * lnhr) / sum(w)
  Q  <- sum(w * (lnhr - b)^2)
  df <- length(lnhr) - 1
  tibble(Q = Q, df = df, p_heterogeneity = pchisq(Q, df, lower.tail = FALSE))
}

subgroup_map <- tribble(
  ~name, ~subgroup_var,
  "sg_60orabove", "age",
  "sg_below60", "age",
  "sg_female", "sex",
  "sg_male", "sex",
  "sg_nonwhite", "ethnicity",
  "sg_white", "ethnicity", 
  "sg_imd1", "deprivation",
  "sg_nonimd1", "deprivation",
  "sg_obese", "bmi",
  "sg_overweight", "bmi",
  "sg_normlowmissing", "bmi",
  "sg_belowHbA1c59", "hba1c",
  "sg_HbA1c59orabove", "hba1c"
)

het_results <- df_results |>
  filter(term == "days0_730",        # the treatment/exposure effect
         model == "mdl_max_adj") |>   # pick the adjustment level you want to report
  inner_join(subgroup_map, by = "name") |>
  group_by(subgroup_var) |>
  group_modify(~ cochran_Q(.x$lnhr, .x$se_lnhr)) |>
  ungroup()

phet_lookup <- het_results %>%
  transmute(
    characteristic = recode(subgroup_var,
                            age = "Age", sex = "Sex", ethnicity = "Ethnicity",
                            deprivation = "Deprivation", bmi = "BMI", hba1c = "HbA1c"),
    p_het = ifelse(p_heterogeneity < 0.001, "<0.001",
                   formatC(p_heterogeneity, format = "f", digits = 3))
  )

## split the combined name into characteristic + level
df_results_sg <- df_results %>%
  filter(grepl("sg_", name)) %>%
  filter(!grepl("_pos_control", name)) %>%
  filter(grepl("mdl_max_adj", model))

df_results_sg <- df_results_sg %>%
  mutate(
    rn = row_number(),
    events_c = if_else(rn %% 2 == 1,
                       lead(N_events_midpoint6),
                       NA_real_)
  ) %>%
  filter(rn %% 2 == 1) %>% # keep only odd rows
  select(-rn)
df_results_sg <- df_results_sg %>%
  rename(events_i = N_events_midpoint6)

df_results_sg <- df_results_sg %>%
  rename(tot_i = N_exposed_midpoint6) %>%
  mutate(tot_c = N_total_midpoint6-tot_i) %>%
  mutate(pct_i = round(100*(events_i/tot_i),2)) %>%
  mutate(pct_c = round(100*(events_c/tot_c),2))

df_results_sg <- df_results_sg %>%
  mutate(
    events_i = paste0(events_i, " (", sprintf("%.2f", pct_i), "%)"),
    events_c = paste0(events_c, " (", sprintf("%.2f", pct_c), "%)")
  )

df_results_sg$name <- recode(df_results_sg$name,
                             "sg_60orabove" = "Age: 60 years or above",
                             "sg_below60" = "Age: below 60 years",
                             "sg_female" = "Sex: Female",
                             "sg_male" = "Sex: Male",
                             "sg_nonwhite" = "Ethnicity: Non-White",
                             "sg_white" = "Ethnicity: White",
                             "sg_imd1" = "Deprivation: Most deprived",
                             "sg_nonimd1" = "Deprivation: Other",
                             "sg_obese" = "BMI: Obese",
                             "sg_overweight" = "BMI: Overweight",
                             "sg_normlowmissing" = "BMI: Normal/low weight or missing",
                             "sg_belowHbA1c59" = "HbA1c: below 59 mmol/mol",
                             "sg_HbA1c59orabove" = "HbA1c: 59 mmol/mol or above"
)
desired_order <- c(
  "Age: 60 years or above",
  "Age: below 60 years",
  "Sex: Female",
  "Sex: Male",
  "Ethnicity: Non-White",
  "Ethnicity: White",
  "Deprivation: Most deprived",
  "Deprivation: Other",
  "BMI: Obese",
  "BMI: Overweight",
  "BMI: Normal/low weight or missing",
  "HbA1c: below 59 mmol/mol",
  "HbA1c: 59 mmol/mol or above"
)
df_results_sg <- df_results_sg %>%
  mutate(
    name = factor(name, levels = desired_order)
  ) %>%
  arrange(name)

df_results_sg <- df_results_sg %>%
  mutate(name = as.character(name)) %>%
  separate(name, into = c("characteristic", "level"), sep = ": ", remove = FALSE)

## assemble grouped rows: bold header + levels + trailing spacer, per characteristic
char_order <- c("Age", "Sex", "Ethnicity", "Deprivation", "BMI", "HbA1c")

make_block <- function(chr) {
  lv <- df_results_sg %>% filter(characteristic == chr)
  ph <- phet_lookup %>% filter(characteristic == chr) %>% pull(p_het)
  
  header <- tibble(
    characteristic = chr, level = "",
    events_i = "", events_c = "", estimates = "", p_het = ph,
    mean = NA_real_, lower = NA_real_, upper = NA_real_, is_summary = TRUE
  )
  rows <- tibble(
    characteristic = "", level = lv$level,
    events_i = lv$events_i, events_c = lv$events_c,
    estimates = paste0(formatC(lv$hr, format = "f", digits = 2), " (",
                       formatC(lv$conf_low,  format = "f", digits = 2), " - ",
                       formatC(lv$conf_high, format = "f", digits = 2), ")"),
    p_het = "",
    mean = lv$hr, lower = lv$conf_low, upper = lv$conf_high, is_summary = FALSE
  )
  spacer <- tibble(
    characteristic = "", level = "", events_i = "", events_c = "",
    estimates = "", p_het = "", mean = NA_real_, lower = NA_real_,
    upper = NA_real_, is_summary = FALSE
  )
  bind_rows(header, rows, spacer)
}

body <- map_dfr(char_order, make_block)
body <- body %>% slice(1:(n() - 1))

## column-title row
title <- tibble(
  characteristic = "Characteristic", level = "Group",
  events_i  = "Events Metformin\n n (%)",
  events_c  = "Events No Antidiabetic\n n (%)",
  estimates = "aHR (95% CI)",
  p_het     = "P-value for\nheterogeneity",
  mean = NA_real_, lower = NA_real_, upper = NA_real_, is_summary = TRUE
)

fp <- bind_rows(title, body)
font <- "sans"

## --- forest plot ---
pdf(file = here::here("output", "figures_tables", "Fp_results_sg.pdf"), width = 24, height = 18)

fp %>%
  forestplot(
    labeltext = c(characteristic, level, events_i, events_c, estimates, p_het),
    mean = mean, lower = lower, upper = upper,
    is.summary = fp$is_summary,                 # bolds header + title rows, no box where mean is NA
    align = c("l", "l", "c", "c", "c", "c"),
    txt_gp = fpTxtGp(
      label   = gpar(fontfamily = font, cex = 1.4),   # row labels (was 1)
      summary = gpar(fontfamily = font, cex = 1.5),   # bold header rows (title + characteristics)
      ticks   = gpar(cex = 1.2),                      # x-axis numbers (was 0.88)
      xlab    = gpar(cex = 1.2)                        # x-axis title (was 0.88)
    ),
    graph.pos = 5,
    graphwidth = unit(8, "cm"),
    lineheight = unit(1.4, "cm"), # <- taller rows; fixes the 2-line header overlap
    hrzl_lines = list("2" = gpar(lty = 1, col = "gray63")),
    xlog = TRUE,
    lty.ci = 1,
    col = fpColors(
      box = "maroon4",
      line = "maroon1",
      summary = "magenta4",
      hrz_lines = "gray63"
    ),
    vertices = TRUE,
    xlab = NULL,
    zero = 1
  )

grid.text(
  "Favours Metformin group < > Favours No Antidiabetic group",
  x = unit(0.652, "npc"),
  y = unit(0.008, "npc"),
  gp = gpar(cex = 1.2) 
)

dev.off()


## KM curves (for those not provided through direct output) ----
exposure  <- "exp_bin_treat" 
subgroups <- character() # no subgroups
exposure_syms  <- syms(exposure)
subgroup_syms  <- syms(subgroups)

# --- annotation text (end-of-follow-up results) ---
annotation_text <- paste(
  "At end of follow-up:",
  "Cumulative risk, metformin: 2.23% (2.02 to 2.45)",
  "Cumulative risk, no antidiabetic: 1.87% (1.72 to 2.03)",
  "Adjusted risk ratio: 1.19 (1.05 to 1.35)",
  "Adjusted risk difference: 0.36% (0.09 to 0.62)",
  sep = "\n"
)

# time = 0 origin row is added once per curve.
data_surv_rounded <-
  df_results_km_primary_weighted |>
  mutate(
    # 0 = No Antidiabetic (reference), 1 = Metformin
    exp_bin_treat = factor(exp_bin_treat,
                           levels = c(0, 1),
                           labels = c("No Antidiabetic", "Metformin"))
  ) |>
  group_by(!!!subgroup_syms, !!!exposure_syms)

km_plot <- function(.data) {
  
  data_with_time0 <-
    if(length(exposure)>0L){
      .data |>
        mutate(
          "{exposure}" := as.factor(!!!exposure_syms)
        )
    } else {
      .data
    }
  
  data_with_time0 <-
    data_with_time0 |>
    mutate(
      lagtime = lag(time, 1, 0), # assumes the time-origin is zero
    ) %>%
    group_modify(
      ~ add_row(
        .x,
        time = 0, # assumes time origin is zero
        lagtime = 0,
        cmlinc = 0,
        cmlinc.low = 0,
        cmlinc.high = 0,
        .before = 0
      )
    )
  
  ggplot_init <- if(length(exposure)==0L){
    ggplot(data_with_time0)
  } else {
    exposure_sym <- sym(exposure)
    ggplot(data_with_time0, aes(group = !!exposure_sym, colour = !!exposure_sym, fill = !!exposure_sym))
  }
  ggplot_init +
    geom_step(aes(x = time, y = cmlinc), direction = "vh") +
    geom_step(aes(x = time, y = cmlinc), direction = "vh", linetype = "dashed", alpha = 0.5) +
    geom_rect(aes(xmin = lagtime, xmax = time, ymin = cmlinc.low, ymax = cmlinc.high), alpha = 0.1, colour = "transparent") +
    # --- results annotation, anchored bottom-right (matches the other figure) ---
    annotate("text",
             x = max(df_results_km_primary_weighted$time),
             y = max(df_results_km_primary_weighted$cmlinc.high, na.rm = TRUE) * 0.04,
             hjust = 1, vjust = 0,
             size = 3.8,
             label = annotation_text) +
    facet_grid(rows = vars(!!!subgroup_syms)) +
    scale_color_manual(
      name = NULL,
      values = c("No Antidiabetic" = "#E7B800",
                 "Metformin"       = "#2E9FDF"),
      na.value = "grey"
    ) +
    scale_fill_manual(
      values = c("No Antidiabetic" = "#E7B800",
                 "Metformin"       = "#2E9FDF"),
      na.value = "grey",
      guide = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    coord_cartesian(xlim = c(0, NA)) +
    xlab("Months") +
    ylab("Cumulative Incidence") +
    theme_minimal() +
    theme(axis.text = element_text(size = 14), legend.position = c(0.2, 0.8),
          axis.line = element_line(colour = "black"),
          legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())
}

p <- km_plot(data_surv_rounded)
ggsave(filename = here::here("output", "figures_tables", "KM_curves_weighted.png"), p, width = 20, height = 20, units = "cm")


## PPA curves ----
# --- annotation text (end-of-follow-up results) ---
annotation_text <- paste(
  "At end of follow-up:",
  "Cumulative risk, metformin: 2.12% (1.87 to 2.35)",
  "Cumulative risk, no antidiabetic: 1.92% (1.73 to 2.10)",
  "Adjusted risk ratio: 1.11 (0.95 to 1.27)",
  "Adjusted risk difference: 0.21% (-0.11 to 0.49)",
  sep = "\n"
)

plot_cum_risk_treat_ltfu_ci_severecovid <- ggplot(df_results_ppa_interaction,
                                                  aes(x = time)) +
  geom_line(aes(y = mean.1_midpoint6, # create line for intervention group
                color = "Metformin"), linewidth = 1.5) +
  geom_ribbon(aes(ymin = ll.1_midpoint6, ymax = ul.1_midpoint6, fill = "Metformin"), alpha = 0.4) +
  geom_line(aes(y = mean.0_midpoint6, # create line for control group
                color = "No Antidiabetic"), linewidth = 1.5) +
  geom_ribbon(aes(ymin = ll.0_midpoint6, ymax = ul.0_midpoint6, fill = "No Antidiabetic"), alpha = 0.4) +
  # --- results annotation, anchored bottom-right ---
  annotate("text",
           x = max(df_results_ppa_interaction$time),
           y = 0,
           hjust = 1, vjust = 0,
           size = 3.8,
           label = annotation_text) +
  xlab("Months") +
  ylab("Cumulative Incidence") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14), legend.position = c(0.2, 0.8),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_color_manual(
    name = "",
    values = c("No Antidiabetic" = "#E7B800",
               "Metformin" = "#2E9FDF")
  ) +
  scale_fill_manual(
    name = "",
    values = c("No Antidiabetic" = "#E7B800",
               "Metformin" = "#2E9FDF")
  )

plot_cum_risk_treat_ltfu_ci_severecovid_untrunc <- ggplot(df_results_ppa_interaction_untrunc,
                                                  aes(x = time)) +
  geom_line(aes(y = mean.1_midpoint6, # create line for intervention group
                color = "Metformin"), linewidth = 1.5) +
  geom_ribbon(aes(ymin = ll.1_midpoint6, ymax = ul.1_midpoint6, fill = "Metformin"), alpha = 0.4) +
  geom_line(aes(y = mean.0_midpoint6, # create line for control group
                color = "No Antidiabetic"), linewidth = 1.5) +
  geom_ribbon(aes(ymin = ll.0_midpoint6, ymax = ul.0_midpoint6, fill = "No Antidiabetic"), alpha = 0.4) +
  xlab("Months") +
  ylab("Cumulative Incidence") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14), legend.position = c(0.2, 0.8),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_color_manual(
    name = "",
    values = c("No Antidiabetic" = "#E7B800",
               "Metformin" = "#2E9FDF")
  ) +
  scale_fill_manual(
    name = "",
    values = c("No Antidiabetic" = "#E7B800",
               "Metformin" = "#2E9FDF")
  )
plot_cum_risk_treat_ltfu_ci_severecovid_no_interaction <- ggplot(df_results_ppa_no_interaction,
                                                  aes(x = time)) +
  geom_line(aes(y = mean.1_midpoint6, # create line for intervention group
                color = "Metformin"), linewidth = 1.5) +
  geom_ribbon(aes(ymin = ll.1_midpoint6, ymax = ul.1_midpoint6, fill = "Metformin"), alpha = 0.4) +
  geom_line(aes(y = mean.0_midpoint6, # create line for control group
                color = "No Antidiabetic"), linewidth = 1.5) +
  geom_ribbon(aes(ymin = ll.0_midpoint6, ymax = ul.0_midpoint6, fill = "No Antidiabetic"), alpha = 0.4) +
  xlab("Months") +
  ylab("Cumulative Incidence") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14), legend.position = c(0.2, 0.8),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_color_manual(
    name = "",
    values = c("No Antidiabetic" = "#E7B800",
               "Metformin" = "#2E9FDF")
  ) +
  scale_fill_manual(
    name = "",
    values = c("No Antidiabetic" = "#E7B800",
               "Metformin" = "#2E9FDF")
  )

# save
ggsave(filename = here::here("output", "figures_tables", "plot_cum_risk_treat_ltfu.png"), plot_cum_risk_treat_ltfu_ci_severecovid,
       width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "figures_tables", "plot_cum_risk_treat_ltfu_untrunc.png"), plot_cum_risk_treat_ltfu_ci_severecovid_untrunc,
       width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "figures_tables", "plot_cum_risk_treat_ltfu_no_interaction.png"), plot_cum_risk_treat_ltfu_ci_severecovid_no_interaction,
       width = 20, height = 20, units = "cm")


# Extract numeric results from weighted KM risk estimates ----
## exposure groups: 1 = metformin (exposed), 0 = no antidiabetic (reference)
## end of follow-up = the last available time point
fup_end <- max(df_results_km_primary_weighted$time)

# reconstruct cmlinc.se per group (Greenwood, as in round_km reusable action)
# the concise df_results_km_primary_weighted file drops the standard error, but it can be rebuilt
# from the n.risk / n.event columns, identical to how the original script does it
df_results_km_primary_weighted_se <-
  df_results_km_primary_weighted |>
  arrange(exp_bin_treat, time) |>
  group_by(exp_bin_treat) |>
  mutate(
    summand   = (1 / (n.risk - n.event)) - (1 / n.risk), # greenwood term
    surv      = cumprod(1 - n.event / n.risk),
    surv.se   = surv * sqrt(cumsum(summand)),
    cmlinc.se = surv.se # cmlinc.se = surv.se
  ) |>
  ungroup()

# per-group cumulative risks at end of follow-up
# CIs are taken straight from the file (complementary log-log limits)
risk_end <-
  df_results_km_primary_weighted_se |>
  filter(time == fup_end) |>
  select(exp_bin_treat, cmlinc, cmlinc.low, cmlinc.high, cmlinc.se)

r0 <- risk_end |> filter(exp_bin_treat == 0) # no antidiabetic (reference)
r1 <- risk_end |> filter(exp_bin_treat == 1) # metformin (exposed)

# contrasts: risk ratio and risk difference (delta method)
ci0 <- r0$cmlinc; se0 <- r0$cmlinc.se
ci1 <- r1$cmlinc; se1 <- r1$cmlinc.se

# risk ratio (log scale)
rr       <- ci1 / ci0
rr.ln.se <- sqrt((se1 / ci1)^2 + (se0 / ci0)^2)
rr.ll    <- exp(log(rr) + qnorm(0.025) * rr.ln.se)
rr.ul    <- exp(log(rr) + qnorm(0.975) * rr.ln.se)

# risk difference (linear scale)
rd    <- ci1 - ci0
rd.se <- sqrt(se0^2 + se1^2)
rd.ll <- rd + qnorm(0.025) * rd.se
rd.ul <- rd + qnorm(0.975) * rd.se

# format and print
pct <- function(x) sprintf("%.2f%%", x * 100)

results <- tibble(
  measure = c(
    "Cumulative risk, metformin (exp_bin_treat==1)",
    "Cumulative risk, no antidiabetic (exp_bin_treat==0)",
    "Risk ratio (metformin vs no antidiabetic)",
    "Risk difference (metformin - no antidiabetic)"
  ),
  estimate = c(pct(ci1), pct(ci0), sprintf("%.2f", rr), pct(rd)),
  ci_95 = c(
    sprintf("%s to %s", pct(r1$cmlinc.low), pct(r1$cmlinc.high)),
    sprintf("%s to %s", pct(r0$cmlinc.low), pct(r0$cmlinc.high)),
    sprintf("%.2f to %.2f", rr.ll, rr.ul),
    sprintf("%s to %s", pct(rd.ll), pct(rd.ul))
  )
)

cat("=== End of follow-up: time =", fup_end, "===\n\n")
print(results, n = Inf)

