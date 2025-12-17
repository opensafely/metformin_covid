####
## This script does the following:
# 1. Import released midpoint rounded disclosure-safe output data
# 2. Create various tables and figures 
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


# Import the data -------------------------------------------------
df_elig <- read_csv(here("output", "data_release_20251215", "data_description", "n_elig_excluded_midpoint6.csv"))
df_qa <- read_csv(here("output", "data_release_20251215", "data_description", "n_qa_excluded_midpoint6.csv"))
df_comp <- read_csv(here("output", "data_release_20251215", "data_description", "n_completeness_excluded_midpoint6.csv"))
df_elig_landmark <- read_csv(here("output", "data_release_20251215", "data_description", "n_restricted_midpoint6.csv"))

df_main <- read_csv(here("output", "data_release_20251215", "data_description", "table1_main_midpoint6.csv"))
df_ps <- read_csv(here("output", "data_release_20251215", "ps", "density_plot_untrimmed.csv"))
df_ps_HbA1c59orabove <- read_csv(here("output", "data_release_20251215", "ps", "ps_density_data_untrimmed_HbA1c59orabove.csv"))
df_ps_HbA1c42to58 <- read_csv(here("output", "data_release_20251215", "ps", "ps_density_data_untrimmed_HbA1c42to58.csv"))
df_ps_HbA1cbelow42 <- read_csv(here("output", "data_release_20251215", "ps", "ps_density_data_untrimmed_HbA1cbelow42.csv"))
df_results_cox <- read_csv(here("output", "data_release_20251215", "make_output", "results_cox_midpoint6.csv"))
df_results_km_primary <- read_csv(here("output", "data_release_20251215", "te", "km_primary", "estimates.csv"))
df_results_km_covid_event <- read_csv(here("output", "data_release_20251215", "te", "km_covid_event", "estimates.csv"))
df_results_km_longvirfat <- read_csv(here("output", "data_release_20251215", "te", "km_longvirfat", "estimates.csv"))


# Baseline table, based on output from table1.R ----------------------

# Hide NAs where not needed
hide_NA <- function(df) {
  df %>%
    mutate(`metformin` = replace_na(as.character(`metformin`), "")) %>% 
    mutate(control = replace_na(as.character(control), ""))
}
df_main <- hide_NA(df_main)

# Reduce to main variables
# df_main <- df_main %>% 
#   filter(!grepl("INT: Any", var_label)) %>% 
#   filter(!grepl("Any LTFU", var_label)) %>% 
#   filter(!grepl("Any death", var_label)) %>% 
#   filter(!grepl("Any Viral Fatigue", var_label)) %>% 
#   filter(!grepl("Any Long COVID diagnosis", var_label))

# Create the gt table
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
gtsave(tbl_gt_main, "baseline_characteristics_main.rtf")


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
                               "covid_event" = "Secondary: Any COVID-19-related event",
                               "longvirfat" = "Secondary: Long COVID",
                               "neg_control_pandemic" = "Negative control outcome (bone fracture)",
                               "pos_control_pandemic" = "Positive control outcome (diabetes-related death)",
)
desired_order <- c(
  "Primary: COVID-19-related hospitalization or death",
  "Secondary: Any COVID-19-related event",
  "Secondary: Long COVID",
  "Negative control outcome (bone fracture)",
  "Positive control outcome (diabetes-related death)"
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
pdf(file = "Fp_results_main.pdf",
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
                      "7" = gpar(lty = 2)),
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
  "Higher risk in No Antidiabetic group < > Higher risk in Metformin group",
  x = unit(0.695, "npc"),
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
                              "covid_event" = "Secondary: Any COVID-19-related event",
                              "longvirfat" = "Secondary: Long COVID",
                              "neg_control_pandemic" = "Negative control outcome (bone fracture)",
                              "pos_control_pandemic" = "Positive control outcome (diabetes-related death)",
)
desired_order <- c(
  "Primary: COVID-19-related hospitalization or death",
  "Secondary: Any COVID-19-related event",
  "Secondary: Long COVID",
  "Negative control outcome (bone fracture)",
  "Positive control outcome (diabetes-related death)"
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
pdf(file = "Fp_results_min.pdf",
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
                      "7" = gpar(lty = 2)),
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
  "Higher risk in No Antidiabetic group < > Higher risk in Metformin group",
  x = unit(0.695, "npc"),
  y = unit(0.015, "npc"),
  gp = gpar(cex = 0.9)
)

dev.off()


## Subgroup analyses ----
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

# build forestplot
base_data <- tibble(
  mean = df_results_sg$hr,
  lower = df_results_sg$conf_low,
  upper = df_results_sg$conf_high,
  name = as.character(df_results_sg$name),
  events_i = as.character(df_results_sg$events_i),
  events_c = as.character(df_results_sg$events_c),
  estimates = paste0(formatC(df_results_sg$hr, format = "f", digits = 2), 
                     " (", formatC(df_results_sg$conf_low, format = "f", digits = 2), 
                     " - ", formatC(df_results_sg$conf_high, format = "f", digits = 2), ")"))
header <- tibble(
  name = "Subgroup on primary outcome",
  events_i = "Events Metformin\n n (%)",
  events_c = "Events No Antidiabetic\n n (%)",
  estimates = "aHR (95% CI)",
  mean = NA, lower = NA, upper = NA)

fp <- bind_rows(header, base_data)
font <- "sans"

# Create the forest plot
pdf(file = "Fp_results_sg.pdf",
    width = 14,
    height = 12)

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
                      "8" = gpar(lty = 2),
                      "9" = gpar(lty = 2),
                      "10" = gpar(lty = 2),
                      "11" = gpar(lty = 2),
                      "12" = gpar(lty = 2),
                      "13" = gpar(lty = 2)),
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
  "Higher risk in No Antidiabetic group < > Higher risk in Metformin group",
  x = unit(0.617, "npc"),
  y = unit(0.008, "npc"),
  gp = gpar(cex = 0.9)
)

dev.off()


# Quality assurance flowchart -------------------------------------

# # Create node labels
# node_labels <- paste(df_qa$Variable, '\n', df_qa$Value)
# 
# # Create flowchart
# graph_qa <- grViz(
#   sprintf(
#     "digraph flowchart {
#       graph [layout = dot, rankdir = TB]
#       node [shape = box, style = filled, fontname = Arial, fontsize = 14]
# 
#       A [label = '%s', fillcolor = lightskyblue]
#       B [label = 'Exclusions\n%s\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
#       C [label = '%s', fillcolor = lightskyblue]
# 
#       A -> B
#       B -> C
#     }",
#     node_labels[1], node_labels[2],
#     node_labels[3], node_labels[4], node_labels[5], node_labels[6], node_labels[7],
#     node_labels[8], node_labels[9]
#   )
# )

# Completeness criteria flowchart -----------------------------------

# # Create node labels
# node_labels <- paste(df_comp$Variable, '\n', df_comp$Value)
# 
# # Create flowchart
# graph_completeness <- grViz(
#   sprintf(
#     "digraph flowchart {
#       graph [layout = dot, rankdir = TB]
#       node [shape = box, style = filled, fontname = Arial, fontsize = 14]
# 
#       A [label = '%s', fillcolor = lightskyblue]
#       B [label = 'Exclusions\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
#       C [label = '%s', fillcolor = lightskyblue]
# 
#       A -> B
#       B -> C
#     }",
#     node_labels[1], node_labels[2],
#     node_labels[3], node_labels[4], node_labels[5], node_labels[6], node_labels[7],
#     node_labels[8], node_labels[9]
#   )
# )


# Eligibility flowchart -----------------------------------

# # Create node labels
# node_labels <- paste(df_elig$Variable, '\n', df_elig$Value)
# 
# # Create flowchart
# graph_elig <- grViz(
#   sprintf(
#     "digraph flowchart {
#       graph [layout = dot, rankdir = TB]
#       node [shape = box, style = filled, fontname = Arial, fontsize = 14]
# 
#       A [label = '%s', fillcolor = lightskyblue]
#       B [label = '%s', fillcolor = lightyellow]
#       C [label = 'Exclusions 1\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
#       D [label = 'Exclusions 2: Events between baseline and landmark\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
#       E [label = '%s', fillcolor = lightskyblue]
# 
#       A -> B
#       B -> C
#       C -> D
#       D -> E
#     }",
#     node_labels[1], node_labels[2], 
#     node_labels[3], node_labels[4], node_labels[5], node_labels[6], node_labels[7],
#     node_labels[8], node_labels[9], node_labels[10], node_labels[11], node_labels[12], 
#     node_labels[13], node_labels[14], node_labels[15], 
#     node_labels[16], node_labels[17], node_labels[18], node_labels[19], node_labels[20], node_labels[21],
#     node_labels[22]
#   )
# )