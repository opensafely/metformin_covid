# # # # # # # # # # # # # # # # # # # # #
# Purpose: Build 1 combined plot from all disclosure-safe Kaplan-Meier estimates of antidiabetic treatment regimen patterns among T2DM patients
# # # # # # # # # # # # # # # # # # # # #

################################################################################
# Import libraries
################################################################################
library(tidyverse)
library(here)
library(arrow)
library(plotly)
library(ggplot2)

################################################################################
# Load the data
################################################################################
# Name of the datasets
file_names <- c("metfin", "metfin_mono", "dpp4_mono", "tzd_mono", 
                "sglt2_mono", "sulfo_mono", "glp1_mono", 
                "megli_mono", "agi_mono", "insulin_mono")

# Load them into a list and assign their names
data_list <- lapply(file_names, function(name) {
  read_csv(here("output", name, "km_estimates.csv"))
})
names(data_list) <- file_names

# Combine into long format and add an identifier (-> group)
combined_data <- bind_rows(
  lapply(seq_along(data_list), function(i) {
    data_list[[i]] %>% mutate(group = file_names[i])
  })
)

# Take out the duplicates with 0 risk
combined_data <- combined_data %>%
  dplyr::filter((exp_bin_metfin_anytime == 1) |
                  (exp_bin_metfin_mono_anytime == 1) |
                  (exp_bin_dpp4_mono_anytime == 1) |
                  (exp_bin_tzd_mono_anytime == 1) |
                  (exp_bin_sglt2_mono_anytime == 1) |
                  (exp_bin_sulfo_mono_anytime == 1) |
                  (exp_bin_glp1_mono_anytime == 1) |
                  (exp_bin_megli_mono_anytime == 1) |
                  (exp_bin_agi_mono_anytime == 1) |
                  (exp_bin_insulin_mono_anytime == 1))

################################################################################
# Create directories for output
################################################################################
fs::dir_create(here::here("output", "data_properties"))

################################################################################
# Plot
################################################################################
# # Using plotly
# # Create a custom color scale (rainbow)
# custom_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(combined_data$group)))
# # Column for months (assuming 30.44 days per month)
# combined_data$month <- combined_data$time / 30.44
# # Plotly without confidence intervals and no area filling
# cum_inc_plot_rounded1 <- plot_ly(data = combined_data, x = ~time, y = ~risk, color = ~group, type = 'scatter', mode = 'lines', colors = custom_colors) %>%
#   layout(
#     title = "Cumulative Incidence Curve by Group",
#     xaxis = list(
#       title = "Time (Days)",
#       tickvals = seq(0, max(combined_data$time), by = 30),  # monthly ticks
#       ticktext = seq(0, max(combined_data$time), by = 30)   # monthly ticks
#     ),
#     yaxis = list(title = "Cumulative Incidence"),
#     showlegend = TRUE
#   )

# Using ggplot, exactly the same as Will's plot in km.R
cum_inc_plot_rounded <- combined_data %>%
  group_modify(
    ~ add_row(
      .x,
      time = 0, # assumes time origin is zero
      lagtime = 0,
      surv = 1,
      surv.low = 1,
      surv.high = 1,
      risk = 0,
      risk.low = 0,
      risk.high = 0,
      .before = 0
    )
  ) %>%
  ggplot(aes(group = group, colour = group, fill = group)) +
  geom_step(aes(x = time, y = risk), direction = "vh") +
  geom_step(aes(x = time, y = risk), direction = "vh", linetype = "dashed", alpha = 0.5) +
  # geom_rect(aes(xmin = lagtime, xmax = time, ymin = risk.low, ymax = risk.high), alpha = 0.1, colour = "transparent") +
  scale_color_brewer(type = "qual", palette = "Set3", na.value = "grey") +
  scale_fill_brewer(type = "qual", palette = "Set3", guide = "none", na.value = "grey") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_x_continuous(breaks = seq(0, max(combined_data$time), by = 30)) +  # Add ticks every 30 days
  coord_cartesian(xlim = c(0, NA)) +
  labs(
    x = "Days since origin",
    y = "Kaplan-Meier estimate",
    colour = NULL,
    title = NULL
  ) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(colour = "black"),
    panel.grid.minor.x = element_blank(),
    legend.position = c(.05, .95),
    legend.justification = c(0, 1),
  )

cum_inc_plot_rounded_withCI <- combined_data %>%
  group_modify(
    ~ add_row(
      .x,
      time = 0, # assumes time origin is zero
      lagtime = 0,
      surv = 1,
      surv.low = 1,
      surv.high = 1,
      risk = 0,
      risk.low = 0,
      risk.high = 0,
      .before = 0
    )
  ) %>%
  ggplot(aes(group = group, colour = group, fill = group)) +
  geom_step(aes(x = time, y = risk), direction = "vh") +
  geom_step(aes(x = time, y = risk), direction = "vh", linetype = "dashed", alpha = 0.5) +
  geom_rect(aes(xmin = lagtime, xmax = time, ymin = risk.low, ymax = risk.high), alpha = 0.1, colour = "transparent") +
  scale_color_brewer(type = "qual", palette = "Set3", na.value = "grey") +
  scale_fill_brewer(type = "qual", palette = "Set3", guide = "none", na.value = "grey") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_x_continuous(breaks = seq(0, max(combined_data$time), by = 30)) +  # Add ticks every 30 days
  coord_cartesian(xlim = c(0, NA)) +
  labs(
    x = "Days since origin",
    y = "Kaplan-Meier estimate",
    colour = NULL,
    title = NULL
  ) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(colour = "black"),
    panel.grid.minor.x = element_blank(),
    legend.position = c(.05, .95),
    legend.justification = c(0, 1),
  )

################################################################################
# Save output
################################################################################
ggsave(filename = here::here("output", "data_properties", "cum_inc_plot_rounded.png"), cum_inc_plot_rounded, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "data_properties", "cum_inc_plot_rounded_withCI.png"), cum_inc_plot_rounded_withCI, width = 20, height = 20, units = "cm")
