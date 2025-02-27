# # # # # # # # # # # # # # # # # # # # #
# Purpose: Build 1 combined plot from all disclosure-safe Kaplan-Meier estimates of antidiabetic treatment regimen patterns among T2DM patients
# # # # # # # # # # # # # # # # # # # # #

################################################################################
# Import libraries
################################################################################
library(tidyverse)
library(here)
library(arrow)
library(ggplot2)

################################################################################
# Load the data
################################################################################
# Name of the datasets
file_names <- c("metfin", "metfin_mono", "insulin_mono", "sglt2_mono"
                , "sulfo_mono", "dpp4_mono", "glp1_mono"
                , "megli_mono", "agi_mono", "tzd_mono"
                )

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

################################################################################
# Create directories for output
################################################################################
fs::dir_create(here::here("output", "data_description"))

################################################################################
# Plot
################################################################################
# Using ggplot, slightly adapt from the reusable action km.R plot: https://github.com/opensafely-actions/kaplan-meier-function/blob/main/analysis/km.R
cum_inc_plot_rounded <- combined_data %>%
  mutate(
    lagtime = lag(time, 1, 0), # assumes the time-origin is zero
  ) %>%
  group_modify(
    ~ add_row(
      .x,
      time = 0,  # Ensure the origin is included
      lagtime = 0,
      cmlinc = 0,
      cmlinc.low = 0,
      cmlinc.high = 0,
      .before = 0
    )
  ) %>%
  ggplot(aes(group = group, colour = group, fill = group)) +
  geom_step(aes(x = time, y = cmlinc), direction = "vh") +
  geom_step(aes(x = time, y = cmlinc), direction = "vh", linetype = "dashed", alpha = 0.5) +
  # geom_rect(aes(xmin = lagtime, xmax = time, ymin = cmlinc.low, ymax = cmlinc.high), alpha = 0.1, colour = "transparent") + # Confidence interval
  scale_color_brewer(type = "qual", palette = "Set3", na.value = "grey") +
  scale_fill_brewer(type = "qual", palette = "Set3", guide = "none", na.value = "grey") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_x_continuous(breaks = seq(0, max(combined_data$time, na.rm = TRUE), by = 30)) +  # Add ticks every 30 days
  coord_cartesian(xlim = c(0, NA)) +
  labs(
    x = "Days since T2DM diagnosis",
    y = "Cumulative Incidence",
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
  mutate(
    lagtime = lag(time, 1, 0), # assumes the time-origin is zero
  ) %>%
  group_modify(
    ~ add_row(
      .x,
      time = 0,  # Ensure the origin is included
      lagtime = 0,
      cmlinc = 0,
      cmlinc.low = 0,
      cmlinc.high = 0,
      .before = 0
    )
  ) %>%
  ggplot(aes(group = group, colour = group, fill = group)) +
  geom_step(aes(x = time, y = cmlinc), direction = "vh") +
  geom_step(aes(x = time, y = cmlinc), direction = "vh", linetype = "dashed", alpha = 0.5) +
  geom_rect(aes(xmin = lagtime, xmax = time, ymin = cmlinc.low, ymax = cmlinc.high), alpha = 0.1, colour = "transparent") + # Confidence interval
  scale_color_brewer(type = "qual", palette = "Set3", na.value = "grey") +
  scale_fill_brewer(type = "qual", palette = "Set3", guide = "none", na.value = "grey") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_x_continuous(breaks = seq(0, max(combined_data$time, na.rm = TRUE), by = 30)) +  # Add ticks every 30 days
  coord_cartesian(xlim = c(0, NA)) +
  labs(
    x = "Days since T2DM diagnosis",
    y = "Cumulative Incidence",
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
ggsave(filename = here::here("output", "data_description", "cum_inc_plot_midpoint6.png"), cum_inc_plot_rounded, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "data_description", "cum_inc_plot_withCI_midpoint6.png"), cum_inc_plot_rounded_withCI, width = 20, height = 20, units = "cm")
