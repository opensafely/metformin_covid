
# # # # # # # # # # # # # # # # # # # # #
# Purpose: Build 1 combined plot from all disclosure-safe Kaplan-Meier estimates of antidiabetic treatment regimen patterns among T2DM patients
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library(tidyverse)
library(here)
library(arrow)
library(plotly)
library(ggplot2)

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
cum_inc_plot_rounded1 <- combined_data %>%
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

# without area under the curves and no confidence intervals
cum_inc_plot_rounded2 <- combined_data %>%
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
  ggplot(aes(group = group, colour = group)) +  # Only 'colour' for lines, no fill
  # Kaplan-Meier survival estimates (main line)
  geom_step(aes(x = time, y = risk), direction = "vh", size = 1) +  # Main line (Kaplan-Meier estimate)
  geom_step(aes(x = time, y = risk), direction = "vh", linetype = "dashed", alpha = 0.5, size = 1) +  # Dashed lines for risk
  scale_color_brewer(type = "qual", palette = "Set3", na.value = "grey") +
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
    panel.grid.major = element_blank(),  # Remove major gridlines (if any)
    panel.background = element_blank(),  # Remove any background that might cause fill
    plot.background = element_blank(),   # Ensure plot background is also clear
    legend.position = c(.05, .95),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.background = element_blank(),  # Remove any legend background
    plot.margin = margin(0, 0, 0, 0)  # Remove margins that may add unexpected space
  )

################################################################################
# Save output
################################################################################
ggsave(filename = here::here("output", "data_properties", "cum_inc_plot_rounded1.png"), cum_inc_plot_rounded1, width = 20, height = 20, units = "cm")
ggsave(filename = here::here("output", "data_properties", "cum_inc_plot_rounded2.png"), cum_inc_plot_rounded2, width = 20, height = 20, units = "cm")
