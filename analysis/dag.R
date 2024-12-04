# # # # # # # # # # # # # # # # # # # # #
# This script creates directed acyclic graphs (DAGs) used for the study
# # # # # # # # # # # # # # # # # # # # #

# Import libraries ----
library(dagitty)
library(ggdag)
library(ggplot2)
library(dplyr)

# Create the DAG ----
dag <- dagify(covidoutcome ~ metformin + agesexsmoking + ethnicity + region + cvdcomorb + othercomorb + DMseverity + hsb,
              metformin ~ agesexsmoking + ethnicity + region + cvdcomorb + hsb,
              otherDMdrugs ~ DMseverity + region,
              DMseverity ~ agesexsmoking,
              cvdcomorb ~ hsb + agesexsmoking,
              othercomorb ~ hsb + agesexsmoking,
              region ~ imd,
              ethnicity ~ imd,
              hsb  ~ ethnicity,
              othercomorb ~ ethnicity + region + agesexsmoking,
              labels = c("covidoutcome" = "COVID outcomes",
                         "metformin" = "Metformin",
                         "agesexsmoking" = "Age, Sex, Smoking",
                         "region" = "Region, Practice,\n Rural/urban",
                         "ethnicity" = "Ethnicity",
                         "cvdcomorb" = "CVDs other\n than T2DM*",
                         "othercomorb" = "Non-CVD comorbidities**",
                         "DMseverity" = "Diabetes severity****",
                         "otherDMdrugs" = "Other antidiabetic treatment",
                         "imd" = "Deprivation",
                         "hsb" = "Health seeking behaviour***"
                         ),
              exposure = "metformin",
              outcome = "covidoutcome"
)
ggdag(dag, text = F, use_labels = "label") +
  theme_void() +
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank()
  )
ggdag_adjustment_set(dag, text = F, use_labels = "label") +
  theme_void() +
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank()
  )
