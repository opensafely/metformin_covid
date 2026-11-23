
################################################################################
# Helper functions to apply exclusion criteria SEQUENTIALLY and record a flowchart
#
# Differences to the landmark functions (fn_quality_assurance_midpoint6.R,
# fn_completeness_criteria_midpoint6.R):
# - Criteria are applied one after the other: each person is counted once, at the
#   first criterion they fail. The excluded counts therefore add up to the difference
#   between the start and end counts (in the landmark functions, criteria were counted
#   independently, so one person could be counted under several criteria).
# - Each criterion must evaluate to TRUE/FALSE for every person; a missing value
#   (NA) stops the script, instead of silently keeping the person.
# - Disclosure control: only the remaining counts are midpoint-6 rounded; the excluded
#   counts are the differences between consecutive rounded remaining counts, so the
#   released flowchart adds up.
################################################################################

# Apply exclusion criteria sequentially
# data:       data frame (one row per person)
# exclusions: named list of expressions, created with rlang::exprs(), each evaluated
#             within `data`; TRUE = exclude. The names are the flowchart labels.
# group:      label for this block of criteria in the flowchart
# Returns a list with the remaining data and the flowchart rows of this block
fn_apply_exclusions <- function(data, exclusions, group) {
  flowchart <- tibble(group = character(), criterion = character(),
                      n_excluded = integer(), n_remaining = integer())
  for (label in names(exclusions)) {
    exclude <- rlang::eval_tidy(exclusions[[label]], data = data)
    if (!is.logical(exclude) || length(exclude) != nrow(data)) {
      stop("Exclusion criterion '", label, "' must return one TRUE/FALSE per person.")
    }
    if (anyNA(exclude)) {
      stop("Exclusion criterion '", label, "' is missing (NA) for ", sum(is.na(exclude)),
           " people. Define how missing values are handled within the criterion.")
    }
    data <- data[!exclude, , drop = FALSE]
    flowchart <- add_row(flowchart, group = group, criterion = label,
                         n_excluded = sum(exclude), n_remaining = nrow(data))
  }
  list(data = data, flowchart = flowchart)
}

# Add the midpoint-6 rounded counts to a flowchart (rows in sequential order; the first
# row must be the starting population, with n_excluded = NA)
fn_round_flowchart <- function(flowchart, threshold) {
  flowchart %>%
    mutate(
      n_remaining_midpoint6 = fn_roundmid_any(n_remaining, threshold),
      n_excluded_midpoint6 = lag(n_remaining_midpoint6) - n_remaining_midpoint6
    )
}
