#' Output Structuring and Formatting Engine
#'
#' @description
#' `ag_rg_format_outputs` rigorously assembles and standardizes the computational results
#' from the AG and RG estimation pipelines into a comprehensive, multi-table `list` object.
#' It maps internal programmatic metadata to human-readable methodological descriptors,
#' compiles analytical and resampling warning arrays, and formats group-level
#' variance metrics for public-facing analysis.
#'
#' @noRd
#' @param disadvantaged_value Numeric.
#' @param advantaged_value Numeric.
#' @param grouped_df Data frame.
#' @param resampling_results Data frame.
#' @param n_groups Integer.
#' @param n_original Integer.
#' @param grouping_approach Character.
#' @param territorial_method Character.
#' @param scenario Character.
#' @param use_counts Logical.
#' @param variance_family Character.
#' @param higher_ineq_is_favorable Logical.
#' @param final_multiplier Numeric.
#' @param log_message Character vector.
#' @param log_warning Character vector.
#' @param results_groups Data frame returned from `.ag_rg_group_summary`.
#' @param methodological_notes_specific Character vector for gap-specific notes.
#' @return An S3 object with the Luna outputs.
#' @author Adrián Vásquez
#' @keywords internal
.ag_rg_format_outputs <- function(metric_name, point_gap, se_out, ci_low, ci_high, conf_level, 
                                    ci_method_used, interval_interpretation, lower_group, upper_group, 
                                    disadvantaged_value, advantaged_value, grouped_df, resampling_results, 
                                    n_groups, n_original, grouping_approach, territorial_method, scenario, 
                                    use_counts, variance_family, higher_ineq_is_favorable, final_multiplier, 
                                    log_message, log_warning, results_groups, methodological_notes_specific = NULL) {

  interval_type <- ifelse(
    grepl("^inferential", interval_interpretation),
    "inferential",
    ifelse(grepl("^ecological_sensitivity", interval_interpretation), "ecological_sensitivity", "descriptive_only")
  )

  grouping_label <- ifelse(grouping_approach == "territorial", paste0("territorial_", territorial_method), grouping_approach)

  inputs_used <- switch(
    scenario,
    counts_no_se = "numerator + denominator + equity stratifier",
    counts_with_se = "numerator + denominator + SE + equity stratifier",
    indicator_population_no_se = "health indicator + population weight + equity stratifier",
    indicator_population_with_se = "health indicator + population weight + SE + equity stratifier"
  )

  point_estimator_detail <- switch(
    scenario,
    counts_no_se = "Group indicators are estimated as sum(numerator) / sum(denominator).",
    counts_with_se = "Group indicators are estimated as sum(numerator) / sum(denominator); uncertainty uses the user-supplied ecological-unit SE.",
    indicator_population_no_se = "Group indicators are estimated as population-weighted ecological averages.",
    indicator_population_with_se = "Group indicators are estimated as population-weighted ecological averages; uncertainty uses the user-supplied ecological-unit SE."
  )

  count_variance_family_used <- ifelse(
    scenario == "counts_no_se" && ci_method_used %in% c("delta", "bootstrap"),
    variance_family,
    NA_character_
  )

  uncertainty_method <- ifelse(
    ci_method_used == "delta" && scenario == "counts_no_se",
    paste0("delta_", variance_family, "_from_counts"),
    ifelse(
      ci_method_used == "bootstrap" && scenario == "counts_no_se",
      paste0("parametric_", variance_family, "_bootstrap_from_counts_boot_package"),
      ifelse(
        ci_method_used == "delta" && scenario %in% c("counts_with_se", "indicator_population_with_se"),
        "delta_from_user_supplied_se",
        ifelse(
          ci_method_used %in% c("bootstrap", "jackknife"),
          ifelse(ci_method_used == "bootstrap", "bootstrap_ecological_sensitivity_boot_package", "jackknife_ecological_sensitivity_bootstrap_package"),
          "not_available"
        )
      )
    )
  )

  summary_stats <- data.frame(
    Function = metric_name,
    Scenario = scenario,
    Inputs_Used = inputs_used,
    Point_Estimator = point_estimator_detail,
    CI_Method = ci_method_used,
    Uncertainty_Method = uncertainty_method,
    Interval_Type = interval_type,
    Grouping = grouping_label,
    Number_of_Groups = n_groups,
    Lower_Group = lower_group,
    Upper_Group = upper_group,
    Units_Original = n_original,
    Units_Analyzed = length(unique(grouped_df$unit_id)),
    Confidence_Level = conf_level,
    Scale_Multiplier = final_multiplier,
    stringsAsFactors = FALSE
  )

  results_overall <- data.frame(
    Metric = metric_name,
    Estimate = as.numeric(point_gap),
    Standard_Error = as.numeric(se_out),
    CI_Lower = as.numeric(ci_low),
    CI_Upper = as.numeric(ci_high),
    Confidence_Level = conf_level,
    CI_Method = ci_method_used,
    Interval_Type = interval_type,
    Disadvantaged_Group = lower_group,
    Advantaged_Group = upper_group,
    Disadvantaged_Value = as.numeric(disadvantaged_value),
    Advantaged_Value = as.numeric(advantaged_value),
    stringsAsFactors = FALSE
  )

  # Keep the group output concise and stable across ag_ineqeco() and rg_ineqeco().
  # .ag_rg_group_summary() returns additional technical columns used internally,
  # but the user-facing table should expose only the core group summary fields.
  group_output_cols <- c(
    "group", "group_role", "n_units_contributing", "population_weight",
    "population_share", "numerator", "denominator", "health_estimate",
    "standard_error"
  )

  for (nm in setdiff(group_output_cols, names(results_groups))) {
    results_groups[[nm]] <- NA_real_
  }

  results_groups <- results_groups[, group_output_cols, drop = FALSE]

  names(results_groups) <- c(
    "Group", "Group_Role", "Units_N", "Population_Weight",
    "Population_Share", "Numerator", "Denominator", "Health_Estimate",
    "Standard_Error"
  )

  unit_cols <- c(
    "unit_id", "group", "allocation_fraction", "allocated_population",
    "equity_stratifier", "health_indicator", "health_se",
    "numerator", "denominator"
  )
  for (nm in setdiff(unit_cols, names(grouped_df))) {
    grouped_df[[nm]] <- NA_real_
  }

  results_units <- grouped_df[, unit_cols, drop = FALSE]
  names(results_units) <- c(
    "Unit_ID", "Group", "Allocation_Fraction", "Allocated_Population",
    "Equity_Stratifier", "Health_Indicator", "Health_SE",
    "Numerator", "Denominator"
  )

  if (nrow(resampling_results) > 0) {
    names(resampling_results) <- c("Replicate", "Method", "Omitted_Unit", "Estimate")
  }

  direction_label <- ifelse(higher_ineq_is_favorable,
                            "higher_stratifier_values_are_more_favorable",
                            "higher_stratifier_values_are_less_favorable")

  methodology <- data.frame(
    Item = c(
      "Scenario",
      "Inputs used",
      "Point estimator",
      "Uncertainty method",
      "Count variance family",
      "Interval interpretation",
      "Resampling engine",
      "Precision policy"
    ),
    Description = c(
      scenario,
      inputs_used,
      point_estimator_detail,
      uncertainty_method,
      ifelse(is.na(count_variance_family_used), "not used", count_variance_family_used),
      interval_interpretation,
      ifelse(ci_method_used == "bootstrap", "boot::boot() with percentile limits", ifelse(ci_method_used == "jackknife", "bootstrap::jackknife() delete-one ecological unit", "not used")),
      "No rounding is applied inside the function. Numeric values are returned as full double-precision numeric values; any apparent shortening is only an R console printing setting."
    ),
    stringsAsFactors = FALSE
  )

  grouping_methodology <- data.frame(
    Item = c(
      "Grouping approach",
      "Territorial method",
      "Number of groups",
      "Stratifier direction",
      "Lower group",
      "Upper group"
    ),
    Description = c(
      grouping_approach,
      ifelse(grouping_approach == "territorial", territorial_method, "not applicable"),
      as.character(n_groups),
      direction_label,
      as.character(lower_group),
      as.character(upper_group)
    ),
    stringsAsFactors = FALSE
  )

  methodological_notes <- c(
    "Ecological data are analyzed as aggregate units; no individual-level survey design variance is estimated.",
    ifelse(scenario == "counts_no_se", paste0("Delta uncertainty uses an analytical ", variance_family, " approximation from the supplied counts."), NA_character_),
    ifelse(scenario %in% c("counts_with_se", "indicator_population_with_se"), "Delta uncertainty uses the user-supplied ecological-unit standard errors.", NA_character_),
    ifelse(scenario == "indicator_population_no_se", "No inferential interval is available unless valid counts or ecological-unit standard errors are supplied.", NA_character_),
    ifelse(ci_method_used == "bootstrap", "Bootstrap resampling is performed with boot::boot(); reported limits are percentile limits from the bootstrap replicate distribution.", NA_character_),
    ifelse(ci_method_used == "jackknife", "Jackknife resampling is performed with bootstrap::jackknife() as a delete-one ecological-unit procedure.", NA_character_),
    ifelse(ci_method_used %in% c("bootstrap", "jackknife") && scenario != "counts_no_se", "Resampling intervals are ecological sensitivity diagnostics, not formal sampling confidence intervals.", NA_character_),
    ifelse(grouping_approach == "fractional", "Fractional grouping allocates ecological units across groups and assumes homogeneous health-indicator values within split units.", NA_character_),
    "All numeric outputs are stored as numeric values without rounding. Returned result tables are base data frames to avoid tibble/pillar display truncation; printing is still controlled by the user's R option `digits`.",
    methodological_notes_specific
  )
  
  methodological_notes <- methodological_notes[!is.na(methodological_notes)]

  diagnostics <- list(
    messages = log_message,
    warnings = log_warning,
    methodological_notes = methodological_notes
  )

  out <- list(
    summary_stats = summary_stats,
    diagnostics = diagnostics,
    results_overall = results_overall,
    results_groups = results_groups,
    results_units = results_units,
    results_replicates = resampling_results,
    methodology = methodology,
    grouping_methodology = grouping_methodology
  )

  return(out)
}
