# lunaineq

`lunaineq` estimates health inequality metrics for ecological or grouped data. The package is intended for analyses where each row represents an aggregate unit—such as a region, province, department, municipality, district, country, quintile, or other population group—and the analyst wants to quantify how a health indicator varies across a socioeconomic or equity stratifier.

The package implements six metrics:

| Function | Metric | Question answered |
|---|---|---|
| `ag_ineqeco()` | Absolute Gap (AG) | How large is the absolute difference between the most disadvantaged and most advantaged groups? |
| `rg_ineqeco()` | Relative Gap (RG) | How many times higher/lower is the indicator in the disadvantaged group compared with the advantaged group? |
| `aci_ineqeco()` | Absolute Concentration Index (ACI) | How much absolute inequality is distributed across the full socioeconomic ranking? |
| `rci_ineqeco()` | Relative Concentration Index (RCI) | Is the outcome concentrated among disadvantaged or advantaged units across the full distribution? |
| `sii_ineqeco()` | Slope Index of Inequality (SII) | What is the model-predicted absolute difference between the top and bottom of the social hierarchy? |
| `rii_ineqeco()` | Relative Index of Inequality (RII) | What is the model-predicted ratio between the top and bottom of the social hierarchy? |

## Installation

```r
# install.packages("devtools")
devtools::install_github("adrianvasquezmd/lunaineq", build_vignettes = TRUE)
```

## Core data requirements

Each analysis requires:

1. one row per ecological unit or group;
2. an equity stratifier supplied through `equity_stratifier_var`;
3. either a precomputed health indicator (`health_indicator_var`) or raw numerator/denominator columns (`health_numerator_var`, `health_denominator_var`);
4. a population or exposure weight, supplied explicitly through `population_weights_var` or implicitly through `health_denominator_var` when counts are used.

If the same ecological unit appears in multiple years, subset before analysis:

```r
df_2020 <- subset(paho_data, year == 2020)
```

## Direction of the equity stratifier

The argument `higher_ineq_is_favorable` controls how the equity stratifier is ordered.

Use `higher_ineq_is_favorable = TRUE` when higher values of the stratifier mean greater advantage, such as income, wealth, education, or HDI.

Use `higher_ineq_is_favorable = FALSE` when higher values mean greater disadvantage, such as poverty, deprivation, marginalization, vulnerability, or unsatisfied basic needs.

Example: if `equity_stratifier_var = "ubn"` and higher UBN means greater deprivation, use:

```r
higher_ineq_is_favorable = FALSE
```

## Health indicator input modes

### Count-based input

Preferred when numerator and denominator are available:

```r
health_indicator_type  = "rate"
health_numerator_var   = "maternal_death"
health_denominator_var = "live_births"
rate_scaling_factor    = 100000
```

### Precomputed indicator input

Use when the rate, percentage, proportion, or ratio is already calculated:

```r
health_indicator_var   = "mmr"
population_weights_var = "live_births"
health_indicator_type  = "rate"
```

If available, `health_se_var` can be supplied for AG/RG/ACI/RCI analytic variance workflows. The current SII/RII functions do not accept `health_se_var`; they estimate uncertainty using Wald, bootstrap, or jackknife methods.

## Quick examples

### Absolute gap

```r
ag <- ag_ineqeco(
  data = df_2020,
  health_indicator_type = "rate",
  health_numerator_var = "maternal_death",
  health_denominator_var = "live_births",
  rate_scaling_factor = 100000,
  analysis_unit_var = "state",
  equity_stratifier_var = "ubn",
  higher_ineq_is_favorable = FALSE,
  grouping_approach = "fractional",
  n_groups = 5,
  ci_method = "delta"
)

ag$summary_stats
ag$results_groups
```

Interpretation: if AG is 40 deaths per 100,000 live births, the maternal mortality ratio is 40 deaths per 100,000 live births higher in the most disadvantaged group than in the most advantaged group.

### Relative concentration index

```r
rci <- rci_ineqeco(
  data = df_2020,
  health_indicator_type = "rate",
  health_numerator_var = "maternal_death",
  health_denominator_var = "live_births",
  rate_scaling_factor = 100000,
  analysis_unit_var = "state",
  equity_stratifier_var = "ubn",
  higher_ineq_is_favorable = FALSE,
  bounded_corrections = "none",
  ci_method = "analytic"
)

rci$summary_stats
rci$results_curve
```

Interpretation: the RCI summarizes whether the indicator is concentrated toward the disadvantaged or advantaged end of the socioeconomic distribution. For adverse indicators, concentration among disadvantaged units indicates an inequitable distribution of health burden.

### Slope index of inequality

```r
sii <- sii_ineqeco(
  data = df_2020,
  health_indicator_type = "rate",
  health_numerator_var = "maternal_death",
  health_denominator_var = "live_births",
  rate_scaling_factor = 100000,
  analysis_unit_var = "state",
  equity_stratifier_var = "ubn",
  higher_ineq_is_favorable = FALSE,
  models = NULL,
  model_selection_metric = "RMSE",
  social_position_method = "classic",
  ci_method = "wald",
  vcov_type = "HC1"
)

sii$summary_stats
sii$model_selection_table
```

Interpretation: SII is the model-predicted absolute difference between the advantaged and disadvantaged ends of the social hierarchy. Inspect `model_selection_table` to verify which model was selected and why.

## Choosing among functions

Use AG/RG when you need a clear comparison between extreme groups. Use ACI/RCI when you want the full ranked distribution, not only extremes. Use SII/RII when you want a regression-based social-gradient estimate that predicts values at the bottom and top of the social hierarchy.

For reporting, use absolute and relative measures together. A relative gap can be large even when the absolute difference is small; an absolute gap can be large even when the relative ratio is modest.

## Main methodological choices

### Grouping for AG/RG

`grouping_approach` controls how the extreme groups are constructed:

- `"territorial"`: groups ecological units based on the distribution of units;
- `"weighted_cut"`: groups using population/exposure weighting;
- `"weighted_midpoint"`: groups using weighted midpoints in the ranked distribution;
- `"fractional"`: uses fractional-rank information to form the comparison.

When `grouping_approach = "territorial"`, `territorial_method` can be `"quantile"`, `"kmeans"`, `"fisher"`, `"manual"`, or `"equal"`.

### Bounded corrections for RCI

`bounded_corrections` controls whether Wagstaff and/or Erreygers corrections are returned for bounded outcomes. Use these mainly for proportions or percentages. For unbounded rates such as mortality rates, use `bounded_corrections = "none"` unless a meaningful lower and upper bound is explicitly justified.

### Model selection for SII/RII

SII/RII select among candidate models using `model_selection_metric`. Prediction-error metrics (`MSE`, `RMSE`, `MAE`, `MAPE`) are computed on the final outcome scale and weighted by population/exposure weights. `AIC` and `BIC` are used only for likelihood-comparable model sets.

For most applied workflows, start with:

```r
models = NULL
model_selection_metric = "RMSE"
social_position_method = "classic"
ci_method = "wald"
```

Then inspect diagnostics and consider bootstrap/jackknife as sensitivity analyses.

## Vignette

```r
vignette("lunaineq", package = "lunaineq")
```

The vignette provides function-by-function documentation, argument-level guidance, and epidemiological interpretation templates.

## Authors and license

Authors: Adrian Vasquez-Mejia, Oscar J. Mujica, Antonio Sanhueza.

License: AGPL (>= 3).
