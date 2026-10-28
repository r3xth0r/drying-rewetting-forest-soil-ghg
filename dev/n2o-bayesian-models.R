# Bayesian N2O GAM models ----
# Fits full and reduced models, evaluates thresholds, and fits depth-specific GAMs.
# Inputs: Klausen N2O, soil moisture, soil temperature, precipitation, irrigation.
# Outputs: Tables 5, 6, and S8; fitted models, curves, and numeric diagnostics.
# Figures are made by n2o-bayesian-outputs.R from the saved artifacts.

suppressPackageStartupMessages({
  library("tidyverse")
  library("brms")
  library("loo")
  library("posterior")
  library("bayesplot")
  library("mgcv")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("n2o-bayesian-models")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))

# Settings ====

n2o_file <- "n2o-klausen.csv"
sm_file <- "soil-moisture-klausen.csv"
st_file <- "soil-temperature-klausen.csv"
precipitation_file <- "precipitation-klausen.csv"
irrigation_file <- "irrigation-klausen.csv"

model_seed <- 1234L
mcmc_chains <- 4L
mcmc_iterations <- 2000L
mcmc_warmup <- 1000L
mcmc_cores <- getOption("mc.cores")
mcmc_adapt_delta <- 0.95
mcmc_max_treedepth <- 12L

lag_days <- 7L
precipitation_lag_days <- 3L
full_curve_points <- 400L
reduced_curve_points <- 300L
depth_gam_curve_points <- 300L
depth_gam_basis_k <- 6L
depth_gam_sensitivity_k <- 10L
gam_check_replicates <- 1000L
threshold_inner_fraction <- 0.8
threshold_min_segment_fraction <- 0.1
threshold_min_curve_points <- 10L

model_depths <- c("5", "15", "30")
model_treatments <- c("NP", "RE")
model_depth_labels <- c("5" = "0-10", "15" = "10-20", "30" = "20-40")
model_plot_depth_labels <- c(
  "0-10" = "0–10 cm",
  "10-20" = "10–20 cm",
  "20-40" = "20–40 cm"
)

n2o_interim_dir <- file.path(dir_interim, "n2o")
n2o_cache_dir <- file.path(n2o_interim_dir, "cache")
diagnostics_dir <- file.path(n2o_interim_dir, "diagnostics")

model_data_file <- file.path(n2o_interim_dir, "n2o-model-data.rds")
full_model_cache <- file.path(n2o_cache_dir, "n2o-full-model")
reduced_model_cache <- file.path(n2o_cache_dir, "n2o-reduced-model")
cache_model_names <- c(
  "n2o-full-no-smooth",
  "n2o-full-no-wfps-lags",
  "n2o-full-no-temperature-lags",
  "n2o-full-no-precipitation-lags",
  "n2o-full-no-rewetting-lags"
)

file_table_05 <- "table-05-n2o-loo-comparison.csv"
file_full_reduced_comparison <- "n2o-full-reduced-model-comparison.csv"
file_depth_wfps_table <- "table-06-n2o-depth-wfps-gams.csv"
file_threshold_results_table <- "n2o-threshold-results.csv"

file_full_wfps_curves <- "n2o-full-partial-curves-wfps.csv"
file_full_temp_curves <- "n2o-full-partial-curves-temperature.csv"
file_reduced_wfps_curves <- "n2o-reduced-partial-curves-wfps.csv"
file_reduced_temp_curves <- "n2o-reduced-partial-curves-temperature.csv"
file_threshold_results <- "n2o-threshold-results.csv"
file_treatment_differences <- "n2o-threshold-treatment-differences.csv"
file_depth_shifts <- "n2o-threshold-depth-shifts.csv"
file_depth_gam_models <- file.path(n2o_interim_dir, "n2o-depth-wfps-gams.rds")
file_depth_gam_data <- file.path(n2o_interim_dir, "n2o-depth-wfps-data.rds")
file_full_curves_rds <- file.path(n2o_interim_dir, "n2o-full-curves.rds")
file_reduced_curves_rds <- file.path(
  n2o_interim_dir, "n2o-reduced-curves.rds"
)
file_full_loo_rds <- file.path(n2o_interim_dir, "n2o-full-loo.rds")
file_reduced_loo_rds <- file.path(n2o_interim_dir, "n2o-reduced-loo.rds")
file_full_model_diagnostics <- "table-s08-n2o-full-model-diagnostics.csv"
file_reduced_model_diagnostics <- "n2o-reduced-model-diagnostics.csv"
file_full_reduced_loo_csv <- "n2o-full-vs-reduced-loo-comparison.csv"
file_full_reduced_loo_text <- "n2o-full-vs-reduced-loo-comparison.txt"
file_sampler_comparison <- "n2o-sampler-comparison.csv"
file_model_notes <- "n2o-model-notes.txt"
file_depth_wfps_results <- "n2o-depth-wfps-gam-results.csv"
file_depth_wfps_predictions <- "n2o-depth-wfps-gam-predictions.csv"
file_depth_wfps_basis_checks <- "n2o-depth-wfps-basis-checks.csv"
file_depth_wfps_k_checks <- "n2o-depth-wfps-k-checks.csv"
file_diagnostic_parameter_suffix <- "parameter-diagnostics.csv"
file_diagnostic_sampler_suffix <- "sampler-diagnostics.csv"
file_diagnostic_loo_suffix <- "loo-summary.txt"

full_diagnostic_model_name <- "n2o-full-gam"
reduced_diagnostic_model_name <- "n2o-reduced-gam"
full_diagnostic_treedepth_limit <- 10L
reduced_diagnostic_treedepth_limit <- 12L
depth_gam_diagnostic_seed <- model_seed
depth_gam_min_observations <- 10L

ensure_dir(n2o_interim_dir)
ensure_dir(n2o_cache_dir)
ensure_dir(diagnostics_dir)
ensure_dir(dir_processed)
ensure_dir(dir_reporting)
ensure_dir(dir_plots)

set.seed(model_seed)

# Helpers ====
.zscale_safe <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(ifelse(is.na(x), NA_real_, 0))
  }
  (x - m) / s
}

depth_range_label <- function(depth_vec) {
  d <- as.character(depth_vec)

  unname(model_depth_labels[d])
}

time_category_n2o <- function(time) {
  # Preserve the original fallback, including midnight and missing times.
  case_when(
    time > "00:00:00" & time <= "06:00:00" ~ "AM1",
    time > "06:00:00" & time <= "12:00:00" ~ "AM2",
    time > "12:00:00" & time <= "18:00:00" ~ "PM1",
    TRUE ~ "PM2"
  )
}

get_looic <- function(x) {
  unname(x$estimates["looic", "Estimate"])
}

get_looic_se <- function(x) {
  unname(x$estimates["looic", "SE"])
}

get_r2_mean <- function(fit) {
  mean(bayes_R2(fit, summary = FALSE))
}

# Precipitation and irrigation ====

precipitation_daily <- read_input(precipitation_file) |>
  add_time("Date") |>
  summarise(
    Precip_mm = sum(Precipitation, na.rm = TRUE),
    .by = SampleDate
  ) |>
  filter(
    SampleDate > as.Date("2021-03-23"),
    SampleDate < as.Date("2021-11-25")
  )

irrigation_daily <- read_input(irrigation_file) |>
  mutate(SampleDate = parse_date(SampleDate)) |>
  rename(Rewet_mm = Irrigation) |>
  filter(
    SampleDate >= min(precipitation_daily$SampleDate, na.rm = TRUE),
    SampleDate <= max(precipitation_daily$SampleDate, na.rm = TRUE)
  )

## Flux table ----
build_flux_tbl <- function(filename) {
  flux_long <- read_input(filename) |>
    pivot_longer(
      cols = 2:17,
      names_to = "KLPlot",
      values_to = "N2ON"
    ) |>
    add_time("Date") |>
    mutate(
      PlotNum = suppressWarnings(as.integer(str_extract(KLPlot, "\\d+")))
    )

  if (anyNA(flux_long$PlotNum)) {
    stop(
      "Could not parse plot numbers: ",
      paste(
        unique(flux_long$KLPlot[is.na(flux_long$PlotNum)]),
        collapse = ", "
      )
    )
  }

  flux_long <- flux_long |> filter(PlotNum != 10)

  treatment_four <- case_when(
    flux_long$PlotNum %in% c(1, 5, 13) ~ "C",
    flux_long$PlotNum %in% c(2, 6, 9, 14) ~ "C+N",
    flux_long$PlotNum %in% c(3, 7, 11, 15) ~ "DRW",
    flux_long$PlotNum %in% c(4, 8, 12, 16) ~ "DRW+N",
    TRUE ~ NA_character_
  )

  if (anyNA(treatment_four)) stop("Unmapped plot numbers present.")

  flux_long |>
    mutate(
      Treatment2 = factor(
        ifelse(treatment_four %in% c("C", "C+N"), "NP", "RE"),
        levels = c("NP", "RE")
      ),
      TimeCat = time_category_n2o(Time)
    ) |>
    group_by(Treatment2, SampleDate, TimeCat) |>
    summarise(N2O = mean(as.numeric(N2ON), na.rm = TRUE), .groups = "drop")
}

## Sensors ----
make_depth_sensors <- function(depth, sm_file, st_file) {
  depth <- as.character(depth)
  stopifnot(depth %in% c("5", "15", "30"))

  moisture_raw <- read_input(sm_file)
  moisture_time_column <- if ("X" %in% names(moisture_raw)) {
    "X"
  } else {
    names(moisture_raw)[1]
  }
  moisture_raw <- moisture_raw |> add_time(moisture_time_column)

  moisture_columns <- switch(depth,
    "5" = c("C1_5cm", "C2_5cm", "D1_5cm", "D2_5cm"),
    "15" = c("C1_15cm", "C2_15cm", "D1_15cm", "D2_15cm"),
    "30" = c("C1_30cm", "C2_30cm", "D1_30cm", "D2_30cm")
  )
  stopifnot(all(moisture_columns %in% names(moisture_raw)))

  moisture <- moisture_raw |>
    select(SampleDate, Time, all_of(moisture_columns)) |>
    pivot_longer(
      cols = -c(SampleDate, Time),
      names_to = "Well",
      values_to = "Moisture_pc"
    ) |>
    mutate(
      TimeCat = time_category_n2o(Time),
      Treatment2 = factor(
        ifelse(grepl("^C", Well), "NP", "RE"),
        levels = c("NP", "RE")
      )
    ) |>
    group_by(Treatment2, SampleDate, TimeCat) |>
    summarise(Moisture_pc = mean(Moisture_pc, na.rm = TRUE), .groups = "drop")

  temperature_raw <- read_input(st_file)
  temperature_time_column <- if ("X" %in% names(temperature_raw)) {
    "X"
  } else {
    names(temperature_raw)[1]
  }
  temperature_raw <- temperature_raw |> add_time(temperature_time_column)

  temperature_columns <- switch(depth,
    "5" = c("C1", "C2", "D1", "D2"),
    "15" = c("C1_15", "C2_15", "D1_15", "D2_15"),
    "30" = c("C1_30", "C2_30", "D1_30", "D2_30")
  )
  stopifnot(all(temperature_columns %in% names(temperature_raw)))

  temperature <- temperature_raw |>
    select(SampleDate, Time, all_of(temperature_columns)) |>
    pivot_longer(
      cols = -c(SampleDate, Time),
      names_to = "Well",
      values_to = "Temp_C"
    ) |>
    mutate(
      TimeCat = time_category_n2o(Time),
      Treatment2 = factor(
        ifelse(grepl("^C", Well), "NP", "RE"),
        levels = c("NP", "RE")
      )
    ) |>
    group_by(Treatment2, SampleDate, TimeCat) |>
    summarise(Temp_C = mean(Temp_C, na.rm = TRUE), .groups = "drop")

  list(sm = moisture, st = temperature, depth = depth)
}

## Join sensors and fluxes ----

join_depth <- function(depth_sens, flux_tbl) {
  d <- depth_sens$depth
  por <- unname(porosity$KL[[d]])
  stopifnot(is.finite(por))

  depth_sens$st |>
    full_join(
      depth_sens$sm,
      by = c("Treatment2", "SampleDate", "TimeCat")
    ) |>
    full_join(
      flux_tbl,
      by = c("Treatment2", "SampleDate", "TimeCat")
    ) |>
    mutate(
      WFPS = (Moisture_pc / por) * 100,
      Depth = factor(d, levels = c("5", "15", "30"))
    )
}

## Daily averages and lags ----

build_daily_tbl <- function(df_depth) {
  lag_functions <- list(
    lag1 = \(x) lag(x, 1),
    lag2 = \(x) lag(x, 2),
    lag3 = \(x) lag(x, 3),
    lag4 = \(x) lag(x, 4),
    lag5 = \(x) lag(x, 5),
    lag6 = \(x) lag(x, 6),
    lag7 = \(x) lag(x, 7)
  )

  df_depth |>
    group_by(Depth, Treatment2, SampleDate) |>
    summarise(
      N2O = mean(N2O, na.rm = TRUE),
      WFPS = mean(WFPS, na.rm = TRUE),
      Temp_C = mean(Temp_C, na.rm = TRUE),
      .groups = "drop"
    ) |>
    group_by(Depth, Treatment2) |>
    complete(
      SampleDate = seq.Date(
        min(SampleDate),
        max(SampleDate),
        by = "day"
      )
    ) |>
    ungroup() |>
    mutate(y = log(pmax(N2O, 0) + 1e-6)) |>
    arrange(Depth, Treatment2, SampleDate) |>
    group_by(Treatment2, Depth) |>
    mutate(
      t_index = row_number(),
      across(c(WFPS, Temp_C), lag_functions, .names = "{.col}_{.fn}")
    ) |>
    ungroup()
}

## Pooled table ----
flux_tbl <- build_flux_tbl(n2o_file)

sens_5 <- make_depth_sensors(5, sm_file, st_file)
sens_15 <- make_depth_sensors(15, sm_file, st_file)
sens_30 <- make_depth_sensors(30, sm_file, st_file)

daily_5 <- build_daily_tbl(join_depth(sens_5, flux_tbl))
daily_15 <- build_daily_tbl(join_depth(sens_15, flux_tbl))
daily_30 <- build_daily_tbl(join_depth(sens_30, flux_tbl))

water_lag_functions <- list(
  lag1 = \(x) lag(x, 1),
  lag2 = \(x) lag(x, 2),
  lag3 = \(x) lag(x, 3)
)

df_pool <- bind_rows(daily_5, daily_15, daily_30) |>
  mutate(
    Depth = factor(Depth, levels = c("5", "15", "30")),
    Treatment2 = factor(Treatment2, levels = c("NP", "RE"))
  ) |>
  left_join(precipitation_daily, by = "SampleDate") |>
  left_join(irrigation_daily, by = "SampleDate") |>
  mutate(
    Precip_mm = coalesce(Precip_mm, 0),
    Rewet_mm = coalesce(Rewet_mm, 0),
    Precip_NP = ifelse(Treatment2 == "NP", Precip_mm, 0),
    Rewet_RE = ifelse(Treatment2 == "RE", Rewet_mm, 0)
  ) |>
  arrange(Treatment2, Depth, SampleDate) |>
  group_by(Treatment2, Depth) |>
  mutate(
    across(
      Precip_NP,
      water_lag_functions,
      .names = "Precip_NP_{.fn}"
    ),
    across(
      Rewet_RE,
      water_lag_functions,
      .names = "Rewet_RE_{.fn}"
    )
  ) |>
  ungroup() |>
  select(
    SampleDate, Depth, Treatment2,
    N2O, y, WFPS, Temp_C, Precip_mm, Rewet_mm,
    starts_with("WFPS_lag"),
    starts_with("Temp_C_lag"),
    starts_with("Precip_NP_lag"),
    starts_with("Rewet_RE_lag"),
    t_index
  )

## Model data ----
df_pool_lag <- df_pool |>
  drop_na(
    y, WFPS, Temp_C,
    WFPS_lag1, WFPS_lag2, WFPS_lag3,
    WFPS_lag4, WFPS_lag5, WFPS_lag6, WFPS_lag7,
    Temp_C_lag1, Temp_C_lag2, Temp_C_lag3,
    Temp_C_lag4, Temp_C_lag5, Temp_C_lag6, Temp_C_lag7,
    Precip_NP_lag1, Precip_NP_lag2, Precip_NP_lag3,
    Rewet_RE_lag1, Rewet_RE_lag2, Rewet_RE_lag3
  ) |>
  group_by(Treatment2, Depth) |>
  mutate(
    across(
      matches(
        "^(WFPS|Temp_C)_lag[1-7]$|^(Precip_NP|Rewet_RE)_lag[1-3]$"
      ),
      .zscale_safe,
      .names = "{.col}_z"
    )
  ) |>
  ungroup() |>
  mutate(
    Depth = factor(Depth, levels = c("5", "15", "30")),
    Treatment2 = factor(Treatment2, levels = c("NP", "RE")),
    DepthTrt = interaction(Depth, Treatment2, drop = TRUE)
  )

stopifnot(nrow(df_pool_lag) > 10)

saveRDS(
  df_pool_lag,
  model_data_file
)
## Priors and formula ----
pri_lag <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(1), class = "sigma"),
  prior(gamma(2, 0.1), class = "nu"),
  prior(exponential(1), class = "sds")
)

fml_depthvary <- bf(
  y ~ 0 + DepthTrt +
    t2(
      WFPS, Temp_C,
      by = DepthTrt,
      k = full_smooth_basis,
      bs = c("ts", "ts")
    ) +
    WFPS_lag1_z + WFPS_lag2_z + WFPS_lag3_z + WFPS_lag4_z +
    WFPS_lag5_z + WFPS_lag6_z + WFPS_lag7_z +
    Temp_C_lag1_z + Temp_C_lag2_z + Temp_C_lag3_z + Temp_C_lag4_z +
    Temp_C_lag5_z + Temp_C_lag6_z + Temp_C_lag7_z +
    Precip_NP_lag1_z + Precip_NP_lag2_z + Precip_NP_lag3_z +
    Rewet_RE_lag1_z + Rewet_RE_lag2_z + Rewet_RE_lag3_z
)

## Partial curves ----

.partial_curve_1d <- function(fit, df, depth, treatment,
                              vary = c("WFPS", "Temp_C"),
                              n = 400) {
  vary <- match.arg(vary)

  g <- df |>
    filter(Depth == depth, Treatment2 == treatment)

  stopifnot(nrow(g) > 5)

  hold_W <- mean(g$WFPS, na.rm = TRUE)
  hold_T <- mean(g$Temp_C, na.rm = TRUE)

  dt_vals <- unique(g$DepthTrt)
  dt_vals <- dt_vals[!is.na(dt_vals)]
  stopifnot(length(dt_vals) == 1)
  dt_val <- dt_vals[1]

  if (vary == "WFPS") {
    xs <- seq(
      min(g$WFPS, na.rm = TRUE),
      max(g$WFPS, na.rm = TRUE),
      length.out = n
    )
    nd <- tibble(
      Depth = depth,
      Treatment2 = treatment,
      WFPS = xs,
      Temp_C = hold_T
    )
  } else {
    xs <- seq(
      min(g$Temp_C, na.rm = TRUE),
      max(g$Temp_C, na.rm = TRUE),
      length.out = n
    )
    nd <- tibble(
      Depth = depth,
      Treatment2 = treatment,
      WFPS = hold_W,
      Temp_C = xs
    )
  }

  nd$Depth <- factor(nd$Depth, levels = levels(df$Depth))
  nd$Treatment2 <- factor(nd$Treatment2, levels = levels(df$Treatment2))
  nd$DepthTrt <- factor(dt_val, levels = levels(df$DepthTrt))

  lag_cols <- grep(
    "^(WFPS|Temp_C)_lag[1-7]_z$|^(Precip_NP|Rewet_RE)_lag[1-3]_z$",
    names(df),
    value = TRUE
  )
  nd[lag_cols] <- 0

  if ("t_index" %in% names(df) && !"t_index" %in% names(nd)) nd$t_index <- 1L

  mu <- fitted(fit, newdata = nd, summary = TRUE)

  tibble(
    Treatment2 = treatment,
    Depth      = depth,
    x          = xs,
    mean_log   = mu[, "Estimate"],
    lo_log     = mu[, "Q2.5"],
    hi_log     = mu[, "Q97.5"],
    mean       = pmax(exp(mu[, "Estimate"]) - 1e-6, 0),
    lo         = pmax(exp(mu[, "Q2.5"]) - 1e-6, 0),
    hi         = pmax(exp(mu[, "Q97.5"]) - 1e-6, 0),
    var        = vary
  )
}

## Full model ----

fit_all <- brm(
  formula = fml_depthvary,
  data = df_pool_lag,
  family = student(),
  prior = pri_lag,
  chains = mcmc_chains,
  iter = mcmc_iterations,
  warmup = mcmc_warmup,
  cores = mcmc_cores,
  backend = "cmdstanr",
  save_cmdstan_config = TRUE,
  file_refit = "always",
  file = full_model_cache
)

print(fit_all)

## Diagnostics ----

n_obs <- nrow(df_pool_lag)

r2_draws <- bayes_R2(fit_all, summary = FALSE)
r2_mean <- mean(r2_draws)
r2_ci <- quantile(r2_draws, c(0.025, 0.975), na.rm = TRUE)

fitted_vals <- fitted(fit_all, summary = TRUE)[, "Estimate"]
y_obs <- df_pool_lag$y

RMSE <- sqrt(mean((y_obs - fitted_vals)^2, na.rm = TRUE))
MAE <- mean(abs(y_obs - fitted_vals), na.rm = TRUE)

loo_full <- loo(fit_all)
saveRDS(loo_full, file_full_loo_rds)

post_sum <- summary(fit_all)

sigma_post <- post_sum$spec_pars["sigma", "Estimate"]
nu_post <- post_sum$spec_pars["nu", "Estimate"]
# Fixed coefficients
max_rhat_fixed <- max(post_sum$fixed[, "Rhat"], na.rm = TRUE)
min_ess_fixed <- min(post_sum$fixed[, "Bulk_ESS"], na.rm = TRUE)

# All saved parameters
draw_summary <- summarise_draws(as_draws_array(fit_all))

max_rhat_all <- max(draw_summary$rhat, na.rm = TRUE)
min_ess_all <- min(draw_summary$ess_bulk, na.rm = TRUE)

diagnostic_table <- tibble(
  Observations = n_obs,
  `Bayesian R2 (mean)` = round(r2_mean, 2),
  `Bayesian R2 (low)` = round(r2_ci[1], 3),
  `Bayesian R2 (high)` = round(r2_ci[2], 3),
  `RMSE (log scale)` = round(RMSE, 2),
  `MAE (log scale)` = round(MAE, 2),
  `LOOIC` = round(get_looic(loo_full), 1),
  `LOOIC SE` = round(get_looic_se(loo_full), 1),
  `Residual SD (sigma)` = round(sigma_post, 3),
  `Degrees of freedom (nu)` = round(nu_post, 2),
  `Maximum Rhat (fixed)` = max_rhat_fixed,
  `Minimum Bulk ESS (fixed)` = min_ess_fixed,
  `Maximum Rhat (all saved)` = max_rhat_all,
  `Minimum Bulk ESS (all saved)` = min_ess_all
)

print(diagnostic_table, width = Inf)
write_table(diagnostic_table, file_full_model_diagnostics, dir_reporting)

# LOO comparison ====

## No smooth ----

fml_no_smooth <- bf(
  y ~ 0 + DepthTrt +
    WFPS_lag1_z + WFPS_lag2_z + WFPS_lag3_z + WFPS_lag4_z +
    WFPS_lag5_z + WFPS_lag6_z + WFPS_lag7_z +
    Temp_C_lag1_z + Temp_C_lag2_z + Temp_C_lag3_z + Temp_C_lag4_z +
    Temp_C_lag5_z + Temp_C_lag6_z + Temp_C_lag7_z +
    Precip_NP_lag1_z + Precip_NP_lag2_z + Precip_NP_lag3_z +
    Rewet_RE_lag1_z + Rewet_RE_lag2_z + Rewet_RE_lag3_z
)

## No WFPS lags ----

fml_no_WFPSlags <- update(
  fml_depthvary,
  . ~ . - WFPS_lag1_z - WFPS_lag2_z - WFPS_lag3_z -
    WFPS_lag4_z - WFPS_lag5_z - WFPS_lag6_z - WFPS_lag7_z
)

## No temperature lags ----

fml_no_Templags <- update(
  fml_depthvary,
  . ~ . - Temp_C_lag1_z - Temp_C_lag2_z - Temp_C_lag3_z -
    Temp_C_lag4_z - Temp_C_lag5_z - Temp_C_lag6_z - Temp_C_lag7_z
)

## No precipitation lags ----

fml_no_Precip <- update(
  fml_depthvary,
  . ~ . - Precip_NP_lag1_z - Precip_NP_lag2_z - Precip_NP_lag3_z
)

## No rewetting lags ----

fml_no_Rewet <- update(
  fml_depthvary,
  . ~ . - Rewet_RE_lag1_z - Rewet_RE_lag2_z - Rewet_RE_lag3_z
)

# standard deviations for smooth terms should be removed from smooth term model

pri_no_smooth <- pri_lag[pri_lag$class != "sds", ]

fit_no_smooth <- brm(
  formula = fml_no_smooth,
  data = df_pool_lag,
  family = student(),
  prior = pri_no_smooth,
  chains = mcmc_chains,
  iter = mcmc_iterations,
  warmup = mcmc_warmup,
  cores = mcmc_cores,
  backend = "cmdstanr",
  save_cmdstan_config = TRUE,
  file_refit = "always",
  file = file.path(n2o_cache_dir, cache_model_names[1])
)

fit_no_WFPSlags <- brm(
  formula = fml_no_WFPSlags,
  data = df_pool_lag,
  family = student(),
  prior = pri_lag,
  chains = mcmc_chains,
  iter = mcmc_iterations,
  warmup = mcmc_warmup,
  cores = mcmc_cores,
  backend = "cmdstanr",
  save_cmdstan_config = TRUE,
  file_refit = "always",
  file = file.path(n2o_cache_dir, cache_model_names[2])
)

fit_no_Templags <- brm(
  formula = fml_no_Templags,
  data = df_pool_lag,
  family = student(),
  prior = pri_lag,
  chains = mcmc_chains,
  iter = mcmc_iterations,
  warmup = mcmc_warmup,
  cores = mcmc_cores,
  backend = "cmdstanr",
  save_cmdstan_config = TRUE,
  file_refit = "always",
  file = file.path(n2o_cache_dir, cache_model_names[3])
)

fit_no_Precip <- brm(
  formula = fml_no_Precip,
  data = df_pool_lag,
  family = student(),
  prior = pri_lag,
  chains = mcmc_chains,
  iter = mcmc_iterations,
  warmup = mcmc_warmup,
  cores = mcmc_cores,
  backend = "cmdstanr",
  save_cmdstan_config = TRUE,
  file_refit = "always",
  file = file.path(n2o_cache_dir, cache_model_names[4])
)

fit_no_Rewet <- brm(
  formula = fml_no_Rewet,
  data = df_pool_lag,
  family = student(),
  prior = pri_lag,
  chains = mcmc_chains,
  iter = mcmc_iterations,
  warmup = mcmc_warmup,
  cores = mcmc_cores,
  backend = "cmdstanr",
  save_cmdstan_config = TRUE,
  file_refit = "always",
  file = file.path(n2o_cache_dir, cache_model_names[5])
)

## Compare models ----
loo_no_smooth <- loo(fit_no_smooth)
loo_no_WFPSlags <- loo(fit_no_WFPSlags)
loo_no_Templags <- loo(fit_no_Templags)
loo_no_Precip <- loo(fit_no_Precip)
loo_no_Rewet <- loo(fit_no_Rewet)


loo_tbl <- tibble(
  term = c(
    "smooth_WFPS-Temp",
    "WFPS_lags",
    "Temp_lags",
    "Precip_lags",
    "Rewet_lags"
  ),
  R2_reduced = c(
    get_r2_mean(fit_no_smooth),
    get_r2_mean(fit_no_WFPSlags),
    get_r2_mean(fit_no_Templags),
    get_r2_mean(fit_no_Precip),
    get_r2_mean(fit_no_Rewet)
  ),
  R2_full = get_r2_mean(fit_all),
  LOO_full = get_looic(loo_full),
  LOO_reduced = c(
    get_looic(loo_no_smooth),
    get_looic(loo_no_WFPSlags),
    get_looic(loo_no_Templags),
    get_looic(loo_no_Precip),
    get_looic(loo_no_Rewet)
  )
) |>
  mutate(
    delta_R2    = R2_full - R2_reduced,
    delta_LOOIC = LOO_reduced - LOO_full
  ) |>
  select(
    term,
    R2_reduced,
    R2_full,
    delta_R2,
    LOO_full,
    LOO_reduced,
    delta_LOOIC
  ) |>
  arrange(desc(delta_LOOIC))

print(loo_tbl, width = Inf)


write_table(loo_tbl, file_table_05, dir_reporting)

cat("\nLOO compare:\n")
print(loo_compare(
  loo_full,
  loo_no_smooth,
  loo_no_WFPSlags,
  loo_no_Templags,
  loo_no_Precip,
  loo_no_Rewet
))

## Full-model curves ----

treatments <- levels(df_pool_lag$Treatment2)
depths <- levels(df_pool_lag$Depth)

curves_W <- map_dfr(
  treatments,
  \(trt) {
    map_dfr(
      depths,
      \(depth) .partial_curve_1d(
        fit = fit_all,
        df = df_pool_lag,
        depth = depth,
        treatment = trt,
        vary = "WFPS",
        n = full_curve_points
      )
    )
  }
)

curves_T <- map_dfr(
  treatments,
  \(trt) {
    map_dfr(
      depths,
      \(depth) .partial_curve_1d(
        fit = fit_all,
        df = df_pool_lag,
        depth = depth,
        treatment = trt,
        vary = "Temp_C",
        n = full_curve_points
      )
    )
  }
)

write_table(curves_W, file_full_wfps_curves, dir_processed)
write_table(curves_T, file_full_temp_curves, dir_processed)
saveRDS(list(wfps = curves_W, temp = curves_T), file_full_curves_rds)

## Console summary ----
model_outputs <- list(
  diagnostic_table     = diagnostic_table,
  loo_comparison_table = loo_tbl,
  curves_W             = curves_W,
  curves_T             = curves_T,
  fit_all              = fit_all
)

# Reduced Bayesian Model ====


df <- df_pool_lag

fml <- bf(
  y ~ 0 + DepthTrt +
    t2(
      WFPS, Temp_C,
      by = DepthTrt,
      k = reduced_smooth_basis,
      bs = c("ts", "ts")
    ) +
    Temp_C_lag1_z + Temp_C_lag2_z + Temp_C_lag3_z +
    Temp_C_lag4_z + Temp_C_lag5_z + Temp_C_lag6_z +
    Temp_C_lag7_z +
    Precip_NP_lag1_z + Precip_NP_lag2_z + Precip_NP_lag3_z
)

fit <- brm(
  formula = fml,
  data = df,
  family = student(),
  prior = pri_lag,
  chains = mcmc_chains,
  iter = mcmc_iterations,
  warmup = mcmc_warmup,
  cores = mcmc_cores,
  backend = "cmdstanr",
  save_cmdstan_config = TRUE,
  file = reduced_model_cache,
  file_refit = "always",
  control = list(
    adapt_delta = mcmc_adapt_delta,
    max_treedepth = mcmc_max_treedepth
  )
)


print(fit)
print(bayes_R2(fit))

## Partial curves ----
grid_tbl <- expand.grid(
  Depth = model_depths,
  Treatment2 = model_treatments,
  stringsAsFactors = FALSE
)

curves_W <- pmap_dfr(
  list(grid_tbl$Depth, grid_tbl$Treatment2),
  \(depth, treatment) {
    .partial_curve_1d(
      fit = fit,
      df = df,
      depth = depth,
      treatment = treatment,
      vary = "WFPS",
      n = reduced_curve_points
    ) |>
      select(
        Depth, Treatment2, var, x, mean_log, lo_log, hi_log,
        mean, lo, hi
      )
  }
)

curves_T <- pmap_dfr(
  list(grid_tbl$Depth, grid_tbl$Treatment2),
  \(depth, treatment) {
    .partial_curve_1d(
      fit = fit,
      df = df,
      depth = depth,
      treatment = treatment,
      vary = "Temp_C",
      n = reduced_curve_points
    ) |>
      select(
        Depth, Treatment2, var, x, mean_log, lo_log, hi_log,
        mean, lo, hi
      )
  }
)

## Threshold helpers ----
compute_slopes <- function(curve_df) {
  curve_df |>
    arrange(x) |>
    mutate(
      dx = lead(x) - x,
      slope = (lead(mean) - mean) / dx
    ) |>
    filter(is.finite(dx), is.finite(slope), dx > 0)
}

find_threshold_piecewise <- function(
  curve_df,
  inner_frac = threshold_inner_fraction,
  min_prop_segment = threshold_min_segment_fraction
) {
  d <- curve_df |> arrange(x)
  n <- nrow(d)

  if (n < threshold_min_curve_points) {
    return(NA_real_)
  }

  x <- d$x
  y <- d$mean

  alpha <- (1 - inner_frac) / 2
  x_lo <- quantile(x, alpha, na.rm = TRUE)
  x_hi <- quantile(x, 1 - alpha, na.rm = TRUE)

  min_n_seg <- ceiling(min_prop_segment * n)

  candidates <- sort(unique(x[x >= x_lo & x <= x_hi]))
  candidates <- candidates[
    sapply(candidates, function(cut) {
      sum(x <= cut) >= min_n_seg && sum(x > cut) >= min_n_seg
    })
  ]

  if (!length(candidates)) {
    return(NA_real_)
  }

  sse_fun <- function(cut) {
    idx1 <- x <= cut
    idx2 <- x > cut
    fit1 <- lm(y[idx1] ~ x[idx1])
    fit2 <- lm(y[idx2] ~ x[idx2])
    sum(residuals(fit1)^2) + sum(residuals(fit2)^2)
  }

  sse_vals <- vapply(candidates, sse_fun, numeric(1))
  candidates[which.min(sse_vals)]
}

find_threshold_slope <- function(curve_df) {
  s <- compute_slopes(curve_df)
  if (nrow(s) == 0) {
    return(NA_real_)
  }
  s$x[which.max(s$slope)]
}

find_temp_ddpeak <- function(curve_df) {
  d <- curve_df |>
    arrange(x) |>
    mutate(
      dx = lead(x) - x,
      d1 = (lead(mean) - mean) / dx
    ) |>
    mutate(
      d2 = (lead(d1) - lag(d1)) / (lead(x) - lag(x))
    ) |>
    filter(is.finite(d2))

  if (nrow(d) < 3) {
    return(NA_real_)
  }

  d$x[which.max(d$d2)]
}

## Curves and thresholds ----
wfps_thresholds <- curves_W |>
  group_by(Depth, Treatment2) |>
  group_modify(\(data, key) {
    piece <- find_threshold_piecewise(data)
    slope <- find_threshold_slope(data)
    tibble(
      wfps_thresh_piece = piece,
      wfps_thresh_slope = slope,
      wfps_thresh = coalesce(piece, slope)
    )
  }) |>
  ungroup()

temp_thresholds <- curves_T |>
  group_by(Depth, Treatment2) |>
  group_modify(\(data, key) {
    tibble(
      temp_thresh = find_temp_ddpeak(data)
    )
  }) |>
  ungroup()

res <- wfps_thresholds |>
  left_join(temp_thresholds, by = c("Depth", "Treatment2")) |>
  mutate(
    Depth_layer = depth_range_label(Depth)
  ) |>
  select(Treatment2, Depth, Depth_layer, wfps_thresh, temp_thresh)

## Comparisons ----
trt_diff <- res |>
  select(Depth, Treatment2, wfps_thresh) |>
  pivot_wider(names_from = Treatment2, values_from = wfps_thresh) |>
  mutate(diff_RE_minus_NP = RE - NP)

depth_shift <- res |>
  select(Treatment2, Depth, wfps_thresh) |>
  pivot_wider(names_from = Depth, values_from = wfps_thresh) |>
  mutate(
    shift_5_15  = `15` - `5`,
    shift_15_30 = `30` - `15`
  )

## Diagnostics ----
n_obs <- nrow(df)

r2_draws <- bayes_R2(fit, summary = FALSE)
r2_mean <- mean(r2_draws)
r2_ci <- quantile(r2_draws, c(0.025, 0.975))

fitted_vals <- fitted(fit)[, "Estimate"]
y_obs <- df$y

RMSE <- sqrt(mean((y_obs - fitted_vals)^2))
MAE <- mean(abs(y_obs - fitted_vals))

loo_fit <- loo(fit)
saveRDS(loo_fit, file_reduced_loo_rds)
looic <- loo_fit$estimates["looic", "Estimate"]
looic_se <- loo_fit$estimates["looic", "SE"]

loo_full_vs_reduced <- loo_compare(
  list(
    Full_GAM = loo_full,
    Reduced_GAM = loo_fit
  )
)

print(loo_full_vs_reduced)

write.csv(
  as.data.frame(loo_full_vs_reduced),
  file.path(diagnostics_dir, file_full_reduced_loo_csv),
  row.names = TRUE
)

post_sum <- summary(fit)

sigma_post <- post_sum$spec_pars["sigma", "Estimate"]
nu_post <- post_sum$spec_pars["nu", "Estimate"]

max_rhat_fixed <- max(post_sum$fixed[, "Rhat"], na.rm = TRUE)
min_ess_fixed <- min(post_sum$fixed[, "Bulk_ESS"], na.rm = TRUE)

# All saved parameters
draw_summary <- summarise_draws(as_draws_array(fit))

max_rhat_all <- max(draw_summary$rhat, na.rm = TRUE)
min_ess_all <- min(draw_summary$ess_bulk, na.rm = TRUE)

diagnostic_table <- tibble(
  Observations = n_obs,
  `Bayesian R2 (mean)` = round(r2_mean, 2),
  `Bayesian R2 (low)` = round(r2_ci[1], 3),
  `Bayesian R2 (high)` = round(r2_ci[2], 3),
  `RMSE (log scale)` = round(RMSE, 2),
  `MAE (log scale)` = round(MAE, 2),
  `LOOIC` = round(looic, 1),
  `LOOIC SE` = round(looic_se, 1),
  `Residual SD (sigma)` = round(sigma_post, 3),
  `Degrees of freedom (nu)` = round(nu_post, 2),
  `Maximum Rhat (fixed)` = max_rhat_fixed,
  `Minimum Bulk ESS (fixed)` = min_ess_fixed,
  `Maximum Rhat (all saved)` = max_rhat_all,
  `Minimum Bulk ESS (all saved)` = min_ess_all
)

print(diagnostic_table, width = Inf)

## Output ----
print(res)
print(trt_diff)
print(depth_shift)

write_table(res, file_threshold_results, dir_processed)
write_table(res, file_threshold_results_table, dir_processed)
write_table(trt_diff, file_treatment_differences, dir_processed)
write_table(depth_shift, file_depth_shifts, dir_processed)
write_table(
  diagnostic_table,
  file_reduced_model_diagnostics,
  n2o_interim_dir
)

model_outputs_reduced <- list(
  diagnostic_table = diagnostic_table,
  curves_W = curves_W,
  curves_T = curves_T,
  fit = fit,
  loo = loo_fit,
  thresholds = res,
  treatment_differences = trt_diff,
  depth_shifts = depth_shift
)

# Full vs reduced GAM comparison ====

full_diag <- model_outputs$diagnostic_table
reduced_diag <- model_outputs_reduced$diagnostic_table

format_diagnostics <- function(d) {
  c(
    sprintf("%.2f", d[["Bayesian R2 (mean)"]]),
    sprintf(
      "%.3f – %.3f",
      d[["Bayesian R2 (low)"]],
      d[["Bayesian R2 (high)"]]
    ),
    sprintf("%.2f", d[["RMSE (log scale)"]]),
    sprintf("%.2f", d[["MAE (log scale)"]]),
    sprintf("%.1f", d[["LOOIC"]]),
    sprintf("%.1f", d[["LOOIC SE"]]),
    sprintf("%.3f", d[["Residual SD (sigma)"]]),
    sprintf("%.2f", d[["Degrees of freedom (nu)"]]),
    sprintf("%.3f", d[["Maximum Rhat (fixed)"]]),
    sprintf("%.0f", d[["Minimum Bulk ESS (fixed)"]]),
    sprintf("%.3f", d[["Maximum Rhat (all saved)"]]),
    sprintf("%.0f", d[["Minimum Bulk ESS (all saved)"]]),
    sprintf("%.0f", d[["Observations"]])
  )
}

comparison_table <- tibble(
  Metric = c(
    "Bayesian R² (mean)",
    "Bayesian R² (95% CI)",
    "RMSE (log scale)",
    "MAE (log scale)",
    "LOOIC",
    "LOOIC SE",
    "Residual SD (σ)",
    "Degrees of freedom (ν)",
    "Maximum R̂ (fixed coefficients)",
    "Minimum Bulk ESS (fixed coefficients)",
    "Maximum R̂ (all saved parameters)",
    "Minimum Bulk ESS (all saved parameters)",
    "Observations"
  ),
  `Full GAM` = format_diagnostics(full_diag),
  `Reduced GAM` = format_diagnostics(reduced_diag)
)

print(comparison_table, n = Inf, width = Inf)

write_table(comparison_table, file_full_reduced_comparison, dir_reporting)


## Partial curves ----
curves_W <- curves_W |>
  mutate(
    Depth_lab = factor(
      depth_range_label(Depth),
      levels = names(model_plot_depth_labels),
      labels = unname(model_plot_depth_labels)
    )
  )

curves_T <- curves_T |>
  mutate(
    Depth_lab = factor(
      depth_range_label(Depth),
      levels = names(model_plot_depth_labels),
      labels = unname(model_plot_depth_labels)
    )
  )

write_table(curves_W, file_reduced_wfps_curves, dir_processed)
write_table(curves_T, file_reduced_temp_curves, dir_processed)
saveRDS(list(wfps = curves_W, temp = curves_T), file_reduced_curves_rds)

write_session_info(dir_interim)


# Diagnostics ====

export_model_diagnostics <- function(
  model_fit, loo_result, model_name, treedepth_limit
) {
  ## Convergence ----
  draws <- as_draws_array(model_fit)
  parameter_diag <- summarise_draws(draws)

  write_table(
    as.data.frame(parameter_diag),
    paste0(model_name, "-", file_diagnostic_parameter_suffix),
    diagnostics_dir
  )

  ## Sampler ----
  nuts <- nuts_params(model_fit)

  sampler_diag <- nuts |>
    group_by(Chain) |>
    summarise(
      Divergences = sum(Value[Parameter == "divergent__"]),
      Treedepth_limit_hits = sum(
        Value[Parameter == "treedepth__"] >= treedepth_limit
      ),
      Maximum_treedepth = max(Value[Parameter == "treedepth__"]),
      .groups = "drop"
    )

  # bfmi
  bfmi_diag <- nuts |>
    filter(Parameter == "energy__") |>
    arrange(Chain, Iteration) |>
    group_by(Chain) |>
    summarise(
      BFMI = mean(diff(Value)^2) / var(Value),
      .groups = "drop"
    )

  sampler_diag <- left_join(sampler_diag, bfmi_diag, by = "Chain")

  write_table(
    sampler_diag,
    paste0(model_name, "-", file_diagnostic_sampler_suffix),
    diagnostics_dir
  )

  ## LOO summary ----
  writeLines(
    capture.output({
      print(loo_result)
      print(pareto_k_table(loo_result))
    }),
    file.path(
      diagnostics_dir,
      paste0(model_name, "-", file_diagnostic_loo_suffix)
    )
  )

  print(sampler_diag)

  invisible(sampler_diag)
}


sampler_full <- export_model_diagnostics(
  fit_all,
  loo_full,
  full_diagnostic_model_name,
  treedepth_limit = full_diagnostic_treedepth_limit
)

sampler_reduced <- export_model_diagnostics(
  fit,
  loo_fit,
  reduced_diagnostic_model_name,
  treedepth_limit = reduced_diagnostic_treedepth_limit
)

## Combined sampler table ----
sampler_comparison <- bind_rows(
  mutate(sampler_full, Model = "Full GAM", .before = 1),
  mutate(sampler_reduced, Model = "Reduced GAM", .before = 1)
)

write_table(sampler_comparison, file_sampler_comparison, diagnostics_dir)

## Model notes ----
writeLines(
  capture.output(print(loo_full_vs_reduced)),
  file.path(diagnostics_dir, file_full_reduced_loo_text)
)

writeLines(
  c(
    paste("Seed:", model_seed),
    paste("Observations per model:", nrow(df_pool_lag)),
    "Both models use the same prepared rows and standardised predictors.",
    "Reduced model removes WFPS lags 1-7 and rewetting lags 1-3.",
    "Response: log(pmax(N2O, 0) + 1e-6).",
    "Family: Student-t.",
    "Predictive checks are shown on the log-response scale."
  ),
  file.path(diagnostics_dir, file_model_notes)
)


# Depth-specific GAMs ====

depth_wfps_data <- df_pool |>
  transmute(
    Treatment = as.character(Treatment2),
    Depth = as.character(Depth),
    N2O,
    WFPS
  ) |>
  filter(is.finite(N2O), is.finite(WFPS))

groups <- depth_wfps_data |>
  distinct(Treatment, Depth) |>
  arrange(Treatment, as.numeric(Depth))

depth_wfps_gams <- list()
summary_rows <- list()
prediction_rows <- list()

for (i in seq_len(nrow(groups))) {
  trt <- groups$Treatment[i]
  depth <- groups$Depth[i]
  model_id <- paste(trt, depth, sep = "_")

  dat <- depth_wfps_data |>
    filter(Treatment == trt, Depth == depth)

  if (
    nrow(dat) < depth_gam_min_observations ||
      n_distinct(dat$WFPS) < depth_gam_basis_k
  ) {
    stop("Insufficient observations or distinct WFPS values: ", model_id)
  }

  # Same basis size and fitting method across groups.
  # "ts" permits the entire smooth to shrink towards zero.
  model <- gam(
    N2O ~ s(WFPS, bs = "ts", k = depth_gam_basis_k),
    data = dat,
    family = gaussian(link = "identity"),
    method = "GCV.Cp"
  )

  depth_wfps_gams[[model_id]] <- model
  sm <- summary(model)

  summary_rows[[i]] <- tibble(
    Treatment = trt,
    Depth = depth,
    n = nobs(model),
    R2_adjusted = sm$r.sq,
    Deviance_explained = sm$dev.expl,
    AIC = AIC(model),
    Smooth_EDF = unname(sm$s.table[1, "edf"])
  )

  prediction_grid <- data.frame(
    WFPS = seq(
      min(dat$WFPS),
      max(dat$WFPS),
      length.out = depth_gam_curve_points
    )
  )

  pred <- predict(
    model,
    newdata = prediction_grid,
    type = "response",
    se.fit = TRUE
  )

  prediction_rows[[i]] <- tibble(
    Treatment = trt,
    Depth = depth,
    WFPS = prediction_grid$WFPS,
    estimate = as.numeric(pred$fit),
    lower = as.numeric(pred$fit - 1.96 * pred$se.fit),
    upper = as.numeric(pred$fit + 1.96 * pred$se.fit)
  )
}

depth_wfps_summary <- bind_rows(summary_rows)
depth_wfps_predictions <- bind_rows(prediction_rows)

# Table matching your requested layout
depth_wfps_table <- depth_wfps_summary |>
  transmute(
    Treatment,
    `Depth (cm)` = unname(model_depth_labels[Depth]),
    n,
    `Adjusted R²` = round(R2_adjusted, 3),
    `Deviance explained` = round(Deviance_explained, 3),
    AIC = round(AIC, 1)
  )

print(depth_wfps_table, n = Inf, width = Inf)

write_table(depth_wfps_table, file_depth_wfps_table, dir_reporting)

# Also retain unrounded results and fitted models
saveRDS(
  depth_wfps_gams,
  file_depth_gam_models
)

write_table(depth_wfps_summary, file_depth_wfps_results, dir_processed)
write_table(
  depth_wfps_predictions,
  file_depth_wfps_predictions,
  dir_processed
)

saveRDS(
  list(data = depth_wfps_data, predictions = depth_wfps_predictions),
  file_depth_gam_data
)

## Depth GAM diagnostics ----

basis_checks <- imap_dfr(
  depth_wfps_gams,
  \(original_model, id) {
    dat <- model.frame(original_model)

    models <- list(
      ts_k6 = original_model,
      ts_k10 = gam(
        N2O ~ s(WFPS, bs = "ts", k = depth_gam_sensitivity_k),
        data = dat,
        family = gaussian(),
        method = "GCV.Cp"
      ),
      tp_k6 = gam(
        N2O ~ s(WFPS, bs = "tp", k = depth_gam_basis_k),
        data = dat,
        family = gaussian(),
        method = "GCV.Cp"
      )
    )

    imap_dfr(models, \(m, setting) {
      sm <- summary(m)

      tibble(
        Group = id,
        Setting = setting,
        EDF = unname(sm$s.table[1, "edf"]),
        Adjusted_R2 = sm$r.sq,
        Deviance_explained = sm$dev.expl,
        AIC = AIC(m),
        GCV = unname(m$gcv.ubre)
      )
    })
  }
)

print(basis_checks, n = Inf, width = Inf)
write_table(basis_checks, file_depth_wfps_basis_checks, diagnostics_dir)

set.seed(depth_gam_diagnostic_seed)

k_check_table <- imap_dfr(
  depth_wfps_gams,
  \(model, id) {
    kc <- k.check(model, n.rep = gam_check_replicates)

    tibble(
      Group = id,
      Smooth = rownames(kc),
      k_available = kc[, "k'"],
      EDF = kc[, "edf"],
      k_index = kc[, "k-index"],
      p_value = kc[, "p-value"]
    )
  }
)

print(k_check_table, n = Inf, width = Inf)
write_table(k_check_table, file_depth_wfps_k_checks, diagnostics_dir)
log_script_end("n2o-bayesian-models")
