# Repeated-measures correlations: CH4 uptake and soil properties ----
# Inputs: CH4 fluxes and soil moisture/temperature at KL and ZB.
# Output: Table 4 CSV files in dat/reporting.

suppressPackageStartupMessages({
  library("tidyverse")
  library("rmcorr")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("table-04-ch4-correlations")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))

# Settings ====

kl_excluded_plot_numbers <- 10L
kl_start_date <- NULL
kl_end_date <- NULL
zb_start_date <- NULL
zb_end_date <- NULL

table_depth_labels <- c(
  `5` = "0 – 10 cm",
  `15` = "10 – 20 cm",
  `30` = "20 – 40 cm"
)

# ZB bulk density from UBA data
ch4_bulk_density_zb <- 0.308477319997918

sensor_map_zb_moisture_ch4 <- tibble(
  Sensor = c("C1", "S1", "S2", "S3", "S4"),
  Depth_cm = 5,
  Treatment = c("NP", "RE", "RE", "RE", "RE")
)

file_table_04 <- "table-04-ch4-correlations"

summary_site_depths <- list(KL = c(5, 15, 30), ZB = 5)

correlation_groups <- imap(summary_site_depths, \(depths, site) {
  expand_grid(
    Site = site,
    Depth_cm = depths,
    Treatment = c("NP", "RE")
  )
}) |>
  list_rbind()

predictors <- c(
  WFPS = "WFPS %",
  Tsoil = "Soil Temperature"
)

minimum_observations_per_plot <- 2L
minimum_plots_per_correlation <- 2L
p_value_thresholds <- c(0.001, 0.01)
p_value_labels <- paste0("<", c("0.001", "0.01"))
p_value_digits <- 2L

# Data ====

restrict_dates <- function(dat, start_date, end_date) {
  if (!is.null(start_date)) {
    dat <- dat |> filter(SampleDate >= .env$start_date)
  }

  if (!is.null(end_date)) {
    dat <- dat |> filter(SampleDate <= .env$end_date)
  }

  dat
}

# CH4 flux preparation ----

prepare_ch4 <- function(filename, site, excluded = integer()) {
  dat <- read_input(filename)

  pattern <- if (site == "KL") {
    "^KL[0-9]+($|\\.)"
  } else {
    "^ZB[0-9]+($|\\.)"
  }

  plot_columns <- names(dat)[str_detect(names(dat), pattern)]

  if (!length(plot_columns)) {
    stop(filename, ": no recognised CH4 plot columns.")
  }

  plot_numbers <- as.integer(str_replace(
    plot_columns, "^[A-Za-z]+0*([0-9]+).*$", "\\1"
  ))

  if (
    anyNA(plot_numbers) ||
      any(!plot_numbers %in% 1:16) ||
      anyDuplicated(plot_numbers)
  ) {
    stop(filename, ": unrecognised or duplicate plot identifiers.")
  }

  if (!all(map_lgl(dat[plot_columns], is.numeric))) {
    stop(filename, ": CH4 plot columns must be numeric.")
  }

  plot_map <- tibble(
    Source_plot = plot_columns,
    Plot_number = plot_numbers
  )

  dat <- dat |>
    add_time("Date") |>
    select(SampleDate, Time, TimeCat, all_of(plot_columns)) |>
    pivot_longer(
      all_of(plot_columns),
      names_to = "Source_plot",
      values_to = "CH4_flux"
    ) |>
    left_join(plot_map, by = "Source_plot") |>
    mutate(
      Plot = sprintf("%s%02d", site, Plot_number),
      Treatment = if_else(
        Plot_number %in% c(1, 2, 5, 6, 9, 10, 13, 14),
        "NP",
        "RE"
      ),
      CH4_uptake = -CH4_flux
    ) |>
    filter(
      !Plot_number %in% excluded,
      !is.na(TimeCat),
      is.finite(CH4_uptake)
    )

  if (site == "KL") {
    dat <- dat |>
      group_by(Plot, Treatment, SampleDate, Time, TimeCat) |>
      summarise(CH4_uptake = mean(CH4_uptake), .groups = "drop")
  }

  dat |>
    group_by(Plot, Treatment, SampleDate, TimeCat) |>
    summarise(CH4_uptake = mean(CH4_uptake), .groups = "drop")
}

# KL flux and soil data ----

kl_flux <- prepare_ch4(
  "ch4-klausen.csv",
  site = "KL",
  excluded = kl_excluded_plot_numbers
) |>
  restrict_dates(kl_start_date, kl_end_date)

kl_moisture <- prepare_sensors(
  "soil-moisture-klausen.csv",
  "X",
  sensor_map_kl_moisture,
  "Moisture_pc",
  require_numeric = TRUE
)

kl_temperature <- prepare_sensors(
  "soil-temperature-klausen.csv",
  "X",
  sensor_map_kl_temperature,
  "Tsoil",
  equal_sensor_weight = TRUE,
  require_numeric = TRUE
)

sensor_keys <- c("SampleDate", "TimeCat", "Depth_cm", "Treatment")
flux_keys <- c("SampleDate", "TimeCat", "Treatment")

kl_soil <- full_join(
  kl_moisture, kl_temperature,
  by = sensor_keys
)

kl_data <- inner_join(
  kl_soil, kl_flux,
  by = flux_keys,
  relationship = "many-to-many"
) |>
  mutate(
    Site = "KL",
    Porosity_pc = unname(porosity$KL[as.character(Depth_cm)]),
    WFPS = 100 * Moisture_pc / Porosity_pc
  )

# ZB flux and soil data ----

zb_flux <- prepare_ch4(
  "ch4-zoe.csv",
  site = "ZB"
) |>
  restrict_dates(zb_start_date, zb_end_date)

zb_moisture <- prepare_sensors(
  "soil-moisture-zb.csv",
  "Date",
  sensor_map_zb_moisture_ch4,
  "Moisture_pc",
  require_numeric = TRUE
)

zb_temperature <- prepare_sensors(
  "soil-temperature-zb.csv",
  "SampleDate",
  sensor_map_zb_temperature,
  "Tsoil",
  require_numeric = TRUE
)

zb_soil <- full_join(
  zb_moisture, zb_temperature,
  by = sensor_keys
)

zb_data <- inner_join(
  zb_soil, zb_flux,
  by = flux_keys,
  relationship = "one-to-many"
) |>
  mutate(
    Site = "ZB",
    Porosity_pc = 100 * (1 - ch4_bulk_density_zb / particle_density),
    WFPS = 100 * Moisture_pc / Porosity_pc
  )

correlation_data <- bind_rows(kl_data, zb_data) |>
  mutate(Depth = unname(table_depth_labels[as.character(Depth_cm)])) |>
  select(
    Site, Depth_cm, Depth, Treatment, Plot, SampleDate, TimeCat,
    CH4_uptake, WFPS, Tsoil, Moisture_pc, Porosity_pc
  ) |>
  arrange(Site, Depth_cm, Treatment, Plot, SampleDate, TimeCat)

# Repeated-measures correlations ====

## Fit models ----

run_rmcorr <- function(dat, predictor, predictor_label) {
  dat$Predictor <- dat[[predictor]]

  dat <- filter(
    dat,
    !is.na(Plot),
    is.finite(CH4_uptake),
    is.finite(Predictor)
  )

  dat <- dat |>
    group_by(Plot) |>
    filter(n() >= minimum_observations_per_plot) |>
    ungroup()
  dat$Plot <- factor(dat$Plot)

  number_plots <- n_distinct(dat$Plot)
  residual_df <- nrow(dat) - number_plots - 1

  if (number_plots < minimum_plots_per_correlation || residual_df <= 0) {
    stop("Insufficient plots or repeated observations.")
  }

  # Check that both measures vary within plots.
  within_flux <- dat$CH4_uptake -
    ave(dat$CH4_uptake, dat$Plot, FUN = mean)

  within_predictor <- dat$Predictor -
    ave(dat$Predictor, dat$Plot, FUN = mean)

  if (
    sum(within_flux^2) <= 0 ||
      sum(within_predictor^2) <= 0
  ) {
    stop("One or both measures have no within-plot variation.")
  }

  fit <- rmcorr(
    participant = Plot,
    measure1 = CH4_uptake,
    measure2 = Predictor,
    dataset = dat
  )

  result <- data.frame(
    Site = dat$Site[1],
    Depth_cm = dat$Depth_cm[1],
    Depth = dat$Depth[1],
    Treatment = dat$Treatment[1],
    Predictor = predictor_label,
    r = as.numeric(fit$r),
    p_value = as.numeric(fit$p),
    CI_lower = as.numeric(fit$CI[1]),
    CI_upper = as.numeric(fit$CI[2]),
    df = as.numeric(fit$df),
    N = nrow(dat),
    N_plots = number_plots,
    stringsAsFactors = FALSE
  )

  list(
    model = fit,
    result = result,
    data = dat
  )
}

rmcorr_models <- list()
result_rows <- list()
used_data_rows <- list()

for (i in seq_len(nrow(correlation_groups))) {
  site_i <- correlation_groups$Site[i]
  depth_i <- correlation_groups$Depth_cm[i]
  treatment_i <- correlation_groups$Treatment[i]

  dat <- filter(
    correlation_data,
    .data$Site == .env$site_i,
    .data$Depth_cm == .env$depth_i,
    .data$Treatment == .env$treatment_i
  )

  for (predictor_i in names(predictors)) {
    label_i <- unname(predictors[predictor_i])

    key <- paste(
      site_i, depth_i, treatment_i, predictor_i,
      sep = "_"
    )

    logger::log_info(paste("Calculating correlation:", key))

    output <- tryCatch(
      run_rmcorr(dat, predictor_i, label_i),
      error = \(err) {
        stop(
          key, ": ", conditionMessage(err),
          call. = FALSE
        )
      }
    )

    rmcorr_models[[key]] <- output$model
    result_rows[[key]] <- output$result

    output$data$Correlation <- label_i
    used_data_rows[[key]] <- output$data
  }
}

rmcorr_results <- bind_rows(result_rows)
used_data <- bind_rows(used_data_rows)

## Export ----

date_description <- function(start, end) {
  paste0(
    "Start: ",
    if (is.null(start)) "unrestricted" else as.character(start),
    "; end: ",
    if (is.null(end)) "unrestricted" else as.character(end)
  )
}

settings <- data.frame(
  Setting = c(
    "Response",
    "Participant",
    "Treatment grouping",
    "KL excluded plots",
    "KL date restrictions",
    "ZB date restrictions",
    "Depth mapping",
    "KL bulk densities: 5 / 15 / 30 cm",
    "ZB bulk density",
    "Particle density",
    "ZB NP moisture sensors",
    "ZB RE moisture sensors",
    "Soil predictor assignment",
    "Missing observations",
    "Minimum observations per plot",
    "Correlation display",
    "P-value display",
    "P-values"
  ),
  Value = c(
    "CH4 uptake = -1 * original CH4 flux; no unit conversion",
    "Plot",
    "NP and RE, pooling nitrogen treatments",
    paste0(
      "KL",
      sprintf("%02d", kl_excluded_plot_numbers),
      collapse = ", "
    ),
    date_description(kl_start_date, kl_end_date),
    date_description(zb_start_date, zb_end_date),
    "5 cm -> 0–10 cm; 15 cm -> 10–20 cm; 30 cm -> 20–40 cm",
    paste(bulk_density$KL, collapse = " / "),
    as.character(ch4_bulk_density_zb),
    as.character(particle_density),
    paste(
      sensor_map_zb_moisture_ch4$Sensor[
        sensor_map_zb_moisture_ch4$Treatment == "NP"
      ],
      collapse = ", "
    ),
    paste(
      sensor_map_zb_moisture_ch4$Sensor[
        sensor_map_zb_moisture_ch4$Treatment == "RE"
      ],
      collapse = ", "
    ),
    "Treatment-average soil values assigned to each corresponding plot",
    "Complete cases selected separately for each predictor",
    paste0(
      minimum_observations_per_plot,
      " usable observations; at least ",
      minimum_plots_per_correlation,
      " plots per correlation"
    ),
    paste(p_value_digits, "decimal places"),
    paste0(
      "<", p_value_thresholds[1], "; <", p_value_thresholds[2],
      "; otherwise ", p_value_digits, " decimals"
    ),
    "Calculated by rmcorr; no multiple-testing adjustment"
  )
)

write_csv_sheets(
  list(
    results = rmcorr_results,
    data_used = used_data,
    settings = settings
  ),
  file_table_04
)
log_script_end("table-04-ch4-correlations")
