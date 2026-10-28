# Cumulative fluxes and CO2 equivalents ----
# Calculates cumulative gas fluxes and plot-level CO2 equivalents.
# Inputs: five gas-flux datasets; outputs: Table S7 CSV files in dat/reporting.

suppressPackageStartupMessages({
  library("tidyverse")
  library("car")
  library("multcompView")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("table-s07-co2-equivalents")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))

# Settings ====

cumulative_file <- "table-s07-cumulative-co2-equivalents"
statistics_file <- "table-s07-ghg-statistics"


site_periods <- tribble(
  ~Site, ~PeriodStart,          ~PeriodEnd,
  "KL",  as.Date("2021-05-04"), as.Date("2021-11-10"),
  "ZB",  as.Date("2022-04-27"), as.Date("2022-11-03")
)

# GWPs
GWP_CH4 <- 29.8
GWP_N2O <- 273

# Molecular masses, g mol-1.
MOLAR_C <- 12.011
MOLAR_N2 <- 28.014
MOLAR_CH4 <- 16.043
MOLAR_N2O <- 44.013
MOLAR_CO2 <- 44.0095

# CH4 uptake
KL_CH4_SIGN <- -1
ZB_CH4_SIGN <- -1


MAXIMUM_BOUNDARY_DISTANCE_DAYS <- 3
alpha <- 0.05
treatments <- c("C", "C + N", "DRW", "DRW + N")
outcomes <- c("CO2", "CH4", "N2O", "NonCO2", "Total")

weights <- list(
  "Precipitation: RE - NP" = c(-0.5, -0.5, 0.5, 0.5),
  "Nitrogen: AddedN - NoN" = c(-0.5, 0.5, -0.5, 0.5),
  "Interaction: nitrogen effect RE - NP" = c(1, -1, -1, 1),
  "Precipitation without N: DRW - C" = c(-1, 0, 1, 0),
  "Precipitation with N: DRW + N - C + N" = c(0, -1, 0, 1)
)

treatment_kl <- c(
  KL01 = "C", KL05 = "C", KL13 = "C",
  KL02 = "C + N", KL06 = "C + N", KL09 = "C + N",
  KL14 = "C + N", KL03 = "DRW", KL07 = "DRW",
  KL11 = "DRW", KL15 = "DRW", KL04 = "DRW + N",
  KL08 = "DRW + N", KL12 = "DRW + N", KL16 = "DRW + N"
)

treatment_kl_co2 <- c(
  KL01 = "C", KL05 = "C", KL13 = "C",
  KL06 = "C + N", KL09 = "C + N", KL14 = "C + N",
  KL07 = "DRW", KL11 = "DRW", KL15 = "DRW",
  KL04 = "DRW + N", KL08 = "DRW + N",
  KL12 = "DRW + N", KL16 = "DRW + N"
)

treatment_zb <- c(
  ZB01 = "C", ZB05 = "C", ZB10 = "C", ZB13 = "C",
  ZB02 = "C + N", ZB06 = "C + N", ZB09 = "C + N",
  ZB14 = "C + N", ZB03 = "DRW", ZB07 = "DRW",
  ZB11 = "DRW", ZB15 = "DRW", ZB04 = "DRW + N",
  ZB08 = "DRW + N", ZB12 = "DRW + N", ZB16 = "DRW + N"
)

gas_specs <- list(
  KL_CH4 = list(
    site = "KL", gas = "CH4", filename = "ch4-klausen.csv",
    flux_column = "CH4", cumulative_column = "CH4_uptake_cum_ug_C_m2",
    prefix = "KL", sign = KL_CH4_SIGN, exclusions = "KL10",
    check_names = TRUE, treatment_map = treatment_kl
  ),
  KL_N2O = list(
    site = "KL", gas = "N2O", filename = "n2o-klausen.csv",
    flux_column = "N2O", cumulative_column = "N2O_cum_ug_N_m2",
    prefix = "KL", sign = 1, exclusions = "KL10",
    check_names = FALSE, treatment_map = treatment_kl
  ),
  KL_CO2 = list(
    site = "KL", gas = "CO2", filename = "co2-klausen.csv",
    flux_column = "CO2", cumulative_column = "CO2_cum_mg_C_m2",
    prefix = "KL", sign = 1,
    exclusions = c("KL02", "KL03", "KL10"),
    check_names = FALSE, treatment_map = treatment_kl_co2
  ),
  ZB_CO2 = list(
    site = "ZB", gas = "CO2", filename = "co2-zoe.csv",
    flux_column = "CO2", cumulative_column = "CO2_cum_mg_C_m2",
    prefix = "ZB", sign = 1, exclusions = character(),
    check_names = FALSE, treatment_map = treatment_zb
  ),
  ZB_CH4 = list(
    site = "ZB", gas = "CH4", filename = "ch4-zoe.csv",
    flux_column = "CH4", cumulative_column = "CH4_uptake_cum_ug_C_m2",
    prefix = "ZB", sign = ZB_CH4_SIGN, exclusions = character(),
    check_names = FALSE, treatment_map = treatment_zb
  )
)

site_fluxes <- list(
  KL = c("KL_CO2", "KL_CH4", "KL_N2O"),
  ZB = c("ZB_CO2", "ZB_CH4")
)


# Helpers ====

parse_datetime <- function(x) {
  cleaned <- trimws(as.character(x))
  cleaned <- gsub("/", ".", cleaned, fixed = TRUE)

  as.POSIXct(
    strptime(
      cleaned,
      format = "%d.%m.%Y %H:%M",
      tz = "UTC"
    ),
    tz = "UTC"
  )
}

standardize_plot_id <- function(x, prefix) {
  plot_number <- suppressWarnings(
    as.integer(
      stringr::str_extract(as.character(x), "[0-9]+")
    )
  )

  result <- rep(NA_character_, length(plot_number))
  valid <- !is.na(plot_number)

  result[valid] <- sprintf(
    "%s%02d",
    prefix,
    plot_number[valid]
  )

  result
}

prepare_daily_flux <- function(spec) {
  read.csv(
    file.path(dir_raw, spec$filename),
    check.names = spec$check_names,
    stringsAsFactors = FALSE
  ) |>
    pivot_longer(
      cols = 2:17,
      names_to = "Plot",
      values_to = spec$flux_column
    ) |>
    mutate(
      DateTime = parse_datetime(Date),
      SampleDate = as.Date(DateTime, tz = "UTC"),
      Plot = standardize_plot_id(Plot, spec$prefix),
      "{spec$flux_column}" := spec$sign *
        as.numeric(.data[[spec$flux_column]])
    ) |>
    filter(
      !is.na(DateTime),
      is.finite(.data[[spec$flux_column]]),
      !is.na(Plot),
      !Plot %in% spec$exclusions
    ) |>
    mutate(
      Treatment = unname(spec$treatment_map[Plot]),
      Drought = ifelse(Treatment %in% c("C", "C + N"), "NP", "RE"),
      Site = spec$site
    ) |>
    filter(!is.na(Treatment)) |>
    group_by(Site, Plot, Treatment, Drought, SampleDate) |>
    summarise(
      "{spec$flux_column}" := mean(.data[[spec$flux_column]]),
      .groups = "drop"
    )
}

get_site_period <- function(site) {
  selected <- site_periods |>
    filter(Site == site)

  if (nrow(selected) != 1L) {
    stop("No unique period was defined for site ", site)
  }

  selected
}

# Cumulative integration ====

calc_cumulative_hourly <- function(
  date,
  flux,
  period_start,
  period_end,
  maximum_boundary_distance_days =
    MAXIMUM_BOUNDARY_DISTANCE_DAYS
) {
  dat <- tibble(
    SampleDate = as.Date(date),
    Flux = as.numeric(flux)
  ) |>
    filter(
      !is.na(SampleDate),
      is.finite(Flux)
    ) |>
    group_by(SampleDate) |>
    summarise(
      Flux = mean(Flux),
      .groups = "drop"
    ) |>
    arrange(SampleDate)

  expected_elapsed_days <- as.numeric(
    period_end - period_start
  )

  expected_calendar_days <- expected_elapsed_days + 1

  observed_days_inside_period <- sum(
    dat$SampleDate >= period_start &
      dat$SampleDate <= period_end
  )

  coverage_percent <- (
    100 * observed_days_inside_period /
      expected_calendar_days
  )

  first_raw_date <- if (nrow(dat) > 0L) {
    min(dat$SampleDate)
  } else {
    as.Date(NA)
  }

  last_raw_date <- if (nrow(dat) > 0L) {
    max(dat$SampleDate)
  } else {
    as.Date(NA)
  }

  if (nrow(dat) < 2L) {
    return(tibble(
      Cumulative = NA_real_,
      FirstRawDate = first_raw_date,
      LastRawDate = last_raw_date,
      IntegratedDays = NA_real_,
      ObservedDaysInsidePeriod = observed_days_inside_period,
      ExpectedCalendarDays = expected_calendar_days,
      CoveragePercent = coverage_percent,
      MaximumGapDays = NA_real_,
      StartBoundaryDistanceDays = NA_real_,
      EndBoundaryDistanceDays = NA_real_,
      StartBoundaryMethod = NA_character_,
      EndBoundaryMethod = NA_character_,
      CompleteIntegrationPeriod = FALSE
    ))
  }

  date_numeric <- as.numeric(dat$SampleDate)
  start_numeric <- as.numeric(period_start)
  end_numeric <- as.numeric(period_end)

  nearest_start_distance <- min(
    abs(date_numeric - start_numeric)
  )

  nearest_end_distance <- min(
    abs(date_numeric - end_numeric)
  )

  interpolate_boundary <- function(target_date) {
    target_numeric <- as.numeric(target_date)

    exact_match <- which(
      dat$SampleDate == target_date
    )

    if (length(exact_match) > 0L) {
      return(list(
        value = dat$Flux[exact_match[1]],
        method = "Observed"
      ))
    }

    before <- which(
      dat$SampleDate < target_date
    )

    after <- which(
      dat$SampleDate > target_date
    )

    # Interpolate between the nearest observations
    # before and after the boundary.
    if (length(before) > 0L && length(after) > 0L) {
      before_row <- max(before)
      after_row <- min(after)

      interpolated <- stats::approx(
        x = as.numeric(
          dat$SampleDate[c(before_row, after_row)]
        ),
        y = dat$Flux[c(before_row, after_row)],
        xout = target_numeric
      )$y

      return(list(
        value = interpolated,
        method = "Interpolated"
      ))
    }

    # Outside the observed range -> nearest
    # observation only within the stated tolerance.
    nearest_row <- which.min(
      abs(date_numeric - target_numeric)
    )

    nearest_distance <- abs(
      date_numeric[nearest_row] - target_numeric
    )

    if (
      nearest_distance <= maximum_boundary_distance_days
    ) {
      return(list(
        value = dat$Flux[nearest_row],
        method = paste0(
          "Nearest observation: ",
          nearest_distance,
          " day(s)"
        )
      ))
    }

    list(
      value = NA_real_,
      method = "Unavailable"
    )
  }

  start_boundary <- interpolate_boundary(period_start)
  end_boundary <- interpolate_boundary(period_end)

  internal_data <- dat |>
    filter(
      SampleDate > period_start,
      SampleDate < period_end
    )

  integration_data <- bind_rows(
    tibble(
      SampleDate = period_start,
      Flux = start_boundary$value
    ),
    internal_data,
    tibble(
      SampleDate = period_end,
      Flux = end_boundary$value
    )
  ) |>
    filter(
      !is.na(SampleDate),
      is.finite(Flux)
    ) |>
    group_by(SampleDate) |>
    summarise(
      Flux = mean(Flux),
      .groups = "drop"
    ) |>
    arrange(SampleDate)

  complete_period <- (
    nrow(integration_data) >= 2L &&
      min(integration_data$SampleDate) == period_start &&
      max(integration_data$SampleDate) == period_end
  )

  if (!complete_period) {
    return(tibble(
      Cumulative = NA_real_,
      FirstRawDate = first_raw_date,
      LastRawDate = last_raw_date,
      IntegratedDays = NA_real_,
      ObservedDaysInsidePeriod = observed_days_inside_period,
      ExpectedCalendarDays = expected_calendar_days,
      CoveragePercent = coverage_percent,
      MaximumGapDays = NA_real_,
      StartBoundaryDistanceDays = nearest_start_distance,
      EndBoundaryDistanceDays = nearest_end_distance,
      StartBoundaryMethod = start_boundary$method,
      EndBoundaryMethod = end_boundary$method,
      CompleteIntegrationPeriod = FALSE
    ))
  }

  elapsed_days <- as.numeric(
    diff(integration_data$SampleDate)
  )

  interval_mean_flux <- (
    head(integration_data$Flux, -1) +
      tail(integration_data$Flux, -1)
  ) / 2

  cumulative <- sum(
    interval_mean_flux * elapsed_days * 24
  )

  tibble(
    Cumulative = cumulative,
    FirstRawDate = first_raw_date,
    LastRawDate = last_raw_date,
    IntegratedDays = as.numeric(
      max(integration_data$SampleDate) -
        min(integration_data$SampleDate)
    ),
    ObservedDaysInsidePeriod = observed_days_inside_period,
    ExpectedCalendarDays = expected_calendar_days,
    CoveragePercent = coverage_percent,
    MaximumGapDays = max(elapsed_days),
    StartBoundaryDistanceDays = nearest_start_distance,
    EndBoundaryDistanceDays = nearest_end_distance,
    StartBoundaryMethod = start_boundary$method,
    EndBoundaryMethod = end_boundary$method,
    CompleteIntegrationPeriod = TRUE
  )
}


# calculate cumulative total for each plot

integrate_plot_data <- function(
  data,
  flux_column,
  output_column,
  site
) {
  selected <- get_site_period(site)

  data |>
    group_by(
      Site,
      Plot,
      Treatment,
      Drought
    ) |>
    nest() |>
    mutate(
      Integration = map(
        data,
        function(plot_data) {
          calc_cumulative_hourly(
            date = plot_data$SampleDate,
            flux = plot_data[[flux_column]],
            period_start = selected$PeriodStart,
            period_end = selected$PeriodEnd
          )
        }
      )
    ) |>
    select(-data) |>
    unnest(Integration) |>
    ungroup() |>
    rename(
      !!output_column := Cumulative
    ) |>
    mutate(
      PrescribedStart = selected$PeriodStart,
      PrescribedEnd = selected$PeriodEnd,
      PrescribedElapsedDays = as.numeric(
        selected$PeriodEnd - selected$PeriodStart
      )
    )
}


# Gas preparation and integration ====

daily_fluxes <- map(gas_specs, prepare_daily_flux)

cumulative_fluxes <- imap(
  gas_specs,
  \(spec, key) {
    integrate_plot_data(
      daily_fluxes[[key]],
      spec$flux_column,
      spec$cumulative_column,
      spec$site
    )
  }
)


# Coverage audit ====

coverage_columns <- c(
  "Site",
  "Plot",
  "Treatment",
  "Drought",
  "FirstRawDate",
  "LastRawDate",
  "IntegratedDays",
  "ObservedDaysInsidePeriod",
  "ExpectedCalendarDays",
  "CoveragePercent",
  "MaximumGapDays",
  "StartBoundaryDistanceDays",
  "EndBoundaryDistanceDays",
  "StartBoundaryMethod",
  "EndBoundaryMethod",
  "CompleteIntegrationPeriod",
  "PrescribedStart",
  "PrescribedEnd",
  "PrescribedElapsedDays"
)

coverage_audit <- imap_dfr(
  cumulative_fluxes,
  \(dat, key) {
    dat |>
      select(all_of(coverage_columns)) |>
      mutate(Gas = gas_specs[[key]]$gas, .before = 1)
  }
) |>
  arrange(Site, Gas, Treatment, Plot)


# Merge plot-level results ====

join_columns <- c(
  "Site", "Plot", "Treatment", "Drought"
)

merge_site_fluxes <- function(site) {
  keys <- site_fluxes[[site]]
  selected <- map2(
    cumulative_fluxes[keys],
    gas_specs[keys],
    \(dat, spec) {
      dat |>
        select(all_of(join_columns), all_of(spec$cumulative_column))
    }
  )

  merged <- reduce(selected, full_join, by = join_columns)
  if (site == "ZB") {
    merged <- merged |> mutate(N2O_cum_ug_N_m2 = NA_real_)
  }
  merged
}

site_plots <- map(names(site_fluxes), merge_site_fluxes) |>
  set_names(names(site_fluxes))


# CO2-equivalent fluxes ====

PLOT_FINAL <- bind_rows(site_plots) |>
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("C", "C + N", "DRW", "DRW + N")
    ),

    # CO2:
    # mg CO2-C -> mg CO2 -> t CO2 equivalent ha-1.
    # CO2 has a GWP of 1.
    #
    # 1 mg m-2 = 1e-5 t ha-1.
    CO2_t_CO2eq_ha =
      CO2_cum_mg_C_m2 *
        (MOLAR_CO2 / MOLAR_C) /
        1e5,

    # CH4:
    # µg CH4-C -> µg CH4 -> µg CO2 equivalent
    # -> t CO2 equivalent ha-1.
    #
    # Positive standardized uptake becomes a negative
    # climate contribution.
    #
    # 1 µg m-2 = 1e-8 t ha-1.
    CH4_t_CO2eq_ha =
      -CH4_uptake_cum_ug_C_m2 *
        (MOLAR_CH4 / MOLAR_C) *
        GWP_CH4 /
        1e8,

    # N2O:
    # µg N2O-N -> µg N2O -> µg CO2 equivalent
    # -> t CO2 equivalent ha-1.
    #
    # MOLAR_N2 represents the mass of both nitrogen atoms.
    N2O_t_CO2eq_ha =
      N2O_cum_ug_N_m2 *
        (MOLAR_N2O / MOLAR_N2) *
        GWP_N2O /
        1e8,

    # KL includes CO2, CH4 and N2O.
    # ZB includes CO2 and CH4 because N2O was not measured.
    #
    # Missing required contributions propagate to the total.
    Total_t_CO2eq_ha = case_when(
      Site == "KL" ~
        CO2_t_CO2eq_ha +
        CH4_t_CO2eq_ha +
        N2O_t_CO2eq_ha,
      Site == "ZB" ~
        CO2_t_CO2eq_ha +
        CH4_t_CO2eq_ha,
      TRUE ~ NA_real_
    ),
    GasesIncluded = case_when(
      Site == "KL" ~ "CO2 + CH4 + N2O",
      Site == "ZB" ~ "CO2 + CH4",
      TRUE ~ NA_character_
    )
  ) |>
  arrange(Site, Treatment, Plot)


# Treatment summaries ====

Treatment_summary <- PLOT_FINAL |>
  group_by(
    Site,
    Treatment,
    GasesIncluded
  ) |>
  summarise(
    Nplots_CO2 = sum(is.finite(CO2_t_CO2eq_ha)),
    CO2_mean_t_CO2eq_ha =
      finite_mean(CO2_t_CO2eq_ha),
    CO2_se_t_CO2eq_ha =
      standard_error(CO2_t_CO2eq_ha),
    Nplots_CH4 = sum(is.finite(CH4_t_CO2eq_ha)),
    CH4_mean_t_CO2eq_ha =
      finite_mean(CH4_t_CO2eq_ha),
    CH4_se_t_CO2eq_ha =
      standard_error(CH4_t_CO2eq_ha),
    Nplots_N2O = sum(is.finite(N2O_t_CO2eq_ha)),
    N2O_mean_t_CO2eq_ha =
      finite_mean(N2O_t_CO2eq_ha),
    N2O_se_t_CO2eq_ha =
      standard_error(N2O_t_CO2eq_ha),
    Nplots_Total = sum(is.finite(Total_t_CO2eq_ha)),
    Total_mean_t_CO2eq_ha =
      finite_mean(Total_t_CO2eq_ha),
    Total_se_t_CO2eq_ha =
      standard_error(Total_t_CO2eq_ha),
    .groups = "drop"
  )

Treatment_summary_kg <- Treatment_summary |>
  mutate(
    across(
      ends_with("_t_CO2eq_ha"),
      ~ .x * 1000
    )
  ) |>
  rename_with(
    ~ sub("_t_CO2eq_ha$", "_kg_CO2eq_ha", .x),
    ends_with("_t_CO2eq_ha")
  )

# Pooled water summaries ====

Drought_summary <- PLOT_FINAL |>
  group_by(
    Site,
    Drought,
    GasesIncluded
  ) |>
  summarise(
    Nplots_CO2 = sum(is.finite(CO2_t_CO2eq_ha)),
    CO2_mean_t_CO2eq_ha =
      finite_mean(CO2_t_CO2eq_ha),
    CO2_se_t_CO2eq_ha =
      standard_error(CO2_t_CO2eq_ha),
    Nplots_CH4 = sum(is.finite(CH4_t_CO2eq_ha)),
    CH4_mean_t_CO2eq_ha =
      finite_mean(CH4_t_CO2eq_ha),
    CH4_se_t_CO2eq_ha =
      standard_error(CH4_t_CO2eq_ha),
    Nplots_N2O = sum(is.finite(N2O_t_CO2eq_ha)),
    N2O_mean_t_CO2eq_ha =
      finite_mean(N2O_t_CO2eq_ha),
    N2O_se_t_CO2eq_ha =
      standard_error(N2O_t_CO2eq_ha),
    Nplots_Total = sum(is.finite(Total_t_CO2eq_ha)),
    Total_mean_t_CO2eq_ha =
      finite_mean(Total_t_CO2eq_ha),
    Total_se_t_CO2eq_ha =
      standard_error(Total_t_CO2eq_ha),
    .groups = "drop"
  )

Drought_summary_kg <- Drought_summary |>
  mutate(
    across(
      ends_with("_t_CO2eq_ha"),
      ~ .x * 1000
    )
  ) |>
  rename_with(
    ~ sub("_t_CO2eq_ha$", "_kg_CO2eq_ha", .x),
    ends_with("_t_CO2eq_ha")
  )
# Period and boundary checks ====

period_check <- site_periods |>
  mutate(
    ElapsedDays = as.numeric(
      PeriodEnd - PeriodStart
    ),
    IncludedCalendarDays = ElapsedDays + 1
  )

if (any(
  !is.finite(period_check$ElapsedDays) |
    period_check$ElapsedDays != 190
)) {
  stop("Each site integration period must span 190 elapsed days.")
}

incomplete_boundaries <- coverage_audit |>
  filter(!CompleteIntegrationPeriod)

if (nrow(incomplete_boundaries) > 0L) {
  warning(
    nrow(incomplete_boundaries),
    " plot-gas series could not cover the complete 190-day period. ",
    "Their cumulative values were returned as NA. ",
    "Inspect the incomplete-boundaries CSV."
  )
}


# Export cumulative results ====

ensure_dir(dir_reporting)
write_csv_sheets(
  list(
    "Plot_level" = PLOT_FINAL,
    "Treatment_summary" = Treatment_summary,
    "Treatment_summary_kg" = Treatment_summary_kg,
    "Drought_summary" = Drought_summary,
    "Drought_summary_kg" = Drought_summary_kg,
    "Period_check" = period_check,
    "Coverage_audit" = coverage_audit,
    "Incomplete_boundaries" = incomplete_boundaries
  ),
  cumulative_file
)

logger::log_info(paste("Cumulative CO2-equivalent CSV files saved under:", dir_reporting))


# Statistical summaries ====

decision <- function(p) {
  case_when(
    !is.finite(p) ~ "Not tested",
    p < alpha ~ "Different",
    TRUE ~ "No difference detected"
  )
}

## Input and summaries ----

required <- c(
  "Site", "Plot", "Treatment", "Drought",
  "CO2_t_CO2eq_ha", "CH4_t_CO2eq_ha",
  "N2O_t_CO2eq_ha", "Total_t_CO2eq_ha"
)
if (!exists("PLOT_FINAL")) stop("PLOT_FINAL is required.")
if (!all(required %in% names(PLOT_FINAL))) {
  stop("Required columns are missing.")
}

plot_data <- PLOT_FINAL |>
  transmute(
    Site = as.character(Site),
    Plot = as.character(Plot),
    Treatment = factor(as.character(Treatment), levels = treatments),
    Drought = as.character(Drought),
    CO2 = CO2_t_CO2eq_ha * 1000,
    CH4 = CH4_t_CO2eq_ha * 1000,
    N2O = N2O_t_CO2eq_ha * 1000,
    NonCO2 = (CH4_t_CO2eq_ha + N2O_t_CO2eq_ha) * 1000,
    Total = Total_t_CO2eq_ha * 1000
  )

if (anyNA(plot_data[c("Site", "Plot", "Treatment", "Drought")])) {
  stop("Missing or unrecognised identifiers.")
}
if (anyDuplicated(plot_data[c("Site", "Plot")])) {
  stop("Duplicate Site/Plot observations.")
}
if (any(
  plot_data$Drought !=
    ifelse(plot_data$Treatment %in% c("C", "C + N"), "NP", "RE")
)) {
  stop("Treatment and precipitation labels disagree.")
}


long_data <- plot_data |>
  pivot_longer(all_of(outcomes), names_to = "Outcome", values_to = "Value")

cell_summary <- long_data |>
  group_by(Site, Outcome, Treatment) |>
  summarise(
    n = sum(is.finite(Value)),
    mean = {
      x <- Value[is.finite(Value)]
      if (length(x)) mean(x) else NA_real_
    },
    sd = {
      x <- Value[is.finite(Value)]
      if (length(x) >= 2L) sd(x) else NA_real_
    },
    .groups = "drop"
  ) |>
  mutate(
    Treatment = as.character(Treatment),
    se = sd / sqrt(n)
  )

status <- cell_summary |>
  group_by(Site, Outcome) |>
  summarise(
    Eligible = n() == 4L &&
      all(n >= 3L) && all(is.finite(sd)) && all(sd > 0),
    .groups = "drop"
  )

eligible <- status |> filter(Eligible)
if (!nrow(eligible)) stop("No eligible outcomes.")

## Diagnostics ----

analysis_data <- long_data |>
  semi_join(eligible, by = c("Site", "Outcome")) |>
  filter(is.finite(Value))

normality <- analysis_data |>
  group_by(Site, Outcome, Treatment) |>
  summarise(
    n = n(),
    Shapiro_p = if (n() <= 5000L) {
      shapiro.test(Value)$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  )

outlier_review <- analysis_data |>
  group_by(Site, Outcome, Treatment) |>
  mutate(
    lower = quantile(Value, 0.25) - 1.5 * IQR(Value),
    upper = quantile(Value, 0.75) + 1.5 * IQR(Value),
    Flagged = Value < lower | Value > upper
  ) |>
  ungroup()

variance_checks <- map_dfr(seq_len(nrow(eligible)), \(i) {
  site <- eligible$Site[i]
  outcome <- eligible$Outcome[i]
  data <- analysis_data |> filter(Site == site, Outcome == outcome)
  levene <- leveneTest(
    Value ~ Treatment,
    data = data,
    center = median
  )

  tibble(
    Site = site,
    Outcome = outcome,
    Levene_p = as.numeric(levene[1, "Pr(>F)"])
  )
})
## Independent-cell contrasts ----

# For weights c_i and v_i = s_i^2/n_i:
# Var(contrast) = sum(c_i^2 v_i)
# df = Var(contrast)^2 / sum[(c_i^2 v_i)^2/(n_i-1)]
contrast_test <- function(s, w, name) {
  parts <- w^2 * s$sd^2 / s$n
  variance <- sum(parts)
  se <- sqrt(variance)
  df <- variance^2 / sum(parts^2 / (s$n - 1))
  estimate <- sum(w * s$mean)
  t <- estimate / se
  p <- 2 * pt(abs(t), df, lower.tail = FALSE)
  margin <- qt(1 - alpha / 2, df) * se

  tibble(
    Comparison = name, estimate = estimate, se = se, df = df,
    t = t, lower = estimate - margin, upper = estimate + margin,
    p = p
  )
}

welch_list <- gh_list <- contrast_list <- letter_list <- list()

for (i in seq_len(nrow(eligible))) {
  site <- eligible$Site[i]
  outcome <- eligible$Outcome[i]
  d <- analysis_data |> filter(Site == site, Outcome == outcome)
  s <- cell_summary |>
    filter(Site == site, Outcome == outcome) |>
    arrange(match(Treatment, treatments))

  tag <- function(x) {
    mutate(x, Site = site, Outcome = outcome, .before = 1)
  }

  omnibus <- oneway.test(Value ~ Treatment, data = d, var.equal = FALSE)
  welch_list[[i]] <- tag(tibble(
    F = unname(omnibus$statistic),
    df1 = unname(omnibus$parameter[1]),
    df2 = unname(omnibus$parameter[2]),
    p = omnibus$p.value
  ))

  # Games-Howell: q = sqrt(2)*abs(difference)/SE.
  # Studentized-range probabilities account for four groups.
  pairs <- combn(4, 2)
  gh <- map_dfr(seq_len(ncol(pairs)), function(j) {
    a <- pairs[1, j]
    b <- pairs[2, j]
    w <- numeric(4)
    w[a] <- 1
    w[b] <- -1

    x <- contrast_test(s, w, paste(treatments[a], "-", treatments[b]))
    q <- sqrt(2) * abs(x$estimate) / x$se
    margin <- qtukey(1 - alpha, 4, x$df) * x$se / sqrt(2)

    tibble(
      Treatment_1 = treatments[a], Treatment_2 = treatments[b],
      estimate = x$estimate, se = x$se, df = x$df,
      lower = x$estimate - margin, upper = x$estimate + margin,
      p_adjusted = ptukey(q, 4, x$df, lower.tail = FALSE)
    )
  })
  gh_list[[i]] <- tag(gh)

  p_matrix <- matrix(
    1, 4, 4,
    dimnames = list(treatments, treatments)
  )
  for (j in seq_len(nrow(gh))) {
    p_matrix[gh$Treatment_1[j], gh$Treatment_2[j]] <- gh$p_adjusted[j]
    p_matrix[gh$Treatment_2[j], gh$Treatment_1[j]] <- gh$p_adjusted[j]
  }
  cld <- if (all(is.finite(gh$p_adjusted))) {
    multcompLetters(
      p_matrix,
      threshold = alpha, Letters = letters
    )$Letters[treatments]
  } else {
    rep(NA_character_, 4)
  }

  letter_list[[i]] <- tag(tibble(
    Treatment = treatments, Letters = unname(cld)
  ))

  contrast_list[[i]] <- tag(
    imap_dfr(weights, ~ contrast_test(s, .x, .y))
  )
}

welch_results <- bind_rows(welch_list) |>
  mutate(Decision = decision(p))
games_howell <- bind_rows(gh_list) |>
  mutate(Decision = decision(p_adjusted))
treatment_letters <- bind_rows(letter_list)
contrasts <- bind_rows(contrast_list)

## Precipitation results ----

precipitation_results <- contrasts |>
  filter(Comparison == "Precipitation: RE - NP") |>
  mutate(Decision = decision(p))

within_n <- contrasts |>
  filter(Comparison %in% names(weights)[4:5]) |>
  group_by(Site, Outcome) |>
  mutate(
    p_adjusted = p.adjust(p, "holm"),
    Decision = decision(p_adjusted)
  ) |>
  ungroup()
# lower/upper above are POINTWISE intervals, not Holm-adjusted.

factorial_results <- contrasts |>
  filter(Comparison %in% names(weights)[1:3]) |>
  mutate(Decision = decision(p))

precipitation_summary <- cell_summary |>
  mutate(Group = if_else(Treatment %in% treatments[1:2], "NP", "RE")) |>
  group_by(Site, Outcome, Group) |>
  summarise(
    n = sum(n),
    mean = if (all(is.finite(mean))) mean(mean) else NA_real_,
    se = if (all(is.finite(se))) sqrt(sum(se^2)) / 2 else NA_real_,
    .groups = "drop"
  ) |>
  left_join(
    precipitation_results |> select(Site, Outcome, p),
    by = c("Site", "Outcome")
  ) |>
  mutate(
    Letters = case_when(
      !is.finite(p) ~ NA_character_,
      p >= alpha ~ "A",
      Group == "NP" ~ "A",
      TRUE ~ "B"
    )
  )

## Presentation table and export ----

treatment_summary <- cell_summary |>
  left_join(treatment_letters, by = c("Site", "Outcome", "Treatment"))

presentation <- bind_rows(
  treatment_summary |>
    transmute(
      Summary = "Treatment", Site, Group = Treatment,
      Outcome, mean, se, Letters
    ),
  precipitation_summary |>
    transmute(
      Summary = "Precipitation: equal N weighting",
      Site, Group, Outcome, mean, se, Letters
    )
)

format_two <- function(x) {
  y <- rep("—", length(x))
  ok <- is.finite(x)
  z <- round(x[ok], 2)
  z[z == 0] <- 0
  y[ok] <- sprintf("%.2f", z)
  y
}

full_table <- presentation |>
  mutate(
    Mean = format_two(mean),
    SE = format_two(se),
    Mean = if_else(
      is.finite(mean) & !is.na(Letters),
      paste(Mean, Letters), Mean
    )
  ) |>
  select(Summary, Site, Group, Outcome, Mean, SE) |>
  pivot_wider(
    names_from = Outcome, values_from = c(Mean, SE),
    names_glue = "{Outcome}_{.value}"
  ) |>
  arrange(
    match(Summary, c("Treatment", "Precipitation: equal N weighting")),
    Site, match(Group, c(treatments, "NP", "RE"))
  ) |>
  select(
    Summary, Site, Group,
    all_of(as.vector(rbind(
      paste0(outcomes, "_Mean"),
      paste0(outcomes, "_SE")
    )))
  )

notes <- tribble(
  ~Topic, ~Explanation,
  "Units", "Estimates, SEs and intervals are in kg CO2-eq ha-1.",
  "Design", "Analysis assumes independent plot-level observations.",
  "Placement", "Shelter locations were chosen for practical reasons, not randomly assigned.",
  "Inference", "Welch-Satterthwaite and Games-Howell are approximate small-sample procedures.",
  "Diagnostics", "Shapiro and Levene are diagnostic only; they do not certify assumptions or trigger test selection.",
  "Outliers", "Flags are descriptive; all QC-approved finite observations are retained.",
  "Eligibility", "At least three observations and positive SD per treatment are required as an operational safeguard.",
  "Treatment letters", "Lowercase letters use Games-Howell adjusted p-values within each site/outcome.",
  "Precipitation", "NP and RE means weight nitrogen levels equally; SEs are propagated from independent cell means.",
  "Precipitation letters", "Uppercase letters use the unadjusted equal-weighted precipitation contrast.",
  "Multiplicity", "Games-Howell covers six treatment pairs; Holm covers two within-nitrogen precipitation contrasts. Other contrasts and cross-outcome/site testing are unadjusted.",
  "Intervals", "Games-Howell intervals use studentized-range adjustment. Other intervals are pointwise 95% intervals, including those displayed beside Holm-adjusted p-values.",
  "Missingness", "Inference uses available finite plot totals and does not adjust for informative exclusions or missingness.",
  "Integration", "Uncertainty from flux interpolation/integration is not propagated.",
  "NonCO2", "CH4 plus N2O is calculated within plots; unavailable at ZB.",
  "Total", "KL includes CO2, CH4 and N2O; ZB includes CO2 and CH4.",
  "Letters", "Shared letters indicate no detected difference, not equivalence. Compare only within site, outcome and comparison family.",
  "Formatting", "Full_table uses two-decimal presentation text; other sheets retain numeric precision."
)

ensure_dir(dir_reporting)
write_csv_sheets(
  list(
    Full_table = full_table,
    Treatment_summary = treatment_summary,
    Precipitation_summary = precipitation_summary,
    Welch_ANOVA = welch_results,
    Games_Howell = games_howell,
    Precipitation_primary = precipitation_results,
    Factorial_contrasts = factorial_results,
    Precipitation_within_N = within_n,
    Normality = normality,
    Variance_checks = variance_checks,
    Outlier_review = outlier_review,
    Eligibility = status,
    Plot_values = plot_data,
    Notes = notes
  ),
  statistics_file
)

logger::log_info(paste("Statistical results CSV files saved under:", dir_reporting))
log_script_end("table-s07-co2-equivalents")
