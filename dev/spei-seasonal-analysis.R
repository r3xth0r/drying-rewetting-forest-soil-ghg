# SPEI seasonal analysis ----
# Extract May–October SPEI at KL and ZB and compare study years to baselines.
# Inputs: annual SPEI30 NetCDF files.
# Outputs: processed summaries and a manuscript report table.

suppressPackageStartupMessages({
  library("tidyverse")
  library("ncdf4")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("spei-seasonal-analysis")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))

# Settings ====

years <- 2016:2022
season_months <- 5:10
input_prefix <- "winfore2-spei30-"
spei_variable <- "SPEI30"
expected_dimensions <- c("x", "y", "time")
time_origin <- as.Date("1961-01-01")
expected_time_units <- paste0(
  "^days since 1961-01-01(?:[ T]00:00:00(?:\\.0+)?)?",
  "(?: ?(?:UTC|Z))?$"
)
supported_calendars <- c("standard", "gregorian", "proleptic_gregorian")

# Coordinates use the same projected CRS as the NetCDF x/y grid.
locations <- tibble(
  Site = c("KL", "ZB"),
  Location = c("Klausen-Leopoldsdorf", "Zöbelboden"),
  x = c(601893.2749875868, 483116.63519076846),
  y = c(472482.32273104903, 438633.5729666162),
  Study_year = c(2021L, 2022L),
  Baseline_start = c(2016L, 2017L),
  Baseline_end = c(2020L, 2021L)
)

output_monthly <- "spei-monthly-summary.csv"
output_seasonal <- "spei-seasonal-means.csv"
output_comparison <- "spei-seasonal-comparison.csv"
output_manuscript <- "spei-seasonal-comparison-manuscript.csv"

# NetCDF reading ====

read_spei_year <- function(year) {
  filename <- file.path(
    dir_raw,
    paste0(input_prefix, year, ".nc")
  )

  if (!file.exists(filename)) {
    stop("File not found: ", filename)
  }

  nc <- nc_open(filename)
  on.exit(nc_close(nc), add = TRUE)

  if (!spei_variable %in% names(nc$var)) {
    stop(spei_variable, " is missing from: ", filename)
  }

  dimension_names <- vapply(
    nc$var[[spei_variable]]$dim,
    \(dimension) dimension$name,
    character(1)
  )

  if (!identical(dimension_names, expected_dimensions)) {
    stop(
      "Expected SPEI30 dimensions x, y, time; found: ",
      paste(dimension_names, collapse = ", "),
      " in ", filename
    )
  }

  x_grid <- ncvar_get(nc, "x")
  y_grid <- ncvar_get(nc, "y")
  time_values <- ncvar_get(nc, "time")

  time_units <- ncatt_get(nc, "time", "units")
  calendar <- ncatt_get(nc, "time", "calendar")

  if (
    !time_units$hasatt ||
      !grepl(
        expected_time_units,
        trimws(time_units$value),
        perl = TRUE
      )
  ) {
    stop(
      "Unexpected time units in ", filename, ": ",
      if (time_units$hasatt) time_units$value else "missing",
      ". Check the time conversion before proceeding."
    )
  }

  if (calendar$hasatt && !calendar$value %in% supported_calendars) {
    stop("Unsupported calendar in ", filename, ": ", calendar$value)
  }

  dates <- time_origin + time_values
  months <- as.integer(format(dates, "%m"))
  date_years <- as.integer(format(dates, "%Y"))
  keep <- date_years == year & months %in% season_months

  if (!any(keep)) {
    stop("No May–October dates found for ", year)
  }

  map_dfr(seq_len(nrow(locations)), \(i) {
    location <- locations[i, ]
    x_index <- which.min(abs(x_grid - location$x))
    y_index <- which.min(abs(y_grid - location$y))

    values <- as.numeric(ncvar_get(
      nc,
      spei_variable,
      start = c(x_index, y_index, 1),
      count = c(1, 1, -1),
      collapse_degen = TRUE
    ))

    if (length(values) != length(dates)) {
      stop("SPEI and time lengths differ in ", filename)
    }

    tibble(
      Site = location$Site,
      Location = location$Location,
      Date = dates[keep],
      Year = date_years[keep],
      Month = months[keep],
      SPEI = values[keep],
      Grid_x = x_grid[x_index],
      Grid_y = y_grid[y_index]
    )
  })
}

# Extract May–October SPEI ====

spei_data <- map_dfr(years, read_spei_year)

if (anyDuplicated(spei_data[c("Site", "Date")])) {
  stop("Duplicate dates found within a site.")
}

cat("\nSelected grid cells:\n")
print(distinct(spei_data, Site, Location, Grid_x, Grid_y))

# Monthly summaries ====

monthly_summary <- spei_data |>
  group_by(Site, Location, Year, Month) |>
  summarise(
    Time_steps = n(),
    Valid_values = sum(is.finite(SPEI)),
    Mean_SPEI = finite_mean(SPEI),
    SD_SPEI = finite_sd(SPEI),
    .groups = "drop"
  )

print(monthly_summary, n = Inf)

# Seasonal means ====

seasonal_spei <- monthly_summary |>
  group_by(Site, Location, Year) |>
  summarise(
    Months_available = sum(is.finite(Mean_SPEI)),
    Seasonal_mean = if (
      n() == length(season_months) && all(is.finite(Mean_SPEI))
    ) {
      mean(Mean_SPEI)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

print(seasonal_spei, n = Inf)

# Check required seasons ====

required_seasons <- map_dfr(seq_len(nrow(locations)), \(i) {
  tibble(
    Site = locations$Site[i],
    Year = c(
      seq.int(locations$Baseline_start[i], locations$Baseline_end[i]),
      locations$Study_year[i]
    )
  )
})

checked_seasons <- required_seasons |>
  left_join(seasonal_spei, by = c("Site", "Year"))

if (any(!is.finite(checked_seasons$Seasonal_mean))) {
  checked_seasons |>
    filter(!is.finite(Seasonal_mean)) |>
    print()

  stop("Missing or incomplete May–October seasons.")
}

if (any(monthly_summary$Valid_values < monthly_summary$Time_steps)) {
  warning(
    "Some months contain missing SPEI values. ",
    "Inspect Valid_values and Time_steps in the monthly summary."
  )
}

# Study years versus preceding five years ====

seasonal_comparison <- map_dfr(seq_len(nrow(locations)), \(i) {
  location <- locations[i, ]
  site_data <- seasonal_spei |>
    filter(Site == location$Site)

  baseline_values <- site_data |>
    filter(
      Year >= location$Baseline_start,
      Year <= location$Baseline_end
    ) |>
    pull(Seasonal_mean)

  study_value <- site_data |>
    filter(Year == location$Study_year) |>
    pull(Seasonal_mean)

  tibble(
    Site = location$Site,
    Location = location$Location,
    Measurement_year = location$Study_year,
    Measurement_SPEI = study_value,
    Baseline_period = paste(
      location$Baseline_start,
      location$Baseline_end,
      sep = "–"
    ),
    Baseline_years = length(baseline_values),
    Baseline_mean = mean(baseline_values),
    Baseline_SD = sd(baseline_values),
    Difference_from_baseline = study_value - mean(baseline_values)
  )
})

cat("\nDescriptive seasonal comparisons:\n")
print(seasonal_comparison, width = Inf)

# Manuscript table ====

spei_manuscript_table <- seasonal_comparison |>
  transmute(
    Site,
    Season = "May–October",
    Measurement_year,
    Measurement_SPEI = sprintf("%.2f", Measurement_SPEI),
    Previous_five_years = Baseline_period,
    Baseline_mean_plus_minus_SD = sprintf(
      "%.2f ± %.2f",
      Baseline_mean,
      Baseline_SD
    ),
    Difference_from_baseline = sprintf("%.2f", Difference_from_baseline)
  )

print(as.data.frame(spei_manuscript_table), row.names = FALSE)

# Export results ====

write_table(monthly_summary, output_monthly, dir_processed)
write_table(seasonal_spei, output_seasonal, dir_processed)
write_table(seasonal_comparison, output_comparison, dir_processed)
write_table(spei_manuscript_table, output_manuscript, dir_reporting)

# Manuscript summaries: mean ± descriptive SE across monthly means ====

spei_reporting <- map_dfr(seq_len(nrow(locations)), \(i) {
  settings <- locations[i, ]

  monthly_summary |>
    filter(
      Site == settings$Site,
      Year == settings$Study_year |
        between(Year, settings$Baseline_start, settings$Baseline_end)
    ) |>
    mutate(
      Period = if_else(
        Year == settings$Study_year,
        as.character(settings$Study_year),
        paste(settings$Baseline_start, settings$Baseline_end, sep = "–")
      )
    ) |>
    group_by(Site, Period) |>
    summarise(
      n_months = sum(is.finite(Mean_SPEI)),
      Mean = finite_mean(Mean_SPEI),
      SE = finite_sd(Mean_SPEI) / sqrt(n_months),
      .groups = "drop"
    ) |>
    mutate(Mean_SE = sprintf("%.2f ± %.2f", Mean, SE))
})

print(spei_reporting)
log_script_end("spei-seasonal-analysis")
