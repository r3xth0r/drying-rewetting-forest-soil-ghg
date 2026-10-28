# Project helpers ----
# Functions shared by all scripts. Requires config.R, dplyr, and ggplot2.

# Input ====

read_input <- function(filename, check.names = TRUE) {
  path <- file.path(dir_raw, filename)
  if (!file.exists(path)) stop("Input file not found: ", path)
  read.csv(
    path,
    check.names = check.names,
    stringsAsFactors = FALSE
  )
}

require_columns <- function(dat, columns, filename) {
  absent <- setdiff(columns, names(dat))
  if (length(absent)) {
    stop(filename, ": missing columns: ", paste(absent, collapse = ", "))
  }
  invisible(dat)
}

# Accepts dd.mm.yyyy HH:MM or dd/mm/yyyy HH:MM
parse_timestamp <- function(x) {
  timestamp <- as.POSIXct(
    gsub("/", ".", trimws(as.character(x))),
    format = "%d.%m.%Y %H:%M",
    tz = data_timezone
  )
  if (anyNA(timestamp)) {
    stop("Unparseable timestamps. Expected dd.mm.yyyy HH:MM.")
  }
  timestamp
}

# Accepts dd.mm.yyyy or dd/mm/yyyy, with optional time part
parse_date <- function(x) {
  dates <- as.Date(
    substr(gsub("/", ".", trimws(as.character(x))), 1, 10),
    format = "%d.%m.%Y"
  )
  if (anyNA(dates)) {
    stop("Unparseable dates: expected dd/mm/yyyy or dd.mm.yyyy.")
  }
  dates
}

# Time of day bin; midnight itself is excluded
time_category <- function(time) {
  dplyr::case_when(
    time > "00:00:00" & time <= "06:00:00" ~ "AM1",
    time > "06:00:00" & time <= "12:00:00" ~ "AM2",
    time > "12:00:00" & time <= "18:00:00" ~ "PM1",
    time > "18:00:00" & time <= "23:59:59" ~ "PM2"
  )
}

# Adds SampleDate, Time (HH:MM:SS) and TimeCat from a timestamp column
add_time <- function(dat, column) {
  require_columns(dat, column, "add_time")
  timestamp <- parse_timestamp(dat[[column]])
  dat$SampleDate <- as.Date(timestamp, tz = data_timezone)
  dat$Time <- format(timestamp, "%H:%M:%S", tz = data_timezone)
  dat$TimeCat <- time_category(dat$Time)
  dat
}

prepare_sensors <- function(filename, date_column, sensor_map, value_name,
                            equal_sensor_weight = FALSE,
                            require_numeric = FALSE) {
  dat <- read_input(filename)
  require_columns(dat, sensor_map$Sensor, filename)

  if (
    require_numeric &&
      !all(vapply(dat[sensor_map$Sensor], is.numeric, logical(1)))
  ) {
    stop(filename, ": selected sensor columns must be numeric.")
  }

  dat <- dat |>
    add_time(date_column) |>
    dplyr::select(SampleDate, TimeCat, dplyr::all_of(sensor_map$Sensor)) |>
    tidyr::pivot_longer(
      dplyr::all_of(sensor_map$Sensor),
      names_to = "Sensor",
      values_to = "Value"
    ) |>
    dplyr::left_join(sensor_map, by = "Sensor") |>
    dplyr::filter(!is.na(TimeCat), is.finite(Value))

  if (equal_sensor_weight) {
    dat <- dat |>
      dplyr::group_by(
        SampleDate, TimeCat, Depth_cm, Treatment, Sensor
      ) |>
      dplyr::summarise(Value = mean(Value), .groups = "drop")
  }

  dat |>
    dplyr::group_by(SampleDate, TimeCat, Depth_cm, Treatment) |>
    dplyr::summarise(Value = mean(Value), .groups = "drop") |>
    dplyr::rename(dplyr::all_of(setNames("Value", value_name)))
}

# Summaries ====

finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

finite_sum <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) sum(x) else NA_real_
}

finite_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) sd(x) else NA_real_
}

standard_error <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_
}

# Output ====

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

mm_to_inches <- function(mm) mm / 25.4

# Writes PNG by default; set save_pdf_figures to also write PDF.
save_figure <- function(
  plot, name, width_mm, height_mm, dir = dir_plots,
  include_pdf = save_pdf_figures
) {
  ensure_dir(dir)
  ggplot2::ggsave(
    file.path(dir, paste0(name, ".png")),
    plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = fig_dpi, bg = "white"
  )
  if (include_pdf) {
    ggplot2::ggsave(
      file.path(dir, paste0(name, ".pdf")),
      plot,
      width = width_mm, height = height_mm, units = "mm",
      device = if (capabilities("cairo")) {
        grDevices::cairo_pdf
      } else {
        grDevices::pdf
      },
      bg = "white"
    )
  }
  logger::log_info(paste("Figure saved:", file.path(dir, name)))
}

write_table <- function(data, name, dir) {
  ensure_dir(dir)
  path <- file.path(dir, name)
  readr::write_csv(data, path)
  logger::log_info(paste("Table saved:", path))
  invisible(path)
}

write_csv_sheets <- function(tables, prefix, dir = dir_reporting) {
  purrr::iwalk(tables, \(data, name) {
    suffix <- tolower(gsub("[^A-Za-z0-9]+", "-", name))
    write_table(data, paste0(prefix, "-", suffix, ".csv"), dir)
  })
}

write_session_info <- function(dir = dir_interim, name = "session-info.txt") {
  ensure_dir(dir)
  path <- file.path(dir, name)
  writeLines(capture.output(sessionInfo()), path)
  logger::log_info(paste("Session information saved:", path))
}
