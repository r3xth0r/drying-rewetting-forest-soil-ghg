# Figure 3: Daily time series ----
# KL and ZB daily precipitation, WFPS, soil temperature and GHG fluxes.
# Inputs: dat/raw flux, moisture, temperature and precipitation files.
# Output: plt/figure-03-gas-fluxes.{png,pdf}

suppressPackageStartupMessages({
  library("tidyverse")
  library("patchwork")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("figure-03-gas-fluxes")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))

# Settings ====

figure_name <- "figure-03-gas-fluxes"
figure_width_mm <- 260
figure_height_mm <- 220

# Plots excluded from daily means
excluded_plots <- list(
  KL_CO2 = c(2L, 3L, 10L),
  KL_CH4 = 10L,
  KL_N2O = 10L,
  ZB_CO2 = integer(),
  ZB_CH4 = integer()
)

# CH4 uptake is plotted as positive values
ch4_multiplier <- -1

# ZB CH4 records up to this date are not used
zb_ch4_start <- as.Date("2022-04-26")

# KL moisture records from this date on, and ZB up to this date, are not used
kl_moisture_end <- as.Date("2021-11-22")
zb_moisture_start <- as.Date("2022-04-26")

# Sensor to environment (NP/RE) allocation
moisture_maps <- list(
  KL = c(C1_5cm = "NP", C2_5cm = "NP", D1_5cm = "RE", D2_5cm = "RE"),
  ZB = c(C1 = "NP", S1 = "RE", S2 = "RE", S4 = "RE")
)

temperature_maps <- list(
  KL = c(C1 = "NP", C2 = "NP", D1 = "RE", D2 = "RE"),
  ZB = c(ZOE.5 = "NP", ZOE.6 = "RE", ZOE.7 = "NP", ZOE.8 = "RE")
)

# Event types and appearance
event_levels <- c(
  "Roof Installation",
  "Irrigation",
  "Fertilisation"
)

event_colours <- c(
  "Roof Installation" = "#000000",
  "Irrigation" = "#444444",
  "Fertilisation" = "#777777"
)

event_linetypes <- c(
  "Roof Installation" = "11",
  "Irrigation" = "66",
  "Fertilisation" = "6414"
)

# Event dates
events <- bind_rows(
  tibble(
    Site = "KL",
    SampleDate = as.Date(c(
      "2021-05-04",
      "2021-06-30", "2021-08-25", "2021-10-20",
      "2021-05-06", "2021-07-06", "2021-08-31", "2021-10-26"
    )),
    Event = c(
      "Roof Installation",
      rep("Irrigation", 3),
      rep("Fertilisation", 4)
    )
  ),
  tibble(
    Site = "ZB",
    SampleDate = as.Date(c(
      "2022-06-20", "2022-08-16", "2022-10-10",
      "2022-06-27", "2022-08-22", "2022-10-17"
    )),
    Event = c(
      rep("Irrigation", 3),
      rep("Fertilisation", 3)
    )
  )
) |>
  mutate(Event = factor(Event, levels = event_levels))

y_limits <- list(
  P = c(0, 85),
  WFPS = c(0, 72),
  ST = c(4, 21),
  CO2 = c(0, 360),
  CH4 = c(0, 145),
  N2O = c(0, 52)
)

y_breaks <- list(
  P = seq(0, 80, 20),
  WFPS = seq(0, 60, 20),
  ST = seq(5, 20, 5),
  CO2 = seq(0, 350, 50),
  CH4 = seq(0, 140, 20),
  N2O = seq(0, 50, 10)
)

y_labels <- list(
  P = "P (mm)",
  WFPS = "WFPS (%)",
  ST = "ST (°C)",
  CO2 = expression(mg ~ CO[2] * "-C" ~ m^-2 ~ h^-1),
  CH4 = expression("μ" * g ~ CH[4] * "-C" ~ m^-2 ~ h^-1),
  N2O = expression("μ" * g ~ N[2] * "O-N" ~ m^-2 ~ h^-1)
)

# Fluxes ====

read_flux <- function(site, gas, filename) {
  raw <- read_input(filename)
  require_columns(raw, "Date", filename)

  # Accept KL1, KL1.C, ZOE1, ZB01, and names containing spaces
  clean_names <- toupper(gsub("\\s+", "", names(raw)))
  plot_columns <- which(
    str_detect(clean_names, "^(KL|ZOE|ZB)[0-9]{1,2}(\\..*)?$")
  )

  if (!length(plot_columns)) {
    stop(filename, ": no plot columns recognised.")
  }

  plot_numbers <- as.integer(
    str_replace(clean_names[plot_columns], "^(KL|ZOE|ZB)([0-9]+).*$", "\\2")
  )

  if (any(!plot_numbers %in% 1:16) || anyDuplicated(plot_numbers)) {
    stop(filename, ": invalid or duplicate plot numbers.")
  }

  if (site == "ZB" && !setequal(plot_numbers, 1:16)) {
    stop(filename, ": expected all 16 ZB plot columns.")
  }

  plot_map <- tibble(Column = names(raw)[plot_columns], Plot = plot_numbers)
  omit <- excluded_plots[[str_c(site, gas, sep = "_")]]
  multiplier <- if (gas == "CH4") ch4_multiplier else 1

  raw |>
    mutate(SampleDate = as.Date(parse_timestamp(Date), tz = data_timezone)) |>
    select(SampleDate, all_of(plot_map$Column)) |>
    pivot_longer(-SampleDate, names_to = "Column", values_to = "Value") |>
    left_join(plot_map, by = "Column") |>
    filter(!Plot %in% omit) |>
    mutate(Value = Value * multiplier) |>
    left_join(select(plot_treatments, Plot, Treatment), by = "Plot") |>
    summarise(
      Value = finite_mean(Value),
      .by = c(Plot, Treatment, SampleDate)
    ) |>
    summarise(
      Mean = finite_mean(Value),
      SE = standard_error(Value),
      N = sum(is.finite(Value)),
      .by = c(Treatment, SampleDate)
    ) |>
    mutate(Site = site, Variable = gas)
}

flux_daily <- bind_rows(
  read_flux("KL", "CO2", "co2-klausen.csv"),
  read_flux("KL", "CH4", "ch4-klausen.csv"),
  read_flux("KL", "N2O", "n2o-klausen.csv"),
  read_flux("ZB", "CO2", "co2-zoe.csv"),
  read_flux("ZB", "CH4", "ch4-zoe.csv")
) |>
  filter(!(Site == "ZB" & Variable == "CH4" & SampleDate <= zb_ch4_start))

# Environment ====

read_sensors <- function(site, variable, filename, date_column, sensor_map,
                         equal_sensor_weight = TRUE) {
  raw <- read_input(filename)
  require_columns(raw, c(date_column, names(sensor_map)), filename)

  timestamp <- parse_timestamp(raw[[date_column]])
  raw$SampleDate <- as.Date(timestamp, tz = data_timezone)

  # Midnight records are dropped for KL moisture only
  if (site == "KL" && variable == "WFPS") {
    time_of_day <- format(timestamp, "%H:%M:%S", tz = data_timezone)
    raw <- raw[time_of_day != "00:00:00", ]
  }

  dat <- raw |>
    select(SampleDate, all_of(names(sensor_map))) |>
    pivot_longer(-SampleDate, names_to = "Sensor", values_to = "Value") |>
    mutate(Treatment = unname(sensor_map[Sensor]))

  if (equal_sensor_weight) {
    dat <- summarise(
      dat,
      Value = finite_mean(Value),
      .by = c(SampleDate, Treatment, Sensor)
    )
  }

  dat <- summarise(
    dat,
    Mean = finite_mean(Value),
    .by = c(SampleDate, Treatment)
  )

  if (variable == "WFPS") {
    dat <- dat |>
      mutate(Mean = 100 * Mean / porosity[[site]][["5"]]) |>
      filter(
        if (site == "KL") {
          SampleDate < kl_moisture_end
        } else {
          SampleDate > zb_moisture_start
        }
      )
  }

  # Complete the date sequence so data gaps stay visible
  dat |>
    group_by(Treatment) |>
    complete(SampleDate = seq(min(SampleDate), max(SampleDate), by = "day")) |>
    ungroup() |>
    mutate(
      Site = site,
      Variable = variable,
      Treatment = factor(Treatment, levels = environment_levels)
    )
}

read_precipitation <- function(site, filename) {
  raw <- read_input(filename)
  require_columns(raw, c("Date", "Precipitation"), filename)

  raw |>
    mutate(SampleDate = as.Date(parse_timestamp(Date), tz = data_timezone)) |>
    summarise(Mean = finite_sum(Precipitation), .by = SampleDate) |>
    mutate(Site = site, Variable = "P")
}

environment_daily <- bind_rows(
  read_sensors(
    "KL", "WFPS", "soil-moisture-klausen.csv", "X",
    moisture_maps$KL,
    equal_sensor_weight = FALSE
  ),
  read_sensors(
    "ZB", "WFPS", "soil-moisture-zb.csv", "Date",
    moisture_maps$ZB
  ),
  read_sensors(
    "KL", "ST", "soil-temperature-klausen.csv", "X",
    temperature_maps$KL
  ),
  read_sensors(
    "ZB", "ST", "soil-temperature-zb.csv", "SampleDate",
    temperature_maps$ZB
  )
)

precipitation_daily <- bind_rows(
  read_precipitation("KL", "precipitation-klausen.csv"),
  read_precipitation("ZB", "precipitation-zb.csv")
)

# Panels ====

make_panel <- function(site, variable, bottom = FALSE) {
  is_flux <- variable %in% c("CO2", "CH4", "N2O")

  source_data <- if (is_flux) {
    flux_daily
  } else if (variable == "P") {
    precipitation_daily
  } else {
    environment_daily
  }

  dat <- filter(
    source_data,
    .data$Site == .env$site,
    .data$Variable == .env$variable
  )

  limits <- study_period[[site]]
  yrange <- y_limits[[variable]]

  month_breaks <- seq(
    as.Date(format(limits[1], "%Y-%m-01")),
    limits[2],
    by = "month"
  )

  p <- ggplot(dat, aes(x = SampleDate))

  if (variable != "P") {
    site_events <- events |>
      filter(.data$Site == .env$site)
    show_event_legend <- site == "KL" && variable == "WFPS"

    p <- p +
      geom_vline(
        data = site_events,
        aes(
          xintercept = SampleDate,
          linetype = Event
        ),
        colour = unname(event_colours[as.character(site_events$Event)]),
        linewidth = 0.65,
        show.legend = show_event_legend,
        key_glyph = "path"
      ) +
      scale_linetype_manual(
        values = event_linetypes,
        limits = event_levels,
        drop = FALSE,
        name = NULL,
        guide = if (show_event_legend) {
          guide_legend(
            order = 3,
            nrow = 1,
            override.aes = list(
              colour = unname(event_colours[event_levels]),
              linewidth = 0.65
            )
          )
        } else {
          "none"
        }
      )
  }

  if (variable == "P") {
    p <- p + geom_col(
      aes(y = Mean),
      width = 0.8,
      fill = colour_precipitation,
      na.rm = TRUE
    )
  } else if (is_flux) {
    show_flux_legend <- site == "ZB" && variable == "CH4"

    p <- p +
      geom_errorbar(
        aes(ymin = Mean - SE, ymax = Mean + SE, colour = Treatment),
        width = 0,
        linewidth = 0.30,
        na.rm = TRUE
      ) +
      geom_point(
        aes(y = Mean, fill = Treatment),
        shape = 21,
        colour = "black",
        size = 1.65,
        stroke = 0.18,
        na.rm = TRUE
      ) +
      scale_colour_manual(
        values = colours_treatment,
        limits = treatment_levels,
        guide = "none"
      ) +
      scale_fill_manual(
        values = colours_treatment,
        limits = treatment_levels,
        drop = FALSE,
        name = NULL,
        guide = if (show_flux_legend) {
          guide_legend(order = 2, nrow = 1, override.aes = list(size = 2.3))
        } else {
          "none"
        }
      )
  } else {
    show_environment_legend <- site == "ZB" && variable == "ST"

    p <- p +
      geom_line(
        aes(y = Mean, colour = Treatment, group = Treatment),
        linewidth = 0.65,
        na.rm = TRUE
      ) +
      scale_colour_manual(
        values = colours_environment,
        limits = environment_levels,
        name = NULL,
        guide = if (show_environment_legend) {
          guide_legend(order = 1, nrow = 1)
        } else {
          "none"
        }
      )
  }

  if (variable == "P") {
    p <- p + labs(
      title = if (site == "KL") "Klausen-Leopoldsdorf" else "Zöbelboden"
    )
  }

  p <- p +
    scale_x_date(
      breaks = month_breaks,
      date_labels = "%b",
      minor_breaks = NULL,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = y_breaks[[variable]],
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(xlim = limits, ylim = yrange, expand = FALSE) +
    labs(
      x = if (bottom) {
        str_c("Date (", if (site == "KL") "2021" else "2022", ")")
      } else {
        NULL
      },
      y = if (site == "KL") y_labels[[variable]] else NULL
    ) +
    theme_exafor() +
    theme(
      axis.title = element_text(size = 11),
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      plot.title.position = "panel",
      plot.margin = margin(
        t = 5,
        r = if (site == "KL") 0 else 6,
        b = 5,
        l = if (site == "ZB") 0 else 6
      ),
      legend.title = element_blank()
    )

  if (!bottom) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }

  if (site == "ZB") {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }

  p
}

row_label <- function(letter) {
  wrap_elements(
    panel = grid::textGrob(
      label = letter,
      x = grid::unit(0.8, "npc"),
      gp = grid::gpar(
        fontsize = 10,
        fontface = "bold",
        fontfamily = "Source Sans 3"
      )
    )
  ) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

panels <- list(
  make_panel("KL", "P"),
  row_label("A"),
  make_panel("ZB", "P"),
  make_panel("KL", "WFPS"),
  row_label("B"),
  make_panel("ZB", "WFPS"),
  make_panel("KL", "ST"),
  row_label("C"),
  make_panel("ZB", "ST"),
  make_panel("KL", "CO2"),
  row_label("D"),
  make_panel("ZB", "CO2"),
  make_panel("KL", "CH4"),
  row_label("E"),
  make_panel("ZB", "CH4", bottom = TRUE),
  make_panel("KL", "N2O", bottom = TRUE),
  row_label("F"),
  guide_area()
)

combined_figure <- wrap_plots(
  panels,
  ncol = 3,
  byrow = TRUE,
  widths = grid::unit(c(1, 4, 1), c("null", "mm", "null")),
  heights = c(1, 1, 1, 1.7, 1.7, 1.7),
  guides = "collect"
) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.y = unit(0.5, "cm")
  )

# Export ====

save_figure(combined_figure, figure_name, figure_width_mm, figure_height_mm)
log_script_end("figure-03-gas-fluxes")
