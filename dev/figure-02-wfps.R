# Figure 2: Klausen WFPS ----
# Daily WFPS at three depths with precipitation (NP) and rewetting (RE) bars.
# Inputs: dat/raw soil-moisture, precipitation, irrigation (Klausen).
# Output: plt/figure-02-wfps.{png,pdf}

suppressPackageStartupMessages({
  library("tidyverse")
  library("patchwork")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("figure-02-wfps")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))

# Settings ====

date_limits <- as.Date(c("2021-06-01", "2021-11-16"))
month_breaks <- seq(as.Date("2021-06-01"), as.Date("2021-11-01"), by = "month")

# Moisture records from this date on are not used
moisture_end <- as.Date("2021-11-22")

figure_name <- "figure-02-wfps"
figure_width_mm <- 180
figure_height_mm <- 180

# Moisture ====

sensor_map <- expand_grid(
  Depth_cm = c(5, 15, 30),
  Position = c("C1", "C2", "D1", "D2")
) |>
  mutate(
    Sensor = str_c(Position, "_", Depth_cm, "cm"),
    Treatment = if_else(Position %in% c("C1", "C2"), "NP", "RE"),
    Depth = factor(
      Depth_cm,
      levels = as.numeric(names(depth_labels)),
      labels = depth_labels
    )
  )

filename <- "soil-moisture-klausen.csv"
raw_moisture <- read_input(filename)
require_columns(raw_moisture, c("X", sensor_map$Sensor), filename)

moisture_daily <- raw_moisture |>
  mutate(
    timestamp = parse_timestamp(X),
    SampleDate = as.Date(timestamp, tz = data_timezone)
  ) |>
  filter(format(timestamp, "%H:%M:%S", tz = data_timezone) != "00:00:00") |>
  select(SampleDate, all_of(sensor_map$Sensor)) |>
  pivot_longer(-SampleDate, names_to = "Sensor", values_to = "Moisture_pc") |>
  left_join(
    select(sensor_map, Sensor, Treatment, Depth_cm, Depth),
    by = "Sensor"
  ) |>
  filter(SampleDate < moisture_end) |>
  summarise(
    Moisture_pc = finite_mean(Moisture_pc),
    .by = c(Treatment, Depth_cm, Depth, SampleDate)
  ) |>
  mutate(
    WFPS = 100 * Moisture_pc / unname(porosity$KL[as.character(Depth_cm)])
  ) |>
  group_by(Treatment, Depth_cm, Depth) |>
  complete(SampleDate = seq(date_limits[1], date_limits[2], by = "day")) |>
  ungroup()

# Water inputs ====

precipitation <- read_input("precipitation-klausen.csv") |>
  mutate(SampleDate = parse_date(Date)) |>
  summarise(Water_mm = finite_sum(Precipitation), .by = SampleDate)

irrigation <- read_input("irrigation-klausen.csv") |>
  mutate(SampleDate = parse_date(SampleDate)) |>
  summarise(Water_mm = finite_sum(Irrigation), .by = SampleDate)

# Figure ====

make_moisture_panel <- function(
  treatment, water_data, water_name, water_colour, letter,
  moisture_limits, water_max, bottom = FALSE
) {
  # Conversion between the two displayed y axes
  water_scale <- diff(moisture_limits) / water_max

  water_data <- water_data |>
    filter(SampleDate >= date_limits[1], SampleDate <= date_limits[2]) |>
    mutate(
      Bar_top = moisture_limits[1] + Water_mm * water_scale,
      Water = water_name
    )

  moisture_data <- filter(moisture_daily, Treatment == treatment)

  p <- ggplot() +
    geom_rect(
      data = water_data,
      aes(
        xmin = SampleDate - 0.4,
        xmax = SampleDate + 0.4,
        ymax = Bar_top,
        fill = Water
      ),
      ymin = moisture_limits[1],
      alpha = 0.6,
      colour = NA,
      na.rm = TRUE
    ) +
    geom_line(
      data = moisture_data,
      aes(x = SampleDate, y = WFPS, colour = Depth, group = Depth),
      linewidth = 0.65,
      na.rm = TRUE
    ) +
    annotate(
      "text",
      x = date_limits[1] + 3,
      y = moisture_limits[2] - 0.03 * diff(moisture_limits),
      label = letter,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      family = "Source Sans 3",
      size = 5
    ) +
    scale_colour_manual(
      values = colours_depth[[treatment]],
      limits = depth_labels,
      labels = str_c("WFPS% ", depth_labels),
      name = NULL,
      guide = guide_legend(order = 1)
    ) +
    scale_fill_manual(
      values = set_names(water_colour, water_name),
      name = NULL,
      guide = guide_legend(order = 2)
    ) +
    scale_x_date(
      breaks = month_breaks,
      date_labels = "%b",
      minor_breaks = NULL,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      name = "WFPS (%)",
      breaks = pretty(moisture_limits, n = 6),
      expand = expansion(mult = 0),
      sec.axis = sec_axis(
        \(y) (y - moisture_limits[1]) / water_scale,
        name = str_c(water_name, " (mm)"),
        breaks = seq(0, water_max, by = 20)
      )
    ) +
    coord_cartesian(
      xlim = date_limits,
      ylim = moisture_limits,
      expand = FALSE
    ) +
    labs(x = if (bottom) "2021" else NULL) +
    theme_exafor(base_size = 12) +
    theme(
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = 0.5
      ),
      axis.title = element_text(size = 13),
      axis.title.x = element_text(margin = margin(t = 12)),
      axis.title.y.right = element_text(margin = margin(l = 10)),
      legend.key.width = unit(1.0, "cm"),
      legend.key.height = unit(0.45, "cm"),
      legend.text = element_text(size = 10),
      legend.box = "vertical",
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(8, 8, 8, 8)
    )

  if (treatment == "NP") {
    p <- p + theme(
      legend.position = c(0.98, 0.94),
      legend.justification = c(1, 1)
    )
  } else {
    p <- p + theme(
      legend.position = c(0.23, 0.87),
      legend.justification = c(0, 1)
    )
  }

  if (!bottom) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }

  p
}

panel_a <- make_moisture_panel(
  treatment = "NP",
  water_data = precipitation,
  water_name = "Precipitation",
  water_colour = colour_precipitation,
  letter = "A",
  moisture_limits = c(30, 80),
  water_max = 100
)

panel_b <- make_moisture_panel(
  treatment = "RE",
  water_data = irrigation,
  water_name = "Rewetting",
  water_colour = colour_water_bars,
  letter = "B",
  moisture_limits = c(20, 60),
  water_max = 160,
  bottom = TRUE
)

moisture_figure <- (panel_a / panel_b) + plot_layout(heights = c(1, 1))

# Export ====

save_figure(moisture_figure, figure_name, figure_width_mm, figure_height_mm)
log_script_end("figure-02-wfps")
