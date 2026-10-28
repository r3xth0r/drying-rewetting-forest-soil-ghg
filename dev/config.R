# Project configuration ----
# Paths, constants, colours and the ggplot theme shared by all scripts.
# Source after attaching packages; run scripts from the project root.

# Paths ====

dir_raw <- here::here("dat", "raw")
dir_interim <- here::here("dat", "interim")
dir_processed <- here::here("dat", "processed")
dir_reporting <- here::here("dat", "reporting")
dir_plots <- here::here("plt")

# General ====

data_timezone <- "Europe/Vienna"

options(
  mc.cores = max(1L, parallel::detectCores() - 1L),
  brms.backend = "cmdstanr"
)

# Study design ====

treatment_levels <- c("C", "C + N", "DRW", "DRW + N")
environment_levels <- c("NP", "RE")

# N2O GAM smooth basis dimensions
full_smooth_basis <- c(6L, 5L)
reduced_smooth_basis <- c(6L, 5L)

# Plot numbers (1-16) per treatment; identical at both sites
plot_treatments <- tibble::tibble(
  Plot = 1:16,
  Treatment = factor(
    dplyr::case_when(
      Plot %in% c(1, 5, 10, 13) ~ "C",
      Plot %in% c(2, 6, 9, 14) ~ "C + N",
      Plot %in% c(3, 7, 11, 15) ~ "DRW",
      Plot %in% c(4, 8, 12, 16) ~ "DRW + N"
    ),
    levels = treatment_levels
  ),
  Environment = factor(
    dplyr::if_else(Plot %in% c(1, 2, 5, 6, 9, 10, 13, 14), "NP", "RE"),
    levels = environment_levels
  )
)

# Measurement periods (Klausen = KL, Zöbelboden = ZB)
study_period <- list(
  KL = as.Date(c("2021-03-24", "2021-11-25")),
  ZB = as.Date(c("2022-04-25", "2022-11-25"))
)

# Sensor maps ====

sensor_map_kl_moisture <- tibble::tibble(
  Sensor = paste0(
    rep(c("C1", "C2", "D1", "D2"), 3), "_",
    rep(c(5, 15, 30), each = 4), "cm"
  ),
  Depth_cm = rep(c(5, 15, 30), each = 4),
  Treatment = rep(c("NP", "NP", "RE", "RE"), 3)
)

sensor_map_kl_temperature <- tibble::tibble(
  Sensor = c(
    "C1", "C2", "D1", "D2",
    "C1_15", "C2_15", "D1_15", "D2_15",
    "C1_30", "C2_30", "D1_30", "D2_30"
  ),
  Depth_cm = rep(c(5, 15, 30), each = 4),
  Treatment = rep(c("NP", "NP", "RE", "RE"), 3)
)

sensor_map_zb_temperature <- tibble::tibble(
  Sensor = c("ZOE.5", "ZOE.6", "ZOE.7", "ZOE.8"),
  Depth_cm = 5,
  Treatment = c("NP", "RE", "NP", "RE")
)

# Soil physics ====

particle_density <- 2.65 # g cm-3

# Bulk density (g cm-3) by sensor depth (cm)
bulk_density <- list(
  KL = c(`5` = 0.606650896, `15` = 0.80177125631552, `30` = 0.99101316636197),
  ZB = c(`5` = 0.261)
)

depth_labels <- c(
  `5` = "0–10 cm",
  `15` = "10–20 cm",
  `30` = "20–40 cm"
)

# Porosity (%) by site and sensor depth (cm)
porosity <- purrr::map(bulk_density, \(x) 100 * (1 - x / particle_density))

# Colours ====

colours_treatment <- c(
  "C" = "#0072B2",
  "C + N" = "#009E73",
  "DRW" = "#CC79A7",
  "DRW + N" = "#D55E00"
)

colours_soil_parameter_treatment <- colours_treatment

colours_soil_parameter_treatment_depth <- c(
  "C:1" = "#56B4E9",
  "C:2" = "#0072B2",
  "C + N:1" = "#56C6A4",
  "C + N:2" = "#009E73",
  "DRW:1" = "#E3A2C3",
  "DRW:2" = "#CC79A7",
  "DRW + N:1" = "#E98768",
  "DRW + N:2" = "#D55E00"
)

colours_environment <- c(NP = "#0072B2", RE = "#D55E00")

colours_depth <- list(
  NP = c(
    "0–10 cm" = "#56B4E9",
    "10–20 cm" = "#2C91C7",
    "20–40 cm" = "#0072B2"
  ),
  RE = c(
    "0–10 cm" = "#F3B49E",
    "10–20 cm" = "#E78360",
    "20–40 cm" = "#D55E00"
  )
)

colour_precipitation <- "#407999"
colour_water_bars <- "#B39040"

# Figure export ====

fig_dpi <- 400
save_pdf_figures <- FALSE

# ggplot theme ====

theme_exafor <- function(base_size = 11, base_family = "Source Sans 3") {
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black", fill = NA, linewidth = 0.55
      ),
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 6, 5, 6)
    )
}

ggplot2::theme_set(theme_exafor())
