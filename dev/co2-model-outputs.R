# CO2 model figures and tables ----
# Reads outputs from co2-models.R and creates manuscript figures and tables.

suppressPackageStartupMessages({
  library("tidyverse")
  library("patchwork")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("co2-model-outputs")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))
source(here::here("dev", "co2-model-helpers.R"))

kl_excluded_plots <- c("KL2", "KL3", "KL10")
kl_end_exclusive <- as.Date("2021-11-22")
zb_start_exclusive <- as.Date("2022-04-26")
site_depths <- list(KL = c(5, 15, 30), ZB = 5)
model_labels <- c("Temperature Function", "Temperature Moisture Function")

file_table_03 <- "table-03-co2-models.csv"
file_table_03_settings <- "table-03-co2-model-settings.csv"
file_table_s05 <- "table-s05-co2-temperature-residuals.csv"
file_models <- "co2-models.rds"
file_results <- "co2-model-results.csv"
file_predictions <- "co2-model-predictions.csv"

residual_width_mm <- 200
residual_height_mm <- c(KL = 220, ZB = 80)
mesh_resolution <- 36L
mesh_theta <- -60
mesh_phi <- 25
mesh_width_mm <- 180
mesh_row_height_mm <- 70
mesh_legend_height_mm <- 20
mesh_colour_power <- 1.3
box_grid_colour <- "#C8C8C8"
box_grid_width <- 0.45
mesh_palette <- colorspace::sequential_hcl(
  n = 256, palette = "Batlow", rev = TRUE
)
treatment_names <- c(NP = "Natural Precipitation", RE = "Rainfall Exclusion")
flux_label <- expression(
  paste("CO"[2], "-C flux (mg C ", m^{
    -2
  }, " ", h^{
    -1
  }, ")")
)
manuscript_keys <- c("KL_5_NP", "KL_5_RE", "KL_15_NP", "KL_15_RE")
supplement_keys <- c(
  manuscript_keys, "KL_30_NP", "KL_30_RE", "ZB_5_NP", "ZB_5_RE"
)

dir_models <- file.path(dir_interim, file_models)
dir_results <- file.path(dir_processed, file_results)
dir_predictions <- file.path(dir_processed, file_predictions)

models <- read_rds(dir_models)
results <- readr::read_csv(dir_results, show_col_types = FALSE)
predictions <- readr::read_csv(dir_predictions, show_col_types = FALSE)

groups <- purrr::imap(site_depths, \(depths, site) {
  tidyr::expand_grid(site = site, depth_cm = depths, treatment = c("NP", "RE"))
}) |>
  purrr::list_rbind()

settings <- tribble(
  ~Setting, ~Value,
  "Input directory", "dat/raw",
  "Excluded KL plots", str_c(kl_excluded_plots, collapse = ", "),
  "KL date restriction", str_c("Before ", kl_end_exclusive),
  "ZB moisture date restriction", str_c("After ", zb_start_exclusive),
  "ZB end date restriction", "None",
  "Depth mapping",
  "5 cm -> 0–10 cm; 15 cm -> 10–20 cm; 30 cm -> 20–40 cm",
  "KL bulk density: 5 / 15 / 30 cm", str_c(bulk_density$KL, collapse = " / "),
  "ZB bulk density", as.character(bulk_density$ZB),
  "Particle density", as.character(particle_density),
  "WFPS calculation", "100 * moisture percent / porosity percent",
  "Observations used",
  "Same complete observations for both models within each group",
  "Temperature–moisture a",
  "Fixed to corresponding temperature-model a at full precision",
  "R² definition", "1 - SSE/SST",
  "KGE definition", "Original 2009 KGE",
  "AIC/BIC interpretation",
  "Moisture-model AIC/BIC are conditional on a being fixed",
  "Equation precision",
  "Six significant figures; parameter sheet retains full precision",
  "Fit statistics", "In-sample goodness of fit; no cross-validation"
)

table_03 <- results |>
  select(
    Model, Site, Depth_cm, Depth, Treatment, Equation, N, R2, RMSE,
    AIC, BIC, KGE
  )

write_table(table_03, file_table_03, dir_reporting)
write_table(settings, file_table_03_settings, dir_reporting)

# Residuals against WFPS ====

## Quadratic fits ----

residual_data <- predictions |>
  filter(is.finite(WFPS), is.finite(Temperature_residual))

fit_residuals <- function(site, depth_cm, treatment) {
  key <- group_key(site, depth_cm, treatment)
  dat <- select_group(residual_data, site, depth_cm, treatment)

  if (nrow(dat) <= 2L || n_distinct(dat$WFPS) < 2L) {
    stop(key, ": insufficient observations or WFPS variation.")
  }

  fit <- lm(Temperature_residual ~ 0 + WFPS + I(WFPS^2), data = dat)

  if (fit$rank != 2L || any(!is.finite(coef(fit)))) {
    stop(key, ": quadratic coefficients cannot be estimated.")
  }

  dat$Residual_prediction <- as.numeric(fitted(fit))
  observed <- dat$Temperature_residual
  errors <- observed - dat$Residual_prediction
  sse <- sum(errors^2)
  sst <- sum((observed - mean(observed))^2)

  list(
    data = dat,
    statistics = tibble(
      Site = site,
      Depth_cm = depth_cm,
      Treatment = treatment,
      N = nrow(dat),
      d = unname(coef(fit)[1]),
      e = unname(coef(fit)[2]),
      R2 = if (sst > 0) 1 - sse / sst else NA_real_,
      RMSE = sqrt(mean(errors^2)),
      AIC = AIC(fit),
      BIC = BIC(fit)
    )
  )
}

residual_fits <- pmap(groups, fit_residuals) |>
  set_names(pmap_chr(groups, group_key))

residual_plot_data <- map(residual_fits, "data")
residual_statistics <- residual_fits |>
  map("statistics") |>
  list_rbind()

## Panels ----

make_residual_panel <- function(site, depth, treatment) {
  key <- group_key(site, depth, treatment)
  dat <- residual_plot_data[[key]]
  metrics <- select_group(residual_statistics, site, depth, treatment)

  if (is.null(dat) || nrow(metrics) != 1L) {
    stop("Missing residual-model results for ", key)
  }

  statistic_label <- sprintf(
    "R² = %.2f\nRMSE = %.2f\nAIC = %.0f\nBIC = %.0f",
    metrics$R2, metrics$RMSE, metrics$AIC, metrics$BIC
  )

  ggplot(dat, aes(x = WFPS)) +
    geom_hline(
      yintercept = 0, colour = colours_environment[["RE"]], linewidth = 0.4
    ) +
    geom_point(
      aes(y = Temperature_residual),
      colour = "black", shape = 18, size = 1.4, alpha = 0.7
    ) +
    geom_point(
      aes(y = Residual_prediction),
      colour = colours_environment[["NP"]], shape = 18, size = 1.4
    ) +
    annotation_custom(
      grid::textGrob(
        statistic_label,
        x = grid::unit(0.04, "npc"), y = grid::unit(0.97, "npc"),
        hjust = 0, vjust = 1,
        gp = grid::gpar(
          fontsize = 3.2 * ggplot2::.pt, lineheight = 1.05,
          fontfamily = "Source Sans 3"
        )
      )
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.15))) +
    labs(x = "WFPS (%)", y = "Residuals of temperature model") +
    theme_classic(base_size = 11, base_family = "Source Sans 3") +
    theme(
      plot.tag = element_text(size = 14),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      plot.margin = margin(10, 12, 10, 10)
    )
}

# Panels run NP, RE within each depth
residual_panels <- function(site) {
  expand_grid(depth = site_depths[[site]], treatment = c("NP", "RE")) |>
    pmap(\(depth, treatment) make_residual_panel(site, depth, treatment)) |>
    wrap_plots(ncol = 2, byrow = TRUE) +
    plot_annotation(tag_levels = "A")
}

## Export ----

save_residual_figure <- function(site, name) {
  save_figure(
    residual_panels(site), name,
    width_mm = residual_width_mm,
    height_mm = residual_height_mm[[site]]
  )
}

save_residual_figure("KL", "figure-s03-kl-temperature-residuals")
save_residual_figure("ZB", "figure-s04-zb-temperature-residuals")

write_table(residual_statistics, file_table_s05, dir_reporting)

# Mesh figures ====

## Surfaces ----

make_mesh_object <- function(site, depth_cm, treatment) {
  key <- group_key(site, depth_cm, treatment)
  pair <- models[[key]]

  dat <- select_group(predictions, site, depth_cm, treatment) |>
    filter(is.finite(Tsoil), is.finite(WFPS), is.finite(Flux))

  if (
    nrow(dat) < 2L || diff(range(dat$Tsoil)) <= 0 ||
      diff(range(dat$WFPS)) <= 0
  ) {
    stop("Insufficient plotting data: ", key)
  }

  temperature <- seq(
    min(dat$Tsoil), max(dat$Tsoil),
    length.out = mesh_resolution
  )
  moisture <- seq(
    min(dat$WFPS), max(dat$WFPS),
    length.out = mesh_resolution
  )

  grid <- expand.grid(
    Tsoil = temperature, WFPS = moisture, KEEP.OUT.ATTRS = FALSE
  )
  grid$a_fixed <- unname(coef(pair$Temperature)["a"])

  fitted_flux <- as.numeric(predict(pair$Temperature_moisture, newdata = grid))

  if (any(!is.finite(fitted_flux))) {
    stop("Non-finite surface predictions: ", key)
  }

  metric <- results |>
    select_group(site, depth_cm, treatment) |>
    filter(Model == model_labels[2])

  if (nrow(metric) != 1L) {
    stop("Expected one temperature–moisture result: ", key)
  }

  list(
    site = site,
    depth = depth_cm,
    treatment = treatment,
    observed = dat,
    x = temperature,
    y = moisture,
    z = matrix(
      fitted_flux,
      nrow = length(temperature), ncol = length(moisture)
    ),
    r2 = metric$R2
  )
}

mesh_objects <- pmap(groups, make_mesh_object) |>
  set_names(pmap_chr(groups, group_key))

## Drawing helpers ----

pad_range <- function(x, fraction = 0.04) {
  limits <- range(x, finite = TRUE)
  span <- diff(limits)

  if (span == 0) {
    span <- max(abs(limits), 1)
  }

  limits + c(-1, 1) * span * fraction
}

flux_text <- function(x) {
  format(signif(x, 3), trim = TRUE, scientific = FALSE)
}

# persp() puts moisture on x and temperature on y
project_mesh <- function(temperature, moisture, flux, projection) {
  trans3d(x = moisture, y = temperature, z = flux, pmat = projection)
}

surface_colours <- function(z, limits) {
  nr <- nrow(z)
  nc <- ncol(z)

  facet_height <- (
    z[-nr, -nc, drop = FALSE] + z[-1, -nc, drop = FALSE] +
      z[-nr, -1, drop = FALSE] + z[-1, -1, drop = FALSE]
  ) / 4

  fraction <- mesh_colour_fraction(facet_height, limits)
  index <- 1L + floor(fraction * (length(mesh_palette) - 1L))

  colours <- matrix(NA_character_, nrow = nr - 1L, ncol = nc - 1L)
  finite <- is.finite(index)
  colours[finite] <- mesh_palette[index[finite]]

  as.vector(colours)
}

mesh_colour_fraction <- function(values, limits) {
  fraction <- (values - limits[1]) / diff(limits)
  pmax(0, pmin(1, fraction))^mesh_colour_power
}

draw_mesh_axes <- function(projection, xlim, ylim, zlim, z_ticks, z_labels) {
  project <- \(x, y, z) {
    project_mesh(
      temperature = x, moisture = y, flux = z, projection = projection
    )
  }

  axis_ticks <- function(limits) {
    ticks <- pretty(limits, n = 4)
    ticks[ticks >= limits[1] & ticks <= limits[2]]
  }

  axis_angle <- function(points) {
    screen_x <- grconvertX(points$x, from = "user", to = "inches")
    screen_y <- grconvertY(points$y, from = "user", to = "inches")
    angle <- atan2(diff(screen_y), diff(screen_x)) * 180 / pi

    if (angle > 90) angle <- angle - 180
    if (angle < -90) angle <- angle + 180

    angle
  }

  xt <- axis_ticks(xlim)
  yt <- axis_ticks(ylim)

  # Temperature axis
  temperature_edge <- project(xlim, rep(ylim[1], 2), rep(zlim[1], 2))

  text(
    project(xt, ylim[1], zlim[1]),
    labels = flux_text(xt), pos = 1, offset = 0.45, cex = 0.88
  )
  text(
    project(mean(xlim), ylim[1], zlim[1]),
    labels = "Soil Temperature (°C)", pos = 1, offset = 2.4,
    srt = axis_angle(temperature_edge), cex = 1.15
  )

  # Moisture axis; the title sits outward of the edge in screen space
  moisture_edge <- project(rep(xlim[1], 2), ylim, rep(zlim[1], 2))

  text(
    project(xlim[1], yt, zlim[1]),
    labels = flux_text(yt), pos = 4, offset = 0.45, cex = 0.88
  )

  edge_x <- grconvertX(moisture_edge$x, from = "user", to = "inches")
  edge_y <- grconvertY(moisture_edge$y, from = "user", to = "inches")
  dx <- diff(edge_x)
  dy <- diff(edge_y)
  edge_length <- sqrt(dx^2 + dy^2)

  title_gap <- 0.32
  title_x <- mean(edge_x) + title_gap * dy / edge_length
  title_y <- mean(edge_y) + title_gap * -dx / edge_length

  text(
    x = grconvertX(title_x, from = "inches", to = "user"),
    y = grconvertY(title_y, from = "inches", to = "user"),
    labels = "WFPS (%)", adj = c(0.5, 0.5),
    srt = axis_angle(moisture_edge), cex = 1.15
  )

  # Flux axis: left vertical edge
  text(
    project(xlim[2], ylim[1], z_ticks),
    labels = z_labels, pos = 2, offset = 0.45, cex = 0.88
  )

  flux_title_position <- project(xlim[2], ylim[1], mean(zlim))
  horizontal_span <- diff(par("usr")[1:2])

  text(
    x = flux_title_position$x - 0.13 * horizontal_span,
    y = flux_title_position$y,
    labels = flux_label, srt = 90, cex = 1.10
  )
}

draw_box_grid <- function(projection, temperature_limits, moisture_limits,
                          flux_limits, flux_ticks) {
  grid_ticks <- function(limits) {
    ticks <- pretty(limits, n = 6)
    ticks[ticks > limits[1] & ticks < limits[2]]
  }

  temperature_ticks <- grid_ticks(temperature_limits)
  moisture_ticks <- grid_ticks(moisture_limits)
  flux_ticks <- flux_ticks[
    flux_ticks > flux_limits[1] & flux_ticks < flux_limits[2]
  ]

  grid_line <- function(temperature, moisture, flux) {
    lines(
      project_mesh(temperature, moisture, flux, projection),
      col = box_grid_colour, lwd = box_grid_width
    )
  }

  t_min <- temperature_limits[1]
  t_max <- temperature_limits[2]
  w_min <- moisture_limits[1]
  w_max <- moisture_limits[2]
  z_min <- flux_limits[1]
  z_max <- flux_limits[2]

  for (temperature in temperature_ticks) {
    grid_line(rep(temperature, 2), c(w_min, w_max), c(z_min, z_min))
    grid_line(rep(temperature, 2), c(w_max, w_max), c(z_min, z_max))
  }

  for (moisture in moisture_ticks) {
    grid_line(c(t_min, t_max), rep(moisture, 2), c(z_min, z_min))
    grid_line(c(t_max, t_max), rep(moisture, 2), c(z_min, z_max))
  }

  for (flux in flux_ticks) {
    grid_line(c(t_min, t_max), c(w_max, w_max), rep(flux, 2))
    grid_line(c(t_max, t_max), c(w_min, w_max), rep(flux, 2))
  }
}

paint_surface <- function(obj, z, limits, zlim, col, border, box) {
  persp(
    x = obj$y, y = obj$x, z = t(z),
    xlim = limits$y, ylim = limits$x, zlim = zlim,
    theta = mesh_theta, phi = mesh_phi, r = 100000,
    expand = 0.90, scale = TRUE, axes = FALSE, box = box,
    col = col, border = border, lwd = 0.18, shade = NA,
    xlab = "", ylab = "", zlab = ""
  )
}

draw_mesh_legend <- function(zlim, z_ticks, z_labels, logarithmic) {
  par(mar = c(0.2, 0.8, 0.2, 0.8), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  bar_left <- 0.16
  bar_right <- 0.70
  bar_bottom <- 0.42
  bar_top <- 0.62
  edges <- seq(bar_left, bar_right, length.out = length(mesh_palette) + 1L)

  rect(
    xleft = head(edges, -1), ybottom = bar_bottom,
    xright = tail(edges, -1), ytop = bar_top,
    col = mesh_palette, border = NA
  )
  rect(bar_left, bar_bottom, bar_right, bar_top, border = "black")

  legend_positions <- bar_left +
    mesh_colour_fraction(z_ticks, zlim) * (bar_right - bar_left)

  segments(legend_positions, bar_bottom, legend_positions, 0.34)
  text(legend_positions, 0.29, labels = z_labels, adj = c(0.5, 1), cex = 0.90)
  text(
    (bar_left + bar_right) / 2, 0.82,
    labels = if (logarithmic) "Flux\n(log scale)" else "Flux",
    adj = c(0.5, 0.5), cex = 1.00
  )
  points(0.78, 0.54, pch = 16, cex = 1.00)
  text(
    0.82, 0.54,
    labels = expression("Measured " * CO[2] * "-C flux"),
    adj = c(0, 0.5), cex = 0.90
  )
}

## Mesh figure ----

draw_mesh_figure <- function(keys, logarithmic = FALSE) {
  if (length(keys) %% 2L != 0L || !all(keys %in% names(mesh_objects))) {
    stop("Panel keys must identify complete NP/RE pairs.")
  }

  objects <- mesh_objects[keys]
  number_of_rows <- length(keys) / 2L
  all_flux <- unlist(map(objects, \(obj) c(obj$z, obj$observed$Flux)))

  # Shared flux limits and colour scale within each figure
  if (logarithmic) {
    positive_flux <- all_flux[is.finite(all_flux) & all_flux > 0]

    if (!length(positive_flux)) {
      stop("No positive flux values for logarithmic plotting.")
    }

    log_range <- range(log10(positive_flux))
    zlim <- c(floor(log_range[1]), ceiling(log_range[2]))

    if (diff(zlim) == 0) {
      zlim <- zlim + c(-1, 1)
    }

    z_grid_ticks <- seq(zlim[1], zlim[2], by = 1)
    z_ticks <- z_grid_ticks

    if (length(z_ticks) > 6L) {
      selected_ticks <- unique(round(
        seq(1, length(z_ticks), length.out = 6)
      ))
      z_ticks <- z_ticks[selected_ticks]
    }

    z_labels <- flux_text(10^z_ticks)
  } else {
    flux_range <- range(c(0, all_flux), finite = TRUE)
    z_ticks <- pretty(flux_range, n = 5)
    zlim <- range(z_ticks)

    if (diff(zlim) == 0) {
      zlim <- zlim + c(-1, 1)
      z_ticks <- pretty(zlim, n = 5)
      zlim <- range(z_ticks)
    }

    z_grid_ticks <- z_ticks
    z_labels <- flux_text(z_ticks)
  }

  # Common predictor ranges per site
  predictor_limits <- unique(map_chr(objects, "site")) |>
    set_names() |>
    map(\(site) {
      selected <- keep(objects, \(obj) obj$site == site)
      list(
        x = pad_range(unlist(map(selected, "x"))),
        y = pad_range(unlist(map(selected, "y")))
      )
    })

  old_par <- par(no.readonly = TRUE)
  on.exit(
    {
      layout(1)
      par(old_par)
    },
    add = TRUE
  )

  panel_layout <- matrix(seq_along(keys), ncol = 2, byrow = TRUE)

  layout(
    rbind(panel_layout, rep(length(keys) + 1L, 2)),
    widths = c(1, 1),
    heights = c(
      rep(1, number_of_rows),
      mesh_legend_height_mm / mesh_row_height_mm
    )
  )
  par(family = "Source Sans 3", oma = c(0.2, 0.4, 0.2, 0.4))

  for (i in seq_along(objects)) {
    obj <- objects[[i]]
    limits <- predictor_limits[[obj$site]]
    displayed_z <- obj$z
    dat <- obj$observed

    if (logarithmic) {
      omitted_points <- sum(dat$Flux <= 0)
      omitted_nodes <- sum(displayed_z <= 0)

      if (omitted_points > 0 || omitted_nodes > 0) {
        logger::log_warn(paste0(
          keys[i], ": log display excludes ", omitted_points,
          " non-positive observations and ", omitted_nodes,
          " non-positive surface nodes."
        ))
      }

      displayed_z[displayed_z <= 0] <- NA_real_
      displayed_z <- log10(displayed_z)
      dat <- dat[dat$Flux > 0, , drop = FALSE]
      displayed_observed <- log10(dat$Flux)
    } else {
      displayed_observed <- dat$Flux
    }

    if (all(!is.finite(displayed_z))) {
      stop("No drawable surface for ", keys[i])
    }

    par(mar = c(3.0, 3.4, 3.4, 2.6), xpd = NA)

    # Grid lines are drawn behind the surface
    projection <- paint_surface(
      obj, displayed_z, limits, zlim,
      col = NA, border = NA, box = FALSE
    )

    draw_box_grid(
      projection,
      temperature_limits = limits$x,
      moisture_limits = limits$y,
      flux_limits = zlim,
      flux_ticks = z_grid_ticks
    )

    par(new = TRUE)

    projection <- paint_surface(
      obj, displayed_z, limits, zlim,
      col = surface_colours(t(displayed_z), zlim),
      border = "#666666", box = TRUE
    )

    draw_mesh_axes(projection, limits$x, limits$y, zlim, z_ticks, z_labels)

    points(
      project_mesh(dat$Tsoil, dat$WFPS, displayed_observed, projection),
      pch = 16, cex = 0.55, col = adjustcolor("#202020", alpha.f = 0.45)
    )

    mtext(
      unname(treatment_names[obj$treatment]),
      side = 3, line = 1.7, font = 2, cex = 0.95
    )
    mtext(
      str_c(obj$site, " ", depth_labels[as.character(obj$depth)]),
      side = 3, line = 0.45, cex = 0.85
    )
    mtext(LETTERS[i], side = 3, line = 0.45, adj = 0, font = 2, cex = 1.15)

    r2_text <- if (is.finite(obj$r2)) {
      formatC(obj$r2, digits = 2, format = "f")
    } else {
      "NA"
    }

    mtext(
      bquote(R^2 == .(r2_text)),
      side = 3, line = 0.45, adj = 1, cex = 0.90
    )
  }

  draw_mesh_legend(zlim, z_ticks, z_labels, logarithmic)
}

## Export ----

save_mesh_figure <- function(keys, logarithmic, name) {
  height_mm <- length(keys) / 2L * mesh_row_height_mm + mesh_legend_height_mm
  ensure_dir(dir_plots)

  if (save_pdf_figures) {
    pdf_file <- file.path(dir_plots, str_c(name, ".pdf"))
    if (capabilities("cairo")) {
      cairo_pdf(
        pdf_file,
        width = mm_to_inches(mesh_width_mm),
        height = mm_to_inches(height_mm),
        family = "Source Sans 3"
      )
    } else {
      pdf(
        pdf_file,
        width = mm_to_inches(mesh_width_mm),
        height = mm_to_inches(height_mm),
        family = "Helvetica"
      )
    }
    tryCatch(draw_mesh_figure(keys, logarithmic), finally = dev.off())
  }

  png(
    file.path(dir_plots, str_c(name, ".png")),
    width = mm_to_inches(mesh_width_mm),
    height = mm_to_inches(height_mm), units = "in",
    res = fig_dpi, bg = "white", family = "Source Sans 3"
  )
  tryCatch(draw_mesh_figure(keys, logarithmic), finally = dev.off())

  logger::log_info(paste("Figure saved:", file.path(dir_plots, name)))
}

save_mesh_figure(
  manuscript_keys,
  logarithmic = TRUE, name = "figure-04-co2-mesh"
)
save_mesh_figure(
  supplement_keys,
  logarithmic = FALSE, name = "figure-s06-co2-mesh-linear"
)
log_script_end("co2-model-outputs")
