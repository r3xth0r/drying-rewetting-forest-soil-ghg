# Bayesian N2O figures and diagnostic plots ----
# Inputs: artifacts written by n2o-bayesian-models.R in dat/interim/n2o/.
# Outputs: Figures 5 and S8; trace, predictive, Pareto-k, and GAM check plots.

suppressPackageStartupMessages({
  library("tidyverse")
  library("brms")
  library("bayesplot")
  library("mgcv")
  library("patchwork")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("n2o-bayesian-outputs")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))

# Settings ====

diagnostic_seed <- 2026L
depth_gam_diagnostic_seed <- 1234L
gam_check_replicates <- 1000L

partial_plot_width_mm <- 200
partial_plot_height_mm <- 90
combined_plot_width_mm <- 150
combined_plot_height_mm <- 200
depth_gam_plot_width_mm <- 160
depth_gam_plot_height_mm <- 160
depth_gam_diagnostic_width_mm <- 200
depth_gam_diagnostic_height_mm <- 200
diagnostic_trace_width_mm <- 200
diagnostic_trace_height_mm <- 200
diagnostic_predictive_width_mm <- 200
diagnostic_predictive_height_mm <- 200
diagnostic_pareto_width_mm <- 200
diagnostic_pareto_height_mm <- 2000
n2o_plot_y_padding <- 5

model_plot_depth_labels <- c(
  "0-10" = "0–10 cm",
  "10-20" = "10–20 cm",
  "20-40" = "20–40 cm"
)
model_depth_labels <- c("5" = "0-10", "15" = "10-20", "30" = "20-40")
n2o_environment_colours <- colours_environment

n2o_interim_dir <- file.path(dir_interim, "n2o")
n2o_cache_dir <- file.path(n2o_interim_dir, "cache")
diagnostics_dir <- file.path(n2o_interim_dir, "diagnostics")

model_data_file <- file.path(n2o_interim_dir, "n2o-model-data.rds")
full_model_file <- file.path(n2o_cache_dir, "n2o-full-model.rds")
reduced_model_file <- file.path(n2o_cache_dir, "n2o-reduced-model.rds")
file_full_curves_rds <- file.path(n2o_interim_dir, "n2o-full-curves.rds")
file_reduced_curves_rds <- file.path(n2o_interim_dir, "n2o-reduced-curves.rds")
file_full_loo_rds <- file.path(n2o_interim_dir, "n2o-full-loo.rds")
file_reduced_loo_rds <- file.path(n2o_interim_dir, "n2o-reduced-loo.rds")
file_depth_gam_models <- file.path(n2o_interim_dir, "n2o-depth-wfps-gams.rds")
file_depth_gam_data <- file.path(n2o_interim_dir, "n2o-depth-wfps-data.rds")

file_final_figure <- "figure-05-n2o-partial-curves"
file_partial_wfps_np <- "figure-05-n2o-partial-wfps-np"
file_partial_wfps_re <- "figure-05-n2o-partial-wfps-re"
file_partial_temp_np <- "figure-05-n2o-partial-temperature-np"
file_partial_temp_re <- "figure-05-n2o-partial-temperature-re"
file_depth_wfps_figure <- "figure-s08-n2o-depth-wfps-gams"
file_depth_wfps_gam_checks <- "n2o-depth-wfps-gam-checks.pdf"
file_diagnostic_trace_suffix <- "trace-plots"
file_diagnostic_predictive_suffix <- "predictive-full-range"
file_diagnostic_predictive_zoom_suffix <- "predictive-zoomed"
file_diagnostic_pareto_suffix <- "pareto-k.pdf"

full_diagnostic_model_name <- "n2o-full-gam"
reduced_diagnostic_model_name <- "n2o-reduced-gam"
diagnostic_trace_parameter_count <- 6L
diagnostic_predictive_draws <- 100L
diagnostic_zoom_padding <- c(-1, 1)

ensure_dir(dir_plots)
ensure_dir(diagnostics_dir)

# Inputs ====

df_pool_lag <- readRDS(model_data_file)
full_curves <- readRDS(file_full_curves_rds)
reduced_curves <- readRDS(file_reduced_curves_rds)
depth_gam_models <- readRDS(file_depth_gam_models)
depth_gam_data <- readRDS(file_depth_gam_data)

# Helpers ====

depth_range_label <- function(depth_vec) {
  unname(model_depth_labels[as.character(depth_vec)])
}

# Figure 5: partial effects ====

plot_partial_by_treatment <- function(curves_df, trt, vary) {
  vary <- match.arg(vary, c("WFPS", "Temp_C"))
  x_label <- if (vary == "WFPS") "WFPS (%)" else "Soil temperature (°C)"
  variable_label <- if (vary == "WFPS") "WFPS" else "temperature"

  ggplot(
    curves_df |> filter(Treatment2 == trt),
    aes(x = x, y = mean)
  ) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~Depth, nrow = 1, scales = "free_x") +
    labs(
      x = x_label,
      y = expression(E[N[2] * O] ~ ("μ" * g ~ m^{
        -2
      } ~ h^{
        -1
      })),
      title = paste0(
        "Partial effect of ", variable_label,
        " by depth (Treatment2 = ", trt, ")"
      )
    ) +
    theme_exafor(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

## Full-model curves ----
p_wfps_np <- plot_partial_by_treatment(full_curves$wfps, "NP", "WFPS")
p_wfps_re <- plot_partial_by_treatment(full_curves$wfps, "RE", "WFPS")
p_temp_np <- plot_partial_by_treatment(full_curves$temp, "NP", "Temp_C")
p_temp_re <- plot_partial_by_treatment(full_curves$temp, "RE", "Temp_C")

save_figure(
  p_wfps_np, file_partial_wfps_np, partial_plot_width_mm,
  partial_plot_height_mm
)
save_figure(
  p_wfps_re, file_partial_wfps_re, partial_plot_width_mm,
  partial_plot_height_mm
)
save_figure(
  p_temp_np, file_partial_temp_np, partial_plot_width_mm,
  partial_plot_height_mm
)
save_figure(
  p_temp_re, file_partial_temp_re, partial_plot_width_mm,
  partial_plot_height_mm
)

## Combined figure (reduced model) ----
df_plot <- df_pool_lag |>
  mutate(
    N2O_obs = pmax(exp(y) - 1e-6, 0),
    Depth_lab = factor(
      depth_range_label(Depth),
      levels = names(model_plot_depth_labels),
      labels = unname(model_plot_depth_labels)
    )
  )

ymax <- max(df_plot$N2O_obs, na.rm = TRUE) + n2o_plot_y_padding

plot_final_curve <- function(curves, x_column, x_label) {
  ggplot() +
    geom_ribbon(
      data = curves,
      aes(x = x, ymin = lo, ymax = hi, fill = Treatment2),
      alpha = 0.15
    ) +
    geom_point(
      data = df_plot,
      aes(x = .data[[x_column]], y = N2O_obs, color = Treatment2),
      alpha = 0.5,
      size = 1.8
    ) +
    geom_line(
      data = curves,
      aes(x = x, y = mean, linetype = "Model"),
      color = "black",
      linewidth = 0.9
    ) +
    scale_color_manual(values = n2o_environment_colours) +
    scale_fill_manual(
      values = n2o_environment_colours,
      guide = "none"
    ) +
    facet_grid(Depth_lab ~ Treatment2) +
    theme_exafor(base_size = 12) +
    theme(legend.position = "none") +
    labs(
      x = x_label,
      y = expression(N[2] * O ~ "μ" * g ~ m^{
        -2
      } ~ h^{
        -1
      }),
      title = NULL
    ) +
    coord_cartesian(ylim = c(0, ymax))
}

p_wfps_curves <- plot_final_curve(reduced_curves$wfps, "WFPS", "WFPS (%)")
p_temp_curves <- plot_final_curve(
  reduced_curves$temp, "Temp_C", "Soil temperature (°C)"
)

# Manual legend plot
legend_points <- data.frame(
  x = c(1, 2, 3),
  y = 1,
  lab = c("NP", "RE", "Measured"),
  col = c(unname(n2o_environment_colours), "black")
)

legend_labels <- data.frame(
  x = c(1.12, 2.12, 3.12, 4.47),
  y = 1,
  lab = c("NP", "RE", "Measured", "Model")
)

p_legend <- ggplot() +
  geom_point(
    data = legend_points,
    aes(x = x, y = y),
    color = legend_points$col,
    size = 2
  ) +
  geom_segment(
    data = data.frame(x = 4, xend = 4.35, y = 1, yend = 1),
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_text(
    data = legend_labels,
    aes(x = x, y = y, label = lab),
    hjust = 0,
    vjust = 0.5,
    size = 4,
    family = "Source Sans 3"
  ) +
  coord_cartesian(xlim = c(0.8, 5.1), ylim = c(0.8, 1.2), expand = FALSE) +
  theme_void(base_family = "Source Sans 3")

p_top <- (p_wfps_curves / p_temp_curves) +
  plot_annotation(tag_levels = "a")

p_top_fixed <- wrap_elements(full = p_top)

p_final <- p_top_fixed / p_legend +
  plot_layout(heights = c(1, 0.12))

save_figure(
  p_final,
  file_final_figure,
  combined_plot_width_mm,
  combined_plot_height_mm
)

# Figure S8: depth-specific WFPS GAMs ====

depth_layer_factor <- function(depth) {
  factor(
    depth,
    levels = c("5", "15", "30"),
    labels = c("0–10 cm", "10–20 cm", "20–40 cm")
  )
}

depth_wfps_plot_data <- depth_gam_data$data |>
  mutate(Depth_layer = depth_layer_factor(Depth))

depth_wfps_plot_predictions <- depth_gam_data$predictions |>
  mutate(Depth_layer = depth_layer_factor(Depth))

depth_wfps_plot <- ggplot() +
  geom_ribbon(
    data = depth_wfps_plot_predictions,
    aes(WFPS, ymin = lower, ymax = upper),
    fill = "#0072B2",
    alpha = 0.2
  ) +
  geom_point(
    data = depth_wfps_plot_data,
    aes(WFPS, N2O),
    alpha = 0.4,
    size = 1.4
  ) +
  geom_line(
    data = depth_wfps_plot_predictions,
    aes(WFPS, estimate),
    colour = "#0072B2",
    linewidth = 0.9
  ) +
  facet_grid(Depth_layer ~ Treatment) +
  labs(
    x = "WFPS (%)",
    y = expression(N[2] * O ~ flux ~ ("μ" * g ~ N ~ m^{
      -2
    } ~ h^{
      -1
    }))
  ) +
  theme_exafor(base_size = 12, base_family = "Source Sans 3")

save_figure(
  depth_wfps_plot,
  file_depth_wfps_figure,
  depth_gam_plot_width_mm,
  depth_gam_plot_height_mm
)

# Diagnostic plots ====

export_model_diagnostic_plots <- function(
  model_fit, loo_result, model_name
) {
  draws <- as_draws_array(model_fit)

  trace_parameters <- posterior::summarise_draws(draws) |>
    filter(is.finite(rhat)) |>
    arrange(desc(rhat)) |>
    slice_head(n = diagnostic_trace_parameter_count) |>
    pull(variable)

  trace_plot <- mcmc_trace(draws, pars = trace_parameters) +
    ggtitle(paste(model_name, "trace plots"))

  save_figure(
    trace_plot,
    paste0(model_name, "-", file_diagnostic_trace_suffix),
    diagnostic_trace_width_mm,
    diagnostic_trace_height_mm,
    dir = diagnostics_dir
  )

  set.seed(diagnostic_seed)

  predictive_plot <- pp_check(
    model_fit,
    type = "dens_overlay", ndraws = diagnostic_predictive_draws
  ) +
    labs(
      title = model_name,
      x = "Log N2O response"
    )

  save_figure(
    predictive_plot,
    paste0(model_name, "-", file_diagnostic_predictive_suffix),
    diagnostic_predictive_width_mm,
    diagnostic_predictive_height_mm,
    dir = diagnostics_dir
  )

  zoom_limits <- range(df_pool_lag$y, na.rm = TRUE) +
    diagnostic_zoom_padding

  save_figure(
    predictive_plot + coord_cartesian(xlim = zoom_limits),
    paste0(model_name, "-", file_diagnostic_predictive_zoom_suffix),
    diagnostic_predictive_width_mm,
    diagnostic_predictive_height_mm,
    dir = diagnostics_dir
  )

  if (save_pdf_figures) {
    pdf(
      file.path(
        diagnostics_dir,
        paste0(model_name, "-", file_diagnostic_pareto_suffix)
      ),
      width = mm_to_inches(diagnostic_pareto_width_mm),
      height = mm_to_inches(diagnostic_pareto_height_mm)
    )
    tryCatch(
      plot(loo_result, diagnostic = "k"),
      finally = dev.off()
    )
  }
}

export_model_diagnostic_plots(
  readRDS(full_model_file),
  readRDS(file_full_loo_rds),
  full_diagnostic_model_name
)

export_model_diagnostic_plots(
  readRDS(reduced_model_file),
  readRDS(file_reduced_loo_rds),
  reduced_diagnostic_model_name
)

## Depth GAM checks ----

if (save_pdf_figures) {
  pdf(
    file.path(diagnostics_dir, file_depth_wfps_gam_checks),
    width = mm_to_inches(depth_gam_diagnostic_width_mm),
    height = mm_to_inches(depth_gam_diagnostic_height_mm)
  )
  tryCatch(
    {
      par(
        mfrow = c(length(depth_gam_models), 4),
        mar = c(2.5, 2.5, 2, 0.5),
        mgp = c(1.5, 0.5, 0),
        cex = 0.65
      )

      set.seed(depth_gam_diagnostic_seed)
      for (id in names(depth_gam_models)) {
        cat("\n---", id, "---\n")
        gam.check(
          depth_gam_models[[id]],
          k.rep = gam_check_replicates
        )
      }
    },
    finally = dev.off()
  )
}

log_script_end("n2o-bayesian-outputs")
