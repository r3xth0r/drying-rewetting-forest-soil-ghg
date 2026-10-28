# Figure 6: Soil parameters ----
# Creates depth treatment comparisons and Figure 6.
# Input: soil-depth-sampling.csv in dat/raw.
# Outputs: Figure 6 and statistics in dat/reporting;
# diagnostics in dat/interim.

suppressPackageStartupMessages({
  library("tidyverse")
  library("patchwork")
  library("car")
  library("emmeans")
  library("multcompView")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("figure-06-soil-parameters")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))
source(here::here("dev", "soil-parameter-helpers.R"))

# Settings ====

depth_input_file <- "soil-depth-sampling.csv"

dot_threshold <- 0.10
dot_symbol <- "·"
depth_levels <- c("1", "2")

depth_key <- tribble(
  ~Site, ~Sampling, ~ID,
  "KL",  "02_21",   "KL2",
  "KL",  "05_21",   "KL5",
  "KL",  "06_21",   "KL6",
  "ZOE", "12_22",   "ZB12",
  "ZOE", "13_22",   "ZB13"
)

depth_axis_maxima <- list(
  KL = c(MBC = 1300, MBN = 120),
  ZOE = c(MBC = 4500, MBN = 600)
)

depth_axis_labels <- list(
  MBC = expression("mg C"["mic"] ~ "kg"^-1 ~ "dw"),
  MBN = expression("mg N"["mic"] ~ "kg"^-1 ~ "dw")
)

depth_plot <- list(
  mean_letter_offset = 0.04,
  group_label_offset = 0.13,
  dot_offset = 0.06,
  treatment_letter_offset = 0.06,
  depth_symbol_offset = 0.05,
  dodge_width = 0.6,
  bar_width = 0.6,
  errorbar_width = 0.2,
  treatment_letter_size = 3,
  depth_symbol_size = 6,
  y_title_size = 12,
  axis_linewidth = 0.3,
  panel_margin = 8,
  x_expand = c(0.15, 0),
  y_expand = c(0, 0)
)

figure_06_name <- "figure-06-soil-parameters"
figure_06_width_mm <- 200
figure_06_height_mm <- 220
depth_report_file <- "figure-06-statistical-results.txt"
depth_summary_file <- "figure-06-means-se-letters.csv"
depth_comparisons_file <- "figure-06-emmeans-comparisons.csv"
depth_symbols_file <- "figure-06-depth-comparisons.csv"
depth_diagnostics_file <- "figure-06-model-diagnostics.pdf"

# Figure 6 ====

## Import and prepare data ----

depth_raw <- read_input(depth_input_file)

require_columns(
  depth_raw,
  c("Site", "Sampling", "treatment", "Depth", "MBC", "MBN"),
  depth_input_file
)

depth_data <- depth_raw |>
  mutate(across(c(Site, Sampling, treatment, Depth), clean_text)) |>
  inner_join(depth_key, by = c("Site", "Sampling"))

missing_samplings <- setdiff(depth_key$ID, unique(depth_data$ID))

if (length(missing_samplings)) {
  stop("No data found for: ", paste(missing_samplings, collapse = ", "))
}

check_treatments(depth_data$treatment)

if (anyNA(depth_data$Depth) || any(!depth_data$Depth %in% depth_levels)) {
  stop("Depth must be coded as 1 or 2.")
}

if (!is.numeric(depth_data$MBC) || !is.numeric(depth_data$MBN)) {
  stop("MBC and MBN must be numeric. Check the CSV formatting.")
}

depth_data <- depth_data |>
  mutate(
    Treatment = as_treatment(treatment),
    Depth = factor(Depth, levels = depth_levels)
  )


## emmeans comparisons and annotations ----

make_annotations <- function(model) {
  emm <- emmeans(model, ~ Treatment * Depth)
  emm_table <- as.data.frame(emm)

  treatments <- levels(model.frame(model)$Treatment)
  depths <- levels(model.frame(model)$Depth)

  if (length(depths) != 2L) {
    stop("Exactly two depths are required.")
  }

  # Preserve Tukey adjustment across all eight treatment × depth groups.
  pairwise <- contrast(
    emm,
    method = "pairwise", adjust = "tukey"
  )

  comparison_summary <- summary(
    pairwise,
    infer = c(TRUE, TRUE),
    level = 1 - alpha,
    adjust = "tukey"
  )

  comparisons <- as.data.frame(comparison_summary)

  # Build reliable labels from coefficients rather than parsing names.
  group_labels <- paste(emm_table$Treatment, emm_table$Depth, sep = ":")
  coefficient_table <- as.data.frame(coef(pairwise))

  comparisons$contrast <- vapply(
    seq_len(nrow(comparisons)),
    function(k) {
      weights <- coefficient_table[[paste0("c.", k)]]
      positive <- which(weights > 0)
      negative <- which(weights < 0)

      if (length(positive) != 1L || length(negative) != 1L) {
        stop("Unexpected coefficients in a pairwise contrast.")
      }

      paste(group_labels[positive], group_labels[negative], sep = " - ")
    },
    character(1)
  )

  get_p <- function(treatment1, depth1, treatment2, depth2) {
    group1 <- paste(treatment1, depth1, sep = ":")
    group2 <- paste(treatment2, depth2, sep = ":")

    forward <- paste(group1, group2, sep = " - ")
    reverse <- paste(group2, group1, sep = " - ")

    index <- match(forward, comparisons$contrast)

    if (is.na(index)) {
      index <- match(reverse, comparisons$contrast)
    }

    if (is.na(index)) {
      stop("Missing comparison: ", forward)
    }

    p <- comparisons$p.value[index]

    if (!is.finite(p)) {
      stop("Non-estimable comparison: ", forward)
    }

    p
  }

  treatment_letters <- map_dfr(depths, function(depth) {
    p_matrix <- matrix(
      1,
      nrow = length(treatments),
      ncol = length(treatments),
      dimnames = list(treatments, treatments)
    )

    for (i in seq_len(length(treatments) - 1L)) {
      for (j in seq.int(i + 1L, length(treatments))) {
        p_matrix[i, j] <- p_matrix[j, i] <- get_p(
          treatments[i], depth, treatments[j], depth
        )
      }
    }

    depth_means <- emm_table |>
      filter(as.character(Depth) == as.character(depth))

    letters <- make_letters(
      p_matrix,
      setNames(depth_means$emmean, as.character(depth_means$Treatment)),
      alpha
    )

    tibble(
      Treatment = treatments,
      Depth = depth,
      significance_treat = unname(letters[treatments])
    )
  })

  depth_stars <- map_dfr(treatments, function(treatment) {
    p <- get_p(treatment, depths[1], treatment, depths[2])

    tibble(
      Treatment = treatment,
      p_adj_depth = p,
      significance_depth = case_when(
        p < alpha ~ "*",
        p < dot_threshold ~ dot_symbol,
        TRUE ~ ""
      )
    )
  })

  list(
    estimated_means = as.data.frame(
      summary(emm, infer = c(TRUE, FALSE), level = 1 - alpha)
    ),
    comparison_summary = comparison_summary,
    comparisons = comparisons,
    treatment_letters = treatment_letters,
    depth_stars = depth_stars
  )
}


## ANOVAs and assumption checks ----

analyse_panel <- function(id, response) {
  panel_data <- depth_data |> filter(ID == id)
  dat <- panel_data |> filter(is.finite(.data[[response]]))
  counts <- dat |> count(Treatment, Depth, .drop = FALSE)

  if (any(counts$n < 2L)) {
    stop(
      id, " / ", response,
      ": each Treatment × Depth group needs at least two observations."
    )
  }

  # Sequential / Type I sums of squares.
  model_formula <- reformulate("Treatment * Depth", response = response)
  model <- aov(model_formula, data = dat)
  residual_values <- residuals(model)

  shapiro <- if (
    length(residual_values) >= 3L &&
      length(residual_values) <= 5000L &&
      sd(residual_values) > 0
  ) {
    shapiro.test(residual_values)
  } else {
    NULL
  }

  levene <- leveneTest(
    model_formula,
    data = dat, center = median
  )

  annotations <- make_annotations(model)

  summary_data <- dat |>
    group_by(Treatment, Depth) |>
    summarise(
      n = n(),
      mean_value = mean(.data[[response]]),
      se_value = sd(.data[[response]]) / sqrt(n()),
      .groups = "drop"
    ) |>
    mutate(across(c(Treatment, Depth), as.character)) |>
    left_join(
      annotations$treatment_letters,
      by = c("Treatment", "Depth")
    )

  list(
    id = id,
    response = response,
    n_excluded = nrow(panel_data) - nrow(dat),
    counts = counts,
    model = model,
    shapiro = shapiro,
    levene = levene,
    estimated_means = annotations$estimated_means,
    comparison_summary = annotations$comparison_summary,
    comparisons = annotations$comparisons,
    summary = summary_data,
    depth_stars = annotations$depth_stars
  )
}

depth_results <- list()

for (id in depth_key$ID) {
  for (response in c("MBC", "MBN")) {
    depth_results[[paste(response, id, sep = "_")]] <-
      analyse_panel(id, response)
  }
}


## Export statistical results and diagnostics ----

depth_report <- capture.output({
  cat(
    "Figure 6 statistical results\n",
    "Required input file: ", depth_input_file, " in dat/raw.\n",
    "Model: response ~ Treatment * Depth\n",
    "Separate models for each site/sampling and response.\n",
    "ANOVA sums of squares: sequential / Type I.\n",
    "With unequal cell sizes, term order can affect results.\n",
    "Post hoc: emmeans, native pairwise Tukey comparisons.\n",
    "Adjustment family: eight Treatment × Depth groups per model.\n",
    "No multiplicity adjustment across the ten models.\n",
    "Levene test: median-centred, across eight groups.\n",
    "Bars: observed means ± SE.\n",
    "Letters: shown only where treatments differ within depth.\n",
    "Letters assigned in descending mean order within depth.\n",
    "Depth symbols: * p < 0.05; dot 0.05 <= p < 0.10.\n",
    sep = ""
  )

  for (key in names(depth_results)) {
    result <- depth_results[[key]]
    print_model_heading(key, result$n_excluded)

    print_section("Sample sizes used:", as.data.frame(result$counts))
    print_section("Two-way ANOVA:", summary(result$model))

    if (is.null(result$shapiro)) {
      cat("\nShapiro–Wilk: residuals unsuitable for this test.\n")
    } else {
      print_section("Shapiro–Wilk test of residuals:", result$shapiro)
    }

    print_section("Levene's test:", result$levene)
    print_section("Estimated marginal means:", result$estimated_means)

    # Retain the emmeans footer documenting adjustments.
    print_section("emmeans pairwise comparisons:", result$comparison_summary)

    print_section(
      "Treatment letters within depth:",
      as.data.frame(select(
        result$summary, Treatment, Depth, significance_treat
      ))
    )

    print_section(
      "Depth comparisons within treatment:",
      as.data.frame(result$depth_stars)
    )
  }
})

ensure_dir(dir_reporting)
writeLines(
  depth_report,
  file.path(dir_reporting, depth_report_file)
)

depth_summary <- collect_results(depth_results, "summary")
depth_comparisons <- collect_results(depth_results, "comparisons")
depth_symbols <- collect_results(depth_results, "depth_stars")

write_table(depth_summary, depth_summary_file, dir_reporting)
write_table(depth_comparisons, depth_comparisons_file, dir_reporting)
write_table(depth_symbols, depth_symbols_file, dir_reporting)

save_diagnostic_pdf(
  depth_results,
  depth_diagnostics_file,
  hide_annotations = TRUE
)


## Create Figure 6 panels ----

fill_levels <- unlist(lapply(treatment_levels, function(treatment) {
  paste(treatment, depth_levels, sep = ":")
}))

fill_colours <- colours_soil_parameter_treatment_depth

fill_labels <- unlist(lapply(treatment_levels, function(treatment) {
  paste(treatment, c("0–10 cm", "10–20 cm"))
}))

depth_sampling_numbers <- c(
  KL2 = "0",
  KL5 = "2A",
  KL6 = "2B",
  ZB12 = "2A",
  ZB13 = "2B"
)

make_depth_panel <- function(
  id, response,
  show_y_title = FALSE,
  show_header = FALSE
) {
  result <- depth_results[[paste(response, id, sep = "_")]]

  site <- if (startsWith(id, "KL")) "KL" else "ZOE"
  y_max <- depth_axis_maxima[[site]][[response]]

  plot_data <- result$summary |>
    mutate(
      Treatment = factor(Treatment, levels = treatment_levels),
      Depth = factor(Depth, levels = depth_levels),
      Fill_Group = factor(
        paste(Treatment, Depth, sep = ":"),
        levels = fill_levels
      ),
      letter_y = mean_value + se_value + depth_plot$mean_letter_offset * y_max
    )

  depth_labels <- plot_data |>
    group_by(Treatment) |>
    summarise(
      y = max(mean_value + se_value) + depth_plot$group_label_offset * y_max,
      .groups = "drop"
    ) |>
    mutate(Treatment = as.character(Treatment)) |>
    left_join(result$depth_stars, by = "Treatment") |>
    mutate(
      Treatment = factor(Treatment, levels = treatment_levels),
      y = if_else(
        significance_depth == dot_symbol,
        y - depth_plot$dot_offset * y_max,
        y
      )
    )

  upper_limit <- max(
    y_max,
    plot_data$mean_value + plot_data$se_value,
    plot_data$letter_y[plot_data$significance_treat != ""] +
      depth_plot$treatment_letter_offset * y_max,
    depth_labels$y[depth_labels$significance_depth != ""] +
      depth_plot$depth_symbol_offset * y_max
  )

  dodge <- position_dodge(width = depth_plot$dodge_width)

  ggplot(
    plot_data,
    aes(Treatment, mean_value, fill = Fill_Group, group = Depth)
  ) +
    geom_col(
      position = dodge,
      width = depth_plot$bar_width,
      colour = "black",
      linewidth = depth_plot$axis_linewidth
    ) +
    geom_errorbar(
      aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      position = dodge, width = depth_plot$errorbar_width
    ) +
    geom_text(
      aes(y = letter_y, label = significance_treat),
      position = dodge, vjust = 0, size = depth_plot$treatment_letter_size,
      family = "Source Sans 3"
    ) +
    geom_text(
      data = depth_labels,
      aes(Treatment, y, label = significance_depth),
      inherit.aes = FALSE, size = depth_plot$depth_symbol_size,
      family = "Source Sans 3"
    ) +
    labs(
      title = if (show_header) {
        paste("Sampling", depth_sampling_numbers[[id]])
      } else {
        NULL
      }
    ) +
    scale_fill_manual(
      name = "Treatment and depth",
      values = fill_colours,
      breaks = fill_levels,
      labels = fill_labels,
      drop = FALSE
    ) +
    scale_x_discrete(expand = expansion(mult = depth_plot$x_expand)) +
    scale_y_continuous(expand = expansion(mult = depth_plot$y_expand)) +
    coord_cartesian(ylim = c(0, upper_limit), clip = "on") +
    labs(
      x = NULL,
      y = if (show_y_title) depth_axis_labels[[response]] else NULL
    ) +
    theme_classic(base_family = "Source Sans 3") +
    theme(
      axis.title.y = if (show_y_title) {
        element_text(size = depth_plot$y_title_size, angle = 90)
      } else {
        element_blank()
      },
      plot.margin = margin(depth_plot$panel_margin, 8, 8, 8),
      plot.title = element_text(
        size = 10,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      plot.title.position = "panel"
    )
}


## Assemble and save Figure 6 ----

depth_panels <- list(
  A = make_depth_panel(
    "KL2", "MBC",
    show_y_title = TRUE,
    show_header = TRUE
  ),
  B = make_depth_panel("KL5", "MBC", show_header = TRUE),
  C = make_depth_panel("KL6", "MBC", show_header = TRUE),
  D = make_depth_panel("KL2", "MBN", show_y_title = TRUE),
  E = make_depth_panel("KL5", "MBN"),
  F = make_depth_panel("KL6", "MBN"),
  G = make_depth_panel("ZB12", "MBC", show_y_title = TRUE),
  H = make_depth_panel("ZB13", "MBC"),
  I = make_depth_panel("ZB12", "MBN", show_y_title = TRUE),
  J = make_depth_panel("ZB13", "MBN"),
  K = guide_area()
)

layout_design <- "
ABC
DEF
KGH
KIJ
"

final_plot2 <- wrap_plots(
  depth_panels,
  design = layout_design,
  widths = c(1, 1, 1),
  heights = c(1, 1, 1, 1),
  guides = "collect"
) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = grid::unit(0.45, "cm"),
    legend.key.width = grid::unit(0.45, "cm")
  )

save_figure(
  final_plot2,
  figure_06_name,
  width_mm = figure_06_width_mm,
  height_mm = figure_06_height_mm
)
log_script_end("figure-06-soil-parameters")
