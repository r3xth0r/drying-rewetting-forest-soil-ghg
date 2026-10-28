# Supplementary Figure S9: Soil parameters ----
# Runs Welch ANOVA and Games-Howell comparisons for the soil survey.
# Inputs: soil-parameter-survey.csv in dat/raw.
# Outputs: Figure S9, statistical summaries, and diagnostic plots.

suppressPackageStartupMessages({
  library("tidyverse")
  library("patchwork")
  library("multcompView")
  library("scales")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("figure-s09-soil-parameters")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))
source(here::here("dev", "soil-parameter-helpers.R"))

soil_input_file <- "soil-parameter-survey.csv"
parameter_levels <- c("MBC", "MBN", "NH4N", "NO3N", "DOC", "DON", "GWC")
event_labels <- c(
  "Sampling 0",
  "Sampling 1A", "Sampling 1B",
  "Sampling 2A", "Sampling 2B",
  "Sampling 3A", "Sampling 3B"
)

soil_key <- tibble(
  Site = rep(c("KL", "ZOE"), each = 7),
  Sampling = c(sprintf("%02d_21", 2:8), sprintf("%02d_22", 9:15)),
  Event = rep(event_labels, times = 2)
)

axis_maxima <- list(
  KL = c(
    MBC = 1100, MBN = 110, NH4N = 30, NO3N = 10,
    DOC = 400, DON = 50, GWC = 40
  ),
  ZOE = c(
    MBC = 4500, MBN = 700, NH4N = 80, NO3N = 50,
    DOC = 1250, DON = 160, GWC = 80
  )
)

axis_labels <- list(
  MBC = expression("mg C"["mic"] ~ "kg"^-1 ~ "dw"),
  MBN = expression("mg N"["mic"] ~ "kg"^-1 ~ "dw"),
  NH4N = expression("μ" * "g NH"[4] * "-N" ~ "g"^-1 ~ "dw"),
  NO3N = expression("μ" * "g NO"[3] * "-N" ~ "g"^-1 ~ "dw"),
  DOC = expression("mg DOC" ~ "kg"^-1 ~ "dw"),
  DON = expression("mg DON" ~ "kg"^-1 ~ "dw"),
  GWC = "GWC %"
)

figure_s09_name <- "figure-s09-soil-parameters"
figure_s09_width_mm <- 200
figure_s09_height_mm <- 220
soil_report_file <- "figure-s09-welch-games-howell-statistics.txt"
soil_welch_file <- "figure-s09-welch-anova.csv"
soil_comparisons_file <- "figure-s09-games-howell-comparisons.csv"
soil_summary_file <- "figure-s09-means-se-letters.csv"
soil_diagnostics_file <- "figure-s09-model-diagnostics.pdf"
# Figure S9 ====

## Import and prepare data ----

soil_raw <- read.csv(
  file.path(dir_raw, soil_input_file),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

names(soil_raw) <- make.names(trimws(names(soil_raw)), unique = TRUE)

require_columns(
  soil_raw,
  c(
    "site", "plot", "sampling", "treatment", "pc_TS",
    "NO3N", "NH4N", "DOC", "MBC", "N_notfum", "MBN"
  ),
  soil_input_file
)

numeric_columns <- c(
  "pc_TS", "NO3N", "NH4N", "DOC", "MBC", "N_notfum", "MBN"
)

for (column in numeric_columns) {
  soil_raw[[column]] <- strict_numeric(soil_raw[[column]], column)
}

soil_data <- soil_raw |>
  transmute(
    Site = site,
    Plot = plot,
    Sampling = sampling,
    treatment,
    pc_TS, MBC, MBN, NH4N, NO3N,
    DOC_raw = DOC,
    N_notfum
  ) |>
  mutate(across(c(Site, Plot, Sampling, treatment), clean_text)) |>
  inner_join(soil_key, by = c("Site", "Sampling"))

missing_events <- soil_key |>
  anti_join(distinct(soil_data, Site, Sampling), by = c("Site", "Sampling"))

if (nrow(missing_events)) {
  print(missing_events)
  stop("Some required site/sampling combinations are missing.")
}

check_treatments(soil_data$treatment)

if (
  any(!is.finite(soil_data$pc_TS)) ||
    any(soil_data$pc_TS <= 0 | soil_data$pc_TS > 100)
) {
  stop("pc_TS must contain finite percentages greater than 0 and <= 100.")
}

# Retain the conversions from the supplied script.
soil_data <- soil_data |>
  mutate(
    Treatment = as_treatment(treatment),
    Event = factor(Event, levels = event_labels),
    DOC = (DOC_raw / pc_TS) * 1000,
    DON = (N_notfum / pc_TS) * 1000,
    GWC = 100 - pc_TS
  )

duplicate_records <- soil_data |>
  count(Site, Sampling, Plot) |>
  filter(n > 1L)

if (nrow(duplicate_records)) {
  print(duplicate_records)
  stop(
    "Multiple records per plot and sampling. ",
    "Check whether these are duplicates or subsamples before continuing."
  )
}

long_data <- soil_data |>
  select(Site, Sampling, Event, Plot, Treatment, all_of(parameter_levels)) |>
  pivot_longer(
    all_of(parameter_levels),
    names_to = "Parameter",
    values_to = "Value"
  )


## Welch ANOVA and Games–Howell ----

games_howell <- function(summary_data, alpha = 0.05) {
  pairs <- combn(seq_len(nrow(summary_data)), 2)

  map_dfr(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1, k]
    j <- pairs[2, k]

    n_i <- summary_data$n[i]
    n_j <- summary_data$n[j]
    v_i <- summary_data$SD[i]^2 / n_i
    v_j <- summary_data$SD[j]^2 / n_j

    estimate <- summary_data$Mean[i] - summary_data$Mean[j]
    se_difference <- sqrt(v_i + v_j)

    df <- (v_i + v_j)^2 /
      (v_i^2 / (n_i - 1) + v_j^2 / (n_j - 1))

    q_statistic <- abs(estimate) / sqrt((v_i + v_j) / 2)

    p_adjusted <- ptukey(
      q_statistic,
      nmeans = nrow(summary_data),
      df = df,
      lower.tail = FALSE
    )

    margin <- qtukey(
      1 - alpha,
      nmeans = nrow(summary_data),
      df = df
    ) * sqrt((v_i + v_j) / 2)

    tibble(
      group1 = as.character(summary_data$Treatment[i]),
      group2 = as.character(summary_data$Treatment[j]),
      n1 = n_i,
      n2 = n_j,
      estimate = estimate,
      SE_difference = se_difference,
      df = df,
      q = q_statistic,
      conf.low = estimate - margin,
      conf.high = estimate + margin,
      p.adj = p_adjusted,
      method = "Games-Howell"
    )
  })
}

analyse_group <- function(site, sampling, parameter) {
  model_name <- paste(site, sampling, parameter, sep = "_")

  original_data <- long_data |>
    filter(Site == site, Sampling == sampling, Parameter == parameter)

  dat <- original_data |>
    filter(is.finite(Value)) |>
    mutate(Treatment = factor(Treatment, levels = treatment_levels))

  counts <- dat |> count(Treatment, .drop = FALSE)

  if (any(counts$n < 2L)) {
    print(counts)
    stop(
      model_name,
      ": each treatment requires at least two finite observations."
    )
  }

  summary_data <- dat |>
    group_by(Treatment, .drop = FALSE) |>
    summarise(
      n = n(),
      Mean = mean(Value),
      SD = sd(Value),
      SE = SD / sqrt(n),
      .groups = "drop"
    )

  if (any(!is.finite(summary_data$SD) | summary_data$SD <= 0)) {
    stop(
      model_name,
      ": a treatment has zero or undefined variance; inspect its values."
    )
  }

  welch <- oneway.test(
    Value ~ Treatment,
    data = dat, var.equal = FALSE
  )

  comparisons <- games_howell(summary_data, alpha)

  if (any(!is.finite(comparisons$p.adj))) {
    stop(model_name, ": a Games–Howell p-value is undefined.")
  }

  treatments <- as.character(summary_data$Treatment)
  p_matrix <- matrix(
    1,
    nrow = length(treatments),
    ncol = length(treatments),
    dimnames = list(treatments, treatments)
  )

  for (k in seq_len(nrow(comparisons))) {
    group1 <- comparisons$group1[k]
    group2 <- comparisons$group2[k]
    p_matrix[group1, group2] <- p_matrix[group2, group1] <-
      comparisons$p.adj[k]
  }

  letters <- make_letters(
    p_matrix,
    setNames(summary_data$Mean, treatments),
    alpha
  )

  event <- as.character(original_data$Event[1])

  summary_data <- summary_data |>
    mutate(
      Letter = unname(letters[as.character(Treatment)]),
      Site = site,
      Sampling = sampling,
      Event = event,
      Parameter = parameter
    )

  welch_table <- tibble(
    Site = site,
    Sampling = sampling,
    Event = event,
    Parameter = parameter,
    F = unname(welch$statistic),
    df1 = unname(welch$parameter[1]),
    df2 = unname(welch$parameter[2]),
    p.value = welch$p.value
  )

  list(
    name = model_name,
    n_excluded = nrow(original_data) - nrow(dat),
    counts = counts,
    # This model supplies descriptive diagnostics only.
    model = lm(Value ~ Treatment, data = dat),
    welch = welch,
    welch_table = welch_table,
    comparisons = comparisons,
    summary = summary_data
  )
}


## Run all within-sampling comparisons ----

soil_results <- list()

for (i in seq_len(nrow(soil_key))) {
  for (parameter in parameter_levels) {
    site <- soil_key$Site[i]
    sampling <- soil_key$Sampling[i]
    key <- paste(site, sampling, parameter, sep = "_")

    soil_results[[key]] <- analyse_group(site, sampling, parameter)
  }
}

soil_summary <- map_dfr(soil_results, "summary") |>
  mutate(
    Treatment = factor(Treatment, levels = treatment_levels),
    Event = factor(Event, levels = event_labels),
    Parameter = factor(Parameter, levels = parameter_levels)
  )

welch_table <- map_dfr(soil_results, "welch_table")
soil_comparisons <- collect_results(soil_results, "comparisons")


## Print and export statistical results ----

soil_report <- capture.output({
  cat(
    "EXPLORATORY SOIL PARAMETER COMPARISONS\n\n",
    "Separate Welch one-way ANOVAs within each\n",
    "site × sampling event × parameter.\n",
    "Required input file: ", soil_input_file, " in dat/raw.\n",
    "Pairwise comparisons: Games–Howell among four treatments.\n",
    "Pairwise estimates: mean(group1) minus mean(group2).\n",
    sep = ""
  )

  cat("Confidence level:", 100 * (1 - alpha), "%.\n")
  cat(
    "No adjustment across the ", length(soil_results),
    " overall tests or their families.\n",
    sep = ""
  )

  cat(
    "Bars: observed means ± SE, in original units.\n",
    "Water content: 100 - pc_TS, as % fresh sample mass.\n",
    "DOC/DON conversions retained from the supplied script.\n",
    "Letters assigned in descending mean order within each event.\n",
    "Treatments sharing a letter have no detected difference;\n",
    "shared letters do not establish equivalence.\n",
    "Letters are not gated by the overall Welch p-value.\n",
    "Residual plots are descriptive diagnostics only.\n",
    sep = ""
  )

  cat("Letters shown only when a Games–Howell p.adj < ", alpha, ".\n", sep = "")

  for (key in names(soil_results)) {
    result <- soil_results[[key]]
    print_model_heading(key, result$n_excluded)

    print_section("Welch one-way ANOVA:", result$welch)

    cat("\nGames–Howell pairwise comparisons:\n")
    print(as.data.frame(result$comparisons), row.names = FALSE, digits = 6)

    cat("\nObserved means, SD, SE, and letters:\n")
    print(as.data.frame(result$summary), row.names = FALSE, digits = 6)
  }

  cat("\nAll Welch results:\n")
  print(as.data.frame(welch_table), row.names = FALSE, digits = 6)
})

cat(paste(soil_report, collapse = "\n"), "\n")

ensure_dir(dir_reporting)
writeLines(
  soil_report,
  file.path(dir_reporting, soil_report_file)
)

write_table(welch_table, soil_welch_file, dir_reporting)
write_table(soil_comparisons, soil_comparisons_file, dir_reporting)
write_table(soil_summary, soil_summary_file, dir_reporting)


## Export diagnostic plots ----

save_diagnostic_pdf(soil_results, soil_diagnostics_file)

for (key in names(soil_results)) {
  save_diagnostic_png(soil_results[[key]], key)
}


## Create Figure S9 panels ----

make_soil_panel <- function(site, parameter) {
  base_max <- unname(axis_maxima[[site]][parameter])

  dat <- soil_summary |>
    filter(Site == site, Parameter == parameter) |>
    mutate(Letter_y = Mean + SE + 0.04 * base_max)

  upper_limit <- max(
    base_max,
    dat$Mean + dat$SE + 0.03 * base_max,
    dat$Letter_y[dat$Letter != ""] + 0.08 * base_max
  )

  show_x <- parameter == "GWC"
  show_y_title <- site == "KL"
  dodge <- position_dodge(width = 0.8)

  ggplot(dat, aes(Event, Mean, fill = Treatment, group = Treatment)) +
    geom_col(
      position = dodge, width = 0.8,
      colour = "black", linewidth = 0.25
    ) +
    geom_errorbar(
      aes(ymin = Mean - SE, ymax = Mean + SE),
      position = dodge, width = 0.15, linewidth = 0.3
    ) +
    geom_text(
      aes(y = Letter_y, label = Letter),
      position = dodge, vjust = 0, size = 2.8,
      family = "Source Sans 3"
    ) +
    scale_fill_manual(
      name = "Treatment",
      values = colours_soil_parameter_treatment,
      breaks = treatment_levels,
      drop = FALSE
    ) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(add = 0.5),
      labels = \(x) str_remove(x, "^Sampling ")
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
      breaks = breaks_pretty(n = 4)
    ) +
    coord_cartesian(ylim = c(0, upper_limit), clip = "on") +
    labs(x = NULL, y = if (show_y_title) axis_labels[[parameter]] else NULL) +
    theme_classic(base_size = 10, base_family = "Source Sans 3") +
    theme(
      axis.title.y = if (show_y_title) {
        element_text(size = 10, angle = 90, margin = margin(r = 5))
      } else {
        element_blank()
      },
      axis.text.x = if (show_x) {
        element_text(size = 8, angle = 0)
      } else {
        element_blank()
      },
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      plot.margin = margin(10, 8, 8, 8),
      plot.tag = element_text(size = 13)
    )
}


## Assemble and save Figure S9 ----

# Rows: parameters; columns: KL, ZOE.
soil_panels <- list()

for (parameter in parameter_levels) {
  for (site in c("KL", "ZOE")) {
    soil_panels[[paste(site, parameter, sep = "_")]] <-
      make_soil_panel(site, parameter)
  }
}

soil_figure <- wrap_plots(
  soil_panels,
  ncol = 2, byrow = TRUE, guides = "collect"
) +
  plot_annotation(
    tag_levels = "A",
    caption = "Sampling Event",
    theme = theme(
      plot.caption = element_text(hjust = 0.5, size = 13)
    )
  ) &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

save_figure(
  soil_figure,
  figure_s09_name,
  width_mm = figure_s09_width_mm,
  height_mm = figure_s09_height_mm
)
log_script_end("figure-s09-soil-parameters")
