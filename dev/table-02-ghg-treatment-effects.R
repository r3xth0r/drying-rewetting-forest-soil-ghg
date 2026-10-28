# Table 2: GHG treatment effects ----
# Fits plot-level LMMs with CAR(1) correlation and treatment-specific variance.
# Inputs: five raw GHG flux datasets. Outputs: Tables 2 and S2, diagnostics.

suppressPackageStartupMessages({
  library("tidyverse")
  library("nlme")
  library("emmeans")
  library("multcompView")
  library("here")
})


source(here::here("dev", "logging.R"))
log_script_start("table-02-ghg-treatment-effects")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))

# Settings ====

config <- tibble(
  ID = c("CO2_KL", "CO2_ZB", "CH4_KL", "CH4_ZB", "N2O_KL"),
  GHG = c("CO2", "CO2", "CH4", "CH4", "N2O"),
  Site = c("KL", "ZB", "KL", "ZB", "KL"),
  File = c(
    "co2-klausen.csv",
    "co2-zoe.csv",
    "ch4-klausen.csv",
    "ch4-zoe.csv",
    "n2o-klausen.csv"
  ),
  Year = c(2021L, 2022L, 2021L, 2022L, 2021L),
  EndDate = as.Date(c(
    "2021-11-16", "2022-11-15",
    "2021-11-16", "2022-11-15",
    "2021-11-16"
  )),
  Sign = c(1, 1, -1, -1, 1),
  ExcludePlots = list(
    c(2L, 3L, 10L), # CO2 KL
    integer(0), # CO2 ZB: keep ZB15
    10L, # CH4 KL
    integer(0), # CH4 ZB: keep ZB15
    10L # N2O KL
  ),
  ExpectedPlots = c(13L, 16L, 15L, 16L, 15L)
)

period_levels <- c(
  "Pre-treatment",
  "DRW Cycle 1",
  "DRW Cycle 2",
  "DRW Cycle 3",
  "Final Period",
  "All cycles"
)

period_boundaries <- list(
  KL = as.Date(c(
    "2021-05-04", "2021-07-06", "2021-08-31",
    "2021-10-26", "2021-11-16"
  )),
  ZB = as.Date(c(
    "2022-06-26", "2022-08-22", "2022-10-16", "2022-11-15"
  ))
)

zb_start_date <- as.Date("2022-04-25")
max_calendar_lag <- 14L
diagnostic_pdf_width_mm <- 280
diagnostic_pdf_height_mm <- 200
event_date <- as.Date("2021-10-20")
table_p_column <- "p_raw"
display_digits <- 2L
model_version <- "date-fixed-plot-car1-varident-v1"
force_refit <- FALSE
reuse_previous_n2o_test <- TRUE
verbose_optimizer <- FALSE
confidence_level <- 0.95
comparison_threshold <- 0.05
small_plot_sd_ratio_threshold <- 0.001
correlation_start_value <- 0.8
correlation_boundary_tolerance <- 1e-7
diagnostic_point_alpha <- 0.3
diagnostic_point_shape <- 16
diagnostic_point_size <- 0.5
diagnostic_point_colour <- "black"
diagnostic_qq_line_colour <- "red"

checkpoint_file_prefix <- "table-02-model-checkpoint"
completed_models_file <- "table-02-models.rds"
daily_data_file <- "table-02-daily-plot-data.csv"
date_check_file <- "table-02-date-coverage.csv"
plot_means_file <- "table-02-observed-plot-means.csv"
treatment_descriptives_file <- "table-02-treatment-descriptives.csv"
water_descriptives_file <- "table-02-water-descriptives.csv"
fit_errors_file <- "table-02-fit-errors.csv"
treatment_means_file <- "table-02-adjusted-treatment-means.csv"
pairwise_file <- "table-02-fisher-lsd-comparisons.csv"
water_means_file <- "table-02-adjusted-water-means.csv"
water_comparisons_file <- "table-02-water-comparisons.csv"
model_status_file <- "table-02-model-status.csv"
variance_parameters_file <- "table-02-residual-variance-parameters.csv"
residual_lags_file <- "table-02-residual-calendar-lags.csv"
manuscript_table_file <- "table-02-ghg-treatment-effects.csv"
notes_file <- "table-02-model-notes.txt"
session_info_file <- "table-02-session-info.txt"
lopo_file <- "table-s02-n2o-kl-cycle-1-lopo.csv"
excel_table_file <- "table-02-ghg-treatment-effects-excel.csv"
excel_letters_table_file <-
  "table-02-ghg-treatment-effects-letters-highest-first-excel.csv"
excel_one_decimal_table_file <-
  "table-02-ghg-treatment-effects-one-decimal-excel.csv"
excel_output_directory <- dir_reporting
treatment_columns <- c("Control", "Control + N", "DRW", "DRW + N")
value_columns <- c(
  treatment_columns,
  "Natural Precipitation",
  "Drying–rewetting"
)
event_model_file <- "table-02-rewetting-day-model.rds"
event_comparisons_file <- "table-02-rewetting-day-comparisons.csv"
event_lags_file <- "table-02-rewetting-day-calendar-lags.csv"
event_diagnostics_file <- "table-02-rewetting-day-diagnostics.csv"
progress_file <- file.path(dir_interim, "table-02-progress-log.txt")

previous_test_directory <- file.path(
  dir_interim,
  "table-02-previous-n2o-test"
)

model_control <- lmeControl(
  opt = "nlminb",
  maxIter = 50,
  msMaxIter = 100,
  msMaxEval = 200,
  niterEM = 25,
  msVerbose = verbose_optimizer,
  returnObject = FALSE
)

# These limits bound iterations, not elapsed wall-clock time.

ensure_dir(dir_interim)
ensure_dir(dir_processed)
ensure_dir(dir_reporting)
ensure_dir(file.path(dir_interim, "table-02-diagnostics"))

# Logging and checkpoints ====

log_progress <- function(text) {
  line <- paste(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    text
  )
  logger::log_info(line)
  cat(line, "\n", file = progress_file, append = TRUE)
}

safe_save_rds <- function(object, destination) {
  temporary <- tempfile(
    pattern = "checkpoint_",
    tmpdir = dirname(destination),
    fileext = ".rds"
  )

  saveRDS(object, temporary)

  # Confirm that the newly written file can be read.
  invisible(readRDS(temporary))

  if (!file.copy(temporary, destination, overwrite = TRUE)) {
    stop(
      "Could not save checkpoint. Temporary copy remains at: ",
      temporary
    )
  }

  unlink(temporary)
}


parse_sample_date <- function(x) {
  x <- trimws(as.character(x))

  # Explicit day/month/year parsing; discard time.
  date_text <- sub("[[:space:]T].*$", "", x)
  date_text <- gsub(".", "/", date_text, fixed = TRUE)

  valid_format <- grepl(
    "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}$",
    date_text
  )

  parsed <- as.Date(rep(NA_character_, length(x)))

  use <- !is.na(x) & nzchar(x) & valid_format

  parsed[use] <- as.Date(
    date_text[use],
    format = "%d/%m/%Y"
  )

  bad <- !is.na(x) & nzchar(x) &
    (!valid_format | is.na(parsed))

  if (any(bad)) {
    stop(
      "Invalid date format: ",
      paste(head(unique(x[bad]), 5), collapse = "; ")
    )
  }

  parsed
}

assign_period <- function(date, site) {
  boundaries <- period_boundaries[[site]]

  if (site == "KL") {
    case_when(
      date < boundaries[[1]] ~ "Pre-treatment",
      date <= boundaries[[2]] ~ "DRW Cycle 1",
      date <= boundaries[[3]] ~ "DRW Cycle 2",
      date <= boundaries[[4]] ~ "DRW Cycle 3",
      date <= boundaries[[5]] ~ "Final Period",
      TRUE ~ NA_character_
    )
  } else {
    case_when(
      date <= boundaries[[1]] ~ "DRW Cycle 1",
      date <= boundaries[[2]] ~ "DRW Cycle 2",
      date <= boundaries[[3]] ~ "DRW Cycle 3",
      date <= boundaries[[4]] ~ "Final Period",
      TRUE ~ NA_character_
    )
  }
}

read_dataset <- function(cfg) {
  file <- file.path(dir_raw, cfg$File[[1]])
  raw <- read_input(cfg$File[[1]], check.names = FALSE)

  names(raw) <- trimws(names(raw))

  if (!"Date" %in% names(raw)) {
    stop("No Date column in ", file)
  }

  plot_columns <- grep(
    "^(KL|ZOE|ZB)[ ._-]*[0-9]+",
    names(raw),
    value = TRUE,
    ignore.case = TRUE
  )

  if (!length(plot_columns)) {
    stop("No recognised plot columns in ", file)
  }

  raw[plot_columns] <- lapply(raw[plot_columns], as.character)

  long <- raw |>
    select(Date, all_of(plot_columns)) |>
    pivot_longer(
      -Date,
      names_to = "OriginalPlot",
      values_to = "RawFlux"
    ) |>
    mutate(
      SampleDate = parse_sample_date(Date),
      PlotNumber = as.integer(
        sub("^[^0-9]*([0-9]+).*$", "\\1", OriginalPlot)
      ),
      Flux = suppressWarnings(as.numeric(RawFlux))
    )

  invalid <- !is.na(long$RawFlux) &
    nzchar(trimws(long$RawFlux)) &
    is.na(long$Flux)

  if (any(invalid)) {
    stop(
      cfg$ID[[1]], ": non-numeric flux entries: ",
      paste(head(unique(long$RawFlux[invalid]), 5), collapse = "; ")
    )
  }

  if (any(!is.finite(long$Flux) & !is.na(long$Flux))) {
    stop(cfg$ID[[1]], ": infinite flux values.")
  }

  plot_map <- distinct(long, OriginalPlot, PlotNumber)

  if (anyDuplicated(plot_map$PlotNumber)) {
    stop(cfg$ID[[1]], ": multiple columns map to the same plot.")
  }

  if (any(!plot_map$PlotNumber %in% 1:16)) {
    stop(cfg$ID[[1]], ": unexpected plot numbers.")
  }

  out <- long |>
    filter(
      !is.na(SampleDate),
      format(SampleDate, "%Y") == as.character(cfg$Year[[1]]),
      SampleDate <= cfg$EndDate[[1]],
      !PlotNumber %in% cfg$ExcludePlots[[1]]
    )

  if (cfg$Site[[1]] == "ZB") {
    out <- out |>
      filter(SampleDate >= zb_start_date)
  }

  out |>
    mutate(Flux = cfg$Sign[[1]] * Flux) |>
    group_by(PlotNumber, SampleDate) |>
    summarise(
      Flux = finite_mean(Flux),
      .groups = "drop"
    ) |>
    mutate(
      ID = cfg$ID[[1]],
      GHG = cfg$GHG[[1]],
      Site = cfg$Site[[1]],
      Plot = paste0(cfg$Site[[1]], sprintf("%02d", PlotNumber)),
      Treatment = case_when(
        PlotNumber %in% c(1, 5, 10, 13) ~ "C",
        PlotNumber %in% c(2, 6, 9, 14) ~ "C + N",
        PlotNumber %in% c(3, 7, 11, 15) ~ "DRW",
        PlotNumber %in% c(4, 8, 12, 16) ~ "DRW + N"
      ),
      Treatment = factor(Treatment, levels = treatment_levels),
      Water = factor(
        ifelse(Treatment %in% c("C", "C + N"), "NP", "DRW"),
        levels = c("NP", "DRW")
      ),
      Period = assign_period(SampleDate, cfg$Site[[1]])
    ) |>
    filter(!is.na(Period))
}


# Data preparation ====

daily_data <- map_dfr(
  seq_len(nrow(config)),
  function(i) read_dataset(config[i, ])
)

observed_data <- daily_data |>
  filter(is.finite(Flux))

plot_check <- observed_data |>
  group_by(ID, GHG, Site) |>
  summarise(
    Nplots = n_distinct(Plot),
    Observations = n(),
    FirstDate = min(SampleDate),
    LastDate = max(SampleDate),
    .groups = "drop"
  ) |>
  left_join(
    select(config, ID, ExpectedPlots),
    by = "ID"
  )

print(plot_check, n = Inf)

if (
  nrow(plot_check) != nrow(config) ||
    any(plot_check$Nplots != plot_check$ExpectedPlots)
) {
  stop("Unexpected dataset or plot counts. Inspect plot_check.")
}

date_check <- observed_data |>
  group_by(ID, Period) |>
  summarise(
    Nplots = n_distinct(Plot),
    Observations = n(),
    DistinctDays = n_distinct(SampleDate),
    FirstDate = min(SampleDate),
    LastDate = max(SampleDate),
    .groups = "drop"
  )

print(date_check, n = Inf)

write_csv(
  daily_data,
  file.path(dir_interim, daily_data_file)
)

write_csv(
  date_check,
  file.path(dir_interim, date_check_file)
)

# Add a separate treatment-period dataset for All cycles.
# Each individual model still contains each plot-day only once.
analysis_data <- bind_rows(
  observed_data,
  observed_data |>
    filter(Period != "Pre-treatment") |>
    mutate(Period = "All cycles")
)


# Descriptive summaries ====

plot_means <- analysis_data |>
  group_by(ID, GHG, Site, Period, Treatment, Water, Plot) |>
  summarise(
    PlotMean = mean(Flux),
    ObservedDays = n(),
    .groups = "drop"
  )

treatment_descriptives <- plot_means |>
  group_by(ID, GHG, Site, Period, Treatment) |>
  summarise(
    Mean = mean(PlotMean),
    SE = standard_error(PlotMean),
    Nplots = n(),
    .groups = "drop"
  )

water_descriptives <- plot_means |>
  group_by(ID, GHG, Site, Period, Water) |>
  summarise(
    Mean = mean(PlotMean),
    SE = standard_error(PlotMean),
    Nplots = n(),
    .groups = "drop"
  )

write_csv(
  plot_means,
  file.path(dir_processed, plot_means_file)
)

write_csv(
  treatment_descriptives,
  file.path(dir_processed, treatment_descriptives_file)
)

write_csv(
  water_descriptives,
  file.path(dir_processed, water_descriptives_file)
)


# Model fitting ====

prepare_model_data <- function(dat) {
  dat |>
    mutate(
      Plot = factor(Plot),
      DateGroup = factor(SampleDate),
      Day = as.numeric(SampleDate - min(SampleDate)),
      Treatment = factor(Treatment, levels = treatment_levels),
      Water = factor(Water, levels = c("NP", "DRW"))
    ) |>
    droplevels() |>
    arrange(Plot, Day)
}

# Standardized representation for checking old versus new input data.
canonical_data <- function(dat) {
  dat <- dat[order(as.character(dat$Plot), dat$SampleDate), ]

  data.frame(
    Plot = as.character(dat$Plot),
    SampleDate = as.character(dat$SampleDate),
    Day = as.numeric(dat$Day),
    DateGroup = as.character(dat$DateGroup),
    Treatment = as.character(dat$Treatment),
    Water = as.character(dat$Water),
    Flux = as.numeric(dat$Flux),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

formula_text <- function(x) {
  gsub(
    "[[:space:]]+",
    "",
    paste(deparse(x), collapse = "")
  )
}

fit_one_model <- function(dat, model_type) {
  dat <- prepare_model_data(dat)

  if (anyDuplicated(dat[c("Plot", "Day")])) {
    stop("Duplicate plot-day observations.")
  }

  if (
    nlevels(dat$Treatment) != 4L ||
      nlevels(dat$Water) != 2L
  ) {
    stop("A treatment or water group is missing.")
  }

  if (model_type == "FourTreatments") {
    fixed_formula <- Flux ~ Treatment + DateGroup
    variance_formula <- ~ 1 | Treatment
  } else {
    fixed_formula <- Flux ~ Water + DateGroup
    variance_formula <- ~ 1 | Water
  }

  X <- model.matrix(fixed_formula, data = dat)

  if (qr(X)$rank < ncol(X)) {
    stop("Fixed effects are rank deficient.")
  }

  id <- unique(dat$ID)
  period <- unique(dat$Period)

  stopifnot(length(id) == 1L, length(period) == 1L)

  tag <- str_to_lower(str_replace_all(
    paste(id, period, model_type),
    "[^A-Za-z0-9]+",
    "-"
  ))

  checkpoint_file <- file.path(
    dir_interim,
    paste0(checkpoint_file_prefix, "-", tag, ".rds")
  )

  settings <- list(
    Version = model_version,
    FixedFormula = formula_text(fixed_formula),
    VarianceFormula = formula_text(variance_formula),
    Random = "~1|Plot",
    Correlation = "~Day|Plot",
    Method = "REML",
    Control = model_control,
    RVersion = R.version.string,
    NlmeVersion = as.character(packageVersion("nlme"))
  )

  if (!force_refit && file.exists(checkpoint_file)) {
    saved <- tryCatch(
      readRDS(checkpoint_file),
      error = function(e) NULL
    )

    if (
      is.list(saved) &&
        inherits(saved$model, "lme") &&
        identical(saved$settings, settings) &&
        identical(canonical_data(saved$data), canonical_data(dat))
    ) {
      log_progress(paste(tag, "restored."))
      return(saved)
    }
  }

  # Attempt to import the already-completed N2O Cycle 1 candidate.
  if (
    !force_refit &&
      reuse_previous_n2o_test &&
      id == "N2O_KL" &&
      period == "DRW Cycle 1"
  ) {
    previous_file <- file.path(
      previous_test_directory,
      paste0("N2O_KL_Cycle1_", model_type, ".rds")
    )

    if (file.exists(previous_file)) {
      previous <- tryCatch(
        readRDS(previous_file),
        error = function(e) NULL
      )

      compatible <- tryCatch(
        {
          is.list(previous) &&
            inherits(previous$model, "lme") &&
            identical(
              canonical_data(previous$data),
              canonical_data(dat)
            ) &&
            identical(
              formula_text(formula(previous$model)),
              formula_text(fixed_formula)
            ) &&
            inherits(
              previous$model$modelStruct$corStruct,
              "corCAR1"
            ) &&
            identical(
              formula_text(formula(
                previous$model$modelStruct$corStruct
              )),
              formula_text(~ Day | Plot)
            ) &&
            inherits(
              previous$model$modelStruct$varStruct,
              "varIdent"
            ) &&
            identical(
              formula_text(formula(
                previous$model$modelStruct$varStruct
              )),
              formula_text(variance_formula)
            ) &&
            identical(previous$model$method, "REML")
        },
        error = function(e) FALSE
      )

      if (isTRUE(compatible)) {
        answer <- list(
          model = previous$model,
          data = dat,
          settings = settings,
          ID = id,
          Period = period,
          model_type = model_type,
          warnings = "",
          source = "Imported successful N2O Cycle 1 test"
        )

        safe_save_rds(answer, checkpoint_file)
        log_progress(paste(tag, "imported previous successful test."))
        return(answer)
      }
    }
  }

  log_progress(paste(tag, "started."))
  started <- Sys.time()
  warning_messages <- character()

  model <- withCallingHandlers(
    lme(
      fixed = fixed_formula,
      random = ~ 1 | Plot,
      correlation = corCAR1(
        value = correlation_start_value,
        form = ~ Day | Plot
      ),
      weights = varIdent(form = variance_formula),
      data = dat,
      method = "REML",
      na.action = na.fail,
      control = model_control
    ),
    warning = function(w) {
      warning_messages <<- c(
        warning_messages,
        conditionMessage(w)
      )
    }
  )

  phi <- unname(coef(
    model$modelStruct$corStruct,
    unconstrained = FALSE
  ))

  if (
    !is.finite(phi) ||
      phi >= 1 - correlation_boundary_tolerance
  ) {
    stop("CAR1 estimate reached its upper boundary.")
  }

  answer <- list(
    model = model,
    data = dat,
    settings = settings,
    ID = id,
    Period = period,
    model_type = model_type,
    warnings = paste(unique(warning_messages), collapse = " | "),
    source = "New fit"
  )

  safe_save_rds(answer, checkpoint_file)

  elapsed <- as.numeric(
    difftime(Sys.time(), started, units = "mins")
  )

  log_progress(sprintf(
    "%s saved after %.2f minutes.",
    tag, elapsed
  ))

  answer
}


# Model results ====

jobs <- analysis_data |>
  distinct(ID, GHG, Site, Period) |>
  mutate(
    DatasetOrder = match(ID, config$ID),
    PeriodOrder = match(Period, period_levels)
  ) |>
  arrange(DatasetOrder, PeriodOrder)

four_fits <- list()
water_fits <- list()
fit_errors <- list()

for (i in seq_len(nrow(jobs))) {
  id <- jobs$ID[i]
  period <- jobs$Period[i]
  key <- paste(id, period, sep = "__")

  dat <- analysis_data |>
    filter(ID == id, Period == period)

  for (model_type in c("FourTreatments", "PooledWater")) {
    result <- tryCatch(
      fit_one_model(dat, model_type),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      fit_errors[[paste(key, model_type)]] <- tibble(
        ID = id,
        Period = period,
        Model = model_type,
        Error = conditionMessage(result)
      )

      write_csv(
        bind_rows(fit_errors),
        file.path(dir_interim, fit_errors_file)
      )

      log_progress(paste(
        id, period, model_type,
        "FAILED:", conditionMessage(result)
      ))
    } else if (model_type == "FourTreatments") {
      four_fits[[key]] <- result
    } else {
      water_fits[[key]] <- result
    }
  }
}

safe_save_rds(
  list(
    FourTreatments = four_fits,
    PooledWater = water_fits
  ),
  file.path(dir_interim, completed_models_file)
)

if (length(fit_errors)) {
  print(
    bind_rows(fit_errors),
    n = Inf,
    width = Inf
  )

  stop(
    "Some models failed. Completed fits remain saved. ",
    "Do not increase iteration limits automatically; inspect the failures."
  )
}


# Treatment comparisons ====

extract_four <- function(fit) {
  emm <- emmeans(
    fit$model,
    ~Treatment,
    data = fit$data,
    weights = "equal",
    mode = "containment"
  )

  means <- as.data.frame(summary(
    emm,
    infer = c(TRUE, FALSE)
  ))

  # Fisher's LSD pairwise comparisons:
  # no multiplicity adjustment is applied.
  pairwise <- as.data.frame(summary(
    pairs(
      emm,
      reverse = FALSE,
      adjust = "none"
    ),
    infer = c(TRUE, TRUE),
    adjust = "none"
  ))

  labels <- as.character(means$Treatment)
  combinations <- combn(seq_along(labels), 2)

  stopifnot(nrow(pairwise) == ncol(combinations))

  if (anyNA(pairwise$p.value)) {
    stop("Non-estimable treatment comparison.")
  }

  # pairs() uses the same pair order as combinations.
  p_matrix <- matrix(
    1,
    nrow = length(labels),
    ncol = length(labels),
    dimnames = list(labels, labels)
  )

  for (j in seq_len(ncol(combinations))) {
    a <- combinations[1, j]
    b <- combinations[2, j]

    p_matrix[a, b] <- pairwise$p.value[j]
    p_matrix[b, a] <- pairwise$p.value[j]
  }

  # Compact-letter display based on Fisher LSD p-values.
  letter_values <- multcompLetters(
    p_matrix,
    threshold = comparison_threshold
  )$Letters

  list(
    means = as_tibble(means) |>
      mutate(
        ID = fit$ID,
        Period = fit$Period,
        .before = 1
      ),
    pairwise = as_tibble(pairwise) |>
      mutate(
        ID = fit$ID,
        Period = fit$Period,
        .before = 1
      ),
    letters = tibble(
      ID = fit$ID,
      Period = fit$Period,
      Treatment = names(letter_values),
      Letters = unname(letter_values)
    )
  )
}

four_results <- map(four_fits, extract_four)

adjusted_treatment_means <- map_dfr(
  four_results,
  "means"
)

pairwise_results <- map_dfr(
  four_results,
  "pairwise"
)

letters_results <- map_dfr(
  four_results,
  "letters"
)


# Pooled water comparisons ====

extract_water <- function(fit) {
  emm <- emmeans(
    fit$model,
    ~Water,
    data = fit$data,
    weights = "equal",
    mode = "containment"
  )

  means <- as_tibble(summary(
    emm,
    infer = c(TRUE, FALSE)
  )) |>
    mutate(
      ID = fit$ID,
      Period = fit$Period,
      .before = 1
    )

  comparison <- as_tibble(summary(
    pairs(emm, reverse = TRUE, adjust = "none"),
    infer = c(TRUE, TRUE),
    adjust = "none"
  )) |>
    rename(p_raw = p.value) |>
    mutate(
      ID = fit$ID,
      Period = fit$Period,
      .before = 1
    )

  list(means = means, comparison = comparison)
}

pooled_results <- map(water_fits, extract_water)

adjusted_water_means <- map_dfr(pooled_results, "means")

water_results <- map_dfr(pooled_results, "comparison") |>
  group_by(ID) |>
  mutate(p_holm = p.adjust(p_raw, method = "holm")) |>
  ungroup()

write_csv(
  adjusted_treatment_means,
  file.path(dir_processed, treatment_means_file)
)

write_csv(
  pairwise_results,
  file.path(dir_processed, pairwise_file)
)

write_csv(
  adjusted_water_means,
  file.path(dir_processed, water_means_file)
)

write_csv(
  water_results,
  file.path(dir_processed, water_comparisons_file)
)


# Model diagnostics ====

calendar_lags <- function(dat, max_lag = max_calendar_lag) {
  map_dfr(seq_len(max_lag), function(k) {
    later <- dat |>
      transmute(
        Plot,
        Day = Day - k,
        LaterResidual = NormalizedResidual
      )

    matched <- dat |>
      select(Plot, Day, NormalizedResidual) |>
      inner_join(later, by = c("Plot", "Day")) |>
      filter(
        is.finite(NormalizedResidual),
        is.finite(LaterResidual)
      )

    r <- if (
      nrow(matched) >= 3L &&
        sd(matched$NormalizedResidual) > 0 &&
        sd(matched$LaterResidual) > 0
    ) {
      cor(matched$NormalizedResidual, matched$LaterResidual)
    } else {
      NA_real_
    }

    tibble(
      LagDays = k,
      MatchedPairs = nrow(matched),
      Correlation = r
    )
  })
}

save_diagnostic_pdf <- function(dat, lag_table, file, title) {
  if (!save_pdf_figures) {
    return(invisible(NULL))
  }

  pdf(
    file,
    width = mm_to_inches(diagnostic_pdf_width_mm),
    height = mm_to_inches(diagnostic_pdf_height_mm)
  )
  on.exit(dev.off(), add = TRUE)

  par(mfrow = c(2, 2))

  point_colour <- adjustcolor(
    diagnostic_point_colour,
    alpha.f = diagnostic_point_alpha
  )

  plot(
    dat$FittedFlux,
    dat$NormalizedResidual,
    pch = diagnostic_point_shape,
    cex = diagnostic_point_size,
    col = point_colour,
    xlab = "Conditional fitted flux",
    ylab = "Normalized residual",
    main = title
  )
  abline(h = 0, lty = 2)

  qqnorm(
    dat$NormalizedResidual,
    pch = diagnostic_point_shape,
    cex = diagnostic_point_size,
    col = point_colour,
    main = "Normalized residual QQ"
  )
  qqline(dat$NormalizedResidual, col = diagnostic_qq_line_colour)

  plot(
    dat$SampleDate,
    dat$NormalizedResidual,
    pch = diagnostic_point_shape,
    cex = diagnostic_point_size,
    col = point_colour,
    xlab = "Date",
    ylab = "Normalized residual"
  )
  abline(h = 0, lty = 2)

  plot(
    lag_table$LagDays,
    lag_table$Correlation,
    type = "b",
    xlab = "Elapsed days",
    ylab = "Matched residual correlation"
  )
  abline(h = 0, lty = 2)
}

diagnose_fit <- function(fit) {
  model <- fit$model

  tag <- str_to_lower(str_replace_all(
    paste(fit$ID, fit$Period, fit$model_type),
    "[^A-Za-z0-9]+",
    "-"
  ))

  dat <- fit$data |>
    mutate(
      FittedFlux = as.numeric(fitted(model, level = 1)),
      NormalizedResidual = as.numeric(residuals(
        model,
        type = "normalized"
      )),
      ResponseResidual = as.numeric(residuals(
        model,
        type = "response"
      ))
    )

  lags <- calendar_lags(dat)

  z <- dat$NormalizedResidual
  scaled_z <- (z - mean(z)) / sd(z)

  vc <- VarCorr(model)
  plot_sd <- as.numeric(vc[1, "StdDev"])

  phi <- unname(coef(
    model$modelStruct$corStruct,
    unconstrained = FALSE
  ))

  status <- tibble(
    ID = fit$ID,
    Period = fit$Period,
    Model = fit$model_type,
    Observations = nrow(dat),
    Nplots = n_distinct(dat$Plot),
    Ndates = n_distinct(dat$SampleDate),
    Phi = phi,
    PlotSD = plot_sd,
    BaselineResidualSD = model$sigma,

    # Screening flag only: not a formal singularity test.
    VerySmallPlotSD =
      plot_sd / model$sigma < small_plot_sd_ratio_threshold,
    Lag1 = lags$Correlation[lags$LagDays == 1],
    Lag2 = lags$Correlation[lags$LagDays == 2],
    ResidualMean = mean(z),
    ResidualSD = sd(z),
    ResidualSkewness = mean(scaled_z^3),
    ResidualExcessKurtosis = mean(scaled_z^4) - 3,
    MaximumAbsoluteResidual = max(abs(z)),
    Warnings = fit$warnings,
    Source = fit$source
  )

  ratios <- coef(
    model$modelStruct$varStruct,
    unconstrained = FALSE,
    allCoef = TRUE
  )

  variance_parameters <- tibble(
    ID = fit$ID,
    Period = fit$Period,
    Model = fit$model_type,
    Group = names(ratios),
    ResidualSDMultiplier = as.numeric(ratios),
    ResidualSD = model$sigma * as.numeric(ratios)
  )

  save_diagnostic_pdf(
    dat,
    lags,
    file.path(
      file.path(dir_interim, "table-02-diagnostics"),
      paste0("lmm-diagnostics-", tag, ".pdf")
    ),
    tag
  )

  write_csv(
    dat,
    file.path(
      dir_interim,
      paste0("table-02-diagnostics-", tag, "-residuals.csv")
    )
  )

  capture.output(
    summary(model),
    file = file.path(
      dir_interim,
      paste0("table-02-diagnostics-", tag, "-model-summary.txt")
    )
  )

  list(
    status = status,
    variance = variance_parameters,
    lag = lags |>
      mutate(
        ID = fit$ID,
        Period = fit$Period,
        Model = fit$model_type,
        .before = 1
      )
  )
}

diagnostic_results <- map(
  c(unname(four_fits), unname(water_fits)),
  diagnose_fit
)

model_status <- map_dfr(diagnostic_results, "status")
variance_parameters <- map_dfr(diagnostic_results, "variance")
residual_lag_results <- map_dfr(diagnostic_results, "lag")

write_csv(
  model_status,
  file.path(dir_interim, model_status_file)
)

write_csv(
  variance_parameters,
  file.path(dir_interim, variance_parameters_file)
)

write_csv(
  residual_lag_results,
  file.path(dir_interim, residual_lags_file)
)


# Manuscript tables ====

format_mean_se <- function(mean_value, se_value) {
  sprintf(
    paste0("%.", display_digits, "f ± %.", display_digits, "f"),
    mean_value, se_value
  )
}

format_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.4f", p))
  )
}

treatment_cells <- treatment_descriptives |>
  mutate(Treatment = as.character(Treatment)) |>
  left_join(
    letters_results,
    by = c("ID", "Period", "Treatment")
  ) |>
  mutate(
    Cell = paste(format_mean_se(Mean, SE), Letters),
    Treatment = recode(
      Treatment,
      "C" = "Control",
      "C + N" = "Control + N"
    )
  ) |>
  select(ID, GHG, Site, Period, Treatment, Cell) |>
  pivot_wider(
    names_from = Treatment,
    values_from = Cell
  )

water_cells <- water_descriptives |>
  mutate(
    Cell = format_mean_se(Mean, SE),
    Water = recode(
      as.character(Water),
      "NP" = "Natural Precipitation",
      "DRW" = "Drying–rewetting"
    )
  ) |>
  select(ID, Period, Water, Cell) |>
  pivot_wider(
    names_from = Water,
    values_from = Cell
  )

pooled_p_values <- water_results |>
  transmute(
    ID,
    Period,
    `p-value` = format_p(.data[[table_p_column]])
  )

manuscript_table <- treatment_cells |>
  left_join(water_cells, by = c("ID", "Period")) |>
  left_join(pooled_p_values, by = c("ID", "Period")) |>
  mutate(
    GHG = factor(GHG, levels = c("CO2", "CH4", "N2O")),
    Site = factor(Site, levels = c("KL", "ZB")),
    Period = factor(Period, levels = period_levels)
  ) |>
  arrange(GHG, Site, Period) |>
  rename(`DRW Cycle` = Period) |>
  select(
    GHG, Site, `DRW Cycle`,
    Control, `Control + N`, DRW, `DRW + N`,
    `Natural Precipitation`, `Drying–rewetting`, `p-value`
  )

write_csv(
  manuscript_table,
  file.path(dir_reporting, manuscript_table_file)
)


# Reporting files ====

writeLines(
  c(
    "Values are observed means ± SE across plot-specific means.",
    "Each plot receives equal descriptive weight within its group.",
    "Daily missing observations are not gap-filled or replaced with zeros.",
    "Recorded zero fluxes are retained.",
    "Positive CH4 values represent uptake.",
    "",
    "Separate models were fitted for every gas, site and period.",
    "The four-treatment model uses Treatment + DateGroup as fixed effects.",
    "The pooled model uses Water + DateGroup as fixed effects.",
    "Nitrogen is omitted from the separately fitted pooled model.",
    "Both models include a plot random intercept and within-plot CAR1.",
    "Residual SD differs by Treatment in the four-treatment model",
    "and by Water in the pooled model, using varIdent.",
    "Sampling date is a fixed factor, not a random effect.",
    "Models are fitted using REML.",
    "Contrasts use containment degrees of freedom.",
    "",
    "Letters use Fisher-adjusted pairwise comparisons within each row.",
    "Shared letters do not establish equivalence.",
    "The pooled p-value tests DRW minus NP in the separate pooled model.",
    paste("The manuscript table displays:", table_p_column),
    "Both raw and Holm-adjusted pooled p-values are exported.",
    "Holm adjustment is across period tests, including All cycles,",
    "separately within each gas/site.",
    "Pooled confidence intervals are pointwise.",
    "",
    "All cycles excludes pre-treatment and includes the Final Period.",
    "The All cycles models are fitted directly to all retained",
    "treatment-period observations; they do not average cycle estimates.",
    "",
    "The pooled analysis reflects the included nitrogen-group composition.",
    "It is not a nitrogen-adjusted water effect.",
    "",
    "VerySmallPlotSD is a screening flag based on PlotSD / sigma < 0.001.",
    "It is not a formal singularity test.",
    "Inspect residual distributions, temporal correlations and model",
    "stability before treating the inferential results as final."
  ),
  file.path(dir_reporting, notes_file)
)

write_session_info(dir_reporting, session_info_file)


# Console summary ====

cat("\nMODEL STATUS\n")
print(model_status, n = Inf, width = Inf)

cat("\nMANUSCRIPT TABLE — SUBJECT TO DIAGNOSTIC REVIEW\n")
print(manuscript_table, n = Inf, width = Inf)

cat("\nSEPARATE POOLED WATER-MODEL COMPARISONS\n")
print(water_results, n = Inf, width = Inf)

cat("\nOutputs saved under:\n", dir_reporting, "\n")


# Targeted sensitivity analysis ====


target_key <- "CO2_KL__DRW Cycle 3"

target_fits <- list(
  FourTreatments = four_fits[[target_key]],
  PooledWater = water_fits[[target_key]]
)

stopifnot(all(vapply(
  target_fits,
  function(x) inherits(x$model, "lme"),
  logical(1)
)))

# Full calendar-lag correlations
residual_lag_results |>
  filter(
    ID == "CO2_KL",
    Period == "DRW Cycle 3"
  ) |>
  print(n = Inf, width = Inf)

# Residuals and fitted values from the existing models
target_diagnostics <- imap_dfr(
  target_fits,
  function(fit, model_name) {
    fit$data |>
      mutate(
        Model = model_name,
        Fitted = as.numeric(fitted(fit$model, level = 1)),
        Residual = as.numeric(residuals(
          fit$model,
          type = "normalized"
        ))
      )
  }
)

# Largest normalized residuals, with observed fluxes and dates
target_diagnostics |>
  group_by(Model) |>
  slice_max(
    order_by = abs(Residual),
    n = 10,
    with_ties = FALSE
  ) |>
  ungroup() |>
  select(
    Model, Plot, Treatment, SampleDate,
    Flux, Fitted, Residual
  ) |>
  print(n = Inf, width = Inf)

# Rewetting-day sensitivity ====

# Existing pooled model: no refitting here.
base_fit <- water_fits[["CO2_KL__DRW Cycle 3"]]
base_model <- base_fit$model

event_data <- base_fit$data |>
  mutate(
    # Additional DRW response on the confirmed rewetting date.
    PulseDRW = as.numeric(
      Water == "DRW" &
        SampleDate == event_date
    )
  ) |>
  arrange(Plot, Day)

stopifnot(
  any(event_data$PulseDRW == 1),
  any(
    event_data$Water == "NP" &
      event_data$SampleDate == event_date
  )
)

event_formula <- Flux ~ Water + DateGroup + PulseDRW

X <- model.matrix(event_formula, data = event_data)

if (qr(X)$rank < ncol(X)) {
  stop("Event model is rank deficient.")
}

# Keep the covariance structure and use its estimates as starts.
event_model <- lme(
  fixed = event_formula,
  random = ~ 1 | Plot,
  correlation = corCAR1(
    value = unname(coef(
      base_model$modelStruct$corStruct,
      unconstrained = FALSE
    )),
    form = ~ Day | Plot
  ),
  weights = varIdent(
    value = coef(
      base_model$modelStruct$varStruct,
      unconstrained = FALSE
    ),
    form = ~ 1 | Water
  ),
  data = event_data,
  method = "REML",
  na.action = na.fail,
  control = lmeControl(
    opt = "nlminb",
    maxIter = 50,
    msMaxIter = 100,
    msMaxEval = 200,
    niterEM = 25,
    returnObject = FALSE
  )
)

# Save immediately.
safe_save_rds(
  list(
    model = event_model,
    data = event_data,
    event_date = event_date
  ),
  file.path(dir_interim, event_model_file)
)


# Pooled treatment contrasts ----


base_emm <- emmeans(
  base_model,
  ~Water,
  weights = "equal",
  mode = "containment"
)

base_comparison <- summary(
  pairs(base_emm, reverse = TRUE, adjust = "none"),
  infer = c(TRUE, TRUE)
) |>
  as_tibble() |>
  mutate(
    Analysis = "Original: constant cycle contrast",
    .before = 1
  )

# Construct the event-aware contrasts explicitly.
# For NP, PulseDRW must always be zero.
event_grid <- emmeans(
  event_model,
  ~ Water * PulseDRW,
  at = list(PulseDRW = c(0, 1)),
  weights = "equal",
  mode = "containment"
)

grid_rows <- as.data.frame(event_grid)

# Equal weight for each sampled date in this cycle.
sampled_dates <- unique(event_data$SampleDate)

event_weight <- mean(
  sampled_dates == event_date
)

ordinary_day <- with(
  grid_rows,
  as.numeric(Water == "DRW" & PulseDRW == 0) -
    as.numeric(Water == "NP" & PulseDRW == 0)
)

rewetting_day <- with(
  grid_rows,
  as.numeric(Water == "DRW" & PulseDRW == 1) -
    as.numeric(Water == "NP" & PulseDRW == 0)
)

cycle_average <- (
  (1 - event_weight) * ordinary_day +
    event_weight * rewetting_day
)

event_comparisons <- summary(
  contrast(
    event_grid,
    method = list(
      "Cycle average: DRW - NP" = cycle_average,
      "Other days: DRW - NP" = ordinary_day,
      "20 October: DRW - NP" = rewetting_day,
      "Extra DRW response on 20 October" =
        rewetting_day - ordinary_day
    ),
    adjust = "none"
  ),
  infer = c(TRUE, TRUE),
  adjust = "none"
) |>
  as_tibble() |>
  mutate(
    Analysis = "Rewetting-day sensitivity",
    .before = 1
  )

comparison_table <- bind_rows(
  base_comparison,
  event_comparisons
)

print(comparison_table, n = Inf, width = Inf)

write_csv(
  comparison_table,
  file.path(dir_processed, event_comparisons_file)
)

# we had a natural autumn drought so pulse in DRW and DRw + N exceeded C
# this explains event sensitivity

# Residual lags ----

event_diagnostics <- event_data |>
  mutate(
    FittedFlux = as.numeric(fitted(event_model, level = 1)),
    NormalizedResidual = as.numeric(residuals(
      event_model,
      type = "normalized"
    ))
  )

event_lags <- calendar_lags(event_diagnostics)

print(event_lags, n = Inf)

write_csv(
  event_lags,
  file.path(dir_interim, event_lags_file)
)

write_csv(
  event_diagnostics,
  file.path(dir_interim, event_diagnostics_file)
)
# Read the original file explicitly as UTF-8.
table_clean <- read_csv(
  file.path(excel_output_directory, manuscript_table_file),
  locale = locale(encoding = "UTF-8"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

# Adds a UTF-8 marker so Excel recognises ± and the dash correctly.
write_excel_csv(
  table_clean,
  file.path(excel_output_directory, excel_table_file)
)

table_relabelled <- read_csv(
  file.path(excel_output_directory, manuscript_table_file),
  locale = locale(encoding = "UTF-8"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

for (i in seq_len(nrow(table_relabelled))) {
  cells <- unlist(
    table_relabelled[i, treatment_columns],
    use.names = FALSE
  )

  # Extract observed means and existing letter labels.
  means <- parse_number(cells)

  old_labels <- regmatches(
    cells,
    regexpr("[a-z]+$", cells)
  )

  if (
    anyNA(means) ||
      length(old_labels) != length(treatment_columns)
  ) {
    stop("Could not read means or letters in row ", i)
  }

  # Visit treatments from highest to lowest observed mean.
  descending_order <- order(means, decreasing = TRUE)

  encountered_letters <- unique(unlist(
    strsplit(old_labels[descending_order], split = "")
  ))

  # Rename existing letters without changing group membership.
  letter_map <- setNames(
    letters[seq_along(encountered_letters)],
    encountered_letters
  )

  new_labels <- vapply(
    strsplit(old_labels, split = ""),
    function(x) {
      paste(sort(unname(letter_map[x])), collapse = "")
    },
    character(1)
  )

  new_cells <- paste(
    sub("\\s+[a-z]+$", "", cells),
    new_labels
  )

  for (j in seq_along(treatment_columns)) {
    table_relabelled[[treatment_columns[j]]][i] <- new_cells[j]
  }
}

write_excel_csv(
  table_relabelled,
  file.path(
    excel_output_directory,
    excel_letters_table_file
  )
)

print(table_relabelled, n = Inf, width = Inf)
format_one_decimal <- function(x) {
  parts <- strsplit(x, "±", fixed = TRUE)

  vapply(parts, function(part) {
    if (length(part) != 2L) {
      stop("Could not identify mean and SE.")
    }

    mean_value <- parse_number(part[1])
    se_value <- parse_number(part[2])

    # Preserve significance letters, where present.
    label <- regmatches(
      part[2],
      regexpr("[a-z]+\\s*$", part[2])
    )

    suffix <- if (length(label)) {
      paste0(" ", trimws(label))
    } else {
      ""
    }

    paste0(
      sprintf("%.1f ± %.1f", mean_value, se_value),
      suffix
    )
  }, character(1))
}

table_relabelled <- table_relabelled |>
  mutate(
    across(
      all_of(value_columns),
      format_one_decimal
    )
  )

write_excel_csv(
  table_relabelled,
  file.path(
    excel_output_directory,
    excel_one_decimal_table_file
  )
)

print(table_relabelled, n = Inf, width = Inf)


# Leave-one-plot-out refits ====

s2_base <- water_fits[["N2O_KL__DRW Cycle 1"]]

if (is.null(s2_base) || !inherits(s2_base$model, "lme")) {
  stop("The N2O KL Cycle 1 pooled model is unavailable.")
}

s2_data <- s2_base$data

if (
  !is.data.frame(s2_data) ||
    !all(c("Flux", "Water", "DateGroup", "Plot", "Day") %in%
      names(s2_data))
) {
  stop("The saved model data are missing required columns.")
}

# Extract the DRW minus NP contrast using explicit model data.
s2_extract <- function(model, omitted_plot, model_data) {
  emm <- emmeans(
    model,
    ~Water,
    data = model_data,
    weights = "equal",
    mode = "containment"
  )

  water_labels <- as.character(as.data.frame(emm)$Water)

  stopifnot(setequal(water_labels, c("NP", "DRW")))

  coefficients <- unname(c(NP = -1, DRW = 1)[water_labels])

  result <- as.data.frame(summary(
    contrast(
      emm,
      method = list("DRW - NP" = coefficients),
      adjust = "none"
    ),
    infer = c(TRUE, TRUE),
    level = confidence_level,
    adjust = "none"
  ))

  tibble(
    OmittedPlot = omitted_plot,
    Estimate = result$estimate,
    SE = result$SE,
    DF = result$df,
    LowerCI = result$lower.CL,
    UpperCI = result$upper.CL
  )
}

# Re-estimate all model parameters after removing one plot.
s2_refit <- function(plot_id) {
  warning_messages <- character()

  result <- tryCatch(
    withCallingHandlers(
      {
        dat <- s2_data[
          as.character(s2_data$Plot) != plot_id, ,
          drop = FALSE
        ]

        dat <- droplevels(dat)
        dat <- dat[order(dat$Plot, dat$Day), ]

        if (length(unique(dat$Water)) != 2L) {
          stop("A water group is missing after omission.")
        }

        X <- model.matrix(~ Water + DateGroup, data = dat)

        if (qr(X)$rank < ncol(X)) {
          stop("Fixed effects are rank deficient after omission.")
        }

        model <- lme(
          fixed = Flux ~ Water + DateGroup,
          random = ~ 1 | Plot,
          correlation = corCAR1(
            value = correlation_start_value,
            form = ~ Day | Plot
          ),
          weights = varIdent(
            form = ~ 1 | Water
          ),
          data = dat,
          method = "REML",
          na.action = na.fail,
          control = model_control
        )

        s2_extract(
          model = model,
          omitted_plot = plot_id,
          model_data = dat
        )
      },
      warning = function(w) {
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      tibble(
        OmittedPlot = plot_id,
        Estimate = NA_real_,
        SE = NA_real_,
        DF = NA_real_,
        LowerCI = NA_real_,
        UpperCI = NA_real_,
        Error = conditionMessage(e)
      )
    }
  )

  if (!"Error" %in% names(result)) {
    result$Error <- NA_character_
  }

  result$Warnings <- paste(
    unique(warning_messages),
    collapse = " | "
  )

  result
}

# Full-data reference: use the existing final model.
s2_full <- s2_extract(
  model = s2_base$model,
  omitted_plot = "None - full data",
  model_data = s2_data
)

s2_full$Error <- NA_character_
s2_full$Warnings <- paste(s2_base$warnings, collapse = " | ")

# Leave out each retained plot in turn.
s2_plots <- sort(unique(as.character(s2_data$Plot)))

s2_results <- bind_rows(
  s2_full,
  map_dfr(s2_plots, \(plot_id) {
    logger::log_info(paste("Table S2: omitting plot", plot_id))
    s2_refit(plot_id)
  })
)

# Export full-precision estimates and diagnostic messages.
s2_output_file <- file.path(
  dir_reporting,
  lopo_file
)

write_csv(s2_results, s2_output_file)

# Display rounded estimates.
s2_display <- s2_results |>
  mutate(
    across(
      c(Estimate, SE, DF, LowerCI, UpperCI),
      ~ round(.x, 2)
    )
  )

print(s2_display, n = Inf, width = Inf)

# Check whether the caption's confidence-interval claim holds.
s2_complete <- all(
  is.na(s2_results$Error) &
    is.finite(s2_results$Estimate) &
    is.finite(s2_results$LowerCI) &
    is.finite(s2_results$UpperCI)
)

s2_has_warnings <- any(
  !is.na(s2_results$Warnings) &
    nzchar(s2_results$Warnings)
)

if (!s2_complete) {
  logger::log_warn(
    "Some fits or contrasts failed. Review Error before using Table S2."
  )
} else if (all(s2_results$UpperCI < 0)) {
  logger::log_info("All reported 95% confidence intervals are below zero.")
} else {
  logger::log_warn(
    "At least one 95% confidence interval reaches or exceeds zero; ",
    "revise the caption."
  )
}

if (s2_has_warnings) {
  logger::log_warn("Model warnings were recorded; review the Warnings column.")
}

logger::log_info(paste("Results saved to:", s2_output_file))
log_script_end("table-02-ghg-treatment-effects")
