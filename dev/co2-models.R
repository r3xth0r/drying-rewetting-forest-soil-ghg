# CO2 model fitting: temperature and temperature x moisture ----
# Inputs: co2-{klausen,zoe}, soil-moisture-{klausen,zb}, soil-temperature-*.
# Outputs: fitted models, model data, results, and predictions.

suppressPackageStartupMessages({
  library("tidyverse")
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("co2-models")

source(here::here("dev", "config.R"))
source(here::here("dev", "helpers.R"))
source(here::here("dev", "co2-model-helpers.R"))

# Settings ====

kl_excluded_plots <- c("KL2", "KL3", "KL10")
kl_end_exclusive <- as.Date("2021-11-22")
zb_start_exclusive <- as.Date("2022-04-26")

# Differs from depth_labels (spaces around the dash); used in model outputs
model_depth_labels <- c(
  `5` = "0 – 10 cm",
  `15` = "10 – 20 cm",
  `30` = "20 – 40 cm"
)

zb_moisture_map <- tibble(
  Sensor = c("C1", "S1", "S2", "S4"),
  Depth_cm = 5,
  Treatment = c("NP", "RE", "RE", "RE")
)

# Start values of the temperature model
temperature_starts <- list(
  list(a = 6.35, b = 0.38, c = -0.005),
  list(a = 15, b = 0.25, c = -0.006),
  list(a = 10, b = 0.30, c = -0.010)
)

# Multipliers of the temperature-model b and c for the moisture-model starts
moisture_start_factors <- c(0.8, 1, 1.2)

site_depths <- list(KL = c(5, 15, 30), ZB = 5)
model_labels <- c("Temperature Function", "Temperature Moisture Function")

# Model artifact names
file_model_data <- "co2-model-data.rds"
file_models <- "co2-models.rds"
file_results <- "co2-model-results.csv"
file_predictions <- "co2-model-predictions.csv"

# Data ====

## Preparation ----

prepare_flux <- function(filename, site, prefix, excluded = character()) {
  plot_columns <- str_c(prefix, 1:16)
  dat <- read_input(filename)
  require_columns(dat, plot_columns, filename)

  dat <- dat |>
    add_time("Date") |>
    select(SampleDate, Time, TimeCat, all_of(plot_columns)) |>
    pivot_longer(
      all_of(plot_columns),
      names_to = "Plot",
      values_to = "Flux"
    ) |>
    mutate(
      Plot_number = match(Plot, plot_columns),
      Treatment4 = case_when(
        Plot_number %in% c(1, 5, 10, 13) ~ "C",
        Plot_number %in% c(2, 6, 9, 14) ~ "C + N",
        Plot_number %in% c(3, 7, 11, 15) ~ "RE",
        Plot_number %in% c(4, 8, 12, 16) ~ "RE + N"
      ),
      Treatment = if_else(Treatment4 %in% c("C", "C + N"), "NP", "RE")
    ) |>
    filter(!Plot %in% excluded, !is.na(TimeCat), is.finite(Flux))

  # KL: plot values -> treatment and timestamp -> treatment and time bin
  if (site == "KL") {
    dat <- dat |>
      group_by(Treatment4, Treatment, SampleDate, Time, TimeCat) |>
      summarise(Flux = mean(Flux), .groups = "drop")
  }

  # Then average the nitrogen treatments within NP and RE
  dat |>
    group_by(Treatment4, Treatment, SampleDate, TimeCat) |>
    summarise(Flux = mean(Flux), .groups = "drop") |>
    group_by(Treatment, SampleDate, TimeCat) |>
    summarise(Flux = mean(Flux), .groups = "drop")
}

## Model data ----

sensor_keys <- c("SampleDate", "TimeCat", "Depth_cm", "Treatment")
flux_keys <- c("SampleDate", "TimeCat", "Treatment")

kl_flux <- prepare_flux(
  "co2-klausen.csv",
  site = "KL",
  prefix = "KL",
  excluded = kl_excluded_plots
) |>
  filter(SampleDate < kl_end_exclusive)

kl_moisture <- prepare_sensors(
  "soil-moisture-klausen.csv", "X", sensor_map_kl_moisture, "Moisture_pc"
)

kl_temperature <- prepare_sensors(
  "soil-temperature-klausen.csv", "X", sensor_map_kl_temperature, "Tsoil",
  equal_sensor_weight = TRUE
)

kl_data <- kl_moisture |>
  inner_join(kl_temperature, by = sensor_keys) |>
  inner_join(kl_flux, by = flux_keys) |>
  mutate(
    Site = "KL",
    Porosity_pc = unname(porosity$KL[as.character(Depth_cm)]),
    WFPS = 100 * Moisture_pc / Porosity_pc
  )

zb_flux <- prepare_flux("co2-zoe.csv", site = "ZB", prefix = "ZOE")

zb_moisture <- prepare_sensors(
  "soil-moisture-zb.csv", "Date", zb_moisture_map, "Moisture_pc"
) |>
  filter(SampleDate > zb_start_exclusive)

zb_temperature <- prepare_sensors(
  "soil-temperature-zb.csv", "SampleDate", sensor_map_zb_temperature, "Tsoil"
)

zb_data <- zb_moisture |>
  inner_join(zb_temperature, by = sensor_keys) |>
  inner_join(zb_flux, by = flux_keys) |>
  mutate(
    Site = "ZB",
    Porosity_pc = unname(porosity$ZB["5"]),
    WFPS = 100 * Moisture_pc / Porosity_pc
  )

model_data <- bind_rows(kl_data, zb_data) |>
  mutate(Depth = unname(model_depth_labels[as.character(Depth_cm)])) |>
  filter(is.finite(Flux), is.finite(Tsoil), is.finite(WFPS)) |>
  select(
    Site, Depth_cm, Depth, Treatment, SampleDate, TimeCat,
    Flux, Tsoil, WFPS, Moisture_pc, Porosity_pc
  ) |>
  arrange(Site, Depth_cm, Treatment, SampleDate, TimeCat)

ensure_dir(dir_interim)
write_rds(model_data, file.path(dir_interim, file_model_data))

# Models ====

## Fitting ----

best_fit <- function(candidates) {
  candidates <- keep(
    candidates,
    \(x) !is.null(x) && isTRUE(x$convInfo$isConv)
  )

  if (!length(candidates)) {
    stop("No starting values produced a converged NLS fit.")
  }

  rss <- map_dbl(candidates, \(x) sum(residuals(x)^2))
  candidates[[which.min(rss)]]
}

fit_temperature <- function(dat) {
  candidates <- map(temperature_starts, \(start) {
    tryCatch(
      nls(
        Flux ~ a * exp(b * Tsoil + c * Tsoil^2),
        data = dat,
        start = start,
        algorithm = "port",
        lower = c(a = 1e-10, b = -Inf, c = -Inf),
        control = nls.control(maxiter = 1000)
      ),
      error = \(err) NULL
    )
  })

  best_fit(candidates)
}

# The temperature-model a is held fixed; d and e start from a linear fit
fit_temperature_moisture <- function(dat, temp_model) {
  temp_coef <- coef(temp_model)
  fixed_a <- unname(temp_coef["a"])

  starts <- expand.grid(
    b = unname(temp_coef["b"]) * moisture_start_factors,
    c = unname(temp_coef["c"]) * moisture_start_factors
  )

  dat$a_fixed <- fixed_a

  candidates <- map(seq_len(nrow(starts)), \(i) {
    tryCatch(
      {
        b_start <- starts$b[i]
        c_start <- starts$c[i]
        temp_component <- fixed_a *
          exp(b_start * dat$Tsoil + c_start * dat$Tsoil^2)

        de_model <- lm(
          Flux ~ 0 + linear + quadratic,
          data = data.frame(
            Flux = dat$Flux,
            linear = temp_component * dat$WFPS,
            quadratic = temp_component * dat$WFPS^2
          )
        )
        de_start <- unname(coef(de_model))

        if (any(!is.finite(de_start))) {
          return(NULL)
        }

        nls(
          Flux ~ a_fixed * exp(b * Tsoil + c * Tsoil^2) *
            (d * WFPS + e * WFPS^2),
          data = dat,
          start = list(
            b = b_start, c = c_start, d = de_start[1], e = de_start[2]
          ),
          algorithm = "port",
          control = nls.control(maxiter = 1000)
        )
      },
      error = \(err) NULL
    )
  })

  best_fit(candidates)
}

## Results ----

# Six significant figures, no scientific notation (small moisture terms)
coefficient_text <- function(x) {
  format(signif(unname(x), 6), scientific = FALSE, trim = TRUE)
}

signed_term <- function(value, variable) {
  str_c(
    if (value < 0) " - " else " + ",
    coefficient_text(abs(value)),
    " * ",
    variable
  )
}

make_equation <- function(p, moisture = FALSE) {
  equation <- str_c(
    coefficient_text(p["a"]), " * exp(",
    coefficient_text(p["b"]), " * Tsoil",
    signed_term(p["c"], "Tsoil²"), ")"
  )

  if (moisture) {
    equation <- str_c(
      equation, " * (",
      coefficient_text(p["d"]), " * WFPS",
      signed_term(p["e"], "WFPS²"), ")"
    )
  }

  equation
}

result_row <- function(model, dat, label, fixed_a = NULL) {
  p <- set_names(rep(NA_real_, 5), c("a", "b", "c", "d", "e"))
  estimates <- coef(model)
  p[names(estimates)] <- estimates

  if (!is.null(fixed_a)) {
    p["a"] <- fixed_a
  }

  observed <- dat$Flux
  predicted <- as.numeric(fitted(model))
  errors <- observed - predicted
  sse <- sum(errors^2)
  sst <- sum((observed - mean(observed))^2)

  # Kling-Gupta efficiency (Gupta et al., 2009)
  kge <- NA_real_

  if (sd(observed) > 0 && sd(predicted) > 0 && mean(observed) != 0) {
    r <- cor(observed, predicted)
    alpha <- sd(predicted) / sd(observed)
    beta <- mean(predicted) / mean(observed)
    kge <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
  }

  data.frame(
    Model = label,
    Site = dat$Site[1],
    Depth_cm = dat$Depth_cm[1],
    Depth = dat$Depth[1],
    Treatment = dat$Treatment[1],
    Equation = make_equation(p, moisture = !is.null(fixed_a)),
    N = nrow(dat),
    R2 = if (sst > 0) 1 - sse / sst else NA_real_,
    RMSE = sqrt(mean(errors^2)),
    AIC = AIC(model),
    BIC = BIC(model),
    KGE = kge,
    as.list(p),
    a_status = if (is.null(fixed_a)) "Estimated" else "Fixed",
    row.names = NULL,
    check.names = FALSE
  )
}

## Fit all groups ----

fit_group <- function(site, depth_cm, treatment) {
  key <- group_key(site, depth_cm, treatment)
  dat <- select_group(model_data, site, depth_cm, treatment)

  if (
    nrow(dat) <= 5 || n_distinct(dat$Tsoil) < 3 || n_distinct(dat$WFPS) < 3
  ) {
    stop(key, ": insufficient complete observations or predictor variation.")
  }

  logger::log_info(paste("Fitting model group:", key))

  fits <- tryCatch(
    {
      temperature <- fit_temperature(dat)
      list(
        Temperature = temperature,
        Temperature_moisture = fit_temperature_moisture(dat, temperature)
      )
    },
    error = \(err) stop(key, ": ", conditionMessage(err), call. = FALSE)
  )

  fixed_a <- unname(coef(fits$Temperature)["a"])

  results <- bind_rows(
    result_row(fits$Temperature, dat, model_labels[1]),
    result_row(
      fits$Temperature_moisture, dat, model_labels[2],
      fixed_a = fixed_a
    )
  )

  predictions <- dat |>
    mutate(
      Temperature_prediction = as.numeric(fitted(fits$Temperature)),
      Temperature_moisture_prediction =
        as.numeric(fitted(fits$Temperature_moisture)),
      Temperature_residual = Flux - Temperature_prediction,
      Temperature_moisture_residual =
        Flux - Temperature_moisture_prediction
    )

  list(models = fits, results = results, predictions = predictions)
}

groups <- imap(site_depths, \(depths, site) {
  expand_grid(site = site, depth_cm = depths, treatment = c("NP", "RE"))
}) |>
  list_rbind()

group_fits <- pmap(groups, fit_group) |>
  set_names(pmap_chr(groups, group_key))

models <- map(group_fits, "models")
results <- group_fits |>
  map("results") |>
  list_rbind()
predictions <- group_fits |>
  map("predictions") |>
  list_rbind()

write_rds(models, file.path(dir_interim, file_models))

ensure_dir(dir_processed)
readr::write_csv(results, file.path(dir_processed, file_results))
readr::write_csv(predictions, file.path(dir_processed, file_predictions))
log_script_end("co2-models")
