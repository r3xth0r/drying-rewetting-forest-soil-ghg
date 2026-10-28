# Project setup ----
# One-time setup: creates data directories and installs CmdStan for brms.
# Usage: source("dev/setup.R") from the project root.

suppressPackageStartupMessages({
  library("here")
})

source(here::here("dev", "logging.R"))
log_script_start("setup")

# Directories ====

for (d in c("dat/interim", "dat/processed", "dat/reporting", "plt")) {
  dir.create(here(d), recursive = TRUE, showWarnings = FALSE)
}

# CmdStan ====

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages(
    "cmdstanr",
    repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
  )
}

cmdstanr::check_cmdstan_toolchain(fix = TRUE)

required_cmdstan_version <- "2.40.0"
cmdstan_install_dir <- file.path(path.expand("~"), ".cmdstan")
cmdstan_path <- file.path(
  cmdstan_install_dir,
  paste0("cmdstan-", required_cmdstan_version)
)

if (!dir.exists(cmdstan_path)) {
  cmdstanr::install_cmdstan(
    version = required_cmdstan_version,
    dir = cmdstan_install_dir,
    cores = max(1L, parallel::detectCores() - 1L)
  )
}

cmdstanr::set_cmdstan_path(cmdstan_path)
installed_cmdstan_version <- cmdstanr::cmdstan_version()
if (installed_cmdstan_version != required_cmdstan_version) {
  stop(
    "Expected CmdStan ", required_cmdstan_version,
    "; found ", installed_cmdstan_version, "."
  )
}

logger::log_info(paste("CmdStan version:", installed_cmdstan_version))
log_script_end("setup")
