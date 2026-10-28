# Analysis

This directory contains all R scripts to reproduce the results presented in the manuscript.

Run analysis scripts from the project root (open `exafor.Rproj`). They attach their packages, source shared configuration/helpers, and write to the shared output directories:

| Output                                                 | Directory        |
| ------------------------------------------------------ | ---------------- |
| figures                                                | `plt/`           |
| interim files (prepared data, model fits, diagnostics) | `dat/interim/`   |
| final model results and derived datasets               | `dat/processed/` |
| tables and statistics reported in the manuscript       | `dat/reporting/` |

Run `make all` to clean previous outputs and reproduce the full pipeline, including session information. Use `make figures` or `make tables` to run either group without cleaning first. Individual scripts can also be run with `Rscript dev/<script>.R`.

The CO2 pipeline has separate `make co2-models` and `make co2-model-outputs` stages; `make co2-outputs` runs both. The N2O and Table 2 scripts combine model fitting with their dependent outputs and are available as `make n2o-analysis` and `make table-02-analysis`.

Figure exports are PNG-only by default. Set `save_pdf_figures <- TRUE` in `config.R` to also create PDFs.

## Shared files

- **[config.R](config.R)**: paths, study design, soil physics, colours and the ggplot theme (`theme_exafor()`)
- **[helpers.R](helpers.R)**: input parsing, summary statistics and output helpers shared by the scripts
- **[logging.R](logging.R)**: concise start, completion, and artifact INFO logs for runnable scripts
- **[soil-parameter-helpers.R](soil-parameter-helpers.R)**: shared validation, statistical-label, and diagnostic helpers for soil analyses
- **[co2-model-helpers.R](co2-model-helpers.R)**: group-selection helpers shared by CO2 modeling and output scripts
- **[setup.R](setup.R)**: one-time setup (creates output directories, installs CmdStan for the Bayesian models)
- **[sessioninfo.qmd](sessioninfo.qmd)**: renders R/platform details and versions for all packages in `renv.lock` to `dat/reporting/session-info.html`.

## Scripts

- **[figure-02-wfps.R](figure-02-wfps.R)**
  - creates Figure 2: panelled WFPS for NP and RE treatments at KL, with precipitation and irrigation events
  - included in `make figures`
  - requires datasets `soil-moisture-klausen.csv`, `precipitation-klausen.csv`, and `irrigation-klausen.csv`

- **[figure-03-gas-fluxes.R](figure-03-gas-fluxes.R)**
  - creates Figure 3: panelled gas fluxes, soil moisture, soil temperature, and precipitation for both sites
  - included in `make figures`
  - requires datasets `co2-klausen.csv`, `ch4-klausen.csv`, `n2o-klausen.csv`, `co2-zoe.csv`, `ch4-zoe.csv`, `soil-moisture-klausen.csv`, `soil-moisture-zb.csv`, `soil-temperature-klausen.csv`, `soil-temperature-zb.csv`, `precipitation-klausen.csv`, and `precipitation-zb.csv`

- **[co2-models.R](co2-models.R)**
  - prepares data and fits temperature and temperature × moisture CO2 models for both sites
  - writes fitted models and typed result/prediction artifacts to `dat/interim/`
  - run with `make co2-models` or as part of `make all`
  - requires datasets `co2-klausen.csv`, `co2-zoe.csv`, `soil-moisture-klausen.csv`, `soil-moisture-zb.csv`, `soil-temperature-klausen.csv`, and `soil-temperature-zb.csv`

- **[co2-model-outputs.R](co2-model-outputs.R)**
  - creates Figures 4, S3, S4, and S6, plus Tables 3 and S5 from saved CO2 model artifacts
  - run `make co2-model-outputs` to render existing artifacts or `make co2-outputs` to fit and render
  - requires artifacts from `co2-models.R` in `dat/interim/`

- **[n2o-bayesian-models.R](n2o-bayesian-models.R)**
  - fits the full, reduced, and depth-specific N2O models; writes Tables 5, 6, and S8, numeric diagnostics, and plot inputs to `dat/interim/n2o/`
  - requires CmdStan; run [setup.R](setup.R) once first
  - run with `make n2o-models` or as part of `make n2o-analysis`
  - requires datasets `n2o-klausen.csv`, `soil-moisture-klausen.csv`, `soil-temperature-klausen.csv`, `precipitation-klausen.csv`, and `irrigation-klausen.csv`

- **[n2o-bayesian-outputs.R](n2o-bayesian-outputs.R)**
  - creates Figures 5 and S8 and the trace, predictive, Pareto-k, and GAM check plots from saved N2O model artifacts
  - run `make n2o-model-outputs` to re-render without refitting, or `make n2o-analysis` to fit and render
  - included in `make figures` and `make tables`
  - requires artifacts from `n2o-bayesian-models.R` in `dat/interim/n2o/`

- **[spei-seasonal-analysis.R](spei-seasonal-analysis.R)**
  - creates the seasonal SPEI comparison table for KL and ZBSM
  - included in `make tables`
  - requires `winfore2-spei30-<YYYY>.nc` for each year from 2016 through 2024

- **[download-spei-data.sh](download-spei-data.sh)**
  - downloads annual SPEI NetCDF files for 2016-2025 into `dat/raw/`, skipping files that already exist
  - run from the project root with `bash dev/download-spei-data.sh`

- **[figure-06-soil-parameters.R](figure-06-soil-parameters.R)**
  - creates Figure 6 from soil-depth samples, with significance annotations and diagnostics
  - included in `make figures`
  - requires `soil-depth-sampling.csv` in `dat/raw`

- **[figure-s09-soil-parameters.R](figure-s09-soil-parameters.R)**
  - creates Figure S9 and statistical summaries from survey soil parameters using Welch ANOVA and Games–Howell comparisons
  - included in `make figures`
  - requires `soil-parameter-survey.csv` in `dat/raw`

- **[table-02-ghg-treatment-effects.R](table-02-ghg-treatment-effects.R)**
  - runs LMMs to test GHG differences between treatments, allowing for autocorrelation and differing variances per treatment
  - creates Tables 2 and S2 (LMM leave-one-plot sensitivity)
  - included in `make tables`; can also be run with `make table-02-analysis`
  - requires datasets `co2-klausen.csv`, `co2-zoe.csv`, `ch4-klausen.csv`, `ch4-zoe.csv`, `n2o-klausen.csv`

- **[table-04-ch4-correlations.R](table-04-ch4-correlations.R)**
  - runs repeated measures correlations between CH4 fluxes and soil moisture and temperature in KL and ZB
  - creates Table 4, including results, data-used, and settings CSVs
  - included in `make tables`
  - requires datasets `ch4-klausen.csv`, `ch4-zoe.csv`, `soil-moisture-klausen.csv`, `soil-moisture-zb.csv`, `soil-temperature-klausen.csv`, and `soil-temperature-zb.csv`

- **[table-s07-co2-equivalents.R](table-s07-co2-equivalents.R)**
  - calculates CO2 equivalents for CH4 and N2O and exports cumulative/statistical tables as CSV files
  - creates Table S7
  - included in `make tables`
  - requires datasets `co2-klausen.csv`, `ch4-klausen.csv`, `n2o-klausen.csv`, `co2-zoe.csv`, and `ch4-zoe.csv`
