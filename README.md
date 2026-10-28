# Analysis of Greenhouse Gas Dynamics in Austrian Forest Soils Under Drying-Rewetting Cycles

<p align="center">
    <a href="https://style.tidyverse.org">
        <img alt="Code style: tidyverse" src="https://img.shields.io/badge/codestyle-tidyverse-1a162d?style=flat-square&logo=r&logoColor=white"></a>
</p>

This repository supplements the manuscript by
Dylan Goff<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0009-0009-2821-1620)</sup>,
Thomas Dirnböck<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-8294-0690)</sup>,
Ika Djukic<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-5144-9321)</sup>,
Markus Gorfer<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-3010-6044)</sup>,
Barbara Kitzler<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-1617-9186)</sup>,
Johannes Kobler<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-0052-4245)</sup>,
Matthias Schlögl<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-4357-523X)</sup>, and
Eugenio Diaz-Pinés<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-9935-106X)</sup>
(2026):
**Drying-Rewetting cycles reduce soil CO<sub>2</sub> and N<sub>2</sub>O emissions while increasing CH<sub>4</sub> uptake in two Austrian broadleaved forests**.
<!-- *Journal*, Issue, Pages. [doi:](https://doi.org/). -->

## Repo structure

The repository follows a very simple cookicutter structure for data science project. The top-level directories `dat`, `dev` and `plt` contain data, code and plots, respectively.

The data directory is further subdivided into:

```sh
dat
├── raw         # the raw, immutable data dump
├── interim     # temporary and interim data files
├── processed   # data files of modelling results
└── reporting   # final results used for reporting in the paper
```

## Reproduce the analysis

From the repository root, with R, `make`, Quarto, and a C++ toolchain available, restore the R packages and set up CmdStan:

```sh
Rscript -e 'renv::restore()'
Rscript -e 'source("dev/setup.R")'
```

Then run the full analysis:

```sh
make all
```

`make all` clears existing generated files in `plt/`, `dat/interim/`, `dat/processed/`, and `dat/reporting/` (preserving `.gitkeep` files), then regenerates the figures, tables, and session information. Use `make figures` or `make tables` to run only that group without first cleaning outputs. Individual targets are listed in [dev/README.md](dev/README.md).

> [!NOTE]
> Running the whole pipeline takes some time, up to 30 minutes depending on the configuration.
