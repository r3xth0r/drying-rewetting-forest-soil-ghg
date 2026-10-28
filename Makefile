RSCRIPT ?= Rscript
QUARTO ?= quarto

LOG ?= make.log

.PHONY: all clean figures tables all-run figures-run tables-run \
	figure-02 figure-03 figure-06 figure-s09 \
	co2-models co2-model-outputs co2-outputs n2o-models n2o-model-outputs n2o-analysis table-02-analysis \
	spei-table table-04 table-s07 session-info

all:
	$(MAKE) all-run > $(LOG) 2>&1

all-run:
	$(MAKE) clean
	$(MAKE) figures-run tables-run session-info

clean:
	find plt dat/interim dat/processed dat/reporting -mindepth 1 -maxdepth 1 ! -name .gitkeep -exec rm -rf -- {} +

figures:
	$(MAKE) figures-run > $(LOG) 2>&1

tables:
	$(MAKE) tables-run > $(LOG) 2>&1

figures-run: figure-02 figure-03 co2-outputs n2o-analysis figure-06 figure-s09

tables-run: co2-outputs n2o-analysis table-02-analysis spei-table table-04 table-s07

figure-02:
	$(RSCRIPT) dev/figure-02-wfps.R

figure-03:
	$(RSCRIPT) dev/figure-03-gas-fluxes.R

co2-models:
	$(RSCRIPT) dev/co2-models.R

co2-model-outputs:
	$(RSCRIPT) dev/co2-model-outputs.R

co2-outputs: co2-models
	$(MAKE) co2-model-outputs

n2o-models:
	$(RSCRIPT) dev/n2o-bayesian-models.R

n2o-model-outputs:
	$(RSCRIPT) dev/n2o-bayesian-outputs.R

n2o-analysis: n2o-models
	$(MAKE) n2o-model-outputs

figure-06:
	$(RSCRIPT) dev/figure-06-soil-parameters.R

figure-s09:
	$(RSCRIPT) dev/figure-s09-soil-parameters.R

table-02-analysis:
	$(RSCRIPT) dev/table-02-ghg-treatment-effects.R

spei-table:
	$(RSCRIPT) dev/spei-seasonal-analysis.R

table-04:
	$(RSCRIPT) dev/table-04-ch4-correlations.R

table-s07:
	$(RSCRIPT) dev/table-s07-co2-equivalents.R

session-info:
	$(QUARTO) render dev/sessioninfo.qmd
	mv dev/sessioninfo.html dat/reporting/session-info.html
