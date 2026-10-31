# Deal or No Deal: applied sampling submission

A region-stratified survey of 10,000 **simulated** contestants, with SRSWOR as a benchmark. The primary outcome is mean winnings in USD; the secondary outcome is the proportion below the regional board's opening expected value.

## Files

All project files are in the `R Code/` folder of this repository.

- `report.pdf`: report, questionnaire and AI declaration.
- `report.Rmd`: editable report source.
- `data/population.csv`: complete finite population (10,000 rows).
- `data/sample.csv`: all 2,670 stratified selections, including nonrespondents.
- `data/sample_srswor.csv`: all 2,670 benchmark selections.
- R scripts, this README, and generated results in `output/` and `figures/`.

## Reproduce

Open `R Code/R Code.Rproj` in RStudio, or open a terminal in the `R Code/` folder. Required packages are `survey`, `ggplot2`, `scales`, `knitr`, and `rmarkdown`. Install missing packages explicitly, for example:

```r
install.packages(c("survey", "ggplot2", "scales", "knitr", "rmarkdown"))
```

Run the complete pipeline:

```sh
Rscript run_all.R
```

PDF rendering requires Pandoc and XeLaTeX (a TeX distribution). If these are unavailable, reproduce the analysis without PDF rendering:

```sh
Rscript run_all.R --no-render
```

To update the report after editing the declaration:

```r
rmarkdown::render("report.Rmd", output_format = "pdf_document")
```

Each script reads saved inputs rather than depending on objects from preceding scripts. `run_all.R` runs them in separate environments. Existing generated CSVs and figures in this submission folder are overwritten. Seeds and software versions are documented in the report and `output/sessionInfo.txt`. Seeds support reproduction with the recorded software; exact cross-version equality is not guaranteed.

## What the code does

1. `01_simulate_population.R` creates the stylised games and saves the population (seed 20261019).
2. `02_sample_size.R` plans sample sizes and non-response allowance, saves allocation tables and sensitivity plots. For tight margins, strata are capped at their population counts and the remaining allocation is recalculated. Integer stratum counts are rounded up conservatively.
3. `03_draw_sample.R` samples independently without replacement within each region and simulates independent 90% response (seed 20261020).
4. `04_estimation.R` uses `survey` for estimates, variance, SE, normal 95% intervals, CV and design effect. It creates an equal-selected-size SRSWOR benchmark (seed 20261021).
5. `05_simulation.R` repeats selection and response 1,000 times per design (seed 20261022), verifies manual formulas against `survey`, and saves replicate and aggregate results.
6. `06_extra_figures.R` produces supplementary figures; only the precision histogram is retained in the shortened report.
7. `validate.R` checks population and sample integrity, weights and replicate completeness.

The respondent-level finite population correction is justified **conditional on response counts**, because response is independent with equal probability within each stratum. This justification does not extend automatically to informative non-response.

## Data dictionary

Population variables:

| Field | Meaning |
|---|---|
| `id` | Unique simulated contestant ID |
| `region` | UK, US or AU; sampling stratum |
| `age`, `gender` | Simulated frame demographics |
| `risk_aversion`, `wealth_usd` | Latent simulation inputs, rounded for storage; questionnaire proxies do not measure them exactly |
| `deal` | 1 if a banker offer was accepted; 0 otherwise |
| `round_ended` | Final round reached |
| `winnings_usd` | Actual simulated payout, in USD at fixed modelling conversion factors |
| `opening_ev_usd` | Mean regional-board prize in USD |
| `below_ev` | 1 if winnings are below opening EV; 0 otherwise |

Additional stratified-sample fields:

| Field | Meaning |
|---|---|
| `stratum_id` | Numeric region identifier |
| `N_h`, `n_h` | Population count and selected count in the region |
| `pi_i` | First-stage selection probability `n_h/N_h` |
| `w_design` | Inverse selection probability |
| `responded` | Simulated response flag |
| `r_h` | Respondent count in the region |
| `w_final` | `N_h/r_h` for respondents; zero otherwise |

Nonrespondents' survey outcomes and latent behavioural variables are `NA`; frame information remains available. The benchmark also records selection probabilities and weights, with respondent weight `N/r`.

`output/D_estimates.csv` reports `survey` design effects; `D_precision.csv` gives ratios of two realised samples' variances. These have different denominators. `E_performance.csv` gives empirical design effects from the 1,000 repetitions, with coverage MCSE, bias MCSE and average estimated variance. `E_replicates.csv` has 4,000 rows (two outcomes by two designs by 1,000 repetitions).

## Assumptions

The regional boards, especially the Australian high-stakes board, are stylised model inputs. They are not asserted to replicate an identified series. Currency factors are modelling constants. Ages are clipped, not generated from a truncated-normal distribution. The offer plot's line is the median fraction, not its expectation. The optional proportion sensitivity plot uses a fixed design-effect approximation and does not size the primary survey.
