# TRACE

An R package for **stable-isotope CN labeling analysis** of LC–MS metabolomics data. TRACE identifies carbon–nitrogen labeling patterns in multi-isotope experiments, builds feature networks (CN labels, adducts, isotopes, fragments), assigns metabolite seeds, and annotates features against a compound database. It is a refactored R implementation of the PAVE-style workflow, built on top of [MSdev](https://github.com/) objects and xcms peak tables.

## Overview

TRACE analyzes LC–MS data from experiments with four isotope-labeled sample groups:

| Sample type | Labeling pattern | Role                                     |
|-------------|------------------|------------------------------------------|
| `S12C14N`   | C0N0 (unlabeled) | Reference / light carbon, light nitrogen |
| `S12C15N`   | C0Ny             | Light carbon, heavy nitrogen             |
| `S13C14N`   | CxN0             | Heavy carbon, light nitrogen             |
| `S13C15N`   | CxNy             | Heavy carbon, heavy nitrogen             |

For each xcms feature, TRACE:

1.  **Detects CN labeling networks** — matches m/z differences to theoretical C/N isotope shifts and scores patterns against ideal labeling ratios (`TRACE_cor`).
2.  **Applies dynamic filtering** — estimates instrument-specific m/z and RT error distributions and removes unreliable edges.
3.  **Assigns seed networks** — partitions CN seeds into compatible subnetworks using adduct, isotope, and fragment connectivity (PAVE-style global assignment).
4.  **Annotates metabolites** — matches seeds to a compound database and scores adduct / RT consistency.
5.  **Adjusts labeling ratios** — evaluates group-level ratio corrections from high-confidence hits and optionally reconstructs CN networks.

Results are stored on the MSdev object at `object@advancedAna$TRACE` (final tables) and `object@advancedAna$TRACE_temp` (intermediate networks, filters, and ratios).

## Dependencies

Install MSdev and its dependencies first, then install TRACE from source.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("xcms", "MSCC"))
install.packages(c("data.table", "openxlsx", "dplyr", "ggplot2", "igraph",
                 "BiocParallel", "patchwork", "ggrastr", "tidyr", "devtools"))

devtools::install("path/to/MSdev")
devtools::install("path/to/TRACE")
```

## Usage

Load a pre-processed MSdev object (with xcms peak data and `sample.type` set to the four isotope groups), run TRACE, and export results:

``` r
library(MSdev)
library(TRACE)

obj <- MSdev_load("path/to/MSdev_processed.Rdata")

obj <- TRACE_workflow(
  obj,
  rt.tol = 10,
  ppm = 5,
  cpdb = "path/to/trace.cp.db.xlsx",
  eval_top = 0.2,
  ratio.plot = FALSE,
  ratio.reconstruct = TRUE
)

TRACE_export(obj, file = "TRACE_results.xlsx")
```

For step-by-step control per polarity (`i.pol`: `0` = negative, `1` = positive):

``` r
obj <- TRACE_get_CN_net(obj, i.pol = 0, rt.tol = 10, ppm = 5)
obj <- TRACE_CN_labelling_ratio_adjust(obj, eval_top = 0.2, reconstruct = TRUE)
obj <- TRACE_dynamic_filter(obj, i.pol = 0)
obj <- TRACE_network_assignment(obj, i.pol = 0)
obj <- TRACE_annotate(obj, i.pol = 0, cpdb = "path/to/trace.cp.db.xlsx")
# repeat for i.pol = 1

TRACE_export(obj, file = "TRACE_results.xlsx")
```

## Main functions

| Function | Description |
|----|----|
| `TRACE_workflow()` | End-to-end TRACE pipeline |
| `TRACE_get_CN_net()` | Build CN labeling candidate networks |
| `TRACE_dynamic_filter()` | Dynamic m/z / RT error filtering |
| `TRACE_network_assignment()` | PAVE-style seed network assignment |
| `TRACE_annotate()` | Compound database annotation |
| `TRACE_CN_labelling_ratio_adjust()` | Group ratio correction and optional CN net rebuild |
| `get_TRACE_CN_labelling_ratio()` | Per-seed labeling ratios across four isotope groups |
| `TRACE_export()` | Export results to Excel |

## Vignettes

| Vignette                | Topic                                      |
|-------------------------|--------------------------------------------|
| `TRACE_Demo`            | End-to-end demo with parameter sensitivity |
| `Compare_PAVE`          | Side-by-side comparison with PAVE-Matlab   |
| `instrument_evaluation` | Multi-instrument TRACE performance         |
| `cor_sensitivity`       | Sensitivity to `TRACE_cor` cutoff          |
| `noise_evaluation`      | Noise and blank evaluation                 |

## License

See `DESCRIPTION` for license information.
