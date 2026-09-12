# Enhanced Elastic Net (EEN) for High-Dimensional Regression

[![R version](https://img.shields.io/badge/R-%3E%3D%204.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Code to reproduce all results in:

> **"An Enhanced Elastic Net for High-Dimensional Regression with Structured Coefficient Patterns"**
> E. N. N. Nortey & E. Somua-Wiafe
> Computational Statistics & Data Analysis (Submitted)

## Overview

EEN is a penalized regression estimator for high-dimensional problems where predictors have a natural ordering (spectroscopy wavelengths, spatial or temporal indices, etc.). It combines an adaptive $\ell_1$ penalty, a ridge initialization, and a second-difference smoothness penalty, fit by coordinate descent.

The repository reproduces the 500-replication simulation study (smooth two-bump signal, AR(1)-correlated predictors, compared against LASSO, Elastic Net, and Adaptive LASSO) and the corn NIR spectroscopy application.
## Repository Structure

| File / Directory | Description |
| :--- | :--- |
| `een_simulation_study.R` | Simulation study. $n = 150$, $p = 500$, $\rho = 0.7$, 500 replications run in four phases (50/150/300/500) to check Monte Carlo convergence. |
| `corn_empirical.R` | Corn NIR application (`pcv::corn`, $n = 80$, $p = 700$): 20 repetitions of 10-fold CV (200 train/test splits), plus the full-sample coefficient fit. |
<!--| `simulation_results/` | Simulation figures and tables (L2 error, TPR, FDR, MCC; cross-phase table). |
| `corn_pcv_moisture_results/` | Corn NIR figures and tables (RMSE comparison, coefficient profiles, roughness/sparsity). |
| `checkpoints/` | Checkpoint files written after each simulation phase. |
-->

## Requirements

R >= 4.0, with:

```r
install.packages(c("MASS", "glmnet", "Matrix", "parallel", "doParallel",
                    "foreach", "ggplot2", "dplyr", "tidyr", "pcv", "gridExtra"))
```

## Running the code

```r
source("een_simulation_study.R")   # writes to simulation_results/
source("corn_empirical.R")         # writes to corn_pcv_moisture_results/
```

The simulation script checkpoints after each phase (50/150/300/500 reps) to `checkpoints/`, so an interrupted run can be resumed rather than restarted from scratch.

## Notes on reproducing exact numbers

Replication seeds follow `seed = 2026 + r`, matching the paper. `sessionInfo.txt` records the package versions the paper's numbers were generated under. `glmnet` in particular has changed its internal CV behavior across versions, so a different version may shift the reported numbers slightly, though not the ranking between methods. Parallelization (via `parallel`/`doParallel`/`foreach`) affects runtime, not the results.

## Citation

```bibtex
@article{nortey2026een,
  title   = {An Enhanced Elastic Net for High-Dimensional Regression with Structured Coefficient Patterns},
  author  = {Nortey, E. N. N. and Somua-Wiafe, E.},
  journal = {Computational Statistics \& Data Analysis},
  year    = {2026},
  note    = {Submitted}
}
```

## Data

The corn NIR data are public, distributed via the `pcv` package for R (originally from Eigenvector Research, Inc.). Nothing proprietary is used.

## License

MIT. See `LICENSE`.

## Contact

E. N. N. N. Nortey, ennnortey@ug.edu.gh;
E. Somua-Wiafe, esomua-wiafe@ug.edu.gh
