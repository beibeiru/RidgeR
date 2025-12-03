
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SecAct: Secreted Protein Activity Inference

## Installation

To install `RidgeR`, we recommend using `devtools`:

``` r
# install.packages("devtools")
devtools::install_github("beibeiru/RidgeR")
```

The package has been installed successfully on Operating Systems:

- Red Hat Enterprise Linux 8.10 (Ootpa)
- macOS Sequoia 15.3.1

## Example

``` r
library(RidgeR)

# Load example data
dataPath <- file.path(system.file(package = "RidgeR"), "extdata/")
expr.diff <- read.table(paste0(dataPath, "Ly86-Fc_vs_Vehicle_logFC.txt"))

# ---- Compare execution time ----
t_old <- system.time({res.old <- SecAct.inference.gsl.old(expr.diff)})
t_new <- system.time({res.new <- SecAct.inference.gsl.new(expr.diff)})
print(t_old); print(t_new)

# ---- Compare output: head of z-scores ----
print(head(res.old$zscore))
print(head(res.new$zscore))

# ---- summarize differences ----
diff <- res.new$zscore - res.old$zscore
print(summary(diff))

# ---- correlation ----
corr <- cor(res.new$zscore, res.old$zscore)
print(corr)
```

## Benchmark Results

| nthreads | nsamples | old_user | old_system | old_elapsed | new_user | new_system | new_elapsed |
|----------|----------|----------|------------|--------------|----------|------------|--------------|
| 4 | 1    | 27.625 | 0.266 | 27.995 | 14.627 | 0.058 | 14.720 |
| 4 | 10   | 31.233 | 0.104 | 31.427 | 20.782 | 0.175 | 21.029 |
| 4 | 100  | 64.277 | 0.118 | 64.568 | 49.821 | 0.221 | 38.042 |
| 4 | 1000 | 366.639 | 0.205 | 367.763 | 673.281 | 0.719 | 182.062 |
| 4 | 10000 |  — | — | — | — | — | — |
| 8 | 1    | 27.663 | 0.294 | 28.045 | 14.990 | 0.099 | 15.143 |
| 8 | 10   | 31.335 | 0.150 | 31.553 | 19.757 | 0.155 | 19.973 |
| 8 | 100  | 64.742 | 0.118 | 65.025 | 49.349 | 0.190 | 37.471 |
| 8 | 1000 | 364.175 | 0.230 | 365.210 | 703.888 | 2.050 | 105.972 |
| 8 | 10000 | — | — | — | — | — | — |
| 12 | 1    | 41.010 | 0.480 | 41.644 | 20.404 | 0.115 | 20.558 |
| 12 | 10   | 39.211 | 0.156 | 39.474 | 25.719 | 0.115 | 25.903 |
| 12 | 100  | 74.971 | 0.161 | 75.402 | 55.432 | 0.236 | 43.653 |
| 12 | 1000 | 388.420 | 0.220 | 389.433 | 690.070 | 2.571 | 86.006 |
| 12 | 10000 | — | — | — | 7002.968 | 99.050 | 629.987 |

