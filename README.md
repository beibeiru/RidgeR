
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

| Threads | Samples | Old (user/sys/elapsed) | New (user/sys/elapsed) |
|:-------:|:-------:|:-----------------------:|:-----------------------:|
| 4  | 1      | 27.63 / 0.27 / 27.99 | 14.63 / 0.06 / 14.72 |
| 4  | 10     | 31.23 / 0.10 / 31.43 | 20.78 / 0.18 / 21.03 |
| 4  | 100    | 64.28 / 0.12 / 64.57 | 49.82 / 0.22 / 38.04 |
| 4  | 1000   | 366.64 / 0.21 / 367.76 | 673.28 / 0.72 / 182.06 |
| 4  | 10000  | 3336.45	/ 1.92 / 3346.15 | 6566.84 / 8.02 / 1668.69 |
| 8  | 1      | 27.66 / 0.29 / 28.05 | 14.99 / 0.10 / 15.14 |
| 8  | 10     | 31.34 / 0.15 / 31.55 | 19.76 / 0.16 / 19.97 |
| 8  | 100    | 64.74 / 0.12 / 65.03 | 49.35 / 0.19 / 37.47 |
| 8  | 1000   | 364.18 / 0.23 / 365.21 | 703.89 / 2.05 / 105.97 |
| 8  | 10000  | — | 6948.65 / 115.15 / 929.69 | 
| 12 | 1      | 41.01 / 0.48 / 41.64 | 20.40 / 0.12 / 20.56 |
| 12 | 10     | 39.21 / 0.16 / 39.47 | 25.72 / 0.12 / 25.90 |
| 12 | 100    | 74.97 / 0.16 / 75.40 | 55.43 / 0.24 / 43.65 |
| 12 | 1000   | 388.42 / 0.22 / 389.43 | 690.07 / 2.57 / 86.01 |
| 12 | 10000  | — | 7002.97 / 99.05 / 629.99 |


