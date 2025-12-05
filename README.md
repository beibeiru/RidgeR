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
- Windows 10

## Example

``` r
library(RidgeR)

# ---- Load example data ----
dataPath <- file.path(system.file(package = "RidgeR"), "extdata")
expr.diff <- read.table(paste0(dataPath, "/Ly86-Fc_vs_Vehicle_logFC.txt"))

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

### Elapsed Time

| Threads | Samples | Old (user/sys/elapsed) | New (user/sys/elapsed) |
|:--:|:------:|:----------------------------:|:--------------------------:|
| 4  | 1      | 27.28 / 0.21 / 27.62         | 14.45 / 0.02 / 14.49       |
| 4  | 10     | 30.55 / 0.07 / 30.71         | 19.31 / 0.09 / 19.44       |
| 4  | 100    | 63.97 / 0.05 / 64.21         | 49.13 / 0.14 / 37.19       |
| 4  | 1000   | 361.94 / 0.19 / 363.19       | 658.43 / 0.36 / 182.86     |
| 4  | 10000  | 3320.56 / 0.83 / 3332.15     | 6517.87 / 2.01 / 1659.43   |
| 4 | 100000  | 32727.57 / 15.36 / 32818.94  | 64459.96 / 6.38 / 16226.34 |
| 8  | 1      | 28.26 / 0.26 / 28.78         | 15.96 / 0.03 / 15.83       |
| 8  | 10     | 33.01 / 0.19 / 33.27         | 21.81 / 0.14 / 21.96       |
| 8  | 100    | 68.16 / 0.07 / 68.36         | 52.54 / 0.23 / 40.56       |
| 8  | 1000   | 373.93 / 0.23 / 374.84       | 679.82 / 1.46 / 106.45     |
| 8  | 10000  | 3418.41 / 14.74 / 3440.61    | 6675.16 / 10.38 / 879.37   |
| 8  | 100000 | 34457.42 / 1264.86 / 35800.8 | 64582.53 / 26.61 / 8190.44 |
| 12 | 1      | 27.69 / 0.22 / 28.02         | 14.75 / 0.03 / 14.59       |
| 12 | 10     | 30.72 / 0.09 / 30.91         | 19.28 / 0.14 / 19.42       |
| 12 | 100    | 64.25 / 0.02 / 64.48         | 49.28 / 0.14 / 37.31       |
| 12 | 1000   | 365.91 / 0.25 / 367.37       | 621.72 / 0.71 / 84.35      |
| 12 | 10000  | 3329.04 / 0.86 / 3342.09     | 6575.52 / 1.6 / 573.09     |
| 12 | 100000 | 32962.06 / 15.71 / 33074.79  | 65136.61 / 5.31 / 5542.53  |
| 24 | 1      | 27.78 / 0.23 / 28.13         | 15.12 / 0.03 / 14.69       |
| 24 | 10     | 30.74 / 0.1 / 30.93          | 19.4 / 0.11 / 19.45        |
| 24 | 100    | 64.39 / 0.02 / 64.56         | 49.5 / 0.15 / 37.49        |
| 24 | 1000   | 368.54 / 0.21 / 369.69       | 473.6 / 0.85 / 53.03       |
| 24 | 10000  | 3309.37 / 0.93 / 3318.24     | 7236.37 / 2.47 / 338.5     |
| 24 | 100000 | 32750.12 / 13.66 / 32835.15  | 73382.74 / 6.3 / 3162.55   |


### Memory (MaxRSS, GB)

| Threads | Samples | Old  | New   |
|:--:|:------:|:-------:|:--------:|
| 4  | 1      | 1.103   | 1.103    |
| 4  | 10     | 1.103   | 1.104    |
| 4  | 100    | 1.103   | 1.309    |
| 4  | 1000   | 1.103   | 2.213    |
| 8  | 1      | 1.103   | 1.110    |
| 8  | 10     | 1.103   | 1.103    |
| 8  | 100    | 1.103   | 1.389    |
| 8  | 1000   | 1.103   | 3.364    |
| 12 | 1      | 1.103   | 1.103    |
| 12 | 10     | 1.103   | 1.103    |
| 12 | 100    | 1.103   | 1.392    |
| 12 | 1000   | 1.103   | 4.516    |
| 24 | 1      | 1.103   | 1.113    |
| 24 | 10     | 1.103   | 1.125    |
| 24 | 100    | 1.103   | 1.404    |
| 24 | 1000   | 1.103   | 5.678    |


