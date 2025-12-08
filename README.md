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

### Time (Elapsed Time)
- **User time** is the CPU time spent executing the function’s own computations. Multithreaded execution accumulates CPU time across cores, so it typically reports higher user time than a single-threaded run.
- **System time** is the CPU time the operating system uses to handle tasks on behalf of the process, such as memory allocation or file I/O.
- **Elapsed time** (wall-clock time) is the total real-world time from start to finish as experienced by the user.

| Threads | Samples | Old (user/sys/elapsed) | New (user/sys/elapsed) |
|:--:|:-------:|:----------------------------:|:------------------------------:|
| 4  | 1       | 32.11 / 0.31 / 32.66         | 17.33 / 0.1 / 17.52            |
| 4  | 10      | 35.21 / 0.13 / 35.52         | 23.06 / 0.2 / 23.38            |
| 4  | 100     | 75.53 / 0.1 / 76.01          | 59.9 / 0.27 / 46.06            |
| 4  | 1000    | 436.16 / 0.33 / 438.47       | 807.29 / 0.83 / 220.09         |
| 4  | 10000   | 3903.03 / 0.91 / 3924.93     | 7811.86 / 23.22 / 1999.87      |
| 4  | 100000  | 34027.72 / 743.25 / 34882.79 | 65122.48 / 727.49 / 16662.49   |
| 8  | 1       | 27.51 / 0.24 / 27.88         | 15.01 / 0.07 / 15.12           |
| 8  | 10      | 30.73 / 0.1 / 30.9           | 19.75 / 0.16 / 19.96           |
| 8  | 100     | 64.44 / 0.08 / 64.65         | 50.19 / 0.16 / 38.38           |
| 8  | 1000    | 371.26 / 0.21 / 372.28       | 670.87 / 1.06 / 102.79         |
| 8  | 10000   | 3326.75 / 0.82 / 3334.58     | 6656.81 / 8.09 / 855.04        |
| 8  | 100000  | 32705.32 / 13 / 32777.81     | 66085.96 / 28.44 / 8358.31     |
| 12 | 1       | 27.64 / 0.24 / 27.98         | 15.17 / 0.09 / 15.29           |
| 12 | 10      | 31.26 / 0.11 / 31.44         | 19.45 / 0.14 / 19.63           |
| 12 | 100     | 64.47 / 0.05 / 64.67         | 50.09 / 0.2 / 38.33            |
| 12 | 1000    | 370.12 / 0.25 / 371.43       | 821.6 / 1.46 / 89.08           |
| 12 | 10000   | 3616.86 / 0.87 / 3630.07     | 8251.81 / 14.23 / 715.02       |
| 12 | 100000  | 33777.23 / 27.64 / 33929.86  | 69783.84 / 35.43 / 5924.52     |
| 24 | 1       | 27.96 / 0.25 / 28.29         | 15.52 / 0.12 / 15.68           |
| 24 | 10      | 30.31 / 0.41 / 30.79         | 19.49 / 0.21 / 19.72           |
| 24 | 100     | 64.66 / 0.07 / 64.88         | 51.09 / 0.22 / 39.1            |
| 24 | 1000    | 375.39 / 0.22 / 376.48       | 523.67 / 1.14 / 72.75          |
| 24 | 10000   | 3789.96 / 1.47 / 3805.26     | 7246.05 / 14.4 / 362.68        |
| 24 | 100000  | 33074.19 / 23.68 / 33205.89  | 71557.8 / 248.73 / 3121.86     |
| 24 | 1000000 | - / - / -                    | 1395025.23 / 950.12 / 60940.62 |

### Memory (MaxRSS, GB)
- MaxRSS (maximum Resident Set Size) is used because it captures the peak amount of physical memory a process occupies during its execution.

| Threads | Samples | Old  | New   |
|:--:|:-------:|:-------:|:--------:|
| 4  | 1       | 1.103   | 1.103    |
| 4  | 10      | 1.103   | 1.104    |
| 4  | 100     | 1.103   | 1.309    |
| 4  | 1000    | 1.103   | 2.213    |
| 4  | 10000   | 4.116   | 44.17    |
| 4  | 100000  | 40.67   | 61.83    |
| 8  | 1       | 1.103   | 1.110    |
| 8  | 10      | 1.103   | 1.103    |
| 8  | 100     | 1.103   | 1.389    |
| 8  | 1000    | 1.103   | 3.364    |
| 8  | 10000   | 4.116   | 44.46    |
| 8  | 100000  | 40.67   | 62.01    |
| 12 | 1       | 1.103   | 1.103    |
| 12 | 10      | 1.103   | 1.103    |
| 12 | 100     | 1.103   | 1.392    |
| 12 | 1000    | 1.103   | 4.516    |
| 12 | 10000   | 4.116   | 44.46    |
| 12 | 100000  | 40.67   | 62.01    |
| 24 | 1       | 1.103   | 1.113    |
| 24 | 10      | 1.103   | 1.125    |
| 24 | 100     | 1.103   | 1.404    |
| 24 | 1000    | 1.103   | 5.678    |
| 24 | 10000   | 4.116   | 44.18    |
| 24 | 100000  | 40.67   | 61.72    |
| 24 | 1000000 | -       | 302.68   |

