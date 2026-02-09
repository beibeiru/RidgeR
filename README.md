
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RidgeR: Ridge Regression with Significance Testing

## Installation

To install `RidgeR`, we recommend using `devtools`:

``` r
# install.packages("devtools")
devtools::install_github("beibeiru/RidgeR")
```

**System Requirements:** GNU Scientific Library (GSL)

The package has been installed successfully on Operating Systems:

- Red Hat Enterprise Linux 8.10 (Ootpa)
- macOS Sequoia 15.3.1
- Windows 10

## Functions

### Main function

`SecAct.inference.gsl.new` is the primary interface. It dispatches to one of 4 optimized C backends via the `method` parameter:

``` r
SecAct.inference.gsl.new(Y, method = "Tcol.mt")  # default: T-col, multi-threaded
SecAct.inference.gsl.new(Y, method = "Tcol.st")  # T-col, single-threaded
SecAct.inference.gsl.new(Y, method = "Yrow.mt")  # Y-row, multi-threaded
SecAct.inference.gsl.new(Y, method = "Yrow.st")  # Y-row, single-threaded
```

### All variants

| # | Function | Permutation | Threading | Interface | Description |
|---|----------|-------------|-----------|-----------|-------------|
| 1 | `SecAct.inference.naive` | T-col (pure R) | Single (R) | `.Call` (perm table only) | Pure R reference implementation |
| 2 | `SecAct.inference.Yrow.st` | Y-row | Single | `.Call` (ncores=1) | Single-threaded Y-row permutation |
| 3 | `SecAct.inference.Tcol.st` | T-col | Single | `.Call` (ncores=1) | Single-threaded T-col permutation |
| 4 | `SecAct.inference.Yrow.mt` | Y-row | Multi (OMP) | `.Call` | Multi-threaded Y-row permutation |
| 5 | `SecAct.inference.Tcol.mt` | T-col | Multi (OMP) | `.Call` | Multi-threaded T-col permutation |
| — | `SecAct.inference.gsl.old` | Y-row | Single | `.C` (legacy, 32-bit) | Legacy single-threaded (preserved) |
| — | `SecAct.inference.gsl.new` | Dispatches | Dispatches | `.Call` (64-bit) | Dispatcher (default: `method="Tcol.mt"`) |

**Permutation strategies:**

- **Y-row**: Permutes rows of the response matrix Y, then multiplies T * Y_perm. Parallelizes over sample strips.
- **T-col**: Permutes columns of the projection matrix T, then multiplies T_perm * Y. Mathematically equivalent (`T[:, inv_perm] @ Y == T @ Y[perm, :]`). Parallelizes over permutations. Matches [SecActPy](https://github.com/data2intelligence/SecActpy) approach.

**Legacy:** `gsl.old` is preserved for backward compatibility. `gsl.new` now dispatches to the 4 variants (previously it was Y-row multi-threaded only).

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Y` | — | Gene expression matrix (genes x samples) |
| `SigMat` | `"SecAct"` | Signature matrix: `"SecAct"` (bundled) or path to file |
| `lambda` | `5e+05` | Ridge regularization parameter |
| `nrand` | `1000` | Number of permutations |
| `ncores` | `NULL` | Number of CPU cores (`NULL` = auto-detect; multi-threaded variants only) |
| `rng_method` | `"srand"` | RNG backend: `"srand"` (C stdlib, matches R SecAct) or `"gsl"` (cross-platform GSL MT19937) |
| `method` | `"Tcol.mt"` | Backend variant: `"Tcol.mt"`, `"Tcol.st"`, `"Yrow.mt"`, or `"Yrow.st"` (`gsl.new` only) |
| `is.group.sig` | `TRUE` | Group correlated signatures before regression |
| `is.group.cor` | `0.9` | Correlation threshold for signature grouping |

## Example

``` r
library(RidgeR)

dataPath <- file.path(system.file(package = "RidgeR"), "extdata/")
expr.diff <- read.table(paste0(dataPath, "Ly86-Fc_vs_Vehicle_logFC.txt"))

# ---- Compare all variants ----
t_naive   <- system.time({res.naive   <- SecAct.inference.naive(expr.diff)})
t_yrow_st <- system.time({res.yrow.st <- SecAct.inference.Yrow.st(expr.diff)})
t_tcol_st <- system.time({res.tcol.st <- SecAct.inference.Tcol.st(expr.diff)})
t_yrow_mt <- system.time({res.yrow.mt <- SecAct.inference.Yrow.mt(expr.diff)})
t_tcol_mt <- system.time({res.tcol.mt <- SecAct.inference.Tcol.mt(expr.diff)})

# ---- Verify equivalence (all should be ~0) ----
cat("naive  vs Tcol.st:", max(abs(res.naive$zscore - res.tcol.st$zscore)), "\n")
cat("Yrow.st vs Yrow.mt:", max(abs(res.yrow.st$zscore - res.yrow.mt$zscore)), "\n")
cat("Tcol.st vs Tcol.mt:", max(abs(res.tcol.st$zscore - res.tcol.mt$zscore)), "\n")

# ---- Elapsed times ----
cat("naive:", t_naive[3], "Yrow.st:", t_yrow_st[3], "Tcol.st:", t_tcol_st[3],
    "Yrow.mt:", t_yrow_mt[3], "Tcol.mt:", t_tcol_mt[3], "\n")
```

## Reproducibility

RidgeR supports two RNG backends via the `rng_method` parameter:

| `rng_method` | Description | Use case |
|---|---|---|
| `"srand"` (default) | C stdlib `srand()`/`rand()` | Match original R SecAct results **on the same platform** |
| `"gsl"` | GSL Mersenne Twister (MT19937) | **Cross-platform** reproducibility; matches [SecActPy](https://github.com/data2intelligence/SecActpy) `rng_method='gsl'` |

``` r
# Default: match R SecAct on same platform (C stdlib rand, platform-dependent)
res <- SecAct.inference.Yrow.mt(expr.diff, rng_method = "srand")

# Cross-platform: matches SecActPy rng_method='gsl' on any OS
res <- SecAct.inference.Yrow.mt(expr.diff, rng_method = "gsl")
```

> **Note:** C `rand()` implementations differ across operating systems, so
> `rng_method="srand"` produces platform-dependent results. Use
> `rng_method="gsl"` when results must be reproducible across Linux, macOS,
> and Windows, or when comparing with
> [SecActPy](https://github.com/data2intelligence/SecActpy).

## Benchmark Results
- A 32-bit signed integer can index up to 2^31 − 1 = 2,147,483,647 (2.147 billion elements).
- A 64-bit signed integer can index up to 2^63 − 1 = 9.22 × 10^18 elements (about 9.22 quintillion).
- A dataset of 16,325 genes × 100,000 samples contains 1.63 billion elements, which is below the 32-bit limit, so it can be indexed with 32-bit integers (works with the .C interface in the old function).
- A dataset of 16,325 genes × 1,000,000 samples contains 16.3 billion elements, which exceeds the 32-bit limit and therefore cannot be indexed with 32-bit integers (fails under the .C interface in the old function).
- The same dataset (16,325 × 1,000,000) is well within the 64-bit indexing range, so it can be indexed with 64-bit integers (works with the .Call interface in the new function).

### Time (Elapsed Time)
- **User time** is the CPU time spent executing the function's own computations. Multithreaded execution accumulates CPU time across cores, so it typically reports higher user time than a single-threaded run.
- **System time** is the CPU time the operating system uses to handle tasks on behalf of the process, such as memory allocation or file I/O.
- **Elapsed time** (wall-clock time) is the total real-world time from start to finish as experienced by the user.

| Threads | Samples | Old (user/sys/elapsed) | New (user/sys/elapsed) |
|:-:|:-------:|:-----------------------------:|:-----------------------------:|
| 2   | 1       | 44.74 / 0.38 / 45.33        | 22.17 / 0.09 / 22.4          |
| 2   | 10      | 44.38 / 0.26 / 44.83        | 30.11 / 0.23 / 30.51         |
| 2   | 100     | 98.87 / 0.17 / 99.42        | 111.44 / 0.38 / 76.45        |
| 2   | 1000    | 536.43 / 0.45 / 538.98      | 1059.65 / 2.46 / 549.39      |
| 2   | 10000   | 4933.87 / 45.55 / 4999.32   | 10342.19 / 190.66 / 5344.48  |
| 2   | 100000  | 50265.06 / 3047.56 / 53531.81 | 98307.27 / 3574.56 / 51537.88 |
| 4   | 1       | 32.11 / 0.31 / 32.66        | 17.33 / 0.1 / 17.52          |
| 4   | 10      | 35.21 / 0.13 / 35.52        | 23.06 / 0.2 / 23.38          |
| 4   | 100     | 75.53 / 0.1 / 76.01         | 59.9 / 0.27 / 46.06          |
| 4   | 1000    | 436.16 / 0.33 / 438.47      | 807.29 / 0.83 / 220.09       |
| 4   | 10000   | 3903.03 / 0.91 / 3924.93    | 7811.86 / 23.22 / 1999.87    |
| 4   | 100000  | 34027.72 / 743.25 / 34882.79 | 65122.48 / 727.49 / 16662.49 |
| 4   | 1000000 | - / - / -                   | 895279.15 / 12060.68 / 231333.51 |
| 8   | 1       | 27.51 / 0.24 / 27.88        | 15.01 / 0.07 / 15.12         |
| 8   | 10      | 30.73 / 0.1 / 30.9          | 19.75 / 0.16 / 19.96         |
| 8   | 100     | 64.44 / 0.08 / 64.65        | 50.19 / 0.16 / 38.38         |
| 8   | 1000    | 371.26 / 0.21 / 372.28      | 670.87 / 1.06 / 102.79       |
| 8   | 10000   | 3326.75 / 0.82 / 3334.58    | 6656.81 / 8.09 / 855.04      |
| 8   | 100000  | 32705.32 / 13 / 32777.81    | 66085.96 / 28.44 / 8358.31   |
| 8   | 1000000 | - / - / -                   | 1026510.18 / 21167.11 / 133237.09 |
| 12  | 1       | 27.64 / 0.24 / 27.98        | 15.17 / 0.09 / 15.29         |
| 12  | 10      | 31.26 / 0.11 / 31.44        | 19.45 / 0.14 / 19.63         |
| 12  | 100     | 64.47 / 0.05 / 64.67        | 50.09 / 0.2 / 38.33          |
| 12  | 1000    | 370.12 / 0.25 / 371.43      | 821.6 / 1.46 / 89.08         |
| 12  | 10000   | 3616.86 / 0.87 / 3630.07    | 8251.81 / 14.23 / 715.02     |
| 12  | 100000  | 33777.23 / 27.64 / 33929.86 | 69783.84 / 35.43 / 5924.52   |
| 12  | 1000000 | - / - / -                   | 1058682.50 / 6238.49 / 91299.53 |
| 24  | 1       | 27.96 / 0.25 / 28.29        | 15.52 / 0.12 / 15.68         |
| 24  | 10      | 30.31 / 0.41 / 30.79        | 19.49 / 0.21 / 19.72         |
| 24  | 100     | 64.66 / 0.07 / 64.88        | 51.09 / 0.22 / 39.1          |
| 24  | 1000    | 375.39 / 0.22 / 376.48      | 523.67 / 1.14 / 72.75        |
| 24  | 10000   | 3789.96 / 1.47 / 3805.26    | 7246.05 / 14.4 / 362.68      |
| 24  | 100000  | 33074.19 / 23.68 / 33205.89 | 71557.8 / 248.73 / 3121.86   |
| 24  | 1000000 | - / - / -                   | 1067790.89 / 509.82 / 46415.26 |


### Memory (MaxRSS, GB)
- MaxRSS (maximum Resident Set Size) is used because it captures the peak amount of physical memory a process occupies during its execution.

| Threads | Samples | Old   | New    |
|:-:|:-------:|:------:|:------:|
| 4   | 1       | 1.103  | 1.103  |
| 4   | 10      | 1.103  | 1.104  |
| 4   | 100     | 1.103  | 1.309  |
| 4   | 1000    | 1.103  | 2.213  |
| 4   | 10000   | 4.116  | 44.17  |
| 4   | 100000  | 40.67  | 61.83  |
| 4   | 1000000 | -      | 302.68 |
| 8   | 1       | 1.103  | 1.110  |
| 8   | 10      | 1.103  | 1.103  |
| 8   | 100     | 1.103  | 1.389  |
| 8   | 1000    | 1.103  | 3.364  |
| 8   | 10000   | 4.116  | 44.46  |
| 8   | 100000  | 40.67  | 62.01  |
| 8   | 1000000 | -      | 302.68 |
| 12  | 1       | 1.103  | 1.103  |
| 12  | 10      | 1.103  | 1.103  |
| 12  | 100     | 1.103  | 1.392  |
| 12  | 1000    | 1.103  | 4.516  |
| 12  | 10000   | 4.116  | 44.46  |
| 12  | 100000  | 40.67  | 62.01  |
| 12  | 1000000 | -      | 302.68 |
| 24  | 1       | 1.103  | 1.113  |
| 24  | 10      | 1.103  | 1.125  |
| 24  | 100     | 1.103  | 1.404  |
| 24  | 1000    | 1.103  | 5.678  |
| 24  | 10000   | 4.116  | 44.18  |
| 24  | 100000  | 40.67  | 61.72  |
| 24  | 1000000 | -      | 302.68 |
