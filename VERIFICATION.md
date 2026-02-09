# Verification: RNG Paths and Cross-Package Equivalence

This document describes how to verify correctness of the 7 RidgeR variants
and their equivalence with R SecAct and SecActPy.

---

## 1. RNG Backends: `srand` vs `gsl`

RidgeR supports two RNG backends for the permutation shuffle:

| Backend | Seed | Shuffle Formula | Bit Width | Platform |
|---------|------|-----------------|-----------|----------|
| `srand` | `srand(0)` | `j = i + rand() / (RAND_MAX / (n-i) + 1)` | 31-bit (RAND_MAX = 2^31-1) | Platform-dependent |
| `gsl` | `gsl_rng_set(rng, 0)` → internal seed 4357 | `j = i + gsl_rng_get() / (2^32-1 / (n-i) + 1)` | 32-bit (MT19937) | Cross-platform deterministic |

### Key differences

1. **Seed mapping**: C stdlib `srand(0)` uses seed 0 directly. GSL MT19937 with
   `gsl_rng_set(rng, 0)` internally maps seed 0 → 4357.

2. **RNG range**: C `rand()` returns values in `[0, RAND_MAX]` where RAND_MAX
   is typically 2^31-1 on Linux/macOS. GSL MT19937 returns values in
   `[0, 2^32-1]`.

3. **Platform behavior**: `srand`/`rand()` produces different sequences on
   different platforms (glibc vs macOS vs Windows). GSL MT19937 is identical
   on all platforms.

4. **Cumulative shuffle**: Both backends use Fisher-Yates with cumulative state
   (the array is NOT reset between permutations — each shuffle continues from
   the previous state). This is critical for reproducibility.

### Verification: srand path

```r
# On Linux (glibc), srand results match R SecAct on the same Linux platform
res_ridger <- SecAct.inference.Yrow.st(expr.diff, rng_method = "srand")
res_secact <- SecAct.inference(expr.diff)  # from original R SecAct package
max(abs(res_ridger$zscore - res_secact$zscore))
# Expected: 0 (exact match on same platform)
```

### Verification: gsl path

```r
# GSL results are identical across all platforms and all RidgeR variants
res_gsl <- SecAct.inference.Yrow.st(expr.diff, rng_method = "gsl")
# This should produce the same result on Linux, macOS, and Windows
```

---

## 2. Permutation Strategies: Y-row vs T-col

### Mathematical equivalence

Given projection matrix T (p x n) and response Y (n x m):

- **Y-row permutation**: `beta_perm = T @ Y[perm, :]`
- **T-col permutation**: `beta_perm = T[:, inv_perm] @ Y`

These are mathematically identical because:
```
T[:, inv_perm] @ Y = T @ P^T @ Y = T @ Y[perm, :]
```
where P is the permutation matrix corresponding to `perm`, and
`inv_perm[perm[j]] = j`.

### Forward vs inverse permutation

The **same** Fisher-Yates shuffle sequence produces `perm` (forward permutation).
The T-col variants compute the **inverse**: `inv_perm[perm[j]] = j`.

Both in RidgeR C code:
```c
// Y-row: store forward perm directly
memcpy(&perm_table[i_rand * n], temp_idx, n * sizeof(int));

// T-col: compute inverse
for (j = 0; j < n; j++) inv_row[temp_idx[j]] = j;
```

And in SecActPy Python code:
```python
# Forward shuffle
self.shuffle_inplace(arr)
# Compute inverse
inv_perm[arr] = np.arange(n)
```

### Verification: Y-row vs T-col equivalence

```r
# Single-threaded: exact match (deterministic accumulation order)
res_yrow <- SecAct.inference.Yrow.st(expr.diff, rng_method = "gsl")
res_tcol <- SecAct.inference.Tcol.st(expr.diff, rng_method = "gsl")
max(abs(res_yrow$zscore - res_tcol$zscore))
# Expected: 0 (exact match, single-threaded)

# Multi-threaded: near-exact (floating-point accumulation order may differ)
res_yrow_mt <- SecAct.inference.Yrow.mt(expr.diff, rng_method = "gsl")
res_tcol_mt <- SecAct.inference.Tcol.mt(expr.diff, rng_method = "gsl")
max(abs(res_yrow_mt$zscore - res_tcol_mt$zscore))
# Expected: < 1e-10 (floating-point tolerance from different accumulation order)
```

---

## 3. Naive Variant Verification

The naive variant uses pure R linear algebra with a C-generated permutation table.

```r
res_naive <- SecAct.inference.naive(expr.diff, rng_method = "gsl")
res_tcol  <- SecAct.inference.Tcol.st(expr.diff, rng_method = "gsl")
max(abs(res_naive$zscore - res_tcol$zscore))
# Expected: < 1e-10 (R LAPACK vs GSL Cholesky may differ slightly)
```

The naive variant gets its permutation table from the same C RNG code
(`generate_perm_table`), so permutation sequences are identical. Small
differences arise from R's `solve()` (LAPACK) vs C's GSL Cholesky for
computing the projection matrix T.

---

## 4. Cross-Package Equivalence

### 4.1 RidgeR ↔ R SecAct (original)

| RidgeR Setting | R SecAct Equivalent | Match Condition |
|---|---|---|
| `rng_method = "srand"` | Default behavior | Same platform (glibc) |
| Any Y-row variant, single-threaded | `SecAct.inference()` | Same RNG + same platform |

```r
# On Linux with glibc:
res_ridger <- SecAct.inference.gsl.old(expr.diff, rng_method = "srand")
res_secact <- SecAct.inference(expr.diff)  # original R SecAct
max(abs(res_ridger$zscore - res_secact$zscore))
# Expected: 0
```

### 4.2 RidgeR ↔ SecActPy

| RidgeR Setting | SecActPy Equivalent | Match Condition |
|---|---|---|
| `rng_method = "gsl"` | `rng_method = "gsl"` | Cross-platform (any OS) |
| Any variant | `ridge()` or `_ridge_permutation_numpy()` | Same rng_method = "gsl" |

SecActPy uses T-column permutation with inverse permutation tables.
RidgeR's T-col variants use the identical algorithm. RidgeR's Y-row variants
produce mathematically identical results.

```r
# RidgeR (R)
res_r <- SecAct.inference.gsl.new(expr.diff, rng_method = "gsl")
```

```python
# SecActPy (Python)
import secactpy
res_py = secactpy.ridge(Y, X, rng_method='gsl')
```

```
# Compare (in R after loading Python results):
max(abs(res_r$zscore - res_py_zscore))
# Expected: < 1e-10
```

### 4.3 RNG Sequence Equivalence Chain

All three packages produce the **same permutation sequence** with `rng_method = "gsl"`:

```
R SecAct (C srand path)  ──── platform-dependent ────  RidgeR (srand path)
                                                            │
                                                        same platform
                                                            │
                                                       R SecAct original

RidgeR (gsl path)  ──── cross-platform identical ────  SecActPy (gsl path)
       │                                                    │
  GSL MT19937                                          GSL MT19937
  seed 0 → 4357                                        seed 0 → 4357
  Fisher-Yates                                         Fisher-Yates
  cumulative state                                     cumulative state
```

---

## 5. All 7 Variants: Expected Equivalence Matrix

With `rng_method = "gsl"`, all variants should produce equivalent results:

| Comparison | Expected Max Diff | Reason |
|---|---|---|
| `naive` vs `Tcol.st` | < 1e-10 | R LAPACK vs GSL Cholesky |
| `Yrow.st` vs `Tcol.st` | 0 | Exact (same accumulation order) |
| `Yrow.st` vs `Yrow.mt` | < 1e-10 | MT floating-point accumulation |
| `Tcol.st` vs `Tcol.mt` | < 1e-10 | MT floating-point accumulation + atomics |
| `gsl.old` vs `Yrow.st` | < 1e-10 | Different C interface, same math |
| `gsl.new(method="Tcol.mt")` vs `Tcol.mt` | 0 | Dispatcher, same function |

### Full verification script

```r
library(RidgeR)

dataPath <- file.path(system.file(package = "RidgeR"), "extdata/")
expr.diff <- read.table(paste0(dataPath, "Ly86-Fc_vs_Vehicle_logFC.txt"))

# Run all 7 variants with gsl RNG
res <- list(
  naive   = SecAct.inference.naive(expr.diff, rng_method = "gsl"),
  yrow.st = SecAct.inference.Yrow.st(expr.diff, rng_method = "gsl"),
  tcol.st = SecAct.inference.Tcol.st(expr.diff, rng_method = "gsl"),
  yrow.mt = SecAct.inference.Yrow.mt(expr.diff, rng_method = "gsl"),
  tcol.mt = SecAct.inference.Tcol.mt(expr.diff, rng_method = "gsl"),
  gsl.old = SecAct.inference.gsl.old(expr.diff, rng_method = "gsl"),
  gsl.new = SecAct.inference.gsl.new(expr.diff, rng_method = "gsl")
)

# Pairwise max absolute z-score differences
names_v <- names(res)
for (i in seq_along(names_v)) {
  for (j in seq(i + 1, length(names_v))) {
    if (j > length(names_v)) break
    d <- max(abs(res[[i]]$zscore - res[[j]]$zscore))
    cat(sprintf("%-10s vs %-10s: %.2e\n", names_v[i], names_v[j], d))
  }
}
```

---

## 6. Threading Model Differences

| Variant | Parallelization Axis | Atomic Writes | Race Conditions |
|---|---|---|---|
| `Yrow.st` | None (single-threaded) | No | None |
| `Yrow.mt` | Sample strips (`#pragma omp for`) | No (each thread owns its strip) | None |
| `Tcol.st` | None (single-threaded) | No | None |
| `Tcol.mt` | Permutations (`#pragma omp for`) | Yes (`#pragma omp atomic`) | None (atomics prevent races) |

**Why Yrow.mt doesn't need atomics**: Each thread processes a unique strip of
samples, so output elements are never shared between threads.

**Why Tcol.mt needs atomics**: Multiple threads process different permutations
but accumulate into the same output elements (pvalue, zscore/aver, se/aver_sq).
The `#pragma omp atomic` ensures correctness.

**Floating-point implications**: Atomic additions in `Tcol.mt` may occur in
different orders across runs, causing tiny floating-point differences (~1e-15
per addition). Over 1000 permutations, these can accumulate to ~1e-12.
Single-threaded variants are fully deterministic.
