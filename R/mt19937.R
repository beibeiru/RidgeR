# =========================================================
# Pure R implementation of GSL's MT19937 (Mersenne Twister)
#
# Produces output identical to gsl_rng_mt19937 with seed 0.
# Algorithm: Matsumoto & Nishimura (1998) with 2002 init.
# No external dependency — uses doubles for 32-bit unsigned
# arithmetic. Bitwise ops decompose into 16-bit halves
# (each in [0,65535], safe as R signed int32).
# =========================================================

# --- 32-bit unsigned helpers (operate on doubles in [0, 2^32-1]) ---

#' @keywords internal
.xor32 <- function(a, b) {
  bitwXor(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 +
    bitwXor(as.integer(a %% 65536), as.integer(b %% 65536))
}

#' @keywords internal
.and32 <- function(a, b) {
  bitwAnd(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 +
    bitwAnd(as.integer(a %% 65536), as.integer(b %% 65536))
}

#' @keywords internal
.or32 <- function(a, b) {
  bitwOr(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 +
    bitwOr(as.integer(a %% 65536), as.integer(b %% 65536))
}

#' @keywords internal
.shr32 <- function(a, k) a %/% (2^k)

#' @keywords internal
.shl32 <- function(a, k) (a * (2^k)) %% 4294967296

#' @keywords internal
.mul32 <- function(a, b) {
  a_lo <- a %% 65536;  a_hi <- a %/% 65536
  b_lo <- b %% 65536;  b_hi <- b %/% 65536
  (a_lo * b_lo + (a_lo * b_hi + a_hi * b_lo) * 65536) %% 4294967296
}

# --- MT19937 state ---

#' @keywords internal
.mt19937_new <- function(seed = 0) {
  N <- 624L
  if (seed == 0) seed <- 4357  # GSL default
  mt <- double(N)
  mt[1] <- seed %% 4294967296
  for (i in 2:N) {
    prev <- mt[i - 1L]
    mt[i] <- (.mul32(1812433253, .xor32(prev, .shr32(prev, 30))) + (i - 1L)) %% 4294967296
  }
  list(mt = mt, mti = N)  # mti = N forces twist on first get
}

#' @keywords internal
.mt19937_twist <- function(mt) {
  N <- 624L; M <- 397L
  UPPER <- 2147483648  # 0x80000000
  LOWER <- 2147483647  # 0x7FFFFFFF
  MAGIC <- 2567483615  # 0x9908b0df
  mt_old <- mt

  # Pre-compute y for kk = 1..N-1 (all use original state:
  # within each loop, mt[kk] and mt[kk+1] are never previously
  # overwritten — writes go to mt[kk], reads are mt[kk] and mt[kk+1],
  # and kk advances sequentially)
  y <- .or32(.and32(mt_old[1:(N - 1L)], UPPER), .and32(mt_old[2:N], LOWER))
  mag <- ifelse(y %% 2 == 1, MAGIC, 0)
  sy <- .xor32(.shr32(y, 1), mag)

  # Loop 1 (vectorized): kk = 1..227
  # XOR targets: mt_old[kk+M] (398..624) — all original
  idx1 <- 1L:(N - M)
  mt[idx1] <- .xor32(mt_old[idx1 + M], sy[idx1])

  # Loop 2a (vectorized): kk = 228..454
  # XOR targets: mt[kk+M-N] = mt[1..227] — all from loop 1
  idx2a <- (N - M + 1L):(2L * (N - M))
  mt[idx2a] <- .xor32(mt[idx2a + M - N], sy[idx2a])

  # Loop 2b (vectorized): kk = 455..623
  # XOR targets: mt[kk+M-N] = mt[228..396] — all from loop 2a
  idx2b <- (2L * (N - M) + 1L):(N - 1L)
  mt[idx2b] <- .xor32(mt[idx2b + M - N], sy[idx2b])

  # Final: kk = N — uses UPDATED mt[1] (from loop 1), not mt_old[1]
  yN <- .or32(.and32(mt_old[N], UPPER), .and32(mt[1L], LOWER))
  magN <- if (yN %% 2 == 1) MAGIC else 0
  mt[N] <- .xor32(mt[M], .xor32(.shr32(yN, 1), magN))
  mt
}

# --- Batch generation with vectorized tempering ---

#' @keywords internal
.mt19937_generate <- function(state, count) {
  values <- double(count)
  pos <- 1L
  while (pos <= count) {
    if (state$mti >= 624L) {
      state$mt <- .mt19937_twist(state$mt)
      state$mti <- 0L
    }
    avail <- 624L - state$mti
    take  <- min(avail, count - pos + 1L)
    raw   <- state$mt[(state$mti + 1L):(state$mti + take)]

    # Tempering (vectorized over batch)
    k <- raw
    k <- .xor32(k, .shr32(k, 11))
    k <- .xor32(k, .and32(.shl32(k, 7), 2636928640))   # 0x9d2c5680
    k <- .xor32(k, .and32(.shl32(k, 15), 4022730752))  # 0xefc60000
    k <- .xor32(k, .shr32(k, 18))

    values[pos:(pos + take - 1L)] <- k
    state$mti <- state$mti + take
    pos <- pos + take
  }
  list(values = values, state = state)
}

# --- Permutation table (matches C generate_perm_table with rng_method=GSL) ---

#' Pure R permutation table using GSL-compatible MT19937
#'
#' Generates the same permutation table as the C \code{generate_perm_table}
#' function with \code{rng_method = "gsl"}, using a pure R implementation
#' of the MT19937 algorithm. No GSL library required.
#'
#' @param n Number of samples (array length to shuffle).
#' @param nrand Number of permutations.
#' @return Integer matrix (nrand x n), 0-indexed (same format as C version).
#' @keywords internal
.gsl_mt19937_perm_table <- function(n, nrand) {
  MT_MAX <- 4294967295  # 2^32 - 1

  # Pre-generate all random values (vectorized RNG + tempering)
  total <- as.integer(n - 1L) * nrand
  state <- .mt19937_new(seed = 0)
  gen   <- .mt19937_generate(state, total)
  vals  <- gen$values

  # Pre-compute swap divisors (same for every permutation)
  ranges  <- n:2L                     # remaining elements per step
  divisors <- MT_MAX %/% ranges + 1   # (n-1) values

  # Pre-compute ALL swap targets at once (vectorized floor division)
  vals_mat <- matrix(vals, nrow = n - 1L, ncol = nrand)
  j_mat    <- floor(vals_mat / divisors) + (0:(n - 2L)) + 1L  # 1-indexed

  # Fisher-Yates shuffles with cumulative state
  arr   <- seq_len(n) - 1L  # 0-indexed: 0, 1, ..., n-1
  table <- matrix(0L, nrow = nrand, ncol = n)

  for (perm in seq_len(nrand)) {
    j_vec <- j_mat[, perm]
    for (i in seq_len(n - 1L)) {
      j <- j_vec[i]
      tmp <- arr[j]; arr[j] <- arr[i]; arr[i] <- tmp
    }
    table[perm, ] <- arr  # cumulative: array NOT reset between permutations
  }
  table
}
