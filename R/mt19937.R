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
  mt_old <- mt  # snapshot for loop 1 reads

  # Loop 1 (vectorized): kk = 1..227 — all reads from old state
  idx1 <- 1L:(N - M)
  y1 <- .or32(.and32(mt_old[idx1], UPPER), .and32(mt_old[idx1 + 1L], LOWER))
  mag1 <- ifelse(y1 %% 2 == 1, MAGIC, 0)
  mt[idx1] <- .xor32(mt_old[idx1 + M], .xor32(.shr32(y1, 1), mag1))

  # Loop 2 (scalar): kk = 228..623 — reads loop-1-updated values
  for (kk in (N - M + 1L):(N - 1L)) {
    y <- .or32(.and32(mt[kk], UPPER), .and32(mt[kk + 1L], LOWER))
    mag <- if (y %% 2 == 1) MAGIC else 0
    mt[kk] <- .xor32(mt[kk + M - N], .xor32(.shr32(y, 1), mag))
  }

  # Final: kk = N (wrap)
  y <- .or32(.and32(mt[N], UPPER), .and32(mt[1L], LOWER))
  mag <- if (y %% 2 == 1) MAGIC else 0
  mt[N] <- .xor32(mt[M], .xor32(.shr32(y, 1), mag))
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

  # Fisher-Yates shuffles with cumulative state
  arr   <- seq_len(n) - 1L  # 0-indexed: 0, 1, ..., n-1
  table <- matrix(0L, nrow = nrand, ncol = n)
  vi    <- 1L  # index into pre-generated values

  for (perm in seq_len(nrand)) {
    for (i in seq_len(n - 1L)) {
      range <- n - i + 1L  # remaining elements (0-indexed count)
      j <- (i - 1L) + floor(vals[vi] / (MT_MAX %/% range + 1))
      vi <- vi + 1L
      j <- j + 1L  # to 1-indexed
      tmp <- arr[j]; arr[j] <- arr[i]; arr[i] <- tmp
    }
    table[perm, ] <- arr  # cumulative: array NOT reset between permutations
  }
  table
}
