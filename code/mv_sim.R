# mv_sim.R
# Multivariate production EEG onset simulation using Hotelling T-squared across nine electrodes.
# Runs 10,000 iterations at 50 trials per condition with 1,000 MAX permutations.

library(tidyverse)
library(changepoint)

source("code/eeg_noise.R")
source("code/erp_template.R")
meanpower <- scan("code/meanpower.txt")

set.seed(26)

Nt   <- 50
w    <- c(0, 0, 0.3, 0.75, 1, 0.75, 0.3, 0, 0)
elec <- 1:9
Nf   <- length(Xf)

signal_9 <- outer(temp2, w)
colnames(signal_9) <- paste0("E", elec)

nperm       <- 1000
cpd_penalty <- "MBIC"
sig_level   <- 0.05


hotelling_t2 <- function(A, B){
  nf <- dim(A)[1]
  nA <- dim(A)[2]
  nB <- dim(B)[2]
  p  <- dim(A)[3]

  if(p >= nA + nB - 1){
    stop(
      "too few trials for ",
      p,
      " electrodes: the pooled covariance is singular"
    )
  }

  const <- (nA * nB) / (nA + nB)
  t2 <- numeric(nf)

  for(k in seq_len(nf)){
    Xa <- matrix(A[k, , ], nrow = nA)
    Xb <- matrix(B[k, , ], nrow = nB)

    d <- colMeans(Xb) - colMeans(Xa)

    S <- (
      (nA - 1) * cov(Xa) +
        (nB - 1) * cov(Xb)
    ) / (nA + nB - 2)

    t2[k] <- tryCatch(
      const * as.numeric(crossprod(d, solve(S, d))),
      error = function(e) NA_real_
    )
  }

  t2
}

hotelling_p <- function(t2, nA, nB, p){
  f <- t2 *
    (nA + nB - p - 1) /
    ((nA + nB - 2) * p)

  pf(
    f,
    df1 = p,
    df2 = nA + nB - p - 1,
    lower.tail = FALSE
  )
}

make_experiment <- function(){
  A <- array(0, dim = c(Nf, Nt, 9))
  B <- array(0, dim = c(Nf, Nt, 9))

  for(e in 1:9){
    A[, , e] <- replicate(
      Nt,
      eeg_noise(Nf, 500, 1, meanpower)
    )

    B[, , e] <- replicate(
      Nt,
      signal_9[, e] + eeg_noise(Nf, 500, 1, meanpower)
    )
  }

  list(A = A, B = B)
}

cp_onset <- function(series, penalty = cpd_penalty){
  tryCatch({
    cp <- cpt.meanvar(
      series,
      method = "BinSeg",
      Q = 2,
      penalty = penalty,
      minseglen = 2
    )

    pts <- cpts(cp)

    if(length(pts) == 0){
      NA_real_
    } else {
      Xf[min(pts)]
    }
  }, error = function(e) NA_real_)
}

first_cross <- function(series, thr){
  s <- which(series > thr)
  if(length(s) == 0) NA_real_ else Xf[s[1]]
}

first_sig <- function(pcol, level = sig_level){
  s <- which(pcol < level)
  if(length(s) == 0) NA_real_ else Xf[s[1]]
}

mv_onsets <- function(A, B, nperm = 1000){
  nf <- dim(A)[1]
  nA <- dim(A)[2]
  nB <- dim(B)[2]
  p  <- dim(A)[3]

  t2 <- hotelling_t2(A, B)
  pv <- hotelling_p(t2, nA, nB, p)

  pool <- array(
    0,
    dim = c(nf, nA + nB, p)
  )

  pool[, 1:nA, ] <- A
  pool[, (nA + 1):(nA + nB), ] <- B

  null_max <- numeric(nperm)

  for(b in seq_len(nperm)){
    perm <- sample(nA + nB)

    null_max[b] <- max(
      hotelling_t2(
        pool[, perm[1:nA], , drop = FALSE],
        pool[, perm[(nA + 1):(nA + nB)], , drop = FALSE]
      ),
      na.rm = TRUE
    )
  }

  thr <- quantile(
    null_max,
    0.95,
    names = FALSE
  )

  tibble(
    method = c(
      "Hotelling FDR BH",
      "Hotelling FDR BY",
      "Hotelling MAX",
      "Hotelling CPD"
    ),
    onset = c(
      first_sig(p.adjust(pv, "BH")),
      first_sig(p.adjust(pv, "BY")),
      first_cross(t2, thr),
      cp_onset(t2)
    )
  )
}


dir.create("data", showWarnings = FALSE)

# ============================================================================
# Monte Carlo simulation: 10,000 iterations, 50 trials per condition
# ============================================================================

niter <- 10000
notify <- 500

mv_methods <- c(
  "Hotelling FDR BH",
  "Hotelling FDR BY",
  "Hotelling MAX",
  "Hotelling CPD"
)

sim_mv <- bind_rows(
  lapply(seq_len(niter), function(it){

    if(it %% notify == 0){
      message("iteration ", it, " / ", niter)
    }

    ex_it <- make_experiment()

    mv_onsets(
      ex_it$A,
      ex_it$B,
      nperm = nperm
    ) |>
      mutate(iter = it)
  })
) |>
  mutate(
    method = factor(
      method,
      levels = mv_methods
    )
  )

method_summary_mv <- sim_mv |>
  group_by(method) |>
  summarise(
    median_onset = median(onset, na.rm = TRUE),
    sd_onset = round(sd(onset, na.rm = TRUE), 1),
    mae = round(
      mean(abs(onset - true_onset), na.rm = TRUE),
      1
    ),
    underestimation = round(
      100 * mean(onset < true_onset, na.rm = TRUE),
      1
    ),
    n_missing = sum(is.na(onset)),
    .groups = "drop"
  ) |>
  mutate(
    bias = round(
      median_onset - true_onset,
      1
    ),
    .after = median_onset
  )

save(
  mv_methods,
  sim_mv,
  method_summary_mv,
  file = "data/simres_n50_mult.RData"
)

message("Saved data/simres_n50_mult.RData")