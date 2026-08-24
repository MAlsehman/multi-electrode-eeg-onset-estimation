# functions.R
# Core functions for the multi-electrode EEG onset simulation.
# Contains the statistical, change-point and onset-estimation functions used by the univariate simulations.

nperm       <- 1000     # permutations for the MAX threshold
cpd_penalty <- "MBIC"   # change point penalty
sig_level   <- 0.05     # significance level for the FDR methods

# variance of each row (time point) of a time x trials matrix
rowVars <- function(X){
  n <- ncol(X); m <- rowMeans(X)
  (rowSums(X^2) - n * m^2) / (n - 1)
}

# Welch squared-t and p at every time point, one electrode
welch_t2p <- function(Ae, Be){
  nA <- ncol(Ae); nB <- ncol(Be)
  mA <- rowMeans(Ae); mB <- rowMeans(Be)
  vA <- rowVars(Ae);  vB <- rowVars(Be)
  se2  <- vA / nA + vB / nB
  tval <- (mB - mA) / sqrt(se2)
  df   <- se2^2 / ((vA / nA)^2 / (nA - 1) + (vB / nB)^2 / (nB - 1))
  list(t2 = tval^2, p = 2 * pt(-abs(tval), df))
}

# helpers that turn a time course into an onset in ms
cp_onset <- function(series, penalty = cpd_penalty){
  tryCatch({
    cp  <- cpt.meanvar(series, method = "BinSeg", Q = 2, penalty = penalty, minseglen = 2)
    pts <- cpts(cp)
    # earliest of the two change points, whatever order cpts() returns them in
    if(length(pts) == 0) NA_real_ else Xf[min(pts)]
  }, error = function(e) NA_real_)
}
first_cross <- function(series, thr){
  s <- which(series > thr); if(length(s) == 0) NA_real_ else Xf[s[1]]
}
first_sig <- function(pcol, level = sig_level){
  s <- which(pcol < level); if(length(s) == 0) NA_real_ else Xf[s[1]]
}

# summaries across electrodes
safe_min <- function(x){ x <- x[!is.na(x)]; if(length(x) == 0) NA_real_ else min(x) }
safe_quantile <- function(x, p){
  x <- x[!is.na(x)]; if(length(x) == 0) NA_real_ else quantile(x, p, names = FALSE)
}

# group onset per method, from the saved per-electrode onsets
group_from_perelec <- function(perelec, virtual){
  bind_rows(
    perelec |>
      group_by(method) |>
      summarise(min    = safe_min(onset),
                q10    = safe_quantile(onset, 0.10),
                q25    = safe_quantile(onset, 0.25),
                median = safe_quantile(onset, 0.50),
                .groups = "drop") |>
      pivot_longer(-method, names_to = "summary", values_to = "onset"),
    virtual |>
      transmute(method, summary = "single", onset)
  )
}

# one experiment: per-electrode onsets, plus one onset per virtual method
experiment_onsets <- function(){

  # trials: time x trial x electrode
  A <- array(0, dim = c(Nf, Nt, 9))   # noise only
  B <- array(0, dim = c(Nf, Nt, 9))   # signal + noise
  for(e in 1:9){
    A[, , e] <- replicate(Nt, eeg_noise(Nf, 500, 1, meanpower))
    B[, , e] <- replicate(Nt, signal_9[, e] + eeg_noise(Nf, 500, 1, meanpower))
  }

  # observed squared-t and p, per electrode and time
  T2obs <- matrix(NA_real_, Nf, 9)
  Pobs  <- matrix(NA_real_, Nf, 9)
  for(e in 1:9){
    tp <- welch_t2p(A[, , e], B[, , e])
    T2obs[, e] <- tp$t2
    Pobs[, e]  <- tp$p
  }
  vobs <- apply(T2obs, 1, max)        # virtual electrode: max squared-t across electrodes

  # permutation null for the MAX threshold: one shared shuffle of the 2*Nt trial
  # labels, applied to every electrode, so the structure across electrodes is kept
  pool <- lapply(1:9, function(e) cbind(A[, , e], B[, , e]))   # time x 2Nt, built once

  null_max <- numeric(nperm)
  for(b in 1:nperm){
    perm <- sample(2 * Nt)
    g1 <- perm[1:Nt]; g2 <- perm[(Nt + 1):(2 * Nt)]
    t2p <- matrix(NA_real_, Nf, 9)
    for(e in 1:9){
      t2p[, e] <- welch_t2p(pool[[e]][, g1, drop = FALSE],
                            pool[[e]][, g2, drop = FALSE])$t2
    }
    null_max[b] <- max(t2p)                              # max across the whole space
  }
  thr <- quantile(null_max, 0.95, names = FALSE)

  # FDR across the whole electrode-by-time search space
  bh <- matrix(p.adjust(as.vector(Pobs), "BH"), nrow = Nf)
  by <- matrix(p.adjust(as.vector(Pobs), "BY"), nrow = Nf)

  # change point everywhere, then the same onsets masked by MAX significance
  cpd_all  <- sapply(1:9, function(e) cp_onset(T2obs[, e]))
  elec_sig <- apply(T2obs, 2, max) > thr
  cpd_sig  <- if_else(elec_sig, cpd_all, NA_real_)

  # onsets at all nine electrodes, no selection on weight
  perelec <- tibble(
    electrode = rep(paste0("E", elec), 5),
    weight    = rep(w, 5),
    method    = rep(c("FDR BH", "FDR BY", "MAX across",
                      "CPD per electrode", "CPD if MAX significant"), each = 9),
    onset     = c(sapply(1:9, function(e) first_sig(bh[, e])),
                  sapply(1:9, function(e) first_sig(by[, e])),
                  sapply(1:9, function(e) first_cross(T2obs[, e], thr)),
                  cpd_all,
                  cpd_sig)
  )

  # virtual electrode methods: one onset each
  virtual <- tibble(
    method = c("MAX virtual", "CPD virtual"),
    onset  = c(first_cross(vobs, thr),
               cp_onset(vobs))
  )

  # raw onsets only: all summarising happens offline
  list(perelec = perelec, virtual = virtual)
}
