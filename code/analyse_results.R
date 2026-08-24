# analyse_results.R
# Post-generation analysis of the saved EEG onset simulation outputs.
# Loads the completed main, diagnostic, trial-size and optional multivariate results; it does not regenerate simulations.

library(tidyverse)

true_onset <- 160

# ============================================================================
# Helpers
# ============================================================================

resolve_file <- function(...){
  candidates <- c(...)
  hit <- candidates[file.exists(candidates)]

  if(length(hit) == 0){
    stop(
      "Required saved result file not found: ",
      paste(candidates, collapse = " | ")
    )
  }

  hit[[1]]
}

safe_min <- function(x){
  x <- x[!is.na(x)]
  if(length(x) == 0) NA_real_ else min(x)
}

safe_quantile <- function(x, p){
  x <- x[!is.na(x)]
  if(length(x) == 0) NA_real_ else quantile(x, p, names = FALSE)
}

load_result_pair <- function(path){
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  objects <- as.list(env)

  pick_by_columns <- function(required, forbidden = character()){
    hits <- objects[vapply(
      objects,
      function(x){
        is.data.frame(x) &&
          all(required %in% names(x)) &&
          !any(forbidden %in% names(x))
      },
      logical(1)
    )]

    if(length(hits) != 1){
      stop(
        "Could not identify exactly one saved table with columns: ",
        paste(required, collapse = ", "),
        " in ",
        path
      )
    }

    hits[[1]]
  }

  pe <- if("perelec_all" %in% names(objects)){
    objects[["perelec_all"]]
  } else {
    pick_by_columns(
      c("iter", "electrode", "weight", "method", "onset")
    )
  }

  vi <- if("virtual_all" %in% names(objects)){
    objects[["virtual_all"]]
  } else {
    pick_by_columns(
      c("iter", "method", "onset"),
      forbidden = c("electrode", "weight")
    )
  }

  list(
    perelec = as_tibble(pe),
    virtual = as_tibble(vi)
  )
}

make_group_onsets <- function(pe, vi, include_Nt = FALSE){
  grouping <- if(include_Nt){
    c("Nt", "iter", "method")
  } else {
    c("iter", "method")
  }

  pe_group <- pe |>
    group_by(across(all_of(grouping))) |>
    summarise(
      min = safe_min(onset),
      q10 = safe_quantile(onset, 0.10),
      q25 = safe_quantile(onset, 0.25),
      median = safe_quantile(onset, 0.50),
      .groups = "drop"
    ) |>
    pivot_longer(
      c(min, q10, q25, median),
      names_to = "summary",
      values_to = "onset"
    )

  if(include_Nt){
    vi_group <- vi |>
      transmute(
        Nt,
        iter,
        method,
        summary = "single",
        onset
      )
  } else {
    vi_group <- vi |>
      transmute(
        iter,
        method,
        summary = "single",
        onset
      )
  }

  bind_rows(pe_group, vi_group)
}

# ============================================================================
# Load saved univariate outputs
# ============================================================================

MAIN_FILE <- resolve_file("data/simres_n50.RData")

DIAG_FILE <- resolve_file("data/simres.RData")

SIZE_FILE <- resolve_file("data/simres_1000_sizes.RData")

main_saved <- load_result_pair(MAIN_FILE)
diag_saved <- load_result_pair(DIAG_FILE)
size_saved <- load_result_pair(SIZE_FILE)

main_pe <- main_saved$perelec
main_vi <- main_saved$virtual
diag_pe <- diag_saved$perelec
diag_vi <- diag_saved$virtual
size_pe <- size_saved$perelec
size_vi <- size_saved$virtual

# Older fixed-Nt files may not include Nt explicitly.
if(!"Nt" %in% names(main_pe)) main_pe$Nt <- 50
if(!"Nt" %in% names(main_vi)) main_vi$Nt <- 50
if(!"Nt" %in% names(diag_pe)) diag_pe$Nt <- 50
if(!"Nt" %in% names(diag_vi)) diag_vi$Nt <- 50

stopifnot(max(main_pe$iter, na.rm = TRUE) == 10000)
stopifnot(max(diag_pe$iter, na.rm = TRUE) == 1000)
stopifnot(
  identical(
    as.numeric(sort(unique(size_pe$Nt))),
    c(25, 50, 75, 100)
  )
)

# The diagnostic file is a separate 1,000-iteration Monte Carlo run at 50 trials.
# Structural checks above verify its size; it is not compared row-for-row with the main run.

group_main <- make_group_onsets(
  main_pe,
  main_vi,
  include_Nt = FALSE
)


group_size <- make_group_onsets(
  size_pe,
  size_vi,
  include_Nt = TRUE
)

# ============================================================================
# Main 10,000-run comparison
# ============================================================================

primary_spec <- tribble(
  ~method, ~summary, ~strategy,
  "FDR BH", "min", "Benjamini-Hochberg",
  "FDR BY", "min", "Benjamini-Yekutieli",
  "MAX across", "min", "Maximum statistic",
  "CPD per electrode", "median", "Unrestricted CPD, median",
  "CPD if MAX significant", "min", "MAX-filtered CPD, minimum",
  "CPD virtual", "single", "Pooled CPD"
)

full_summary <- group_main |>
  group_by(method, summary) |>
  summarise(
    median_onset = median(onset, na.rm = TRUE),
    bias = median_onset - true_onset,
    mae = mean(abs(onset - true_onset), na.rm = TRUE),
    sd = sd(onset, na.rm = TRUE),
    underestimation = 100 * mean(onset < true_onset, na.rm = TRUE),
    n_missing = sum(is.na(onset)),
    .groups = "drop"
  ) |>
  mutate(
    across(
      c(median_onset, bias, mae, sd, underestimation),
      ~round(.x, 1)
    )
  )

primary_long <- group_main |>
  inner_join(
    primary_spec,
    by = c("method", "summary")
  ) |>
  mutate(
    strategy = factor(
      strategy,
      levels = primary_spec$strategy
    )
  )

primary_summary <- primary_long |>
  group_by(strategy) |>
  summarise(
    median_onset = median(onset, na.rm = TRUE),
    bias = median_onset - true_onset,
    mae = mean(abs(onset - true_onset), na.rm = TRUE),
    sd = sd(onset, na.rm = TRUE),
    underestimation = 100 * mean(onset < true_onset, na.rm = TRUE),
    n_missing = sum(is.na(onset)),
    .groups = "drop"
  ) |>
  mutate(
    across(
      c(median_onset, bias, mae, sd, underestimation),
      ~round(.x, 1)
    )
  )

# ============================================================================
# CPD aggregation summaries
# ============================================================================

cpd_summaries <- group_main |>
  filter(
    method %in% c(
      "CPD per electrode",
      "CPD if MAX significant"
    ),
    summary %in% c("min", "q10", "q25", "median")
  ) |>
  group_by(method, summary) |>
  summarise(
    median_onset = median(onset, na.rm = TRUE),
    bias = median_onset - true_onset,
    mae = mean(abs(onset - true_onset), na.rm = TRUE),
    sd = sd(onset, na.rm = TRUE),
    underestimation = 100 * mean(onset < true_onset, na.rm = TRUE),
    n_missing = sum(is.na(onset)),
    .groups = "drop"
  ) |>
  mutate(
    across(
      c(median_onset, bias, mae, sd, underestimation),
      ~round(.x, 1)
    )
  )

# ============================================================================
# 1,000-run electrode-wise diagnostic at 50 trials
# ============================================================================

electrode_summary <- diag_pe |>
  group_by(electrode, weight, method) |>
  summarise(
    median_onset = median(onset, na.rm = TRUE),
    detection_rate = 100 * mean(!is.na(onset)),
    no_onset_rate = 100 * mean(is.na(onset)),
    median_bias = if(
      first(weight) > 0
    ){
      median(onset, na.rm = TRUE) - true_onset
    } else {
      NA_real_
    },
    mae = if(
      first(weight) > 0
    ){
      mean(abs(onset - true_onset), na.rm = TRUE)
    } else {
      NA_real_
    },
    sd = if(
      first(weight) > 0
    ){
      sd(onset, na.rm = TRUE)
    } else {
      NA_real_
    },
    underestimation = if(
      first(weight) > 0
    ){
      100 * mean(onset < true_onset, na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) |>
  mutate(
    across(
      c(
        median_onset,
        detection_rate,
        no_onset_rate,
        median_bias,
        mae,
        sd,
        underestimation
      ),
      ~round(.x, 1)
    )
  )

# ============================================================================
# Trial-number sensitivity
# ============================================================================

size_primary <- group_size |>
  inner_join(
    primary_spec,
    by = c("method", "summary")
  ) |>
  mutate(
    strategy = factor(
      strategy,
      levels = primary_spec$strategy
    )
  )

size_summary <- size_primary |>
  group_by(Nt, strategy) |>
  summarise(
    median_onset = median(onset, na.rm = TRUE),
    bias = median_onset - true_onset,
    mae = mean(abs(onset - true_onset), na.rm = TRUE),
    sd = sd(onset, na.rm = TRUE),
    underestimation = 100 * mean(onset < true_onset, na.rm = TRUE),
    no_onset = 100 * mean(is.na(onset)),
    .groups = "drop"
  ) |>
  mutate(
    across(
      c(
        median_onset,
        bias,
        mae,
        sd,
        underestimation,
        no_onset
      ),
      ~round(.x, 1)
    )
  )

# ============================================================================
# Optional multivariate results
# ============================================================================

mv_summary <- NULL

MV_FILE <- "data/simres_n50_mult.RData"
MV_FILE <- MV_FILE[file.exists(MV_FILE)]

if(length(MV_FILE) > 0){
  mv_env <- new.env(parent = emptyenv())
  load(MV_FILE[[1]], envir = mv_env)

  if(exists("method_summary_mv", envir = mv_env, inherits = FALSE)){
    mv_summary <- as_tibble(
      get(
        "method_summary_mv",
        envir = mv_env,
        inherits = FALSE
      )
    )
  } else if(exists("sim_mv", envir = mv_env, inherits = FALSE)){
    sim_mv <- as_tibble(
      get(
        "sim_mv",
        envir = mv_env,
        inherits = FALSE
      )
    )

    mv_summary <- sim_mv |>
      group_by(method) |>
      summarise(
        median_onset = median(onset, na.rm = TRUE),
        bias = median_onset - true_onset,
        mae = mean(abs(onset - true_onset), na.rm = TRUE),
        sd_onset = sd(onset, na.rm = TRUE),
        underestimation = 100 * mean(onset < true_onset, na.rm = TRUE),
        n_missing = sum(is.na(onset)),
        .groups = "drop"
      ) |>
      mutate(
        across(
          c(
            median_onset,
            bias,
            mae,
            sd_onset,
            underestimation
          ),
          ~round(.x, 1)
        )
      )
  }
}

# ============================================================================
# Save derived analysis tables
# ============================================================================

dir.create(
  "analysis_outputs",
  showWarnings = FALSE
)

write_csv(
  full_summary,
  "analysis_outputs/full_summary.csv"
)

write_csv(
  primary_summary,
  "analysis_outputs/primary_summary.csv"
)

write_csv(
  cpd_summaries,
  "analysis_outputs/cpd_summaries.csv"
)

write_csv(
  electrode_summary,
  "analysis_outputs/electrode_summary_1000.csv"
)

write_csv(
  size_summary,
  "analysis_outputs/trial_size_summary.csv"
)

if(!is.null(mv_summary)){
  write_csv(
    mv_summary,
    "analysis_outputs/multivariate_summary.csv"
  )
}

save(
  full_summary,
  primary_summary,
  cpd_summaries,
  electrode_summary,
  size_summary,
  mv_summary,
  file = "analysis_outputs/derived_results.RData"
)

message(
  "Analysis complete. Derived tables saved in analysis_outputs/."
)