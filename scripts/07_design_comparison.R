# =============================================================================
# 07_design_comparison.R
# Does the efficient design (script 06) actually estimate better than the
# random design (script 01)? Test it empirically, using the SAME 300 simulated
# respondents (same true utilities) and the SAME choice process.
#
#   A. Pooled clogit, N_REP replications per design:
#      compare standard errors and the spread of estimates across replications.
#      (Pooled estimates are biased under heterogeneity -- see README Finding 1 --
#      so precision, not bias, is the fair comparison here.)
#   B. Hierarchical Bayes, N_HB_REP replications per design:
#      compare individual-level recovery (correlation with true utilities).
#
# Each replication redraws the random design / block assignment and the Gumbel
# choice noise, so the result is not an artefact of one lucky seed.
#
# Outputs
#   data/simulated_choices_efficient.csv            one efficient-design dataset
#   outputs/efficient_design/clogit_precision.csv
#   outputs/efficient_design/hb_recovery_comparison.csv
# =============================================================================

source("scripts/helpers_design.R")

N_REP    <- 30     # clogit replications per design
N_HB_REP <- 3      # HB replications per design (each ~ tens of seconds)
R_HB     <- 8000   # MCMC draws (same as script 04)

true_indiv <- read.csv("data/respondent_true_utilities.csv")
eff_design <- read.csv("data/efficient_design.csv")
n_resp     <- nrow(true_indiv)
n_blocks   <- length(unique(eff_design$block))

# Attach the efficient design to respondents (each sees ONE block of 10 tasks)
efficient_for_respondents <- function() {
  assign_tbl <- data.frame(resp_id = 1:n_resp,
                           block = sample(rep(1:n_blocks, length.out = n_resp)))
  d <- merge(assign_tbl, eff_design, by = "block")
  d[, c("resp_id", "task_id", "alt_id", "positioning", "ingredient", "format", "price")]
}

make_dataset <- function(kind, seed) {
  set.seed(seed)
  design <- if (kind == "random") random_design(n_resp, 10) else efficient_for_respondents()
  simulate_choice_data(design, true_indiv)
}

# A. Pooled clogit precision
cat("A. Pooled clogit,", N_REP, "replications per design...\n")
fits <- list()
for (rep in 1:N_REP) {
  for (kind in c("random", "efficient")) {
    ds <- make_dataset(kind, seed = 1000 + rep)
    f  <- fit_clogit(add_dummies(ds)); f$design <- kind; f$rep <- rep
    fits[[length(fits) + 1]] <- f
  }
}
fits <- do.call(rbind, fits)

precision <- fits %>%
  group_by(parameter, design) %>%
  summarise(mean_se = mean(se), sd_of_estimates = sd(estimate), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = design, values_from = c(mean_se, sd_of_estimates)) %>%
  mutate(se_ratio = mean_se_efficient / mean_se_random,
         spread_ratio = sd_of_estimates_efficient / sd_of_estimates_random) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  arrange(match(parameter, param_names))
write.csv(precision, "outputs/efficient_design/clogit_precision.csv", row.names = FALSE)
cat("\nPooled clogit: mean standard error by design (lower = more precise)\n")
print(as.data.frame(precision), row.names = FALSE)
cat(sprintf("\nAverage SE ratio (efficient / random): %.3f  ->  ~%.0f%% smaller standard errors\n",
            mean(precision$se_ratio), 100 * (1 - mean(precision$se_ratio))))

# Save one efficient-design dataset for reference / re-use
write.csv(make_dataset("efficient", seed = 1001), "data/simulated_choices_efficient.csv", row.names = FALSE)

# B. HB individual-level recovery
cat("\nB. Hierarchical Bayes,", N_HB_REP, "replications per design (R =", R_HB, ")...\n")
rec <- list()
for (rep in 1:N_HB_REP) {
  for (kind in c("random", "efficient")) {
    t0 <- Sys.time()
    ds <- make_dataset(kind, seed = 2000 + rep)
    ib <- fit_hb(ds, R = R_HB, seed = 42 + rep)
    r  <- hb_recovery(ib, true_indiv); r$design <- kind; r$rep <- rep
    rec[[length(rec) + 1]] <- r
    cat(sprintf("  rep %d %-9s done in %.0fs\n", rep, kind, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
rec <- do.call(rbind, rec)

hb_cmp <- rec %>%
  group_by(parameter, design) %>%
  summarise(correlation = mean(correlation), mean_abs_error = mean(mean_abs_error), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = design, values_from = c(correlation, mean_abs_error)) %>%
  mutate(corr_change = correlation_efficient - correlation_random) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  arrange(match(parameter, param_names))
write.csv(hb_cmp, "outputs/efficient_design/hb_recovery_comparison.csv", row.names = FALSE)
cat("\nHB individual-level recovery, mean over replications\n")
print(as.data.frame(hb_cmp), row.names = FALSE)
