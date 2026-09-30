# =============================================================================
# helpers_design.R
# Shared functions for scripts 06-08: efficient-design search, choice
# simulation, and HB fitting. Sourced, not run directly.
#
# Coding matches scripts 01/02/04 exactly:
#   pos_derm, pos_luxury, ing_vitc, ing_pept, fmt_pump  = 0/1 dummies
#   price_eur = price in EUR (continuous), none_flag = 1 for the "None" option
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ---- Attribute definitions (same levels as script 01) ------------------------
pos_levels   <- c("Mass-Prestige", "Dermatologist-Recommended", "Luxury Heritage")
ing_levels   <- c("Retinol", "Vitamin C", "Peptide Complex")
fmt_levels   <- c("Standard Bottle", "Airless Pump")
price_levels <- c(45, 65, 85, 110)
n_lv         <- c(3, 3, 2, 4)           # levels per attribute (pos, ing, fmt, price)
param_names  <- c("pos_derm", "pos_luxury", "ing_vitc", "ing_pept",
                  "fmt_pump", "price_eur", "none_flag")
K <- length(param_names)

# ---- Prohibitions (implausible product concepts) -----------------------------
# A profile is a length-4 vector of level INDICES: (positioning, ingredient,
# format, price). Prohibited combinations are never shown to respondents:
#   1. Luxury Heritage positioning at the lowest price (EUR 45): not credible.
#   2. Mass-Prestige positioning at the highest price (EUR 110): not credible.
# Edit here to change the constraints; everything downstream adapts.
profile_allowed <- function(lv) {
  !(lv[1] == 3 && lv[4] == 1) && !(lv[1] == 1 && lv[4] == 4)
}

# ---- Design matrix and information matrix ------------------------------------
profile_row <- function(lv) {
  c(as.numeric(lv[1] == 2), as.numeric(lv[1] == 3),
    as.numeric(lv[2] == 2), as.numeric(lv[2] == 3),
    as.numeric(lv[3] == 2), price_levels[lv[4]], 0)
}

# task_lv: 3 x 4 matrix of level indices (one row per product alternative).
# A fourth row (the "None" option) is appended.
task_X <- function(task_lv) {
  rbind(t(apply(task_lv, 1, profile_row)), c(0, 0, 0, 0, 0, 0, 1))
}

# Fisher information of one MNL choice task, evaluated at prior beta.
task_info <- function(task_lv, beta) {
  X <- task_X(task_lv)
  u <- as.vector(X %*% beta)
  p <- exp(u - max(u)); p <- p / sum(p)
  crossprod(X, (diag(p) - outer(p, p)) %*% X)
}

log_det <- function(M) {
  d <- determinant(M, logarithm = TRUE)
  if (d$sign <= 0 || !is.finite(d$modulus)) -Inf else as.numeric(d$modulus)
}

# D-error: (det of the asymptotic covariance matrix)^(1/K). Lower = better.
d_error <- function(info) exp(-log_det(info) / K)

# ---- Random valid task -------------------------------------------------------
random_task <- function() {
  repeat {
    m <- cbind(sample(3, 3, TRUE), sample(3, 3, TRUE),
               sample(2, 3, TRUE), sample(4, 3, TRUE))
    if (all(apply(m, 1, profile_allowed)) && !any(duplicated(m))) return(m)
  }
}

# ---- Level-frequency helper (for the minimum-balance constraint) -------------
level_counts <- function(D) {
  all_rows <- do.call(rbind, D)
  lapply(1:4, function(a) tabulate(all_rows[, a], nbins = n_lv[a]))
}

# ---- Coordinate-exchange search for a Dp-efficient design --------------------
# Repeatedly tries changing one attribute level of one alternative at a time,
# keeping any change that raises log det(information) while respecting:
#   * prohibitions,
#   * no duplicate alternatives within a task,
#   * a minimum level frequency (each level >= min_share of its equal share),
#     so the search cannot collapse onto extreme price levels only.
ce_search <- function(n_tasks, beta, n_starts = 10, max_pass = 30,
                      min_share = 0.7, verbose = TRUE) {
  best <- NULL
  n_alts_total <- n_tasks * 3
  min_count <- floor(min_share * n_alts_total / n_lv)
  for (s in seq_len(n_starts)) {
    D     <- replicate(n_tasks, random_task(), simplify = FALSE)
    infos <- lapply(D, task_info, beta = beta)
    tot   <- Reduce(`+`, infos)
    cur   <- log_det(tot)
    counts <- level_counts(D)
    for (pass in seq_len(max_pass)) {
      changed <- FALSE
      for (t in seq_len(n_tasks)) for (j in 1:3) for (a in 1:4) {
        old_l <- D[[t]][j, a]
        best_l <- old_l; best_ld <- cur; best_info <- NULL
        for (l in setdiff(seq_len(n_lv[a]), old_l)) {
          if (counts[[a]][old_l] - 1 < min_count[a]) next
          Dt <- D[[t]]; Dt[j, a] <- l
          if (!profile_allowed(Dt[j, ]) || any(duplicated(Dt))) next
          info_new <- task_info(Dt, beta)
          ld <- log_det(tot - infos[[t]] + info_new)
          if (ld > best_ld + 1e-9) { best_ld <- ld; best_l <- l; best_info <- info_new }
        }
        if (best_l != old_l) {
          counts[[a]][old_l] <- counts[[a]][old_l] - 1
          counts[[a]][best_l] <- counts[[a]][best_l] + 1
          D[[t]][j, a] <- best_l
          tot <- tot - infos[[t]] + best_info
          infos[[t]] <- best_info
          cur <- best_ld
          changed <- TRUE
        }
      }
      if (!changed) break
    }
    if (verbose) cat(sprintf("  start %2d: log det = %.3f (%d passes)\n", s, cur, pass))
    if (is.null(best) || cur > best$logdet) best <- list(D = D, logdet = cur)
  }
  best
}

# ---- Choice simulation (same data-generating process as script 01) -----------
rgumbel <- function(n) -log(-log(runif(n)))

# design: data.frame with resp_id, task_id, alt_id (1-3), positioning,
#         ingredient, format, price (one row per respondent x task x alternative)
# indiv_utilities: respondent-level true utilities (data/respondent_true_utilities.csv)
simulate_choice_data <- function(design, indiv_utilities) {
  d <- merge(design, indiv_utilities, by = "resp_id")
  d$utility <- ifelse(d$positioning == "Dermatologist-Recommended", d$pos_derm,
               ifelse(d$positioning == "Luxury Heritage", d$pos_luxury, 0)) +
               ifelse(d$ingredient == "Vitamin C", d$ing_vitc,
               ifelse(d$ingredient == "Peptide Complex", d$ing_pept, 0)) +
               ifelse(d$format == "Airless Pump", d$fmt_pump, 0) +
               d$price_coef * d$price

  none_rows <- unique(d[, c("resp_id", "task_id", "segment", "none_asc")])
  none_rows$alt_id <- 0
  none_rows$positioning <- NA; none_rows$ingredient <- NA
  none_rows$format <- NA;      none_rows$price <- NA
  none_rows$utility <- none_rows$none_asc

  cols <- c("resp_id", "task_id", "alt_id", "segment", "positioning",
            "ingredient", "format", "price", "utility")
  full <- rbind(d[, cols], none_rows[, cols])
  full <- full[order(full$resp_id, full$task_id, full$alt_id), ]
  full$is_none <- as.integer(full$alt_id == 0)
  full$rand_utility <- full$utility + rgumbel(nrow(full))
  full$chosen <- ave(full$rand_utility, paste(full$resp_id, full$task_id),
                     FUN = function(x) as.integer(x == max(x)))
  full[, c("resp_id", "task_id", "alt_id", "segment", "positioning",
           "ingredient", "format", "price", "is_none", "chosen")]
}

# Random design exactly as in script 01 (levels drawn independently per alternative)
random_design <- function(n_resp, n_tasks) {
  d <- expand.grid(resp_id = 1:n_resp, task_id = 1:n_tasks, alt_id = 1:3)
  d$positioning <- sample(pos_levels,   nrow(d), replace = TRUE)
  d$ingredient  <- sample(ing_levels,   nrow(d), replace = TRUE)
  d$format      <- sample(fmt_levels,   nrow(d), replace = TRUE)
  d$price       <- sample(price_levels, nrow(d), replace = TRUE)
  d
}

# ---- Estimation helpers -------------------------------------------------------
add_dummies <- function(choices) {
  choices %>%
    mutate(
      positioning = ifelse(is_none == 1, "Mass-Prestige", positioning),
      ingredient  = ifelse(is_none == 1, "Retinol", ingredient),
      format      = ifelse(is_none == 1, "Standard Bottle", format),
      price       = ifelse(is_none == 1, 0, price),
      pos_derm   = as.integer(positioning == "Dermatologist-Recommended" & is_none == 0),
      pos_luxury = as.integer(positioning == "Luxury Heritage" & is_none == 0),
      ing_vitc   = as.integer(ingredient == "Vitamin C" & is_none == 0),
      ing_pept   = as.integer(ingredient == "Peptide Complex" & is_none == 0),
      fmt_pump   = as.integer(format == "Airless Pump" & is_none == 0),
      price_eur  = price,
      none_flag  = is_none,
      stratum    = paste(resp_id, task_id)
    ) %>%
    arrange(resp_id, task_id, alt_id)
}

fit_clogit <- function(choices) {
  library(survival)
  m <- clogit(chosen ~ pos_derm + pos_luxury + ing_vitc + ing_pept + fmt_pump +
                price_eur + none_flag + strata(stratum), data = choices)
  s <- summary(m)$coefficients
  data.frame(parameter = rownames(s), estimate = s[, "coef"], se = s[, "se(coef)"],
             row.names = NULL)
}

# HB-MNL via bayesm, same specification as script 04 (ncomp = 1)
fit_hb <- function(choices, R = 8000, keep = 5, seed = 42) {
  library(bayesm)
  choices <- add_dummies(choices)
  resp_ids <- sort(unique(choices$resp_id))
  lgtdata <- lapply(resp_ids, function(rid) {
    d_r <- choices[choices$resp_id == rid, ]
    tasks <- sort(unique(d_r$task_id))
    y <- vapply(tasks, function(t) which(d_r$chosen[d_r$task_id == t] == 1), integer(1))
    list(y = y, X = as.matrix(d_r[, param_names]))
  })
  set.seed(seed)
  # bayesm prints a long progress log; suppress it (assignment still happens here)
  invisible(utils::capture.output(
    out <- rhierMnlRwMixture(Data = list(lgtdata = lgtdata, p = 4),
                             Prior = list(ncomp = 1), Mcmc = list(R = R, keep = keep))))
  n_kept <- dim(out$betadraw)[3]
  keep_idx <- (floor(n_kept / 2) + 1):n_kept          # discard first half as burn-in
  betas <- apply(out$betadraw[, , keep_idx], c(1, 2), mean)
  colnames(betas) <- param_names
  data.frame(resp_id = resp_ids, betas)
}

# Correlation / MAE between HB individual estimates and true utilities
hb_recovery <- function(indiv_betas, true_indiv) {
  names(true_indiv)[names(true_indiv) == "price_coef"] <- "price_eur"
  names(true_indiv)[names(true_indiv) == "none_asc"]   <- "none_flag"
  cmp <- merge(indiv_betas, true_indiv, by = "resp_id", suffixes = c("_hb", "_true"))
  data.frame(
    parameter   = param_names,
    correlation = sapply(param_names, function(v) cor(cmp[[paste0(v, "_hb")]], cmp[[paste0(v, "_true")]])),
    mean_abs_error = sapply(param_names, function(v) mean(abs(cmp[[paste0(v, "_hb")]] - cmp[[paste0(v, "_true")]])))
  )
}
