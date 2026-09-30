# =============================================================================
# 06_efficient_design.R
# Replace the random choice-set assignment of script 01 with a D-efficient,
# blocked design that respects prohibitions, and quantify the gain.
#
# What this script does
#   1. Sets planning priors (what a researcher would assume from a pilot).
#   2. Searches for a Dp-efficient design (coordinate exchange) of 60 choice
#      tasks x 3 alternatives + None, subject to:
#        - prohibitions (implausible concepts never shown; see helpers_design.R)
#        - no duplicate alternatives within a task
#        - minimum level balance (every level appears often enough)
#   3. Splits the 60 tasks into 6 blocks of 10 (each respondent sees ONE block,
#      so the survey stays at 10 tasks per respondent, as in script 01).
#   4. Compares the new design against the random design actually used in
#      script 01: predicted standard errors, D-error, level balance, and the
#      number of implausible concepts shown.
#
# Outputs (all under outputs/efficient_design/ and data/)
#   data/efficient_design.csv                    the design, ready to field
#   outputs/efficient_design/design_diagnostics.csv
#   outputs/efficient_design/level_balance.csv
# =============================================================================

source("scripts/helpers_design.R")
dir.create("outputs/efficient_design", showWarnings = FALSE, recursive = TRUE)
set.seed(42)

n_resp <- 300
n_blocks <- 6
tasks_per_block <- 10
n_tasks_total <- n_blocks * tasks_per_block

# 1. Planning priors
# Assumed pilot-style guesses:
# modest positive part-worths, mildly negative price. The None constant is
# calibrated so the predicted None share matches the None share observed in the
# existing random-design data, which plays the role of a pilot wave.
pilot <- read.csv("data/simulated_choices.csv")
pilot_none_share <- mean(pilot$chosen[pilot$is_none == 1])

prior <- c(pos_derm = 0.5, pos_luxury = 0.5, ing_vitc = 0.25, ing_pept = 0.25,
           fmt_pump = 0.25, price_eur = -0.02, none_flag = 0)

set.seed(42)
calib_tasks <- replicate(2000, random_task(), simplify = FALSE)
none_share_at <- function(none_b) {
  b <- prior; b["none_flag"] <- none_b
  mean(sapply(calib_tasks, function(tk) {
    u <- as.vector(task_X(tk) %*% b); p <- exp(u - max(u)); (p / sum(p))[4]
  }))
}
prior["none_flag"] <- uniroot(function(z) none_share_at(z) - pilot_none_share,
                              interval = c(-8, 4))$root
cat("Pilot None share:", round(pilot_none_share, 3),
    "| calibrated None prior:", round(prior["none_flag"], 3), "\n\n")

# 2. Search for the efficient design
set.seed(2026)
cat("Searching for a Dp-efficient design (", n_tasks_total, "tasks )...\n")
res <- ce_search(n_tasks_total, beta = prior, n_starts = 10)
D <- res$D

# 3. Blocking: 6 blocks of 10 tasks
# Random partitions; keep the one whose WORST block is most informative.
infos <- lapply(D, task_info, beta = prior)
best_part <- NULL; best_min <- -Inf
for (it in 1:3000) {
  perm <- sample(n_tasks_total)
  blk  <- rep(1:n_blocks, each = tasks_per_block)[order(perm)]
  worst <- min(sapply(1:n_blocks, function(b) log_det(Reduce(`+`, infos[blk == b]))))
  if (worst > best_min) { best_min <- worst; best_part <- blk }
}
block_of_task <- best_part

# 4. Write the fieldable design
rows <- do.call(rbind, lapply(seq_len(n_tasks_total), function(t) {
  m <- D[[t]]
  data.frame(design_task = t, block = block_of_task[t], alt_id = 1:3,
             positioning = pos_levels[m[, 1]], ingredient = ing_levels[m[, 2]],
             format = fmt_levels[m[, 3]], price = price_levels[m[, 4]])
}))
rows <- rows[order(rows$block, rows$design_task, rows$alt_id), ]
rows$task_id <- as.integer(factor(rows$design_task)) # placeholder, reset below
rows$task_id <- ave(rows$design_task, rows$block, FUN = function(x) match(x, sort(unique(x))))
rows <- rows[, c("block", "task_id", "design_task", "alt_id", "positioning",
                 "ingredient", "format", "price")]
write.csv(rows, "data/efficient_design.csv", row.names = FALSE)

# 5. Compare against the random design used in script 01
prod_rows <- pilot[pilot$is_none == 0, ]
prod_rows <- prod_rows[order(prod_rows$resp_id, prod_rows$task_id, prod_rows$alt_id), ]
idx <- cbind(match(prod_rows$positioning, pos_levels), match(prod_rows$ingredient, ing_levels),
             match(prod_rows$format, fmt_levels), match(prod_rows$price, price_levels))
grp <- factor(paste(prod_rows$resp_id, prod_rows$task_id),
              levels = unique(paste(prod_rows$resp_id, prod_rows$task_id)))
random_tasks <- split(seq_len(nrow(prod_rows)), grp)
info_random <- Reduce(`+`, lapply(random_tasks, function(i) task_info(idx[i, , drop = FALSE], prior)))

resp_per_block <- n_resp / n_blocks
info_eff <- Reduce(`+`, lapply(1:n_blocks, function(b)
  resp_per_block * Reduce(`+`, infos[block_of_task == b])))

se_random <- sqrt(diag(solve(info_random)))
se_eff    <- sqrt(diag(solve(info_eff)))
diag_tbl <- data.frame(
  parameter = param_names,
  predicted_se_random = round(se_random, 4),
  predicted_se_efficient = round(se_eff, 4),
  se_ratio_efficient_vs_random = round(se_eff / se_random, 3)
)
write.csv(diag_tbl, "outputs/efficient_design/design_diagnostics.csv", row.names = FALSE)

d_rand <- d_error(info_random); d_eff <- d_error(info_eff)
rel_eff <- d_rand / d_eff       # >1: efficient design is more informative
cat("\n=== Design comparison (same 300 respondents x 10 tasks) ===\n")
print(diag_tbl, row.names = FALSE)
cat(sprintf("\nD-error (lower is better):  random = %.5f | efficient = %.5f\n", d_rand, d_eff))
cat(sprintf("Relative D-efficiency = %.2f  ->  the random design would need ~%.0f%% more respondents\n",
            rel_eff, 100 * (rel_eff - 1)))
cat("to match the efficient design's overall precision.\n")

# Implausible concepts shown
n_bad_random <- sum(!apply(idx, 1, profile_allowed))
eff_idx <- do.call(rbind, D)
n_bad_eff <- sum(!apply(eff_idx, 1, profile_allowed))
cat(sprintf("\nProhibited concepts shown: random design = %d of %d alternatives (%.1f%%) | efficient = %d\n",
            n_bad_random, nrow(idx), 100 * n_bad_random / nrow(idx), n_bad_eff))

# Level balance
bal <- do.call(rbind, lapply(1:4, function(a) {
  lv_names <- list(pos_levels, ing_levels, fmt_levels, price_levels)[[a]]
  data.frame(attribute = c("Positioning", "Ingredient", "Format", "Price")[a],
             level = lv_names,
             share_random = as.numeric(tabulate(idx[, a], n_lv[a]) / nrow(idx)),
             share_efficient = as.numeric(tabulate(eff_idx[, a], n_lv[a]) / nrow(eff_idx)))
}))
bal[, 3:4] <- round(bal[, 3:4], 3)
write.csv(bal, "outputs/efficient_design/level_balance.csv", row.names = FALSE)
cat("\nLevel balance (share of alternatives):\n"); print(bal, row.names = FALSE)

cat("\nDesign written to data/efficient_design.csv (", nrow(rows), "rows )\n")
