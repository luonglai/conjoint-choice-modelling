# =============================================================================
# 04_hb_estimation.R
# Fit a Hierarchical Bayes multinomial logit (HB-MNL) using bayesm's
# rhierMnlRwMixture.
#
# Unlike the pooled model in script 02, HB estimates a separate set of
# part-worths for every respondent (shrunk toward a population distribution
# via a Bayesian hierarchical prior), which is what makes segmentation and
# individual-level preference simulation possible.
# =============================================================================

library(bayesm)
library(dplyr)

choices <- read.csv("data/simulated_choices.csv")
true_indiv   <- read.csv("data/respondent_true_utilities.csv")
true_shares  <- read.csv("data/true_segment_shares.csv")

# Align naming with the model's dummy-variable names (price_coef -> price_eur,
# none_asc -> none_flag) so the two can be merged and compared directly.
names(true_indiv)[names(true_indiv) == "price_coef"] <- "price_eur"
names(true_indiv)[names(true_indiv) == "none_asc"]   <- "none_flag"

# 1. Build dummy-coded predictors (same coding as script 02)
choices <- choices %>%
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
    price_eur  = ifelse(is_none == 1, 0, price),
    none_flag  = is_none
  ) %>%
  arrange(resp_id, task_id, alt_id)  # alt_id: 0=None, 1-3=products -> consistent order per task

var_names <- c("pos_derm","pos_luxury","ing_vitc","ing_pept","fmt_pump","price_eur","none_flag")
nvar <- length(var_names)
p    <- 4  # alternatives per choice task (3 products + None)

# 2. Reshape into bayesm's per-respondent list format
resp_ids <- sort(unique(choices$resp_id))
lgtdata <- vector("list", length(resp_ids))

for (i in seq_along(resp_ids)) {
  rid <- resp_ids[i]
  d_r <- choices[choices$resp_id == rid, ]
  n_tasks_r <- length(unique(d_r$task_id))

  X_r <- as.matrix(d_r[, var_names])              # (n_tasks_r * p) x nvar
  y_r <- integer(n_tasks_r)
  for (t in seq_len(n_tasks_r)) {
    block <- d_r[d_r$task_id == sort(unique(d_r$task_id))[t], ]
    y_r[t] <- which(block$chosen == 1)             # position within the p-row block
  }
  lgtdata[[i]] <- list(y = y_r, X = X_r)
}

cat("Prepared", length(lgtdata), "respondents,", n_tasks_r, "tasks each, p =", p, ", nvar =", nvar, "\n")

# 3. Fit HB-MNL
# ncomp = 1: a single multivariate-normal distribution of individual-level
# betas (standard "random coefficients" HB). Segments are recovered afterward
# via clustering on the individual posterior means. This two-step approach 
# (HB then cluster) mirrors common commercial practice and avoids assuming the 
# "true" number of latent segments inside the model itself, even though we 
# happen to know it here (3) because this is simulated.
set.seed(42)
R_draws <- 8000
Data_hb <- list(lgtdata = lgtdata, p = p)
Mcmc_hb <- list(R = R_draws, keep = 5)
Prior_hb <- list(ncomp = 1)

cat("Running HB-MNL Gibbs sampler (", R_draws, "draws)")
t0 <- Sys.time()
out_hb <- rhierMnlRwMixture(Data = Data_hb, Prior = Prior_hb, Mcmc = Mcmc_hb)
cat("Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "minutes.\n")

saveRDS(out_hb, "outputs/hb_model.rds")

# 4. Extract individual-level posterior means (post-burn-in)
n_kept   <- dim(out_hb$betadraw)[3]
burn_in  <- floor(n_kept / 2)  # discard first half as burn-in
keep_idx <- (burn_in + 1):n_kept

indiv_betas <- apply(out_hb$betadraw[, , keep_idx], c(1, 2), mean)
colnames(indiv_betas) <- var_names
indiv_betas <- data.frame(resp_id = resp_ids, indiv_betas)

write.csv(indiv_betas, "outputs/hb_individual_partworths.csv", row.names = FALSE)

# 5. Validate: HB individual estimates vs. true individual utilities
compare <- merge(indiv_betas, true_indiv, by = "resp_id", suffixes = c("_hb", "_true"))
validation_indiv <- data.frame(
  parameter = var_names,
  correlation = sapply(var_names, function(v) {
    round(cor(compare[[paste0(v, "_hb")]], compare[[paste0(v, "_true")]]), 3)
  }),
  mean_abs_error = sapply(var_names, function(v) {
    round(mean(abs(compare[[paste0(v, "_hb")]] - compare[[paste0(v, "_true")]])), 3)
  })
)
write.csv(validation_indiv, "outputs/hb_individual_validation.csv", row.names = FALSE)
cat("\n HB individual-level recovery (correlation with true individual utilities)\n")
print(validation_indiv)

# 6. Aggregate (population mean) comparison: true vs clogit vs HB
clogit_val <- read.csv("outputs/model_validation.csv")
hb_pop_mean <- colMeans(indiv_betas[, var_names])

three_way <- data.frame(
  parameter        = clogit_val$parameter,
  true_population  = clogit_val$true_value,
  clogit_aggregate = clogit_val$estimated_value,
  hb_population_mean = round(as.numeric(hb_pop_mean[c("pos_derm","pos_luxury","ing_vitc",
                                                        "ing_pept","fmt_pump","price_eur","none_flag")]), 4)
)
write.csv(three_way, "outputs/three_way_comparison.csv", row.names = FALSE)
cat("\n True vs. pooled clogit vs. HB population mean\n")
print(three_way)

# 7. Post-hoc segmentation: k-means on individual HB part-worths
set.seed(42)
km <- kmeans(scale(indiv_betas[, var_names]), centers = 3, nstart = 25)
indiv_betas$cluster <- km$cluster

# Match cluster numbers to true segment labels by best overlap (Hungarian-style
# greedy match, since k-means cluster numbering is arbitrary).
compare$cluster <- indiv_betas$cluster[match(compare$resp_id, indiv_betas$resp_id)]
cross_tab <- table(true_segment = compare$segment, kmeans_cluster = compare$cluster)
write.csv(as.data.frame.matrix(cross_tab), "outputs/segment_recovery_crosstab.csv")
cat("\n Cross-tab: true segment vs. recovered k-means cluster\n")
print(cross_tab)

# Simple recovery accuracy: for each true segment, % falling in its dominant cluster
best_match_acc <- sum(apply(cross_tab, 1, max)) / sum(cross_tab)
cat("\nSegment recovery accuracy (best cluster-to-segment matching):",
    round(100 * best_match_acc, 1), "%\n")

cat("\nAll HB outputs written to outputs/\n")
