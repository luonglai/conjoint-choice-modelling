# =============================================================================
# 01_simulate_data.R
# Design a choice-based conjoint (CBC) experiment for a premium skincare serum
# and simulate respondent-level choice data from three latent preference
# segments, each with known ground-truth utilities plus individual-level
# noise around their segment mean.
#
# Why simulated data? Real conjoint panels are proprietary and not publicly
# available. Simulating from known parameters lets us validate that the
# estimation pipeline (scripts 02 and 04) correctly recovers the true
# preference structure. This check is normally impossible with real data,
# where the "truth" is never known.
#
# Why three segments instead of one shared set of preferences? A single
# global preference vector gives an aggregate model nothing to fail at since
# every respondent is identical, so there is no heterogeneity for Hierarchical 
# Bayes to recover, and no business case for segmentation. Three segments with 
# genuinely different preferences let the project show what an aggregate model 
# (script 02) can and cannot see and what HB (script 04) recovers that it can't.
# =============================================================================

set.seed(42)

# 1. Attributes and levels
positioning_levels <- c("Mass-Prestige", "Dermatologist-Recommended", "Luxury Heritage")
ingredient_levels  <- c("Retinol", "Vitamin C", "Peptide Complex")
format_levels      <- c("Standard Bottle", "Airless Pump")
price_levels       <- c(45, 65, 85, 110)

# 2. Ground-truth segment-level utilities (unknown to the "researcher")
# Three latent segments with distinct preference structures. Coefficients are
# on the same effects-coded / continuous-price scale as before.
segment_names  <- c("Clinical/Dermatologist-Focused", "Luxury Heritage Seekers", "Price-Sensitive Pragmatists")
segment_shares <- c(0.35, 0.30, 0.35)

segment_utilities <- data.frame(
  segment    = segment_names,
  pos_derm   = c(1.60, 0.30, 0.40),
  pos_luxury = c(0.10, 1.30, 0.10),
  ing_vitc   = c(0.20, 0.50, 0.20),
  ing_pept   = c(0.90, 0.30, 0.20),
  fmt_pump   = c(0.30, 0.60, 0.10),
  price_coef = c(-0.020, -0.010, -0.045),
  none_asc   = c(-1.00, -1.20, -0.60)
)

# Population-weighted average 
population_utilities <- list(
  positioning = c("Mass-Prestige" = 0,
                   "Dermatologist-Recommended" = sum(segment_shares * segment_utilities$pos_derm),
                   "Luxury Heritage" = sum(segment_shares * segment_utilities$pos_luxury)),
  ingredient  = c("Retinol" = 0,
                   "Vitamin C" = sum(segment_shares * segment_utilities$ing_vitc),
                   "Peptide Complex" = sum(segment_shares * segment_utilities$ing_pept)),
  format      = c("Standard Bottle" = 0,
                   "Airless Pump" = sum(segment_shares * segment_utilities$fmt_pump)),
  price_coef  = sum(segment_shares * segment_utilities$price_coef),
  none_asc    = sum(segment_shares * segment_utilities$none_asc)
)

# 3. Assign respondents to segments, with individual noise
# Individual heterogeneity = segment mean + small Gaussian noise, matching
# the mixture-of-normals structure that Hierarchical Bayes (script 04) is
# designed to recover.
n_resp  <- 300
n_tasks <- 10
n_alts  <- 3

resp_segment <- sample(segment_names, n_resp, replace = TRUE, prob = segment_shares)

indiv_noise_sd <- 0.15  # within-segment individual-level noise
indiv_utilities <- do.call(rbind, lapply(1:n_resp, function(i) {
  seg <- segment_utilities[segment_utilities$segment == resp_segment[i], ]
  data.frame(
    resp_id    = i,
    segment    = resp_segment[i],
    pos_derm   = rnorm(1, seg$pos_derm, indiv_noise_sd),
    pos_luxury = rnorm(1, seg$pos_luxury, indiv_noise_sd),
    ing_vitc   = rnorm(1, seg$ing_vitc, indiv_noise_sd),
    ing_pept   = rnorm(1, seg$ing_pept, indiv_noise_sd),
    fmt_pump   = rnorm(1, seg$fmt_pump, indiv_noise_sd),
    price_coef = rnorm(1, seg$price_coef, indiv_noise_sd * 0.05),  # smaller noise, different scale
    none_asc   = rnorm(1, seg$none_asc, indiv_noise_sd)
  )
}))

write.csv(indiv_utilities, "data/respondent_true_utilities.csv", row.names = FALSE)

# 4. Build choice tasks
# 300 respondents x 10 choice tasks x 3 product alternatives + a "None" option.
# Levels are assigned at random per alternative (a simple random design; a
# production study would use a D-efficient / orthogonal design -- see README).

design <- expand.grid(resp_id = 1:n_resp, task_id = 1:n_tasks, alt_id = 1:n_alts)
design$positioning <- sample(positioning_levels, nrow(design), replace = TRUE)
design$ingredient  <- sample(ingredient_levels, nrow(design), replace = TRUE)
design$format      <- sample(format_levels, nrow(design), replace = TRUE)
design$price       <- sample(price_levels, nrow(design), replace = TRUE)

# 5. Compute deterministic utility using EACH respondent's OWN utilities
design <- merge(design, indiv_utilities, by = "resp_id")
design$utility <- ifelse(design$positioning == "Dermatologist-Recommended", design$pos_derm,
                   ifelse(design$positioning == "Luxury Heritage", design$pos_luxury, 0)) +
                  ifelse(design$ingredient == "Vitamin C", design$ing_vitc,
                   ifelse(design$ingredient == "Peptide Complex", design$ing_pept, 0)) +
                  ifelse(design$format == "Airless Pump", design$fmt_pump, 0) +
                  design$price_coef * design$price

# 6. Add a "None of these" alternative to every choice task
none_rows <- unique(design[, c("resp_id", "task_id", "segment", "none_asc")])
none_rows$alt_id      <- 0
none_rows$positioning <- NA
none_rows$ingredient  <- NA
none_rows$format      <- NA
none_rows$price       <- NA
none_rows$utility     <- none_rows$none_asc

common_cols <- c("resp_id","task_id","alt_id","segment","positioning","ingredient","format","price","utility")
full <- rbind(design[, common_cols], none_rows[, common_cols])
full <- full[order(full$resp_id, full$task_id, full$alt_id), ]
full$is_none <- as.integer(full$alt_id == 0)

# 7. Simulate choices: utility + Gumbel(0,1) noise, respondent picks max
rgumbel <- function(n) -log(-log(runif(n)))
full$rand_utility <- full$utility + rgumbel(nrow(full))
full$chosen <- ave(full$rand_utility, paste(full$resp_id, full$task_id),
                    FUN = function(x) as.integer(x == max(x)))

# 8. Save design + choice data
write.csv(full[, c("resp_id","task_id","alt_id","segment","positioning","ingredient",
                    "format","price","is_none","chosen")],
          "data/simulated_choices.csv", row.names = FALSE)

saveRDS(population_utilities, "data/true_utilities.rds")
saveRDS(segment_utilities, "data/true_segment_utilities.rds")
write.csv(data.frame(segment = segment_names, share = segment_shares),
          "data/true_segment_shares.csv", row.names = FALSE)

cat("Simulated", n_resp, "respondents across", length(segment_names), "segments,",
    n_tasks, "tasks x", n_alts, "alternatives (+None).\n")
cat("Segment sizes:\n"); print(table(resp_segment))
cat("Rows written:", nrow(full), "\n")
cat("Choice rate check (should be ~1 chosen per resp/task):",
    mean(ave(full$chosen, paste(full$resp_id, full$task_id), FUN = sum)), "\n")
