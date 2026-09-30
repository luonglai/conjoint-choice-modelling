# =============================================================================
# 08_market_simulator.R
# Individual-level market simulator built on the Hierarchical Bayes part-worths
# from script 04 (outputs/hb_individual_partworths.csv).
#
# Script 03 already simulates preference share from the POOLED model. This
# script upgrades that in three ways a client would care about:
#   1. Uses each respondent's own part-worths, then aggregates -- so
#      heterogeneity is not averaged away before shares are computed.
#   2. Supports two standard decision rules: share of preference (logit) and
#      first choice.
#   3. Runs what-if scenarios: a price sweep for one concept (with a revenue
#      index to locate a price optimum) and shares by recovered segment.
#
# Reusable entry point: simulate_shares(concepts, betas, rule)
#
# Outputs (outputs/simulator/)
#   sim_base_market.csv, sim_price_scenario.csv, sim_price_curve.png,
#   sim_shares_by_segment.csv, sim_shares_by_segment.png
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })
dir.create("outputs/simulator", showWarnings = FALSE, recursive = TRUE)

betas <- read.csv("outputs/hb_individual_partworths.csv")
attr_cols <- c("pos_derm", "pos_luxury", "ing_vitc", "ing_pept", "fmt_pump", "price_eur")

# Core simulator 
# concepts: data.frame with columns name, positioning, ingredient, format, price
# rule:     "share_of_preference" (logit share per respondent, averaged) or
#           "first_choice" (each respondent picks their highest-utility option)
# The "None of these" option is always in the competitive set, at each
# respondent's own None constant.
concept_matrix <- function(concepts) {
  cbind(
    pos_derm   = as.numeric(concepts$positioning == "Dermatologist-Recommended"),
    pos_luxury = as.numeric(concepts$positioning == "Luxury Heritage"),
    ing_vitc   = as.numeric(concepts$ingredient == "Vitamin C"),
    ing_pept   = as.numeric(concepts$ingredient == "Peptide Complex"),
    fmt_pump   = as.numeric(concepts$format == "Airless Pump"),
    price_eur  = concepts$price
  )
}

respondent_utilities <- function(concepts, betas) {
  B <- as.matrix(betas[, attr_cols])
  U <- B %*% t(concept_matrix(concepts))              # respondents x concepts
  cbind(U, None = betas$none_flag)                    # add the None option
}

simulate_shares <- function(concepts, betas, rule = c("share_of_preference", "first_choice"),
                            by = NULL) {
  rule <- match.arg(rule)
  U <- respondent_utilities(concepts, betas)
  P <- if (rule == "share_of_preference") {
    E <- exp(U - apply(U, 1, max)); E / rowSums(E)
  } else {
    M <- matrix(0, nrow(U), ncol(U)); M[cbind(seq_len(nrow(U)), max.col(U, ties.method = "first"))] <- 1; M
  }
  colnames(P) <- c(concepts$name, "None of these")
  if (is.null(by)) {
    data.frame(option = colnames(P), share_pct = round(100 * colMeans(P), 1))
  } else {
    as.data.frame(P) %>% mutate(group = by) %>%
      group_by(group) %>% summarise(across(everything(), ~ round(100 * mean(.x), 1)), .groups = "drop") %>%
      pivot_longer(-group, names_to = "option", values_to = "share_pct")
  }
}

# 1. Base market: the same three concepts as script 03
concepts <- data.frame(
  name        = c("A: Derm-Recommended Peptide, Pump, EUR 85",
                  "B: Luxury Heritage Vitamin C, Bottle, EUR 65",
                  "C: Mass-Prestige Retinol, Bottle, EUR 45"),
  positioning = c("Dermatologist-Recommended", "Luxury Heritage", "Mass-Prestige"),
  ingredient  = c("Peptide Complex", "Vitamin C", "Retinol"),
  format      = c("Airless Pump", "Standard Bottle", "Standard Bottle"),
  price       = c(85, 65, 45)
)

base <- data.frame(
  option = c(concepts$name, "None of these"),
  HB_share_of_preference = simulate_shares(concepts, betas, "share_of_preference")$share_pct,
  HB_first_choice        = simulate_shares(concepts, betas, "first_choice")$share_pct
)
pooled_path <- "outputs/preference_share.csv"           # from script 03 (pooled clogit)
if (file.exists(pooled_path)) {
  pooled <- read.csv(pooled_path)
  base$pooled_model_share <- pooled$preference_share_pct[match(gsub("EUR ", "\u20ac", base$option),
                                                                pooled$concept)]
  if (anyNA(base$pooled_model_share)) base$pooled_model_share <- pooled$preference_share_pct
}
write.csv(base, "outputs/simulator/sim_base_market.csv", row.names = FALSE)
cat("Base market -- share of preference (%), by model and decision rule\n")
print(base, row.names = FALSE)

# 2. Price scenario: sweep concept A's price
# Competitors B and C stay fixed. Price utility is linear in the model, so
# intermediate prices (e.g. EUR 50) are interpolated within the tested range
# (EUR 45-110); nothing is extrapolated beyond it.
#
# Three curves are compared:
#   HB      individual-level simulation (this script)
#   Pooled  aggregate clogit model (script 02/03)
#   Truth   the known simulated utilities -- only available because the data
#           are simulated; it is the benchmark a real study never has.
price_grid <- seq(45, 110, by = 5)

sweep_for <- function(betas_df, label) {
  do.call(rbind, lapply(price_grid, function(p) {
    cs <- concepts; cs$price[1] <- p
    s <- simulate_shares(cs, betas_df, "share_of_preference")
    data.frame(model = label, price_A = p, share_A = s$share_pct[1], share_None = s$share_pct[4])
  }))
}

# Truth: rename ground-truth columns to the model's parameter names
true_b <- read.csv("data/respondent_true_utilities.csv")
names(true_b)[names(true_b) == "price_coef"] <- "price_eur"
names(true_b)[names(true_b) == "none_asc"]   <- "none_flag"

# Pooled: one "respondent" carrying the pooled clogit coefficients
pooled_b <- as.data.frame(t(coef(readRDS("outputs/clogit_model.rds"))[c(attr_cols, "none_flag")]))

sweep <- rbind(sweep_for(betas, "HB (individual-level)"),
               sweep_for(pooled_b, "Pooled clogit"),
               sweep_for(true_b, "Truth (simulated)")) %>%
  mutate(revenue_index = share_A * price_A / 100)
write.csv(sweep, "outputs/simulator/sim_price_scenario.csv", row.names = FALSE)

peaks <- sweep %>% group_by(model) %>% slice_max(revenue_index, n = 1, with_ties = FALSE) %>% ungroup()
truth_share <- sweep$share_A[sweep$model == "Truth (simulated)"]
acc <- sweep %>% group_by(model) %>%
  summarise(mean_abs_share_error_pts = round(mean(abs(share_A - truth_share)), 1), .groups = "drop")
cat("\nPrice scenario for concept A (share of preference, %):\n")
print(as.data.frame(sweep %>% filter(price_A %in% c(45, 65, 85, 110)) %>%
                      select(model, price_A, share_A) %>% pivot_wider(names_from = price_A, values_from = share_A)),
      row.names = FALSE)
cat("\nMean absolute error of the share curve vs. truth (percentage points):\n")
print(as.data.frame(acc), row.names = FALSE)
cat("\nPrice at which the revenue index (share x price) peaks:\n")
print(as.data.frame(peaks[, c("model", "price_A", "share_A")]), row.names = FALSE)
if (any(peaks$price_A %in% range(price_grid)))
  cat("NOTE: a peak on the edge of the tested range (EUR 45-110) means the optimum is not\n",
      "identified inside the range. Read it as 'still rising at the top', not as an optimum.\n")

plot_df <- sweep %>%
  transmute(model, price_A, `Share of preference (%)` = share_A,
            `Revenue index (share x price)` = revenue_index) %>%
  pivot_longer(c(`Share of preference (%)`, `Revenue index (share x price)`),
               names_to = "metric", values_to = "value")
p_price <- ggplot(plot_df, aes(price_A, value, colour = model, linetype = model)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("HB (individual-level)" = "#1F3A5F", "Pooled clogit" = "#C0392B",
                                 "Truth (simulated)" = "#7F8C8D")) +
  scale_linetype_manual(values = c("HB (individual-level)" = "solid", "Pooled clogit" = "solid",
                                   "Truth (simulated)" = "dashed")) +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "What-if: price of concept A (Derm-Recommended Peptide, Pump)",
       subtitle = "Competitors B and C held fixed. Truth is known only because the data are simulated",
       x = "Price of concept A (EUR)", y = NULL, colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave("outputs/simulator/sim_price_curve.png", p_price, width = 9.5, height = 4.5, dpi = 150)

# 3. Shares by recovered segment
# Same clustering as script 04 (k-means, k = 3, seed 42) so cluster IDs match.
set.seed(42)
km <- kmeans(scale(betas[, c(attr_cols, "none_flag")]), centers = 3, nstart = 25)
betas$cluster <- km$cluster

# Because the data are simulated, name each cluster after the true segment it
# overlaps most (validation only -- with real data you would name clusters from
# their part-worth profiles).
true_seg <- read.csv("data/respondent_true_utilities.csv")[, c("resp_id", "segment")]
xt <- table(cluster = betas$cluster[match(true_seg$resp_id, betas$resp_id)], segment = true_seg$segment)
cluster_names <- setNames(colnames(xt)[apply(xt, 1, which.max)], rownames(xt))
if (anyDuplicated(cluster_names)) cluster_names <- setNames(paste("Cluster", rownames(xt)), rownames(xt))
seg_label <- paste0(cluster_names[as.character(betas$cluster)], " (n=",
                    table(betas$cluster)[as.character(betas$cluster)], ")")

by_seg <- simulate_shares(concepts, betas, "share_of_preference", by = seg_label)
write.csv(by_seg, "outputs/simulator/sim_shares_by_segment.csv", row.names = FALSE)
cat("\nShare of preference (%) by recovered segment\n")
print(as.data.frame(pivot_wider(by_seg, names_from = group, values_from = share_pct)), row.names = FALSE)

p_seg <- ggplot(by_seg, aes(x = option, y = share_pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_x_discrete(labels = function(x) sub(":.*", "", x)) +
  scale_fill_manual(values = c("#1F3A5F", "#C0392B", "#7F8C8D")) +
  labs(title = "Who prefers which concept?", subtitle = "Share of preference by recovered HB segment",
       x = NULL, y = "Share of preference (%)", fill = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("outputs/simulator/sim_shares_by_segment.png", p_seg, width = 9, height = 4.5, dpi = 150)

cat("\nAll simulator outputs written to outputs/simulator/\n")
