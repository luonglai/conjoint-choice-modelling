# =============================================================================
# 02_estimate_model.R
# Estimate a multinomial logit (MNL) choice model from the simulated CBC data
# using survival::clogit (conditional logistic regression on matched strata.
# =============================================================================

library(survival)
library(dplyr)

choices <- read.csv("data/simulated_choices.csv")
true_utilities <- readRDS("data/true_utilities.rds")

# 1. Prepare effects-coded dummies (reference level = first level)
choices <- choices %>%
  mutate(
    positioning = ifelse(is_none == 1, "Mass-Prestige", positioning), # placeholder for NA rows
    ingredient  = ifelse(is_none == 1, "Retinol", ingredient),
    format      = ifelse(is_none == 1, "Standard Bottle", format),
    price       = ifelse(is_none == 1, 0, price),

    pos_derm   = as.integer(positioning == "Dermatologist-Recommended" & is_none == 0),
    pos_luxury = as.integer(positioning == "Luxury Heritage" & is_none == 0),
    ing_vitc   = as.integer(ingredient == "Vitamin C" & is_none == 0),
    ing_pept   = as.integer(ingredient == "Peptide Complex" & is_none == 0),
    fmt_pump   = as.integer(format == "Airless Pump" & is_none == 0),
    price_eur  = ifelse(is_none == 1, 0, price),
    none_flag  = is_none,
    stratum    = paste(resp_id, task_id)
  )

# 2. Estimate via conditional logistic regression
# clogit models P(chosen | choice set) using a stratified Cox partial-likelihood 
# trick (strata = each respondent-task).
model <- clogit(
  chosen ~ pos_derm + pos_luxury + ing_vitc + ing_pept + fmt_pump +
           price_eur + none_flag + strata(stratum),
  data = choices
)

model_summary <- summary(model)
print(model_summary)

# 3. Compare recovered vs. true parameters (validation)
recovered <- coef(model)
comparison <- data.frame(
  parameter = c("Dermatologist-Recommended", "Luxury Heritage", "Vitamin C",
                "Peptide Complex", "Airless Pump", "Price (per EUR)", "None ASC"),
  true_value = c(
    true_utilities$positioning["Dermatologist-Recommended"],
    true_utilities$positioning["Luxury Heritage"],
    true_utilities$ingredient["Vitamin C"],
    true_utilities$ingredient["Peptide Complex"],
    true_utilities$format["Airless Pump"],
    true_utilities$price_coef,
    true_utilities$none_asc
  ),
  estimated_value = c(
    recovered["pos_derm"], recovered["pos_luxury"], recovered["ing_vitc"],
    recovered["ing_pept"], recovered["fmt_pump"], recovered["price_eur"],
    recovered["none_flag"]
  ),
  std_error = c(
    model_summary$coefficients["pos_derm", "se(coef)"],
    model_summary$coefficients["pos_luxury", "se(coef)"],
    model_summary$coefficients["ing_vitc", "se(coef)"],
    model_summary$coefficients["ing_pept", "se(coef)"],
    model_summary$coefficients["fmt_pump", "se(coef)"],
    model_summary$coefficients["price_eur", "se(coef)"],
    model_summary$coefficients["none_flag", "se(coef)"]
  )
)
comparison$abs_diff <- round(abs(comparison$true_value - comparison$estimated_value), 3)
comparison[, 2:5] <- round(comparison[, 2:5], 4)

write.csv(comparison, "outputs/model_validation.csv", row.names = FALSE)
cat("\n Recovered vs. true parameters\n")
print(comparison)

saveRDS(model, "outputs/clogit_model.rds")
cat("\nModel saved. Validation table saved to outputs/model_validation.csv\n")
