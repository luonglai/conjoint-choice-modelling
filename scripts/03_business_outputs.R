# =============================================================================
# 03_business_outputs.R
# Turn the estimated model into the outputs a client actually asks for:
#   - relative attribute importance
#   - willingness-to-pay (WTP) per attribute level
#   - preference-share simulation for candidate product concepts
# =============================================================================

library(ggplot2)
library(dplyr)
library(scales)

model <- readRDS("outputs/clogit_model.rds")
b <- coef(model)

# 1. Part-worth utilities (reference level = 0)
part_worths <- data.frame(
  attribute = c("Positioning","Positioning","Positioning",
                "Active Ingredient","Active Ingredient","Active Ingredient",
                "Format","Format"),
  level = c("Mass-Prestige","Dermatologist-Recommended","Luxury Heritage",
            "Retinol","Vitamin C","Peptide Complex",
            "Standard Bottle","Airless Pump"),
  utility = c(0, b["pos_derm"], b["pos_luxury"],
              0, b["ing_vitc"], b["ing_pept"],
              0, b["fmt_pump"])
)

# 2. Relative attribute importance
importance <- part_worths %>%
  group_by(attribute) %>%
  summarise(range = max(utility) - min(utility)) %>%
  mutate(importance_pct = round(100 * range / sum(range), 1)) %>%
  arrange(desc(importance_pct))

write.csv(importance, "outputs/attribute_importance.csv", row.names = FALSE)

p1 <- ggplot(importance, aes(x = reorder(attribute, importance_pct), y = importance_pct)) +
  geom_col(fill = "#1F3A5F", width = 0.6) +
  geom_text(aes(label = paste0(importance_pct, "%")), hjust = -0.15, size = 4.2) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(importance$importance_pct) * 1.2), expand = c(0,0)) +
  labs(title = "Relative attribute importance",
       subtitle = "Share of total preference driven by each attribute",
       x = NULL, y = "Importance (%)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggsave("outputs/attribute_importance.png", p1, width = 7, height = 4, dpi = 150)

# 3. Willingness-to-pay (EUR) per level, relative to reference
price_coef <- b["price_eur"]
wtp <- part_worths %>%
  mutate(wtp_eur = round(utility / -price_coef, 1)) %>%
  filter(utility != 0)

write.csv(wtp, "outputs/willingness_to_pay.csv", row.names = FALSE)

p2 <- ggplot(wtp, aes(x = reorder(level, wtp_eur), y = wtp_eur, fill = attribute)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0("+\u20ac", wtp_eur)), hjust = -0.1, size = 4) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(wtp$wtp_eur) * 1.25), expand = c(0,0)) +
  scale_fill_manual(values = c("Positioning" = "#1F3A5F", "Active Ingredient" = "#C0392B", "Format" = "#7F8C8D")) +
  labs(title = "Willingness-to-pay premium vs. reference level",
       subtitle = "EUR a consumer is estimated to pay above the base option",
       x = NULL, y = "WTP premium (\u20ac)", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "bottom")
ggsave("outputs/willingness_to_pay.png", p2, width = 9, height = 4.5, dpi = 150)

# 4. Preference-share simulation for candidate product concepts
# A typical client question: "if we launched these 3 concepts against each
# other, what share of preference would each capture?"
concepts <- data.frame(
  concept = c("A: Derm-Recommended Peptide, Pump, \u20ac85",
              "B: Luxury Heritage Vitamin C, Bottle, \u20ac65",
              "C: Mass-Prestige Retinol, Bottle, \u20ac45"),
  pos_derm   = c(1, 0, 0),
  pos_luxury = c(0, 1, 0),
  ing_vitc   = c(0, 1, 0),
  ing_pept   = c(1, 0, 0),
  fmt_pump   = c(1, 0, 0),
  price_eur  = c(85, 65, 45)
)

concepts$utility <- b["pos_derm"]*concepts$pos_derm + b["pos_luxury"]*concepts$pos_luxury +
                     b["ing_vitc"]*concepts$ing_vitc + b["ing_pept"]*concepts$ing_pept +
                     b["fmt_pump"]*concepts$fmt_pump + b["price_eur"]*concepts$price_eur

# Include the "None" option in the competitive set at its estimated ASC
none_utility <- b["none_flag"]
all_utils <- c(concepts$utility, none_utility)
shares <- exp(all_utils) / sum(exp(all_utils))

share_table <- data.frame(
  concept = c(concepts$concept, "None of these"),
  utility = round(all_utils, 3),
  preference_share_pct = round(100 * shares, 1)
)
write.csv(share_table, "outputs/preference_share.csv", row.names = FALSE)

p3 <- ggplot(share_table, aes(x = reorder(concept, preference_share_pct), y = preference_share_pct)) +
  geom_col(fill = "#1F3A5F", width = 0.55) +
  geom_text(aes(label = paste0(preference_share_pct, "%")), hjust = -0.15, size = 4.2) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(share_table$preference_share_pct) * 1.25), expand = c(0,0)) +
  labs(title = "Simulated preference share",
       subtitle = "Three candidate concepts competing head-to-head, incl. \"None\"",
       x = NULL, y = "Preference share (%)") +
  theme_minimal(base_size = 12.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggsave("outputs/preference_share.png", p3, width = 9, height = 4, dpi = 150)

cat("Attribute importance:\n"); print(importance)
cat("\nWillingness-to-pay:\n"); print(wtp[, c("attribute","level","wtp_eur")])
cat("\nPreference share:\n"); print(share_table)
cat("\nAll outputs written to outputs/\n")
