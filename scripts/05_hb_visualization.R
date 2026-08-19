# =============================================================================
# 05_hb_visualization.R
# Visualize what HB recovers that the pooled model (script 02) cannot:
# individual-level heterogeneity that separates by true latent segment, and
# business-usable segment profiles from post-hoc clustering.
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

indiv_betas <- read.csv("outputs/hb_individual_partworths.csv")
true_indiv  <- read.csv("data/respondent_true_utilities.csv")
d <- merge(indiv_betas, true_indiv[, c("resp_id", "segment")], by = "resp_id")

seg_colors <- c("Clinical/Dermatologist-Focused" = "#1F3A5F",
                 "Luxury Heritage Seekers" = "#C0392B",
                 "Price-Sensitive Pragmatists" = "#7F8C8D")

# 1. Individual heterogeneity recovered by HB, split by TRUE segment
plot_data <- d %>%
  select(segment, pos_derm, price_eur) %>%
  pivot_longer(cols = c(pos_derm, price_eur), names_to = "attribute", values_to = "value") %>%
  mutate(attribute = recode(attribute,
                             pos_derm = "Dermatologist-Recommended part-worth",
                             price_eur = "Price coefficient (per \u20ac)"))

p1 <- ggplot(plot_data, aes(x = segment, y = value, fill = segment)) +
  geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.8) +
  facet_wrap(~attribute, scales = "free_y") +
  scale_fill_manual(values = seg_colors) +
  labs(title = "HB recovers individual-level heterogeneity across true segments",
       subtitle = "Each point is one respondent's HB-estimated part-worth, grouped by their true segment",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("outputs/hb_heterogeneity_by_segment.png", p1, width = 10, height = 5, dpi = 150)

# 2. Recovered segment profiles (k-means clusters) vs. true segments
km <- kmeans(scale(indiv_betas[, c("pos_derm","pos_luxury","ing_vitc","ing_pept",
                                     "fmt_pump","price_eur","none_flag")]), centers = 3, nstart = 25)
indiv_betas$cluster <- km$cluster

# Best-match cluster numbers to true segment names by dominant overlap
d$cluster <- indiv_betas$cluster[match(d$resp_id, indiv_betas$resp_id)]
cross_tab <- table(d$segment, d$cluster)
match_map <- apply(cross_tab, 1, which.max)  # true segment -> dominant cluster number
cluster_label <- setNames(names(match_map), match_map)

cluster_profile <- indiv_betas %>%
  mutate(recovered_segment = cluster_label[as.character(cluster)]) %>%
  group_by(recovered_segment) %>%
  summarise(across(c(pos_derm, pos_luxury, ing_vitc, ing_pept, fmt_pump), mean), n = n()) %>%
  pivot_longer(cols = c(pos_derm, pos_luxury, ing_vitc, ing_pept, fmt_pump),
               names_to = "attribute", values_to = "utility") %>%
  mutate(attribute = recode(attribute,
                             pos_derm = "Derm-Recommended", pos_luxury = "Luxury Heritage",
                             ing_vitc = "Vitamin C", ing_pept = "Peptide Complex",
                             fmt_pump = "Airless Pump"))

write.csv(cluster_profile, "outputs/hb_segment_profiles.csv", row.names = FALSE)

p2 <- ggplot(cluster_profile, aes(x = attribute, y = utility, fill = recovered_segment)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = seg_colors) +
  labs(title = "Recovered segment profiles (k-means on HB part-worths)",
       subtitle = "Mean part-worth utility per attribute, by segment recovered from HB estimates",
       x = NULL, y = "Mean part-worth utility", fill = "Recovered segment") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
ggsave("outputs/hb_segment_profiles.png", p2, width = 9, height = 5.5, dpi = 150)

cat("Cluster to true-segment mapping:\n")
print(cluster_label)
cat("\nHB visualizations written to outputs/\n")
