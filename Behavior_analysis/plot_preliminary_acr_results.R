args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

fly <- read_csv(file.path(out, "fly_level_metrics.csv"), show_col_types = FALSE) %>%
  mutate(group = if_else(experimental, "Experimental ACR", "GFP / empty control"))

long <- fly %>%
  select(cohort, fly, group, response_less_heading_change,
         response_more_straight, response_more_consistent, response_faster_walking) %>%
  pivot_longer(starts_with("response_"), names_to = "outcome", values_to = "change") %>%
  mutate(outcome = recode(outcome,
    response_less_heading_change = "Less heading change (deg)",
    response_more_straight = "More path straightness",
    response_more_consistent = "More heading consistency (rho)",
    response_faster_walking = "Faster while walking"))

p1 <- ggplot(long, aes(group, change, color = group)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.12) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.8, size = 2) +
  facet_wrap(vars(cohort, outcome), scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Experimental ACR" = "#D1495B", "GFP / empty control" = "#777777")) +
  labs(title = "Fly-level optogenetic responses",
       subtitle = "Positive change is in the hypothesized direction",
       x = NULL, y = "Hypothesis-oriented during - pre change") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out, "fly_level_changes.png"), p1, width = 14, height = 7, dpi = 200)

p2 <- ggplot(fly, aes(response_faster_walking, response_less_heading_change, color = group)) +
  geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") +
  geom_point(size = 2.5, alpha = 0.85) +
  facet_wrap(~ cohort, scales = "free") +
  scale_color_manual(values = c("Experimental ACR" = "#D1495B", "GFP / empty control" = "#777777")) +
  labs(title = "Joint response space", subtitle = "Faster and less heading change is upper right",
       x = "Change in forward velocity while walking",
       y = "Negative change in total heading change", color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(out, "joint_response_space.png"), p2, width = 10, height = 5, dpi = 200)
