suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

root <- file.path("Behavior_analysis", "smp544_pivot_results")
d <- read_csv(file.path(root, "presentation_metrics.csv"), show_col_types = FALSE)
traces <- read_csv(file.path(root, "smp544_opto_aligned_traces.csv"), show_col_types = FALSE)

permutation_group_p <- function(x, g, max_exact = 1e6, n_mc = 1e5) {
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]; g <- as.character(g[keep])
  lev <- unique(g); stopifnot(length(lev) == 2)
  n1 <- sum(g == lev[1]); obs <- abs(mean(x[g == lev[1]]) - mean(x[g == lev[2]]))
  if (choose(length(x), n1) <= max_exact) {
    sets <- combn(seq_along(x), n1)
    null <- apply(sets, 2, function(i) abs(mean(x[i]) - mean(x[-i])))
    return(mean(null >= obs - 1e-12))
  }
  set.seed(544 + length(x) + n1)
  null <- replicate(n_mc, {
    i <- sample.int(length(x), n1)
    abs(mean(x[i]) - mean(x[-i]))
  })
  (sum(null >= obs - 1e-12) + 1) / (n_mc + 1)
}

pivot_flag <- function(data, phase, window, degrees, displacement, directionality = .7) {
  data[[paste0("abs_net_heading_deg_", phase, "_", window)]] >= degrees &
    data[[paste0("displacement_", phase, "_", window)]] <= displacement &
    data[[paste0("turn_directionality_", phase, "_", window)]] >= directionality
}

# Primary operational definition: a directed turn of at least 45 degrees with
# no more than 2 mm net displacement in the matched 1 s window.
d <- d %>%
  mutate(
    pivot_pre = pivot_flag(., "pre", "1s", 45, 2),
    pivot_opto = pivot_flag(., "opto", "1s", 45, 2),
    comparison_group = case_when(
      dataset == "smp544" & group == "experimental" ~ "SMP544 experimental",
      dataset == "smp544" & group == "control" ~ "SMP544 control",
      dataset == "ves041_activation" & group == "experimental" ~ "VES041 activation experimental",
      dataset == "ves041_activation" & group == "control" ~ "VES041 activation control",
      TRUE ~ dataset
    )
  )

fly <- d %>%
  group_by(dataset, fly, group, comparison_group) %>%
  summarise(
    n_presentations = n(),
    pivot_rate_pre = mean(pivot_pre),
    pivot_rate_opto = mean(pivot_opto),
    delta_pivot_rate = pivot_rate_opto - pivot_rate_pre,
    mean_delta_abs_net_heading = mean(delta_abs_net_heading_deg_1s),
    mean_delta_p95_angular_velocity = mean(delta_p95_abs_angular_velocity_1s),
    median_pivot_index = median(pivot_index_1s),
    mad_pivot_index = mad(pivot_index_1s, constant = 1),
    sd_delta_abs_net_heading = sd(delta_abs_net_heading_deg_1s),
    mad_delta_abs_net_heading = mad(delta_abs_net_heading_deg_1s, constant = 1),
    .groups = "drop"
  )
write_csv(fly, file.path(root, "fly_level_summary.csv"))

metrics <- c(
  delta_pivot_rate = "Change in pivot-like rate",
  mean_delta_abs_net_heading = "Change in absolute net turn",
  mean_delta_p95_angular_velocity = "Change in p95 angular velocity",
  mad_pivot_index = "Within-fly MAD of pivot index",
  sd_delta_abs_net_heading = "Within-fly SD of turn response",
  mad_delta_abs_net_heading = "Within-fly MAD of turn response"
)

tests <- list(); k <- 1
for (metric in names(metrics)) {
  ves <- filter(fly, dataset == "ves041_activation")
  tests[[k]] <- tibble(
    comparison = "VES041 activation experimental vs control", metric = metric,
    estimate = mean(ves[[metric]][ves$group == "experimental"]) -
      mean(ves[[metric]][ves$group == "control"]),
    p_value = permutation_group_p(ves[[metric]], ves$group)
  ); k <- k + 1
}
tests <- bind_rows(tests) %>%
  group_by(comparison) %>% mutate(bh_q = p.adjust(p_value, "BH")) %>% ungroup()
write_csv(tests, file.path(root, "ves041_activation_tests.csv"))

# Sensitivity across plausible turn and displacement cutoffs.
sensitivity <- list(); k <- 1
for (degrees in c(30, 45, 60)) for (disp in c(1, 2, 3)) {
  temp <- d %>%
    mutate(pre = pivot_flag(., "pre", "1s", degrees, disp),
           opto = pivot_flag(., "opto", "1s", degrees, disp)) %>%
    group_by(dataset, fly, group) %>%
    summarise(pre = mean(pre), opto = mean(opto), delta = opto - pre, .groups = "drop")
  smp <- filter(temp, dataset == "smp544")
  smp_exp <- filter(temp, dataset == "smp544", group == "experimental")
  ves <- filter(temp, dataset == "ves041_activation")
  ves_exp <- filter(ves, group == "experimental")
  sensitivity[[k]] <- tibble(
    turn_threshold_deg = degrees, displacement_threshold_mm = disp,
    smp_exp_pre = mean(smp_exp$pre), smp_exp_opto = mean(smp_exp$opto),
    smp_exp_delta = mean(smp_exp$delta), smp_ctrl_delta = mean(smp$delta[smp$group == "control"]),
    ves_activation_exp_delta = mean(ves_exp$delta),
    ves_activation_ctrl_delta = mean(ves$delta[ves$group == "control"])
  ); k <- k + 1
}
sensitivity <- bind_rows(sensitivity)
write_csv(sensitivity, file.path(root, "pivot_threshold_sensitivity.csv"))

theme_set(theme_bw(base_size = 11))
colors <- c("SMP544 experimental" = "#C73E4D", "SMP544 control" = "#606060",
            "VES041 activation experimental" = "#3A76AF",
            "VES041 activation control" = "#8DABC6")

# Rotation versus translation: upper-left is the pivot-like region.
p_scatter <- ggplot(d, aes(displacement_opto_1s, abs_net_heading_deg_opto_1s,
                           color = comparison_group)) +
  annotate("rect", xmin = 0, xmax = 2, ymin = 45, ymax = Inf,
           fill = "#E8C547", alpha = .12) +
  geom_vline(xintercept = 2, linetype = 2) + geom_hline(yintercept = 45, linetype = 2) +
  geom_point(alpha = .22, size = 1) +
  facet_wrap(~comparison_group, scales = "free", nrow = 1) +
  scale_color_manual(values = colors) +
  coord_cartesian(xlim = c(0, 18), ylim = c(0, 180)) +
  labs(title = "Rotation and translation during the first second of opto",
       subtitle = "Shading marks the rotation/displacement part of the pivot-like definition",
       x = "Net displacement (mm)", y = "Absolute net heading change (degrees)") +
  theme(legend.position = "none")
ggsave(file.path(root, "pivot_rotation_vs_translation.png"), p_scatter,
       width = 13, height = 4.4, dpi = 220)

# Fly-level pre and opto rates.
rate_long <- fly %>%
  filter(dataset == "smp544") %>%
  select(fly, comparison_group, pivot_rate_pre, pivot_rate_opto) %>%
  pivot_longer(starts_with("pivot_rate"), names_to = "period", values_to = "rate") %>%
  mutate(period = factor(period, levels = c("pivot_rate_pre", "pivot_rate_opto"),
                         labels = c("Pre-opto", "Opto")))
p_rate <- ggplot(rate_long, aes(period, rate, group = fly, color = comparison_group)) +
  geom_line(alpha = .7) + geom_point(size = 2.5) +
  facet_wrap(~comparison_group) + scale_color_manual(values = colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Pivot-like presentations by fly",
       subtitle = "At least 45 degrees directed rotation and at most 2 mm displacement",
       x = NULL, y = "Presentations meeting the definition") +
  theme(legend.position = "none")
ggsave(file.path(root, "pivot_rate_pre_vs_opto.png"), p_rate,
       width = 8.5, height = 4.5, dpi = 220)

# Robust and conventional within-fly variability comparisons.
var_long <- fly %>%
  select(fly, comparison_group, mad_pivot_index, sd_delta_abs_net_heading) %>%
  pivot_longer(c(mad_pivot_index, sd_delta_abs_net_heading), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    mad_pivot_index = "MAD of pivot index",
    sd_delta_abs_net_heading = "SD of opto-minus-pre net turn"))
p_var <- ggplot(var_long, aes(comparison_group, value, color = comparison_group)) +
  geom_boxplot(outlier.shape = NA, alpha = .08, width = .55) +
  geom_jitter(width = .09, size = 2.2) +
  facet_wrap(~metric, scales = "free_y") + scale_color_manual(values = colors) +
  labs(title = "Within-fly presentation-to-presentation variability",
       x = NULL, y = "Within-fly variability") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 22, hjust = 1))
ggsave(file.path(root, "pivot_response_variability.png"), p_var,
       width = 10.5, height = 4.7, dpi = 220)

# Opto-aligned fly-weighted profiles.
profile <- traces %>%
  mutate(abs_angular_velocity = abs(angular_velocity),
         comparison_group = if_else(group == "experimental", "SMP544 experimental", "SMP544 control"),
         time_bin = round(time_from_onset_s / .04) * .04) %>%
  group_by(fly, comparison_group, time_bin) %>%
  summarise(abs_angular_velocity = mean(abs_angular_velocity, na.rm = TRUE),
            forward_velocity = mean(forward_velocity, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(c(abs_angular_velocity, forward_velocity), names_to = "metric", values_to = "value") %>%
  group_by(comparison_group, time_bin, metric) %>%
  summarise(mean = mean(value), se = sd(value) / sqrt(n()), .groups = "drop") %>%
  mutate(metric = recode(metric, abs_angular_velocity = "Absolute angular velocity",
                         forward_velocity = "Forward velocity"))
p_profile <- ggplot(profile, aes(time_bin, mean, color = comparison_group, fill = comparison_group)) +
  annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#E8C547", alpha = .08) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), color = NA, alpha = .14) +
  geom_line(linewidth = .8) + geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) + scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(title = "SMP544 opto-aligned kinematics",
       subtitle = "Fly-weighted mean and standard error; shaded interval is the first second of opto",
       x = "Time from opto onset (s)", y = NULL, color = NULL, fill = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(root, "smp544_opto_aligned_kinematics.png"), p_profile,
       width = 8.8, height = 7, dpi = 220)

# Show the highest-scoring SMP events, including both groups, as trajectories.
top_events <- d %>% filter(dataset == "smp544") %>%
  group_by(group) %>% slice_max(pivot_index_1s, n = 6, with_ties = FALSE) %>%
  ungroup() %>% select(fly, group, recording, pulse, pivot_index_1s)
examples <- traces %>%
  inner_join(top_events, by = c("fly", "group", "recording", "pulse")) %>%
  filter(time_from_onset_s >= -1, time_from_onset_s <= 1) %>%
  mutate(period = if_else(time_from_onset_s < 0, "Pre-opto", "Opto"),
         panel = paste(group, fly, "pulse", pulse))
p_examples <- ggplot(examples, aes(x_relative, y_relative, color = period,
                                   group = interaction(panel, period))) +
  geom_path(linewidth = .75) + geom_point(data = filter(examples, abs(time_from_onset_s) < .011),
                                         color = "black", size = 1.5) +
  facet_wrap(~panel, scales = "free", ncol = 6) +
  scale_color_manual(values = c("Pre-opto" = "#9A9A9A", "Opto" = "#C73E4D")) +
  labs(title = "Highest pivot-index SMP544 presentations",
       subtitle = "Paths are centered at opto onset; panels were selected without using genotype-specific thresholds",
       x = "Relative x position (mm)", y = "Relative y position (mm)", color = NULL) +
  theme(legend.position = "bottom", axis.text = element_blank(), axis.ticks = element_blank())
ggsave(file.path(root, "top_pivot_like_trajectories.png"), p_examples,
       width = 14, height = 6.3, dpi = 220)

cat("Wrote SMP544 pivot analysis to", root, "\n")
