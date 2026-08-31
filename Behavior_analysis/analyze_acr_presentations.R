args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else file.path("Behavior_analysis", "preliminary_acr_results")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

options(contrasts = c("contr.treatment", "contr.poly"))

bh <- function(p) p.adjust(p, method = "BH")

pulse <- read_csv(file.path(out, "pulse_level_metrics.csv"), show_col_types = FALSE) %>%
  mutate(
    experimental = as.logical(experimental),
    moving_at_onset = as.logical(moving_at_onset),
    moving_pre_1s = as.logical(moving_pre_1s),
    movement_state = factor(if_else(moving_at_onset, "moving", "not_moving"),
                            levels = c("not_moving", "moving")),
    group = factor(if_else(experimental, "experimental", "control"),
                   levels = c("control", "experimental")),
    pulse_z = as.numeric(scale(pulse)),
    response_less_heading_change = -(hchange_deg_during - hchange_deg_pre),
    response_more_straight = straightness_during - straightness_pre,
    response_more_consistent = rho_during - rho_pre,
    response_faster_walking = fwd_walking_during - fwd_walking_pre,
    response_faster_all = fwd_all_during - fwd_all_pre,
    wiggle_heading_change = hchange_deg_pre,
    wiggle_low_straightness = 1 - straightness_pre,
    wiggle_low_rho = 1 - rho_pre
  )

write_csv(pulse, file.path(out, "presentation_level_metrics.csv"))

movement_counts <- pulse %>%
  group_by(cohort, group, movement_state) %>%
  summarise(n_presentations = n(), n_flies = n_distinct(fly),
            median_presentations_per_fly = median(table(fly)), .groups = "drop") %>%
  group_by(cohort, group) %>%
  mutate(proportion_presentations = n_presentations / sum(n_presentations)) %>%
  ungroup()
write_csv(movement_counts, file.path(out, "movement_state_counts.csv"))

outcomes <- c(
  less_heading_change = "response_less_heading_change",
  more_straightness = "response_more_straight",
  more_heading_consistency = "response_more_consistent",
  faster_while_walking = "response_faster_walking",
  higher_mean_forward_velocity = "response_faster_all"
)

response_rows <- list()
diagnostic_rows <- list()
row_i <- 1
diag_i <- 1

for (coh in unique(pulse$cohort)) {
  cohort_data <- filter(pulse, cohort == coh)
  for (outcome_name in names(outcomes)) {
    outcome <- outcomes[[outcome_name]]
    model_data <- cohort_data %>%
      select(fly, group, movement_state, pulse_z, all_of(outcome)) %>%
      drop_na()
    formula <- as.formula(paste0(outcome, " ~ group * movement_state + pulse_z + (1 + movement_state | fly)"))
    model <- suppressWarnings(lmer(formula, data = model_data, REML = TRUE,
                                   control = lmerControl(check.conv.singular = "ignore")))
    singular <- isSingular(model, tol = 1e-5)

    means <- as.data.frame(emmeans(model, ~ group | movement_state))
    tests0 <- as.data.frame(test(emmeans(model, ~ group | movement_state), null = 0, adjust = "none"))
    contrasts <- as.data.frame(contrast(emmeans(model, ~ group | movement_state),
                                       method = "revpairwise", adjust = "none"))
    movement_contrasts <- as.data.frame(contrast(emmeans(model, ~ movement_state | group),
                                                method = "revpairwise", adjust = "none"))

    for (j in seq_len(nrow(means))) {
      test_row <- tests0[j, ]
      response_rows[[row_i]] <- tibble(
        cohort = coh, outcome = outcome_name,
        movement_state = as.character(means$movement_state[j]),
        comparison = as.character(means$group[j]),
        estimate = means$emmean[j], SE = means$SE[j], df = means$df[j],
        lower_95 = means$emmean[j] - qt(.975, means$df[j]) * means$SE[j],
        upper_95 = means$emmean[j] + qt(.975, means$df[j]) * means$SE[j],
        t_ratio = test_row$t.ratio, p_value = test_row$p.value,
        n_presentations = nrow(model_data), n_flies = n_distinct(model_data$fly),
        singular_fit = singular
      )
      row_i <- row_i + 1
    }
    for (j in seq_len(nrow(contrasts))) {
      response_rows[[row_i]] <- tibble(
        cohort = coh, outcome = outcome_name,
        movement_state = as.character(contrasts$movement_state[j]),
        comparison = "experimental_minus_control",
        estimate = contrasts$estimate[j], SE = contrasts$SE[j], df = contrasts$df[j],
        lower_95 = contrasts$estimate[j] - qt(.975, contrasts$df[j]) * contrasts$SE[j],
        upper_95 = contrasts$estimate[j] + qt(.975, contrasts$df[j]) * contrasts$SE[j],
        t_ratio = contrasts$t.ratio[j], p_value = contrasts$p.value[j],
        n_presentations = nrow(model_data), n_flies = n_distinct(model_data$fly),
        singular_fit = singular
      )
      row_i <- row_i + 1
    }
    for (j in seq_len(nrow(movement_contrasts))) {
      response_rows[[row_i]] <- tibble(
        cohort = coh, outcome = outcome_name,
        movement_state = "moving_minus_not_moving",
        comparison = paste0(as.character(movement_contrasts$group[j]), "_movement_effect"),
        estimate = movement_contrasts$estimate[j], SE = movement_contrasts$SE[j],
        df = movement_contrasts$df[j],
        lower_95 = movement_contrasts$estimate[j] - qt(.975, movement_contrasts$df[j]) * movement_contrasts$SE[j],
        upper_95 = movement_contrasts$estimate[j] + qt(.975, movement_contrasts$df[j]) * movement_contrasts$SE[j],
        t_ratio = movement_contrasts$t.ratio[j], p_value = movement_contrasts$p.value[j],
        n_presentations = nrow(model_data), n_flies = n_distinct(model_data$fly),
        singular_fit = singular
      )
      row_i <- row_i + 1
    }
    coefficient_table <- coef(summary(model))
    interaction_term <- "groupexperimental:movement_statemoving"
    if (interaction_term %in% rownames(coefficient_table)) {
      z <- coefficient_table[interaction_term, ]
      response_rows[[row_i]] <- tibble(
        cohort = coh, outcome = outcome_name,
        movement_state = "moving_minus_not_moving",
        comparison = "group_by_movement_interaction",
        estimate = unname(z["Estimate"]), SE = unname(z["Std. Error"]),
        df = unname(z["df"]),
        lower_95 = estimate - qt(.975, df) * SE,
        upper_95 = estimate + qt(.975, df) * SE,
        t_ratio = unname(z["t value"]), p_value = unname(z["Pr(>|t|)"]),
        n_presentations = nrow(model_data), n_flies = n_distinct(model_data$fly),
        singular_fit = singular
      )
      row_i <- row_i + 1
    }
    diagnostic_rows[[diag_i]] <- tibble(
      cohort = coh, outcome = outcome_name, n_presentations = nrow(model_data),
      n_flies = n_distinct(model_data$fly), singular_fit = singular,
      fly_intercept_variance = as.data.frame(VarCorr(model))$vcov[1],
      residual_variance = sigma(model)^2
    )
    diag_i <- diag_i + 1
  }
}

response_tests <- bind_rows(response_rows) %>%
  group_by(cohort, movement_state, comparison) %>%
  mutate(bh_q_within_cohort_state_comparison = bh(p_value)) %>%
  ungroup()
write_csv(response_tests, file.path(out, "presentation_mixed_model_response_tests.csv"))
write_csv(bind_rows(diagnostic_rows), file.path(out, "presentation_model_diagnostics.csv"))

predictors <- c(
  high_pre_heading_change = "wiggle_heading_change",
  low_pre_straightness = "wiggle_low_straightness",
  low_pre_rho = "wiggle_low_rho"
)

baseline_rows <- list()
base_i <- 1
for (coh in unique(pulse$cohort)) {
  for (state in levels(pulse$movement_state)) {
    state_data <- pulse %>% filter(cohort == coh, movement_state == state)
    for (predictor_name in names(predictors)) {
      predictor <- predictors[[predictor_name]]
      for (outcome_name in names(outcomes)) {
        outcome <- outcomes[[outcome_name]]
        model_data <- state_data %>%
          select(fly, group, pulse_z, all_of(c(predictor, outcome))) %>% drop_na() %>%
          mutate(z_wiggle = as.numeric(scale(.data[[predictor]])))
        if (nrow(model_data) < 20 || n_distinct(model_data$fly) < 4) next
        formula <- as.formula(paste0(outcome, " ~ z_wiggle * group + pulse_z + (1 + z_wiggle | fly)"))
        model <- suppressWarnings(lmer(formula, data = model_data, REML = TRUE,
                                       control = lmerControl(check.conv.singular = "ignore")))
        trends <- as.data.frame(emtrends(model, ~ group, var = "z_wiggle"))
        trend_diff <- as.data.frame(contrast(emtrends(model, ~ group, var = "z_wiggle"),
                                             method = "revpairwise", adjust = "none"))
        for (j in seq_len(nrow(trends))) {
          baseline_rows[[base_i]] <- tibble(
            cohort = coh, movement_state = state, predictor = predictor_name,
            outcome = outcome_name, comparison = as.character(trends$group[j]),
            slope_per_sd_wiggliness = trends$z_wiggle.trend[j], SE = trends$SE[j],
            df = trends$df[j], t_ratio = trends$z_wiggle.trend[j] / trends$SE[j],
            p_value = 2 * pt(-abs(t_ratio), df = trends$df[j]),
            n_presentations = nrow(model_data), n_flies = n_distinct(model_data$fly),
            singular_fit = isSingular(model, tol = 1e-5)
          )
          base_i <- base_i + 1
        }
        baseline_rows[[base_i]] <- tibble(
          cohort = coh, movement_state = state, predictor = predictor_name,
          outcome = outcome_name, comparison = "experimental_minus_control_slope",
          slope_per_sd_wiggliness = trend_diff$estimate[1], SE = trend_diff$SE[1],
          df = trend_diff$df[1], t_ratio = trend_diff$t.ratio[1],
          p_value = trend_diff$p.value[1], n_presentations = nrow(model_data),
          n_flies = n_distinct(model_data$fly), singular_fit = isSingular(model, tol = 1e-5)
        )
        base_i <- base_i + 1
      }
    }
  }
}

baseline_tests <- bind_rows(baseline_rows) %>%
  group_by(cohort, movement_state, comparison) %>%
  mutate(bh_q_within_cohort_state_comparison = bh(p_value)) %>%
  ungroup()
write_csv(baseline_tests, file.path(out, "presentation_baseline_wiggliness_models.csv"))

plot_long <- pulse %>%
  select(cohort, group, movement_state, fly, pulse,
         response_less_heading_change, response_more_straight,
         response_more_consistent, response_faster_all) %>%
  pivot_longer(starts_with("response_"), names_to = "outcome", values_to = "change") %>%
  mutate(outcome = recode(outcome,
    response_less_heading_change = "Less heading change (deg)",
    response_more_straight = "More path straightness",
    response_more_consistent = "More heading consistency (rho)",
    response_faster_all = "Higher mean forward velocity"),
    panel = paste(cohort, movement_state, sep = " | "))

set.seed(20260807)
p1 <- ggplot(plot_long, aes(group, change, color = group)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(position = position_jitter(width = .18, height = 0), alpha = .08, size = .7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  facet_wrap(vars(panel, outcome), scales = "free_y", ncol = 4) +
  scale_color_manual(values = c(control = "#777777", experimental = "#D1495B")) +
  labs(title = "Individual optogenetic presentations",
       subtitle = "Each faint point is one presentation; black diamond is the presentation mean",
       x = NULL, y = "Hypothesis-oriented during - pre change") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out, "presentation_level_changes_by_movement.png"), p1,
       width = 15, height = 11, dpi = 200)

sequence <- plot_long %>%
  mutate(pulse_bin = cut(pulse, breaks = c(0, 10, 20, 30, 40, 50, 60, 70, Inf),
                         labels = c("1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71+"))) %>%
  group_by(cohort, group, movement_state, outcome, fly, pulse_bin) %>%
  summarise(change = mean(change, na.rm = TRUE), .groups = "drop") %>%
  group_by(cohort, group, movement_state, outcome, pulse_bin) %>%
  summarise(mean_change = mean(change, na.rm = TRUE),
            se = sd(change, na.rm = TRUE) / sqrt(sum(is.finite(change))), .groups = "drop") %>%
  mutate(panel = paste(cohort, movement_state, sep = " | "))

p2 <- ggplot(sequence, aes(pulse_bin, mean_change, color = group, group = group)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean_change - se, ymax = mean_change + se), width = .15) +
  facet_wrap(vars(panel, outcome), scales = "free_y", ncol = 4) +
  scale_color_manual(values = c(control = "#777777", experimental = "#D1495B")) +
  labs(title = "Response across repeated optogenetic presentations",
       subtitle = "Fly-weighted mean +/- SE within presentation-number bins",
       x = "Presentation number", y = "Hypothesis-oriented change", color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out, "presentation_sequence_by_movement.png"), p2,
       width = 15, height = 11, dpi = 200)

cat("Presentation-level analysis written to", out, "\n")
