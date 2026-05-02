### FOR THE ANALYSIS OF CSCHRIM BEHAVIOR DATA ###

### BLINDED CsChrimson ###
library(dplyr)
library(readxl)

savepath = "/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/Behavior_analysis/Behavior_Crunch/CsChrim_figures"

# LOAD IN AND LABEL DATA --------------------------------------------------

behavior_kin<- read.csv("/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/Behavior_analysis/Chrimson_Data_final.csv")

# from blinding sheet
# ves041_1_Chrim <- [2, 7, 16, 19, 21, 26, 53, 58, 65, 70, 72, 74, 77, 79, 81, 85]
# ves041_2_Chrim <- [1, 14, 18, 24, 29, 40, 44, 66, 73, 75, 80, 90, 100, 102, 105]
# ves041_1_gfp <- [3, 4, 5, 8, 12, 15, 25, 27, 69, 99, 102, 103]
# ves041_2_gfp <- [35, 52, 59, 63, 89, 93]
# empty_Chrim <- [6, 17, 20, 23, 31, 34, 39, 46, 50, 54, 57, 62, 84, 87, 88, 91, 96]


# right now segmented to be something > 5mm/s, but with 4mm/s we could also have data
behavior_kin <- behavior_kin %>%
  mutate(
    Fly = as.numeric(as.character(Fly)), # Ensure Fly is numeric for matching
    Genotype = case_when(
      Fly %in% c(2, 7, 16, 19, 21, 26, 53, 58, 65, 70, 72, 74, 77, 79, 81, 85) ~ "ves041_1_Chrim",
      Fly %in% c(1, 14, 18, 24, 29, 40, 44, 66, 73, 75, 80, 90, 100, 102, 105) ~ "ves041_2_Chrim",
      Fly %in% c(3, 4, 5, 8, 12, 15, 25, 27, 69, 99, 101, 103) ~ "ves041_1_gfp",
      Fly %in% c(35, 52, 59, 63, 89, 93) ~ "ves041_2_gfp",
      Fly %in% c(6, 17, 20, 23, 31, 34, 39, 46, 50, 54, 57, 62, 84, 87, 88, 91, 96) ~ "empty_Chrim",
      TRUE ~ NA_character_  # Or another default
    )
  )


# FORWARD SPEED -----------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

behavior_long <- behavior_kin %>%
  pivot_longer(
    cols = c("Pre.Opto.Forward.Velocity", "Opto.On.Forward.Velocity"),
    names_to = "Condition",
    values_to = "ForwardVelocity"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre.Opto.Forward.Velocity" = "Pre-Opto",
                       "Opto.On.Forward.Velocity" = "Opto On"),
    # Fix condition order
    Condition = factor(Condition, levels = c("Pre-Opto", "Opto On")),
    # Fix genotype order - replace with your desired order
    Genotype = factor(Genotype, levels = c("ves041_1_Chrim", "ves041_2_Chrim", "ves041_1_gfp", "ves041_2_gfp", "empty_Chrim"))
  )
p_fw_vel <- ggplot(behavior_long, aes(x = Condition, y = ForwardVelocity, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width = 0.6, alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width = 0.2) +
  geom_line(aes(group = Fly), color = "gray40", alpha = 0.5, linewidth = 0.4) +
  geom_point(aes(group = Fly), size = 1.5, alpha = 0.7, shape = 21, color = "black") +
  facet_wrap(~ Genotype, nrow = 1) +
  scale_fill_manual(values = c("Pre-Opto" = "blue", "Opto On" = "red")) +
  theme_minimal() +
  labs(
    title = "Forward Velocity Pre and During Cs Chrimson Opto by Genotype",
    x = "Condition",
    y = "Forward Velocity (mean ± SE)",
    fill = "Condition"
  )
p_fw_vel

ggsave(
  filename = "forward_velocity_plot.png",
  plot = p_fw_vel,
  path = savepath,
  width = 12,
  height = 9,
  dpi = 300
)

# ROTATIONAL SPEED -----------------------------------------------------------

behavior_long <- behavior_kin %>%
  pivot_longer(
    cols = c("Pre.Opto.Rotational.Speed", "Opto.On.Rotational.Speed"),
    names_to = "Condition",
    values_to = "RotationalSpeed"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre.Opto.Rotational.Speed" = "Pre-Opto",
                       "Opto.On.Rotational.Speed" = "Opto On"),
    # Fix condition order
    Condition = factor(Condition, levels = c("Pre-Opto", "Opto On")),
    # Fix genotype order - replace with your desired order
    Genotype = factor(Genotype, levels = c("ves041_1_Chrim", "ves041_2_Chrim", "ves041_1_gfp", "ves041_2_gfp", "empty_Chrim"))
  )
p_rot_sp <- ggplot(behavior_long, aes(x = Condition, y = RotationalSpeed, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width = 0.6, alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width = 0.2) +
  geom_line(aes(group = Fly), color = "gray40", alpha = 0.5, linewidth = 0.4) +
  geom_point(aes(group = Fly), size = 1.5, alpha = 0.7, shape = 21, color = "black") +
  facet_wrap(~ Genotype, nrow = 1) +
  scale_fill_manual(values = c("Pre-Opto" = "blue", "Opto On" = "red")) +
  theme_minimal() +
  labs(
    title = "Rotational Speed Pre and During Cs Chrimson Opto by Genotype",
    x = "Condition",
    y = "Rotational Speed (mean ± SE)",
    fill = "Condition"
  )
p_rot_sp

ggsave(
  filename = "rotational_speed_plot.png",
  plot = p_rot_sp,
  path = savepath,
  width = 12,
  height = 9,
  dpi = 300
)




# FORWARD ACCELERATION -----------------------------------------------------------
behavior_long <- behavior_kin %>%
  pivot_longer(
    cols = c("Pre.Opto.Forward.Acceleration", "Opto.On.Forward.Acceleration"),
    names_to = "Condition",
    values_to = "ForwardAcceleration"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre.Opto.Forward.Acceleration" = "Pre-Opto",
                       "Opto.On.Forward.Acceleration" = "Opto On"),
    # Fix condition order
    Condition = factor(Condition, levels = c("Pre-Opto", "Opto On")),
    # Fix genotype order - replace with your desired order
    Genotype = factor(Genotype, levels = c("ves041_1_Chrim", "ves041_2_Chrim", "ves041_1_gfp", "ves041_2_gfp", "empty_Chrim"))
  )
p_fwd_acc <- ggplot(behavior_long, aes(x = Condition, y = ForwardAcceleration, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width = 0.6, alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width = 0.2) +
  geom_line(aes(group = Fly), color = "gray40", alpha = 0.5, linewidth = 0.4) +
  geom_point(aes(group = Fly), size = 1.5, alpha = 0.7, shape = 21, color = "black") +
  facet_wrap(~ Genotype, nrow = 1) +
  scale_fill_manual(values = c("Pre-Opto" = "blue", "Opto On" = "red")) +
  theme_minimal() +
  labs(
    title = "Forward Acceleration Pre and During Cs Chrimson Opto by Genotype",
    x = "Condition",
    y = "Forward Acceleration (mean ± SE)",
    fill = "Condition"
  )
p_fwd_acc

ggsave(
  filename = "fwd_accel_plot.png",
  plot = p_fwd_acc,
  path = savepath,
  width = 12,
  height = 9,
  dpi = 300
)


# ROTATIONAL ACCELERATION -------------------------------------------------
behavior_long <- behavior_kin %>%
  pivot_longer(
    cols = c("Pre.Opto.Rotational.Acceleration", "Opto.On.Rotational.Acceleration"),
    names_to = "Condition",
    values_to = "RotationalAcceleration"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre.Opto.Rotational.Acceleration" = "Pre-Opto",
                       "Opto.On.Rotational.Acceleration" = "Opto On"),
    # Fix condition order
    Condition = factor(Condition, levels = c("Pre-Opto", "Opto On")),
    # Fix genotype order - replace with your desired order
    Genotype = factor(Genotype, levels = c("ves041_1_Chrim", "ves041_2_Chrim", "ves041_1_gfp", "ves041_2_gfp", "empty_Chrim"))
  )
p_rot_acc <- ggplot(behavior_long, aes(x = Condition, y = RotationalAcceleration, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width = 0.6, alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width = 0.2) +
  geom_line(aes(group = Fly), color = "gray40", alpha = 0.5, linewidth = 0.4) +
  geom_point(aes(group = Fly), size = 1.5, alpha = 0.7, shape = 21, color = "black") +
  facet_wrap(~ Genotype, nrow = 1) +
  scale_fill_manual(values = c("Pre-Opto" = "blue", "Opto On" = "red")) +
  theme_minimal() +
  labs(
    title = "Rotational Acceleration Pre and During Cs Chrimson Opto by Genotype",
    x = "Condition",
    y = "Rotational Acceleration (mean ± SE)",
    fill = "Condition"
  )

p_rot_acc

ggsave(
  filename = "rot_accel_plot.png",
  plot = p_rot_acc,
  path = savepath,
  width = 12,
  height = 9,
  dpi = 300
)

# SANDBOX -----------------------------------------------------------------
# rot speed stats

library(lme4)
library(lmerTest)
library(emmeans)

# Fit model with Condition * genotype interaction, Fly as random effect
model <- lmer(RotationalSpeed ~ Condition * Genotype + (1|Fly), data = behavior_long)
summary(model)

# Get estimated marginal means
emm <- emmeans(model, ~ Condition * Genotype)

# Option 1: Compare opto response (interaction) between all genotype pairs
# This directly asks "do genotypes differ in how they respond to opto?"
contrast(emm, interaction = "pairwise", by = "Condition", adjust = "bonferroni")

# Option 2: Compare genotypes within each condition separately
emmeans(model, pairwise ~ Genotype | Condition, adjust = "bonferroni")

# Option 3: Compare conditions within each genotype (like paired t-tests)
emmeans(model, pairwise ~ Condition | Genotype, adjust = "bonferroni")




# probably not this
library(rstatix)

# Compute delta per fly
delta_data <- behavior_kin %>%
  mutate(
    delta = Opto.On.Rotational.Speed - Pre.Opto.Rotational.Speed,
    genotype = factor(Genotype, levels = c("ves041_1_Chrim", "ves041_2_Chrim", 
                                           "ves041_1_gfp", "ves041_2_gfp", "empty_Chrim"))
  )

# One-way ANOVA
delta_data %>% anova_test(delta ~ genotype)

# Post-hoc pairwise comparisons
delta_tukey <- delta_data %>%
  tukey_hsd(delta ~ genotype)
print(delta_tukey)
