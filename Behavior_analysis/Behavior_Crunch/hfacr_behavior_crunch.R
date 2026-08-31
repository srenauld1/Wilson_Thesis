### FOR THE ANALYSIS OF BEHAVIOR DATA ###

### BLINDED RUBYACR ###
library(dplyr)
library(readxl)

savepath = "/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/Behavior_analysis/Behavior_Crunch/ACR_figures"


# LOAD IN AND LABEL DATA --------------------------------------------------

behavior_kin<- read.csv("/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/Behavior_analysis/rubyACR pilot data(blinded).csv")

# from blinding sheet
# ves041_1_acr <- [3, 5, 8, 15, 22]
# ves041_2_acr <- [2, 11, 17, 38, 40, 42]
# ves041_1_gfp <- [1, 13, 16, 19, 27, 30]
# ves041_2_gfp <- [43, 45]
# empty_acr <- [4, 9, 26, 31, 35, 37]

behavior_kin <- behavior_kin %>%
  mutate(
    Fly = as.numeric(as.character(Fly)), # Ensure Fly is numeric for matching
    genotype = case_when(
      Fly %in% c(3, 5, 8, 15, 22, 32) ~ "ves041_1_acr",
      Fly %in% c(2, 11, 17, 38, 40, 42, 47) ~ "ves041_2_acr",
      Fly %in% c(1, 13, 16, 19, 27, 30) ~ "ves041_1_gfp",
      Fly %in% c(43, 45, 48, 50, 51) ~ "ves041_2_gfp",
      Fly %in% c(4, 9, 26, 31, 35, 37) ~ "empty_acr",
      TRUE ~ NA_character_  # Or another default
    )
  )


# FORWARD SPEED -----------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

behavior_long <- behavior_kin %>%
  pivot_longer(
    cols = c("Pre.Opto.Forward.Velocity.", "Opto.On.Forward.Velocity"),
    names_to = "Condition",
    values_to = "ForwardVelocity"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre.Opto.Forward.Velocity." = "Pre-Opto",
                       "Opto.On.Forward.Velocity" = "Opto On"),
    # Fix condition order
    Condition = factor(Condition, levels = c("Pre-Opto", "Opto On")),
    # Fix genotype order - replace with your desired order
    genotype = factor(genotype, levels = c("ves041_1_acr", "ves041_2_acr", "ves041_1_gfp", "ves041_2_gfp", "empty_acr"))
  )
p_fw_vel <- ggplot(behavior_long, aes(x = Condition, y = ForwardVelocity, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width = 0.6, alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width = 0.2) +
  geom_line(aes(group = Fly), color = "gray40", alpha = 0.5, linewidth = 0.4) +
  geom_point(aes(group = Fly), size = 1.5, alpha = 0.7, shape = 21, color = "black") +
  facet_wrap(~ genotype, nrow = 1) +
  scale_fill_manual(values = c("Pre-Opto" = "blue", "Opto On" = "red")) +
  theme_minimal() +
  labs(
    title = "Forward Velocity Pre and During RubyACR Opto by genotype",
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
    genotype = factor(genotype, levels = c("ves041_1_acr", "ves041_2_acr", "ves041_1_gfp", "ves041_2_gfp", "empty_acr"))
  )
p_rot_sp <- ggplot(behavior_long, aes(x = Condition, y = RotationalSpeed, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width = 0.6, alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width = 0.2) +
  geom_line(aes(group = Fly), color = "gray40", alpha = 0.5, linewidth = 0.4) +
  geom_point(aes(group = Fly), size = 1.5, alpha = 0.7, shape = 21, color = "black") +
  facet_wrap(~ genotype, nrow = 1) +
  scale_fill_manual(values = c("Pre-Opto" = "blue", "Opto On" = "red")) +
  theme_minimal() +
  labs(
    title = "Rotational Speed Pre and During RubyACR Opto by genotype",
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
    genotype = factor(genotype, levels = c("ves041_1_acr", "ves041_2_acr", "ves041_1_gfp", "ves041_2_gfp", "empty_acr"))
  )
p_fwd_acc <- ggplot(behavior_long, aes(x = Condition, y = ForwardAcceleration, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width = 0.6, alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width = 0.2) +
  geom_line(aes(group = Fly), color = "gray40", alpha = 0.5, linewidth = 0.4) +
  geom_point(aes(group = Fly), size = 1.5, alpha = 0.7, shape = 21, color = "black") +
  facet_wrap(~ genotype, nrow = 1) +
  scale_fill_manual(values = c("Pre-Opto" = "blue", "Opto On" = "red")) +
  theme_minimal() +
  labs(
    title = "Forward Acceleration Pre and During RubyACR Opto by genotype",
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
    cols = c("Pre.Opto.Rotational.Acceleration.", "Opto.On.Rotational.Acceleration"),
    names_to = "Condition",
    values_to = "RotationalAcceleration"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre.Opto.Rotational.Acceleration." = "Pre-Opto",
                       "Opto.On.Rotational.Acceleration" = "Opto On"),
    # Fix condition order
    Condition = factor(Condition, levels = c("Pre-Opto", "Opto On")),
    # Fix genotype order - replace with your desired order
    genotype = factor(genotype, levels = c("ves041_1_acr", "ves041_2_acr", "ves041_1_gfp", "ves041_2_gfp", "empty_acr"))
  )
p_rot_acc <- ggplot(behavior_long, aes(x = Condition, y = RotationalAcceleration, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width = 0.6, alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width = 0.2) +
  geom_line(aes(group = Fly), color = "gray40", alpha = 0.5, linewidth = 0.4) +
  geom_point(aes(group = Fly), size = 1.5, alpha = 0.7, shape = 21, color = "black") +
  facet_wrap(~ genotype, nrow = 1) +
  scale_fill_manual(values = c("Pre-Opto" = "blue", "Opto On" = "red")) +
  theme_minimal() +
  labs(
    title = "Rotational Acceleration Pre and During RubyACR Opto by genotype",
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

# Un-blinded Data ----------------------------------------------------------
behavior_kin_unblinded<- read.csv("/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/Behavior_analysis/rubyACR pilot data(2 second stim).csv")

# Gather data into long format
behavior_long_unblinded <- behavior_kin_unblinded %>%
  pivot_longer(
    cols = c("Pre.Opto.Forward.Velocity", "Opto.On..Forward.Velocity"),
    names_to = "Condition",
    values_to = "ForwardVelocity"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre.Opto.Forward.Velocity" = "Pre-Opto",
                       "Opto.On..Forward.Velocity" = "Opto On")
  )

p_fw_vel <- ggplot(behavior_long_unblinded, aes(x = Condition, y = ForwardVelocity, fill = Condition)) +
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
    title = "UNBLINDED PILOT: Forward Speed Pre and During Opto by genotype",
    x = "Genotype",
    y = "Forward Velocity (mean ± SE)",
    fill = "Condition"
  )
p_fw_vel

ggsave(
  filename = "unblinded_forward_velocity_plot.png",
  plot = p_fw_vel,
  path = savepath,
  width = 12,
  height = 9,
  dpi = 300
)


# Option 1: Compare opto response (interaction) between all genotype pairs
# This directly asks "do genotypes differ in how they respond to opto?"
contrast(emm, interaction = "pairwise", by = "Condition", adjust = "bonferroni")


# rotation!! unblinded!!

# Gather data into long format
behavior_long_unblinded <- behavior_kin_unblinded %>%
  pivot_longer(
    cols = c("Pre.Opto.Rotational.Speed", "Opto.On.Rotational.Speed"),
    names_to = "Condition",
    values_to = "RotationalSpeed"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre.Opto.Rotational.Speed" = "Pre-Opto",
                       "Opto.On.Rotational.Speed" = "Opto On")
  )

p_rot_vel <- ggplot(behavior_long_unblinded, aes(x = Condition, y = RotationalSpeed, fill = Condition)) +
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
    title = "UNBLINDED PILOT: Rotational Speed Pre and During Opto by genotype",
    x = "Genotype",
    y = "Rotational Speed (mean ± SE)",
    fill = "Condition"
  )
p_rot_vel



# Pooled Data -------------------------------------------------------------

pooled_behavior <- read_xlsx("/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/Behavior_analysis/pooled_ves041behavior.xlsx")



# Gather data into long format
behavior_long <- pooled_behavior %>%
  pivot_longer(
    cols = c("Pre Opto Forward Velocity", "Opto On Forward Velocity"),
    names_to = "Condition",
    values_to = "ForwardVelocity"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre Opto Forward Velocity" = "Pre-Opto",
                       "Opto On Forward Velocity" = "Opto On")
  )

# Bar plot: mean with error bars (SE) for each genotype and condition
ggplot(behavior_long, aes(x = Genotype, y = ForwardVelocity, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width=0.6, alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width=0.2) +
  theme_minimal() +
  labs(
    title = "Forward Speed Pre and During Opto by genotype",
    x = "genotype",
    y = "Forward Velocity (mean ± SE)",
    fill = "Condition"
  )











