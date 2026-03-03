### FOR THE ANALYSIS OF BEHAVIOR DATA ###

### BLINDED RUBYACR ###
library(dplyr)

# LOAD IN AND LABEL DATA --------------------------------------------------

behavior_kin<- read.csv("/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/Behavior_analysis/rubyACR pilot data(blinded).csv")

# from blinding sheet
# ves041_1_acr <- [3, 5, 8, 15, 22]
# ves041_2_acr <- [2, 11, 17, 38, 40, 42]
# ves041_1_gfp <- [1, 13, 16, 19, 27, 30]
# ves041_2_gfp <- [43, 45]
# empty_acr <- [4, 9, 26, 31, 35, 37]

# behavior_kin <- behavior_kin %>%
#   mutate(
#     Fly = as.numeric(as.character(Fly)), # Ensure Fly is numeric for matching
#     genotype = case_when(
#       Fly %in% c(3, 5, 8, 15, 22, 32) ~ "ves041_1_acr",
#       Fly %in% c(2, 11, 17, 38, 40, 42) ~ "ves041_2_acr",
#       Fly %in% c(1, 13, 16, 19, 27, 30) ~ "ves041_1_gfp",
#       Fly %in% c(43, 45) ~ "ves041_2_gfp",
#       Fly %in% c(4, 9, 26, 31, 35, 37) ~ "empty_acr",
#       TRUE ~ NA_character_  # Or another default
#     )
#   )


# FORWARD SPEED -----------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# Gather data into long format
behavior_long <- behavior_kin %>%
  pivot_longer(
    cols = c("Pre.Opto.Forward.Velocity.", "Opto.On.Forward.Velocity"),
    names_to = "Condition",
    values_to = "ForwardVelocity"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre.Opto.Forward.Velocity." = "Pre-Opto",
                       "Opto.On.Forward.Velocity" = "Opto On")
  )

# Bar plot: mean with error bars (SE) for each genotype and condition
ggplot(behavior_long, aes(x = Genotype, y = ForwardVelocity, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width=0.6, alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width=0.2) +
  theme_minimal() +
  labs(
    title = "Forward Speed Pre and During Opto by Genotype",
    x = "Genotype",
    y = "Forward Velocity (mean ± SE)",
    fill = "Condition"
  )

# ROTATIONAL SPEED -----------------------------------------------------------

# Gather data into long format
behavior_long <- behavior_kin %>%
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

# Bar plot: mean with error bars (SE) for each genotype and condition
ggplot(behavior_long, aes(x = Genotype, y = RotationalSpeed, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width=0.6, alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width=0.2) +
  theme_minimal() +
  labs(
    title = "Rotational Speed Pre and During Opto by Genotype",
    x = "Genotype",
    y = "Forward Velocity (mean ± SE)",
    fill = "Condition"
  )


# FORWARD ACCELERATION -----------------------------------------------------------

# Gather data into long format
behavior_long <- behavior_kin %>%
  pivot_longer(
    cols = c("Pre.Opto.Forward.Acceleration", "Opto.On.Forward.Acceleration"),
    names_to = "Condition",
    values_to = "ForwardAcceleration"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Pre.Opto.Forward.Acceleration" = "Pre-Opto",
                       "Opto.On.Forward.Acceleration" = "Opto On")
  )

# Bar plot: mean with error bars (SE) for each genotype and condition
ggplot(behavior_long, aes(x = Genotype, y = ForwardAcceleration, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width=0.6, alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width=0.2) +
  theme_minimal() +
  labs(
    title = "Forward Acceleration Pre and During Opto by Genotype",
    x = "Genotype",
    y = "Forward Velocity (mean ± SE)",
    fill = "Condition"
  )



# ## percent change
# behavior_pct <- behavior_kin %>%
#   mutate(
#     pct_change = (Opto.On.Forward.Velocity - Pre.Opto.Forward.Velocity.) 
#   )
# # Now you have one row per fly, with its genotype and its % change
# ggplot(behavior_pct, aes(x = Genotype, y = pct_change, fill = Genotype)) +
#   stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6, alpha = 0.8) +
#   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
#   geom_jitter(width = 0.15, alpha = 0.5) +
#   theme_minimal() +
#   labs(
#     title = "Change in Forward Velocity Pre vs Opto by Genotype",
#     x = "Genotype",
#     y = "% Change in Forward Velocity",
#     fill = "Genotype"
#   )

# Unblinded Data ----------------------------------------------------------
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

# Bar plot: mean with error bars (SE) for each genotype and condition
ggplot(behavior_long_unblinded, aes(x = Genotype, y = ForwardVelocity, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width=0.6, alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width=0.2) +
  theme_minimal() +
  labs(
    title = "UNBLINDED PILOT: Forward Speed Pre and During Opto by Genotype",
    x = "Genotype",
    y = "Forward Velocity (mean ± SE)",
    fill = "Condition"
  )

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

# Bar plot: mean with error bars (SE) for each genotype and condition
ggplot(behavior_long_unblinded, aes(x = Genotype, y = RotationalSpeed, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.7), color = "black", width=0.6, alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.7), width=0.2) +
  theme_minimal() +
  labs(
    title = "UNBLINDED PILOT: Rotational Speed Pre and During Opto by Genotype",
    x = "Genotype",
    y = "Rotational Speed (mean ± SE)",
    fill = "Condition"
  )



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
    title = "Forward Speed Pre and During Opto by Genotype",
    x = "Genotype",
    y = "Forward Velocity (mean ± SE)",
    fill = "Condition"
  )











