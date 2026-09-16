library(dplyr)
library(ggplot2)
library(tidyverse)

# loading main dataset
merged_manybabies <- read.csv("final_merged_manybabies.csv")

## some inspection
str(merged_manybabies) # checking variable formats
install.packages("skimr") # overview of dataset
library(skimr)
skim(merged_manybabies_final)

install.packages("showtext") # for fonts
library(showtext)
showtext.auto()

# there are some NA values for age_group_new, checking them 
missing_bastards_df = merged_manybabies_final[is.na(merged_manybabies_final$age_group_new), ]
unique(missing_bastards_df$unique_subject_identifier) # these are 2 babies from 
# MB1, one below 3mos and another over 21mos
sum(is.na(merged_manybabies_final$age_group_new))

## outlier removal 
# highest LTs for 9-12 mos (because this group showed a spike in the plot)
merged_manybabies %>%
  arrange(desc(looking_time)) %>%
  slice_head(n = 50) %>%
  select(
    project_id,
    age_group_new,
    unique_subject_identifier,
    looking_time
  )

# removing outlier datapoints from MB4_babylabmpib_mb4_02, MB4_babylabmpib_mb4_15, 
# MB3_UPFbabylab_MB3_59637, MB3_UPFbabylab_MB3_59649
merged_manybabies_final <- merged_manybabies %>%       # from original dataset
  filter(!looking_time %in% c(41852, 41518, 751, 472))

unique(merged_manybabies_final$age_group_new)

# extracting this new csv file
# this is the final dataset
write.csv(merged_manybabies_final, "C:/Arya/UvA/Thesis/merged_manybabies_final.csv", row.names = FALSE)
unique(merged_manybabies_final$age_group_new)

## some plots for visualisation

# 1. mean LTs grouped my method 
ggplot(merged_manybabies_final,
       aes(x = method, y = looking_time, fill = method)) +
  stat_summary(fun = mean, geom = "col", width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  labs(
    x = NULL,
    y = "Mean LTs",
    title = "Mean Infant LTs by Method"
  ) + 
  scale_fill_brewer(palette = "Reds") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      family = "Times New Roman")
  )

# 2. individual infant LTs (not so informative due to large N)
ggplot(merged_manybabies_final, aes(x = age_group_new, y = looking_time)) +
  geom_violin(width = 0.2, alpha = 0.4) +
  labs(
    x = "Age Group",
    y = "Looking Time",
    title = "Distribution of Infant Looking Times by Age Group"
  ) +
  theme_minimal() + ylim(0, 20)  
  # LTs as a function of age , distribution/density plot

# ordering age_group_new
merged_manybabies_final$age_group_new <- factor(
  merged_manybabies_final$age_group_new,
  levels = c(
    "3–6 months",
    "6–9 months",
    "9–12 months",
    "12–15 months",
    "15–18 months",
    "18-21 months",
  ),
  ordered = TRUE
)

# 3. mean LTs by age group
ggplot(merged_manybabies_final,
       aes(x = age_group_new, y = looking_time, fill = age_group_new)) +
  stat_summary(fun = mean, geom = "col", width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  labs(
    x = NULL,
    y = "Mean Looking Time(s)",
    title = "Mean Infant LTs by Age Group"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# line graph age without method
ggplot(merged_manybabies_final, aes(x = age_group_new, y = looking_time, group = 1)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point") +
  labs(
    x = "Age Group",
    y = "Mean Looking Time",
    title = "Developmental Changes in Looking Time"
  ) +
  theme_minimal() 

# line graph age without method seperated by method
ggplot(merged_manybabies_final, aes(x = age_group_new, y = looking_time, group = method, colour = method)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point") +
  labs(
    x = "Age Group",
    y = "Mean Looking Time",
    colour = "Method",
    title = "Developmental Changes in Looking Time"
  ) +
  theme_minimal() 

# with error bars
ggplot(merged_manybabies_final, aes(x = age_group_new, y = looking_time, group = 1)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  labs(
    x = "Age Group",
    y = "Mean Looking Time"
  ) +
  theme_minimal()

print(merged_manybabies_final$trial_type)
head(merged_manybabies_final %>% dplyr::select(trial_type, trial_number))


## habituation curves
# filtering TRAIN and fam trials
hab_data <- merged_manybabies_final %>%         # habituation dataset
  dplyr::filter(habituation == "TRUE",
                trial_type != "fam")
# hab curves by trials
ggplot(hab_data,
       aes(x = trial_number, y = looking_time, group = 1)) +
  stat_summary(fun = mean, geom = "line", color = "black", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", color = "black", size = 2) +
  labs(
    x = "Habituation/Familiarisation Trials",
    y = "Looking Times (s)"
  ) +
  theme_classic(base_size = 12)


# habuation curves by age group
# ordering the ages 
hab_data$age_group_new <- factor(
  hab_data$age_group_new,
  levels = c(
    "3–6 months",
    "6–9 months",
    "9–12 months",
    "12–15 months",
    "15–18 months",
    "18–21 months"
  ),
  ordered = TRUE
)

unique(hab_data$age_group_new)

# 1. LTs across hab trials 
# removing age group 18-21 (weird spike, for visualisation purposes only)
ggplot(
  hab_data %>% 
         dplyr::filter(
           !is.na(age_group_new),
           age_group_new != "18–21 months"),
       aes(x = trial_number, y = looking_time,
           color = age_group_new, group = age_group_new)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  scale_color_brewer(palette = "Reds") +
  labs(
    x = "Habituation/Familiarisation Trials",
    y = "Looking Times",
    color = "Age Group"
  ) +
  theme_minimal(base_size = 12)

# 2. LTs across hab trials with all age groups
ggplot(
  hab_data %>% filter(
      !is.na(age_group_new)),
  aes(x = trial_number, y = looking_time,
      color = age_group_new, group = age_group_new)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  scale_color_brewer(palette = "Reds") +
  labs(
    x = "Habituation/Familiarisation Trials",
    y = "Looking Times",
    color = "Age Group"
  ) +
  theme_minimal(base_size = 12)

# modelling curves with simulated data 
## dual process theory
library(ggplot2)
library(patchwork)
library(tidyr)

# Trials
trial <- seq(0, 20, length.out = 300)

# ---- Panel A: habituation dominates ----
S1 <- 0.8 * (1 - exp(-trial/2))
H1 <- -2.5 * (1 - exp(-trial/4))
Data1 <- S1 + H1

df1 <- data.frame(trial, Sensitisation = S1, Habituation = H1, Data = Data1) |>
  pivot_longer(-trial, names_to = "Process", values_to = "value")

p1 <- ggplot(df1, aes(trial, value, colour = Process, linetype = Process)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(
    "Sensitisation" = "red",
    "Habituation" = "blue",
    "Data" = "black"
  )) +
  scale_linetype_manual(values = c(
    "Sensitisation" = "dashed",
    "Habituation" = "dashed",
    "Data" = "solid"
  )) +
  scale_y_continuous(limits = c(-2, 4), expand = c(0,0)) +
  labs(
    x = "Trials",
    y = "Change in looking time(s)",
    colour = "",
    linetype = ""
  ) +
  theme_classic()

# ---- Panel B: sensitisation initially dominates ----
S2 <- 2.5 * trial * exp(-trial/4)
H2 <- -1.8 * (1 - exp(-trial/4))
Data2 <- S2 + H2

df2 <- data.frame(trial, Sensitisation = S2, Habituation = H2, Data = Data2) |>
  pivot_longer(-trial, names_to = "Process", values_to = "value")

p2 <- ggplot(df2, aes(trial, value, colour = Process, linetype = Process)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(
    "Sensitisation" = "red",
    "Habituation" = "blue",
    "Data" = "black"
  )) +
  scale_linetype_manual(values = c(
    "Sensitisation" = "dashed",
    "Habituation" = "dashed",
    "Data" = "solid"
  )) +
  scale_y_continuous(limits = c(-2, 4), expand = c(0,0)) +
  labs(
    x = "Trials",
    y = "Change in looking time(s)",
    colour = "",
    linetype = ""
  ) +
  theme_classic()

# Combine panels
p1 + p2

