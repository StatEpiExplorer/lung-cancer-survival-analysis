# Install and load necessary packages
install.packages("survminer")
library(survminer)

install.packages("survival")
library(survival)

install.packages("dplyr")
library(dplyr)

install.packages("ggplot2")
library(ggplot2)

# Create an output folder if it doesn't exist
if (!dir.exists("outputs")) {
  dir.create("outputs")
}

# Loading the dataset
df <- read.csv("C:/Users/Dr.Aditi Sharada/OneDrive/Documents/BIOSTAT COMPUTING/CANCER PROJECT/CANCER/DATA/cancer.csv")
head(df)
str(df)

# Checking for null values
sink("outputs/null_value_summary.txt")
print("Null Value Summary:")
print(colSums(is.na(df)))
sink()

# Dropping unnecessary columns
df <- df[, !names(df) %in% c("Unnamed..0", "inst")]
head(df)

# Convert status column to binary (1 for death, 0 for censored)
df <- df %>% mutate(dead = ifelse(status == 1, 0, 1))

# Define survival object
surv_obj <- Surv(time = df$time, event = df$dead)

# Fit Kaplan-Meier estimator
km_fit <- survfit(Surv(time, dead) ~ 1, data = df)

# Save event table to a text file
sink("outputs/km_event_table.txt")
print("Kaplan-Meier Event Table:")
print(summary(km_fit)$table)
sink()

# Survival probability at specific time points
surv_prob_0 <- km_fit$surv[1]
surv_prob_11 <- summary(km_fit, times = 11)$surv
surv_prob_13 <- summary(km_fit, times = 13)$surv

# Save survival probabilities to a text file
sink("outputs/survival_probabilities.txt")
print("Survival Probabilities at Specific Time Points:")
print(paste("Survival at time t = 0:", surv_prob_0))
print(paste("Survival at time t = 11:", round(surv_prob_11, 3)))
print(paste("Survival at time t = 13:", round(surv_prob_13, 3)))
sink()

# Kaplan-Meier survival plot
ggsave("outputs/kaplan_meier_plot.png",
       plot = ggsurvplot(km_fit, data = df,
                         conf.int = FALSE,
                         title = "Kaplan-Meier Estimation",
                         xlab = "Number of days",
                         ylab = "Probability of survival")$plot,
       width = 8, height = 6, units = "in")

# Median survival time
median_surv_time <- summary(km_fit)$table["median"]
if (is.null(median_surv_time)) {
  median_surv_time <- km_fit$time[which.min(abs(km_fit$surv - 0.5))]
}

# Save median survival time to a text file
sink("outputs/median_survival_time.txt")
print("Median Survival Time:")
print(paste("The median survival time:", median_surv_time, "days."))
sink()

# Confidence interval for survival function plot
ggsave("outputs/survival_with_ci.png",
       plot = plot(km_fit, conf.int = TRUE, main = "Survival Function with Confidence Interval",
                   xlab = "Number of days", ylab = "Survival Probability"),
       width = 8, height = 6, units = "in")

# Cumulative hazard function (Nelson-Aalen estimator) using basehaz()
cox_model <- coxph(Surv(time, dead) ~ 1, data = df)
naf_fit <- basehaz(cox_model, centered = FALSE)

# Create data frame for plotting the cumulative hazard
hazard_df <- data.frame(time = naf_fit$time, cumhaz = naf_fit$hazard)

# Plot cumulative hazard
hazard_plot <- ggplot(hazard_df, aes(x = time, y = cumhaz)) +
  geom_step(color = "blue") +
  labs(title = "Cumulative Hazard Function (Nelson-Aalen Estimator)",
       x = "Number of Days",
       y = "Cumulative Hazard")
ggsave("outputs/cumulative_hazard_plot.png", plot = hazard_plot, width = 8, height = 6, units = "in")

# Load necessary libraries (redundant here but kept for clarity if this section is run independently)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# Create copy of dataframe and rename columns
df2 <- df %>%
  rename(
    meal_cal = meal.cal,
    wt_loss = wt.loss,
    ph_karno = ph.karno,
    ph_ecog = ph.ecog,
    pat_karno = pat_karno
  )

# Survival analysis by sex (separate plots)
males <- df2 %>% filter(sex == 1)
females <- df2 %>% filter(sex == 2)

kmf_males <- survfit(Surv(time, dead) ~ 1, data = males)
kmf_females <- survfit(Surv(time, dead) ~ 1, data = females)

ggsave("outputs/survival_males.png",
       plot = ggsurvplot(kmf_males, data = males, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival Probability of Males", legend = "none")$plot,
       width = 8, height = 6, units = "in")
ggsave("outputs/survival_females.png",
       plot = ggsurvplot(kmf_females, data = females, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival Probability of Females", legend = "none")$plot,
       width = 8, height = 6, units = "in")

# Age category analysis (separate plots)
df_age_cats <- df2 %>%
  mutate(age_cat = ifelse(age >= 70, 1, 0))

old <- df_age_cats %>% filter(age_cat == 1)
young <- df_age_cats %>% filter(age_cat == 0)

kmf_old <- survfit(Surv(time, dead) ~ 1, data = old)
kmf_young <- survfit(Surv(time, dead) ~ 1, data = young)

ggsave("outputs/survival_age_geq_70.png",
       plot = ggsurvplot(kmf_old, data = old, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival Probability (Age ≥ 70)", legend = "none")$plot,
       width = 8, height = 6, units = "in")
ggsave("outputs/survival_age_lt_70.png",
       plot = ggsurvplot(kmf_young, data = young, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival Probability (Age < 70)", legend = "none")$plot,
       width = 8, height = 6, units = "in")

# Weight loss analysis (separate plots)
df_wt_loss <- df2
mean_wt_loss <- mean(df_wt_loss$wt_loss, na.rm = TRUE)

sink("outputs/mean_weight_loss.txt")
print("Mean Weight Loss:")
print(paste("Mean weight loss:", mean_wt_loss))
sink()

df_wt_loss <- df_wt_loss %>%
  mutate(wt_loss_cat = ifelse(wt_loss > mean_wt_loss, 1, 0))

above_avg <- df_wt_loss %>% filter(wt_loss_cat == 1)
below_avg <- df_wt_loss %>% filter(wt_loss_cat == 0)

kmf_above_avg_wt_loss <- survfit(Surv(time, dead) ~ 1, data = above_avg)
kmf_below_avg_wt_loss <- survfit(Surv(time, dead) ~ 1, data = below_avg)

ggsave("outputs/survival_wt_loss_above_mean.png",
       plot = ggsurvplot(kmf_above_avg_wt_loss, data = above_avg, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival (Weight Loss > Mean)", legend = "none")$plot,
       width = 8, height = 6, units = "in")
ggsave("outputs/survival_wt_loss_below_mean.png",
       plot = ggsurvplot(kmf_below_avg_wt_loss, data = below_avg, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival (Weight Loss ≤ Mean)", legend = "none")$plot,
       width = 8, height = 6, units = "in")

# Patient Karnofsky score analysis (separate plots)
df_pat_karno <- df2 %>%
  mutate(pat_karno_cat = ifelse(pat_karno >= 80, 1, 0))

healthy_pat <- df_pat_karno %>% filter(pat_karno_cat == 1)
sick_pat <- df_pat_karno %>% filter(pat_karno_cat == 0)

kmf_pat_karno_healthy <- survfit(Surv(time, dead) ~ 1, data = healthy_pat)
kmf_pat_karno_sick <- survfit(Surv(time, dead) ~ 1, data = sick_pat)

ggsave("outputs/survival_pat_karno_geq_80.png",
       plot = ggsurvplot(kmf_pat_karno_healthy, data = healthy_pat, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival (Pat-Karno ≥ 80)", legend = "none")$plot,
       width = 8, height = 6, units = "in")
ggsave("outputs/survival_pat_karno_lt_80.png",
       plot = ggsurvplot(kmf_pat_karno_sick, data = sick_pat, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival (Pat-Karno < 80)", legend = "none")$plot,
       width = 8, height = 6, units = "in")

# Physician Karnofsky score analysis (separate plots)
df_ph_karno <- df2 %>%
  mutate(ph_karno_cat = ifelse(ph_karno >= 80, 1, 0))

healthy_ph <- df_ph_karno %>% filter(ph_karno_cat == 1)
sick_ph <- df_ph_karno %>% filter(ph_karno_cat == 0)

kmf_ph_karno_healthy <- survfit(Surv(time, dead) ~ 1, data = healthy_ph)
kmf_ph_karno_sick <- survfit(Surv(time, dead) ~ 1, data = sick_ph)

ggsave("outputs/survival_ph_karno_geq_80.png",
       plot = ggsurvplot(kmf_ph_karno_healthy, data = healthy_ph, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival (Ph-Karno ≥ 80)", legend = "none")$plot,
       width = 8, height = 6, units = "in")
ggsave("outputs/survival_ph_karno_lt_80.png",
       plot = ggsurvplot(kmf_ph_karno_sick, data = sick_ph, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival (Ph-Karno < 80)", legend = "none")$plot,
       width = 8, height = 6, units = "in")

# Physician ECOG score analysis (separate plots)
df_ph_ecog <- df2 %>%
  mutate(ph_ecog_cat = ifelse(ph_ecog <= 1, 1, 0))

healthy_ecog <- df_ph_ecog %>% filter(ph_ecog_cat == 1)
sick_ecog <- df_ph_ecog %>% filter(ph_ecog_cat == 0)

kmf_ph_ecog_healthy <- survfit(Surv(time, dead) ~ 1, data = healthy_ecog)
kmf_ph_ecog_sick <- survfit(Surv(time, dead) ~ 1, data = sick_ecog)

ggsave("outputs/survival_ph_ecog_leq_1.png",
       plot = ggsurvplot(kmf_ph_ecog_healthy, data = healthy_ecog, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival (Ph-Ecog ≤ 1)", legend = "none")$plot,
       width = 8, height = 6, units = "in")
ggsave("outputs/survival_ph_ecog_gt_1.png",
       plot = ggsurvplot(kmf_ph_ecog_sick, data = sick_ecog, conf.int = FALSE,
                         xlab = "Days passed", ylab = "Survival Probability",
                         title = "Survival (Ph-Ecog > 1)", legend = "none")$plot,
       width = 8, height = 6, units = "in")


# looking at the null values: (Already saved at the beginning)

# running a base analysis first by removing all the null values without any imputation

# temporary copy of the original dataframe
df_cph <- df

# identify columns to check for NA (all columns in this case)
dropper_subset <- names(df_cph)

# remove rows with NA in the specified columns
df_cph <- df_cph[complete.cases(df_cph[, dropper_subset]), ]

# dropping the status column also as the dead column is sufficient
df_cph <- df_cph %>% select(-status)

# sanity check on whether all null values were dropped (Already saved)

# creating a Cox proportional hazards model
library(survival)

# fit the Cox PH model
cph <- coxph(Surv(time, dead) ~ ., data = df_cph)

# checking whether the proportional hazards assumption is alright
# (This is more involved in R and often requires visual inspection and statistical tests)
# One common approach is to examine Schoenfeld residuals:
test.zph <- cox.zph(cph)

sink("outputs/cox_ph_assumption_test.txt")
print("Cox Proportional Hazards Assumption Test (Schoenfeld Residuals):")
print(test.zph)
sink()

png("outputs/cox_ph_assumption_plots.png", width = 1200, height = 800)
plot(test.zph) # Visual inspection of Schoenfeld residuals
dev.off()

# getting the summary from the Cox PH model
sink("outputs/cox_ph_model_summary.txt")
print("Cox Proportional Hazards Model Summary:")
print(summary(cph))
sink()

# Display the last 50 rows of the dataframe (optional, not typically saved as a primary output)
# tail(df_cph, 50)

# selecting 3 values from the actual observations (adjust indices for R)
df_checker <- df_cph[18:20, ] # R uses 1-based indexing

# Create a survival object for the selected individuals
surv_checker <- survfit(cph, newdata = df_checker)

# Plot the survival functions for selected individuals
ggsave("outputs/survival_predictions.png",
       plot = ggsurvplot(surv_checker,
                         data = df_checker,
                         fun = "pct", # Show survival probability as percentage
                         xlab = "Time",
                         ylab = "Survival Probability",
                         title = "Survival Probability of Subjects",
                         conf.int = FALSE, # Remove confidence intervals for simplicity
                         legend.title = "Subject ID")$plot,
       width = 8, height = 6, units = "in")