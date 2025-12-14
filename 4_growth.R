# ==============================================================================
# Multilevel Growth Curve Models (Chapter 4)
# Dataset: PBC (Primary Biliary Cirrhosis) Longitudinal Data
# Theoretical Framework: Equations 4.1-4.3 (Multilevel Models), 4.16-4.21 (Growth Curves)
# ==============================================================================

library(lme4)
library(lmerTest)
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(lattice)

if(!dir.exists("plots_growth")) dir.create("plots_growth")
cat("Saving plots to: plots_growth/\n\n")

# ---- Plot style (shared) ----
P_COL <- c(
  blue   = "#1F77B4",
  orange = "#FF7F0E",
  green  = "#2CA02C",
  red    = "#D62728",
  purple = "#9467BD",
  gray   = "#7F7F7F"
)

theme_project <- function(base_size = 14, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

set_base_par <- function() {
  par(
    family = "sans",
    cex.main = 1.2,
    cex.lab = 1.15,
    cex.axis = 1.05
  )
}

# ==============================================================================
# Part 1: Data Preparation and Exploration
# ==============================================================================

cat("========================================\n")
cat("Loading PBC longitudinal data\n")
cat("========================================\n")

# Load the longitudinal PBC data (pbcseq from survival package)
data(pbcseq)

# Data preprocessing
df <- pbcseq %>%
  filter(!is.na(trt)) %>%                    # Remove non-randomized patients
  mutate(
    log_bili = log(bili + 0.1),              # Log-transform bilirubin (right-skewed)
    log_protime = log(protime),              # Log-transform prothrombin time
    year = day / 365.25,                     # Convert days to years
    trt = factor(trt, labels = c("D-penicil", "Placebo")),
    id = factor(id),                         # Grouping factor (level 2)
    sex = factor(sex),
    stage = factor(stage)                    # Disease stage
  )

n_obs <- nrow(df)
n_patients <- length(unique(df$id))
obs_per_patient <- df %>% group_by(id) %>% summarise(n = n())

cat("\nData overview:\n")
cat("  Total observations:", n_obs, "\n")
cat("  Patients:", n_patients, "\n")
cat("  Mean visits per patient:", round(mean(obs_per_patient$n), 2), "\n")
cat("  Median visits per patient:", median(obs_per_patient$n), "\n")
cat("  Follow-up time range:", paste(round(range(df$year), 2), collapse = " - "), "years\n\n")

# ==============================================================================
# Part 2: Exploratory Data Analysis
# ==============================================================================

cat("========================================\n")
cat("Exploratory data analysis\n")
cat("========================================\n")

# Figure 1: Distribution of observations per patient
pdf("plots_growth/01_data_overview.pdf", width=12, height=5)
par(mfrow=c(1,3))
set_base_par()

hist(obs_per_patient$n, breaks=15, col="skyblue", border="white",
     main="Number of visits per patient", xlab="Visits", ylab="Patients")

hist(df$year, breaks=30, col="lightcoral", border="white",
     main="Visit time distribution", xlab="Years", ylab="Observations")

boxplot(log_bili ~ trt, data=df, col=c("lightblue", "lightgreen"),
        main="Baseline bilirubin by treatment group", ylab="log(bilirubin + 0.1)", xlab="Treatment")

dev.off()
cat("Saved: plots_growth/01_data_overview.pdf\n")

# Figure 2: Spaghetti plot (individual trajectories)
set.seed(123)
sample_ids <- sample(unique(df$id), 40)
df_sample <- df %>% filter(id %in% sample_ids)

p1 <- ggplot(df_sample, aes(x = year, y = log_bili, group = id, color = trt)) +
  geom_line(alpha = 0.7, size = 0.8) +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_wrap(~ trt) +
  labs(
    title = "Individual trajectories (40 randomly sampled patients)",
    subtitle = "Stratified by treatment group",
    y = "log(bilirubin + 0.1)",
    x = "Years",
    color = "Treatment"
  ) +
  scale_color_manual(values = c("D-penicil" = unname(P_COL["blue"]), "Placebo" = unname(P_COL["orange"]))) +
  theme_project()

ggsave("plots_growth/02_spaghetti_plot.pdf", p1, width=12, height=6)
cat("Saved: plots_growth/02_spaghetti_plot.pdf\n")

# Figure 3: Mean trajectories by treatment
p2 <- ggplot(df, aes(x = year, y = log_bili, color = trt)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.25, aes(fill = trt), color = NA) +
  labs(
    title = "Mean trajectories by treatment group",
    subtitle = "Mean ± SE",
    y = "log(bilirubin + 0.1)",
    x = "Years",
    color = "Treatment",
    fill = "Treatment"
  ) +
  scale_color_manual(values = c("D-penicil" = unname(P_COL["blue"]), "Placebo" = unname(P_COL["orange"]))) +
  scale_fill_manual(values = c("D-penicil" = unname(P_COL["blue"]), "Placebo" = unname(P_COL["orange"]))) +
  theme_project()

ggsave("plots_growth/03_mean_trajectories.pdf", p2, width=10, height=6)
cat("Saved: plots_growth/03_mean_trajectories.pdf\n")

# ==============================================================================
# Part 3: Multilevel Models (Chapter 4.2 - Equations 4.1-4.3)
# ==============================================================================

cat("\n========================================\n")
cat("Fitting multilevel models\n")
cat("========================================\n")

# Model 0: Null Model (Unconditional Means Model)
# Purpose: Calculate ICC (Intraclass Correlation Coefficient)
# Equation: y_ij = β_0 + u_j + ε_ij
cat("\nModel 0: Null model (ICC)\n")

m0 <- lmer(log_bili ~ 1 + (1 | id), data = df, REML = TRUE)
summary(m0)

# Calculate ICC
vars <- as.data.frame(VarCorr(m0))
var_between <- vars$vcov[1]  # Between-patient variance (τ²)
var_within <- vars$vcov[2]   # Within-patient variance (σ²)
icc <- var_between / (var_between + var_within)

cat("\n========================================\n")
cat("Intraclass Correlation (ICC)\n")
cat("========================================\n")
cat("ICC =", round(icc, 4), "\n")
cat("Interpretation:", round(icc*100, 1), "% variance is between patients\n")
cat("               ", round((1-icc)*100, 1), "% variance is within patients\n\n")

# Model 1: Random Intercept Model (Chapter 4.2, Equation 4.1)
# y_ij = x_ij β + z_ij u_j + ε_ij
cat("\nModel 1: Random intercept model\n")

m1 <- lmer(log_bili ~ year + trt + age + sex + (1 | id), 
           data = df, REML = FALSE)  # Use ML for model comparison
summary(m1)

# Model 2: Random Slope Model (Chapter 4.6, Equations 4.16-4.17)
# y_it = b_i1 + b_i2*t + ε_it
# (b_i1, b_i2) ~ N([β₁, β₂], Σ_b)
cat("\nModel 2: Random slope model\n")

m2 <- lmer(log_bili ~ year + trt + age + sex + (1 + year | id), 
           data = df, REML = FALSE)
summary(m2)

# Model comparison
cat("\n========================================\n")
cat("Model comparison (Likelihood Ratio Test)\n")
cat("========================================\n")
print(anova(m0, m1, m2))

# Model 3: Quadratic Growth Curve (Equation 4.18)
# y_it = b_i1 + b_i2*t + b_i3*t² + ε_it
cat("\nModel 3: Quadratic growth curve\n")

df$year2 <- df$year^2
m3 <- lmer(log_bili ~ year + year2 + trt + age + sex + (1 + year | id), 
           data = df, REML = FALSE)
summary(m3)

# Model 4: Cross-Level Interaction (Equation 4.20)
# Does treatment affect growth rate?
# b_ik = X_i γ_k + δ_ik
cat("\nModel 4: Treatment × time interaction\n")

m4 <- lmer(log_bili ~ year * trt + age + sex + (1 + year | id), 
           data = df, REML = FALSE)
summary(m4)

cat("\n========================================\n")
cat("Final model comparison\n")
cat("========================================\n")
print(anova(m2, m3, m4))

# ==============================================================================
# Part 4: Model Diagnostics
# ==============================================================================

cat("\n========================================\n")
cat("Model diagnostics\n")
cat("========================================\n")

# Extract residuals and fitted values
df$fitted_m4 <- fitted(m4)
df$resid_m4 <- residuals(m4)

# Get random effects
random_effects <- ranef(m4)$id
random_effects$id <- rownames(random_effects)
colnames(random_effects)[1:2] <- c("u0_intercept", "u1_slope")

# Figure 4: Residual diagnostics
pdf("plots_growth/04_residual_diagnostics.pdf", width=12, height=8)
par(mfrow=c(2,2))
set_base_par()

# Residuals vs Fitted
plot(df$fitted_m4, df$resid_m4, 
     main="Residuals vs Fitted", 
     xlab="Fitted values", ylab="Residuals",
     pch=19, cex=0.6, col=rgb(0,0,0,0.3))
abline(h=0, col=P_COL["red"], lwd=2, lty=2)
lowess_fit <- lowess(df$fitted_m4, df$resid_m4)
lines(lowess_fit, col=P_COL["blue"], lwd=2)

# Q-Q plot
qqnorm(df$resid_m4, main="Normal Q-Q (residuals)",
       pch=19, cex=0.6, col=rgb(0,0,0,0.3))
qqline(df$resid_m4, col=P_COL["red"], lwd=2)

# Scale-Location
sqrt_abs_resid <- sqrt(abs(df$resid_m4))
plot(df$fitted_m4, sqrt_abs_resid,
     main="Scale-Location",
     xlab="Fitted values", ylab="sqrt(|residuals|)",
     pch=19, cex=0.6, col=rgb(0,0,0,0.3))
lowess_fit2 <- lowess(df$fitted_m4, sqrt_abs_resid)
lines(lowess_fit2, col=P_COL["blue"], lwd=2)

# Histogram of residuals
hist(df$resid_m4, breaks=50, col=P_COL["blue"], border="white",
     main="Residual distribution", xlab="Residuals", freq=FALSE)
curve(dnorm(x, mean=mean(df$resid_m4), sd=sd(df$resid_m4)), 
      add=TRUE, col=P_COL["red"], lwd=2)

dev.off()
cat("Saved: plots_growth/04_residual_diagnostics.pdf\n")

# Figure 5: Random effects diagnostics
pdf("plots_growth/05_random_effects.pdf", width=12, height=8)
par(mfrow=c(2,2))
set_base_par()

# Random intercept distribution
hist(random_effects$u0_intercept, breaks=30, col=P_COL["blue"], border="white",
     main="Random intercepts", xlab="u0 (intercept)", freq=FALSE)
curve(dnorm(x, mean=0, sd=sd(random_effects$u0_intercept)), 
      add=TRUE, col=P_COL["red"], lwd=2)

# Random slope distribution
hist(random_effects$u1_slope, breaks=30, col=P_COL["orange"], border="white",
     main="Random slopes", xlab="u1 (slope)", freq=FALSE)
curve(dnorm(x, mean=0, sd=sd(random_effects$u1_slope)), 
      add=TRUE, col=P_COL["red"], lwd=2)

# Q-Q plots for random effects
qqnorm(random_effects$u0_intercept, main="Normal Q-Q: random intercepts",
       pch=19, col=rgb(0,0,0,0.35))
qqline(random_effects$u0_intercept, col=P_COL["red"], lwd=2)

qqnorm(random_effects$u1_slope, main="Normal Q-Q: random slopes",
       pch=19, col=rgb(0,0,0,0.35))
qqline(random_effects$u1_slope, col=P_COL["red"], lwd=2)

dev.off()
cat("Saved: plots_growth/05_random_effects.pdf\n")

# Figure 6: Correlation between random effects
pdf("plots_growth/06_random_effects_correlation.pdf", width=8, height=8)
set_base_par()
plot(random_effects$u0_intercept, random_effects$u1_slope,
     main="Correlation of random effects",
     xlab="u0 (intercept)",
     ylab="u1 (slope)",
     pch=19, col=rgb(0,0,0,0.35))
abline(lm(u1_slope ~ u0_intercept, data=random_effects), 
       col=P_COL["red"], lwd=2)
cor_value <- cor(random_effects$u0_intercept, random_effects$u1_slope)
text(min(random_effects$u0_intercept), max(random_effects$u1_slope),
     paste0("r = ", round(cor_value, 3)), pos=4, cex=1.3, col=P_COL["red"])
dev.off()
cat("Saved: plots_growth/06_random_effects_correlation.pdf\n")

# ==============================================================================
# Part 5: Visualizing Model Predictions
# ==============================================================================

cat("\n========================================\n")
cat("Visualizing model predictions\n")
cat("========================================\n")

# Figure 7: Fitted trajectories by treatment
df$pred_m4 <- predict(m4)

p3 <- ggplot(df, aes(x = year, y = log_bili)) +
  geom_point(aes(color = trt), alpha = 0.15) +
  geom_line(aes(y = pred_m4, group = id, color = trt), alpha = 0.4) +
  stat_summary(aes(y = pred_m4, color = trt), 
               fun = mean, geom = "line", linewidth = 1.4) +
  facet_wrap(~ trt) +
  labs(
    title = "Fitted trajectories (Model 4: interaction)",
    subtitle = "Thin lines: individual predictions; thick line: mean prediction",
    y = "log(bilirubin + 0.1)",
    x = "Years",
    color = "Treatment"
  ) +
  scale_color_manual(values = c("D-penicil" = unname(P_COL["blue"]), "Placebo" = unname(P_COL["orange"]))) +
  theme_project()

ggsave("plots_growth/07_fitted_trajectories.pdf", p3, width=12, height=6)
cat("Saved: plots_growth/07_fitted_trajectories.pdf\n")

# Figure 8: Treatment effect visualization
# Create prediction data
pred_data <- expand.grid(
  year = seq(0, max(df$year), length.out = 100),
  trt = c("D-penicil", "Placebo"),
  age = mean(df$age),
  sex = levels(df$sex)[1]
)

pred_data$pred <- predict(m4, newdata = pred_data, re.form = NA)

p4 <- ggplot(pred_data, aes(x = year, y = pred, color = trt)) +
  geom_line(linewidth = 1.4) +
  geom_ribbon(aes(ymin = pred - 0.1, ymax = pred + 0.1, fill = trt),
              alpha = 0.2) +
  labs(
    title = "Treatment effect (fixed effects)",
    subtitle = "Trajectory adjusted for age and sex",
    y = "Predicted log(bilirubin + 0.1)",
    x = "Years",
    color = "Treatment",
    fill = "Treatment"
  ) +
  scale_color_manual(values = c("D-penicil" = P_COL["blue"], "Placebo" = P_COL["orange"])) +
  scale_fill_manual(values = c("D-penicil" = P_COL["blue"], "Placebo" = P_COL["orange"])) +
  theme_project()

ggsave("plots_growth/08_treatment_effect.pdf", p4, width=10, height=6)
cat("Saved: plots_growth/08_treatment_effect.pdf\n")

# ==============================================================================
# Part 6: Individual Growth Curves
# ==============================================================================

cat("\n========================================\n")
cat("Individual growth curves\n")
cat("========================================\n")

# Select patients with extreme random effects
top_fast <- random_effects %>% 
  arrange(desc(u1_slope)) %>% 
  slice(1:5)

top_slow <- random_effects %>% 
  arrange(u1_slope) %>% 
  slice(1:5)

extreme_ids <- c(top_fast$id, top_slow$id)
df_extreme <- df %>% filter(id %in% extreme_ids)

# Figure 9: Extreme growth curves
p5 <- ggplot(df_extreme, aes(x = year, y = log_bili, group = id)) +
  geom_point(aes(color = id), size = 2) +
  geom_line(aes(y = pred_m4, color = id), size = 1) +
  facet_wrap(~ id, ncol = 5) +
  labs(
    title = "Example extreme trajectories",
    subtitle = "Top 5 fastest vs top 5 slowest slopes",
    y = "log(bilirubin + 0.1)",
    x = "Years"
  ) +
  theme_project() +
  theme(legend.position = "none")

ggsave("plots_growth/09_extreme_curves.pdf", p5, width=15, height=6)
cat("Saved: plots_growth/09_extreme_curves.pdf\n")

# ==============================================================================
# Part 7: Variance Components Analysis
# ==============================================================================

cat("\n========================================\n")
cat("Variance components\n")
cat("========================================\n")

# Extract variance components from Model 4
vc <- as.data.frame(VarCorr(m4))
cat("\n variance components (Model 4):\n")
print(vc)

# Calculate proportions
total_var <- sum(vc$vcov)
vc$proportion <- vc$vcov / total_var * 100

cat("\nVariance proportions:\n")
cat("Random intercept:", round(vc$proportion[1], 2), "%\n")
cat("Random slope:", round(vc$proportion[2], 2), "%\n")
cat("Residual:", round(vc$proportion[4], 2), "%\n")

# Figure 10: Variance decomposition
pdf("plots_growth/10_variance_decomposition.pdf", width=10, height=6)
par(mfrow=c(1,2))
set_base_par()

# Pie chart
pie(vc$vcov[c(1,2,4)], 
    labels = c("Random intercept", "Random slope", "Residual"),
    col = c(P_COL["blue"], P_COL["orange"], P_COL["green"]),
    main = "Variance decomposition")

# Bar plot
barplot(vc$proportion[c(1,2,4)],
        names.arg = c("Random intercept", "Random slope", "Residual"),
        col = c(P_COL["blue"], P_COL["orange"], P_COL["green"]),
        main = "Variance proportions",
        ylab = "Percent (%)",
        ylim = c(0, max(vc$proportion) * 1.2))
text(x = seq(0.7, by = 1.2, length.out = 3),
     y = vc$proportion[c(1,2,4)] + 2,
     labels = paste0(round(vc$proportion[c(1,2,4)], 1), "%"),
     font = 2)

dev.off()
cat("Saved: plots_growth/10_variance_decomposition.pdf\n")

# ==============================================================================
# Part 8: Summary Table
# ==============================================================================

cat("\n========================================\n")
cat("Model summary\n")
cat("========================================\n")

# Create summary table
model_summary <- data.frame(
  Model = c("M0: Null", "M1: Random Intercept", "M2: Random Slope", 
            "M3: Quadratic", "M4: Interaction"),
  LogLik = c(logLik(m0), logLik(m1), logLik(m2), logLik(m3), logLik(m4)),
  AIC = c(AIC(m0), AIC(m1), AIC(m2), AIC(m3), AIC(m4)),
  BIC = c(BIC(m0), BIC(m1), BIC(m2), BIC(m3), BIC(m4))
)

cat("\nModel comparison table:\n")
print(model_summary)

# Identify best model
best_aic <- which.min(model_summary$AIC)
best_bic <- which.min(model_summary$BIC)
cat("\nBest model (AIC):", model_summary$Model[best_aic], "\n")
cat("Best model (BIC):", model_summary$Model[best_bic], "\n")

# ==============================================================================
# Part 9: Clinical Interpretation
# ==============================================================================

cat("\n========================================\n")
cat("Interpretation\n")
cat("========================================\n")

# Extract key coefficients from Model 4
coef_m4 <- summary(m4)$coefficients

cat("\nKey fixed effects (Model 4):\n")
cat("1. Intercept:", round(coef_m4["(Intercept)", "Estimate"], 3), "\n")
cat("2. Year slope:", round(coef_m4["year", "Estimate"], 3), "\n")

if("year:trtPlacebo" %in% rownames(coef_m4)) {
  interaction_coef <- coef_m4["year:trtPlacebo", "Estimate"]
  interaction_p <- coef_m4["year:trtPlacebo", "Pr(>|t|)"]
  cat("3. Treatment × year interaction:", round(interaction_coef, 3),
      "(p =", signif(interaction_p, 3), ")\n")
}

cat("\nRandom effects:\n")
cat("SD(random intercept):", round(sqrt(vc$vcov[1]), 3), "\n")
cat("SD(random slope):", round(sqrt(vc$vcov[2]), 3), "\n")
cat("Corr(intercept, slope):", round(cor_value, 3), "\n")

# ==============================================================================
# Part 10: Stratified Growth Curve Analysis by GMM Classes
# ==============================================================================
# 理论整合: Chapter 5 (GMM) + Chapter 4 (Growth Curves)
# 先用GMM识别潜在患者亚组，再在每个亚组内拟合生长曲线模型

cat("\n========================================\n")
cat("Part 10: Stratified growth curve analysis (GMM classes)\n")
cat("========================================\n")
cat("Combining Chapter 5 (GMM) and Chapter 4 (Growth Curves)\n\n")

# Step 1: Load and fit GMM model from Chapter 5
# ------------------------------------------------------------------------------
cat("Step 1: Fit GMM (K=2)\n")

if(!require(lcmm)) install.packages("lcmm")
library(lcmm)

# Load longitudinal data for GMM
data_gmm <- read.csv("pbc_longitudinal.csv")

# Fit optimal GMM model (K=2 based on previous analysis)
cat("  Fitting GMM (K=2)...\n")

# Single class baseline
m_gmm1 <- hlme(bili ~ day, 
               random = ~day, 
               subject = 'id', 
               data = data_gmm, 
               ng = 1,
               verbose = FALSE)

# Two-class GMM
m_gmm2 <- hlme(bili ~ day, 
               mixture = ~day,
               random = ~day, 
               subject = 'id', 
               data = data_gmm, 
               ng = 2, 
               B = m_gmm1,
               verbose = FALSE)

cat("  ✓ GMM模型拟合完成\n")
cat("  AIC(K=1):", round(m_gmm1$AIC, 1), "\n")
cat("  AIC(K=2):", round(m_gmm2$AIC, 1), "\n")
cat("  最佳模型: K=2\n\n")

# Extract class membership
posterior <- m_gmm2$pprob
patient_classes <- posterior[, c("id", "class")]

# Convert id to factor to match df
patient_classes$id <- factor(patient_classes$id)

# Display class proportions
class_props <- table(patient_classes$class) / nrow(patient_classes)
cat("GMM class proportions:\n")
for(k in names(class_props)) {
  cat(sprintf("  Class %s: %.1f%% (N=%d)\n", 
              k, class_props[k]*100, 
              sum(patient_classes$class == k)))
}
cat("\n")

# Step 2: Merge GMM classes with Growth Curve data
# ------------------------------------------------------------------------------
cat("Step 2: Merge class labels into growth-curve data...\n")

# Merge class information
df_stratified <- df %>%
  left_join(patient_classes, by = c("id" = "id"))

# Check merge success
n_with_class <- sum(!is.na(df_stratified$class))
cat("  Merged", n_with_class, "/", nrow(df_stratified), "rows\n")
cat("  Patient class counts:\n")
print(table(df_stratified %>% distinct(id, class) %>% pull(class)))
cat("\n")

# Step 3: Fit Growth Curves within each GMM class
# ------------------------------------------------------------------------------
cat("Step 3: Fit growth curves within each class...\n\n")

# Get unique classes (may not be consecutive, e.g., 1 and 3)
unique_classes <- sort(unique(patient_classes$class))
models_by_class <- list()
coefs_by_class <- list()

for(k in unique_classes) {
  cat("--------------------\n")
  cat(sprintf("Class %d: random intercept + slope model\n", k))
  cat("--------------------\n")
  
  # Subset data
  df_k <- df_stratified %>% filter(class == k)
  n_patients_k <- length(unique(df_k$id))
  n_obs_k <- nrow(df_k)
  
  cat(sprintf("  Patients: %d\n", n_patients_k))
  cat(sprintf("  Observations: %d\n", n_obs_k))
  
  # Fit growth curve model
  # Model: y_ij = (β₀ + u₀j) + (β₁ + u₁j)*year + β₂*trt + β₃*age + ε_ij
  model_k <- lmer(log_bili ~ year + trt + age + sex + (1 + year | id), 
                  data = df_k, 
                  REML = FALSE)
  
  # Store model
  models_by_class[[as.character(k)]] <- model_k
  
  # Extract coefficients
  coef_k <- summary(model_k)$coefficients
  coefs_by_class[[as.character(k)]] <- coef_k
  
  # Print key findings
  cat(sprintf("  Intercept: %.3f (SE=%.3f, p%s)\n",
              coef_k["(Intercept)", "Estimate"],
              coef_k["(Intercept)", "Std. Error"],
              ifelse(coef_k["(Intercept)", "Pr(>|t|)"] < 0.001, "<0.001",
                     paste0("=", round(coef_k["(Intercept)", "Pr(>|t|)"], 3)))))
  
  cat(sprintf("  Year slope: %.3f (SE=%.3f, p%s)\n",
              coef_k["year", "Estimate"],
              coef_k["year", "Std. Error"],
              ifelse(coef_k["year", "Pr(>|t|)"] < 0.001, "<0.001",
                     paste0("=", round(coef_k["year", "Pr(>|t|)"], 3)))))
  
  # Variance components
  vc_k <- as.data.frame(VarCorr(model_k))
  cat(sprintf("  SD(random intercept): %.3f\n", sqrt(vc_k$vcov[1])))
  cat(sprintf("  SD(random slope): %.3f\n", sqrt(vc_k$vcov[2])))
  cat(sprintf("  SD(residual): %.3f\n\n", sqrt(vc_k$vcov[nrow(vc_k)])))
}

# Step 4: Compare slopes between classes
# ------------------------------------------------------------------------------
cat("========================================\n")
cat("Step 4: Compare slopes between classes\n")
cat("========================================\n")

# Create comparison table
slope_comparison <- data.frame(
  Class = integer(),
  Intercept = numeric(),
  Slope_Year = numeric(),
  SE_Slope = numeric(),
  P_value = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

for(k in unique_classes) {
  coef_k <- coefs_by_class[[as.character(k)]]
  slope_comparison <- rbind(slope_comparison, data.frame(
    Class = k,
    Intercept = coef_k["(Intercept)", "Estimate"],
    Slope_Year = coef_k["year", "Estimate"],
    SE_Slope = coef_k["year", "Std. Error"],
    P_value = coef_k["year", "Pr(>|t|)"],
    Interpretation = ifelse(coef_k["year", "Estimate"] > 0.05,
                            "Faster progression",
                            ifelse(coef_k["year", "Estimate"] > 0,
                                   "Slower progression",
                                   "Improvement"))
  ))
}

cat("\nSlope comparison:\n")
print(slope_comparison)

# Statistical test for slope difference
if(length(unique_classes) == 2) {
  slope_diff <- slope_comparison$Slope_Year[1] - slope_comparison$Slope_Year[2]
  se_diff <- sqrt(slope_comparison$SE_Slope[1]^2 + slope_comparison$SE_Slope[2]^2)
  z_stat <- slope_diff / se_diff
  p_diff <- 2 * (1 - pnorm(abs(z_stat)))
  
  cat("\nSlope difference test:\n")
  cat(sprintf("  Class %d vs Class %d:\n", unique_classes[1], unique_classes[2]))
  cat(sprintf("  Difference: %.4f (SE=%.4f)\n", slope_diff, se_diff))
  cat(sprintf("  Z: %.3f\n", z_stat))
  cat(sprintf("  P: %.4f %s\n", p_diff,
              ifelse(p_diff < 0.05, "***", "")))
}

# Step 5: Visualizations
# ------------------------------------------------------------------------------
cat("\n========================================\n")
cat("Step 5: Generate plots\n")
cat("========================================\n")

# Figure 11: Stratified trajectories
p_strat1 <- ggplot(df_stratified, aes(x = year, y = log_bili, color = factor(class))) +
  geom_point(alpha = 0.15, size = 0.5) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.4) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.25, 
               aes(fill = factor(class))) +
  facet_wrap(~ class, labeller = labeller(class = function(x) paste0("Class ", x))) +
  scale_color_manual(
    values = c("1" = unname(P_COL["blue"]), "2" = unname(P_COL["orange"]), 
               "3" = unname(P_COL["green"]), "4" = unname(P_COL["purple"]))
  ) +
  scale_fill_manual(
    values = c("1" = unname(P_COL["blue"]), "2" = unname(P_COL["orange"]), 
               "3" = unname(P_COL["green"]), "4" = unname(P_COL["purple"]))
  ) +
  labs(
    title = "Stratified growth curves by GMM class",
    subtitle = "Thick line: class mean; ribbon: SE",
    y = "log(bilirubin + 0.1)",
    x = "Years",
    color = "Class",
    fill = "Class"
  ) +
  theme_project()

ggsave("plots_growth/11_stratified_trajectories.pdf", p_strat1, width=12, height=6)
cat("Saved: plots_growth/11_stratified_trajectories.pdf\n")

# Figure 12: Slope comparison visualization
slope_comparison$Class <- factor(slope_comparison$Class)

p_strat2 <- ggplot(slope_comparison, 
                   aes(x = Class, y = Slope_Year, fill = Class)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = Slope_Year - 1.96*SE_Slope, 
                    ymax = Slope_Year + 1.96*SE_Slope),
                width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = P_COL["red"], size = 1) +
  scale_fill_manual(
    values = c("1" = unname(P_COL["blue"]), "2" = unname(P_COL["orange"]), 
               "3" = unname(P_COL["green"]), "4" = unname(P_COL["purple"]))
  ) +
  labs(
    title = "Progression rate by class",
    subtitle = "Error bars show 95% CI",
    y = "Year slope (log bilirubin / year)",
    x = "Class"
  ) +
  theme_project() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# Add significance annotation if applicable
if(length(unique_classes) == 2 && p_diff < 0.05) {
  p_strat2 <- p_strat2 + 
    annotate("text", x = 1.5, 
             y = max(slope_comparison$Slope_Year + 1.96*slope_comparison$SE_Slope) * 1.1,
             label = sprintf("p = %.4f ***", p_diff), 
             size = 5, fontface = "bold", color = P_COL["red"])
}

ggsave("plots_growth/12_slope_comparison.pdf", p_strat2, width=8, height=6)
cat("Saved: plots_growth/12_slope_comparison.pdf\n")

# Figure 13: Individual trajectories by class with model fits
set.seed(456)
sample_ids_per_class <- df_stratified %>%
  group_by(class) %>%
  distinct(id) %>%
  slice_sample(n = 10) %>%
  pull(id)

df_sample_stratified <- df_stratified %>% 
  filter(id %in% sample_ids_per_class)

# Add predictions
df_sample_stratified$pred_stratified <- NA
for(k in unique_classes) {
  idx_k <- df_sample_stratified$class == k
  df_sample_stratified$pred_stratified[idx_k] <- 
    predict(models_by_class[[as.character(k)]], 
            newdata = df_sample_stratified[idx_k, ])
}

p_strat3 <- ggplot(df_sample_stratified, 
                   aes(x = year, y = log_bili, group = id)) +
  geom_point(aes(color = factor(class)), alpha = 0.6, size = 1.5) +
  geom_line(aes(y = pred_stratified, color = factor(class)), 
            alpha = 0.8, size = 0.8) +
  facet_wrap(~ class, labeller = labeller(class = function(x) paste0("Class ", x))) +
  scale_color_manual(
    values = c(
      "1" = unname(P_COL["blue"]),   # Class 1: blue
      "2" = unname(P_COL["orange"])  # Class 2: orange
    ),
    breaks = c("1","2"),
    labels = c("Class 1","Class 2")
  ) +
  labs(
    title = "Individual fits by class",
    subtitle = "10 randomly sampled patients per class",
    y = "log(bilirubin + 0.1)",
    x = "Years",
    color = "Class"
  ) +
  theme_project()

ggsave("plots_growth/13_individual_fits_by_class.pdf", p_strat3, width=12, height=6)
cat("Saved: plots_growth/13_individual_fits_by_class.pdf\n")

# Figure 14: Fixed effects prediction comparison
pred_data_strat <- expand.grid(
  year = seq(0, max(df$year), length.out = 100),
  trt = levels(df$trt)[1],  # Use first treatment level
  age = mean(df$age),
  sex = levels(df$sex)[1],
  class = unique_classes
)

# Predict for each class
pred_data_strat$pred <- NA
for(k in unique_classes) {
  idx_k <- pred_data_strat$class == k
  pred_data_strat$pred[idx_k] <- 
    predict(models_by_class[[as.character(k)]], 
            newdata = pred_data_strat[idx_k, ], 
            re.form = NA)  # Fixed effects only
}

p_strat4 <- ggplot(pred_data_strat, aes(x = year, y = pred, color = factor(class))) +
  geom_line(linewidth = 1.4) +
  geom_ribbon(aes(ymin = pred - 0.15, ymax = pred + 0.15, fill = factor(class)), 
              alpha = 0.2) +
  scale_color_manual(
    values = c(
      "1" = unname(P_COL["blue"]),
      "2" = unname(P_COL["orange"])
    ),
    name = "Class",
    breaks = c("1","2"),
    labels = c("Class 1","Class 2")
  ) +
  scale_fill_manual(
    values = c(
      "1" = adjustcolor(unname(P_COL["blue"]), alpha.f = 0.35),
      "2" = adjustcolor(unname(P_COL["orange"]), alpha.f = 0.35)
    ),
    name = "Class",
    breaks = c("1","2"),
    labels = c("Class 1","Class 2")
  ) +
  labs(
    title = "Fixed-effects predictions by class",
    subtitle = "Adjusted for treatment, age, and sex",
    y = "Predicted log(bilirubin + 0.1)",
    x = "Years"
  ) +
  theme_project()

ggsave("plots_growth/14_fixed_effects_by_class.pdf", p_strat4, width=10, height=6)
cat("Saved: plots_growth/14_fixed_effects_by_class.pdf\n")

# ==============================================================================
# Final Summary with Stratified Analysis
# ==============================================================================
cat("========================================\n")
cat("Done.\n")
cat("ICC =", round(icc, 3), "\n")
cat("Best model (BIC):", model_summary$Model[best_bic], "\n")
cat("Stratified analysis: K =", length(unique_classes), "\n")
cat("========================================\n")
