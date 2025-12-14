# ============================================
# Missing Data Analysis (Chapter 6)
# Methods: Naive OLS, Heckman Selection, Pattern Mixture (IPW), Multiple Imputation, Common Factor Model (factor-score proxy)
# ============================================

library(sampleSelection)
library(mice)
library(naniar)
library(ggplot2)
library(dplyr)
library(tidyr)
library(psych)

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

if(!dir.exists("plots")) dir.create("plots")
cat("Saving plots to: plots/\n\n")

data <- read.csv("pbc_longitudinal.csv")

# ============================================
# Part 1: Missingness overview
# ============================================

cat("\n=== Missingness summary ===\n")
miss_var_summary(data)

data$chol_observed <- ifelse(is.na(data$chol), 0, 1)

pdf("plots/01_missing_pattern.pdf", width = 10, height = 6)
set_base_par()
vis_miss(data[, c("chol", "age", "sex", "bili", "albumin", "alk.phos")])
dev.off()
cat("Saved: plots/01_missing_pattern.pdf\n")

p1 <- ggplot(data, aes(x = factor(chol_observed), y = log(bili + 0.1), fill = factor(chol_observed))) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(
    values = c("0" = unname(P_COL["red"]), "1" = unname(P_COL["blue"])),
    labels = c("0" = "Missing", "1" = "Observed")
  ) +
  labs(
    title = "Bilirubin by Cholesterol Missingness",
    x = "Cholesterol observed? (0 = missing, 1 = observed)",
    y = "log(bilirubin + 0.1)",
    fill = "Observed"
  ) +
  theme_project()
ggsave("plots/02_bili_vs_missing.pdf", p1, width = 8, height = 6)
cat("Saved: plots/02_bili_vs_missing.pdf\n")

t.test(bili ~ chol_observed, data = data)

# ============================================
# Part 2: Method 1 - Naive (Complete Case)
# ============================================

cat("\n=== Method 1: Naive Complete Case Analysis ===\n")
naive_model <- lm(chol ~ age + sex + bili + albumin, 
                  data = data, na.action = na.omit)
summary(naive_model)

# ============================================
# Part 3: Method 2 - Heckman Selection Model (MNAR)
# ============================================

cat("\n=== Method 2: Heckman Selection Model (Eq. 6.1) ===\n")

selection_eq <- chol_observed ~ age + sex + bili + albumin

outcome_eq <- chol ~ age + sex + bili

heckman_model <- selection(selection_eq, outcome_eq, data = data)
heckman_summary <- summary(heckman_model)
print(heckman_summary)

rho_value <- heckman_summary$estimate["rho", "Estimate"]
rho_pvalue <- heckman_summary$estimate["rho", "Pr(>|t|)"]

cat("\n" , "="  , rep("=", 50), "\n", sep="")
cat("Key parameter rho (error correlation): ", round(rho_value, 4), "\n")
cat("="  , rep("=", 50), "\n", sep="")
if(rho_pvalue < 0.05) {
  cat("rho is significantly different from 0 (p < 0.05): strong evidence for MNAR.\n")
} else {
  cat("rho is not significant: missingness may be close to MAR.\n")
}
cat("="  , rep("=", 50), "\n\n", sep="")

# ============================================
# Part 4: Method 3 - Pattern Mixture Model (Eq. 6.2) via IPW
# ============================================

cat("\n=== Method 3: Pattern Mixture Model (IPW) ===\n")

missing_model <- glm(chol_observed ~ age + sex + bili + albumin, 
                     data = data, family = binomial(link = "logit"))

data$ps <- predict(missing_model, type = "response")

data$ipw_weight <- ifelse(data$chol_observed == 1, 1/data$ps, 0)

pmm_ipw_model <- lm(chol ~ age + sex + bili + albumin, 
                    data = data[data$chol_observed == 1, ],
                    weights = ipw_weight[data$chol_observed == 1])
cat("\nWeighted regression (IPW) results:\n")
summary(pmm_ipw_model)

pdf("plots/07_pattern_mixture_ps.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))
set_base_par()
hist(data$ps[data$chol_observed == 1], main = "Propensity score (Observed)", xlab = "P(observed)", col = P_COL["blue"], breaks = 30)
hist(data$ps[data$chol_observed == 0], main = "Propensity score (Missing)",  xlab = "P(observed)", col = P_COL["red"],  breaks = 30)
dev.off()
cat("Saved: plots/07_pattern_mixture_ps.pdf\n")

p_pattern <- ggplot(data, aes(x = bili, y = ps, color = factor(chol_observed))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "loess") +
  labs(
    title = "Propensity Score vs. Bilirubin",
    subtitle = "Colored by cholesterol missingness",
    x = "Bilirubin",
    y = "P(observed | covariates)",
    color = "Observed"
  ) +
  scale_color_manual(values = c("0" = unname(P_COL["red"]), "1" = unname(P_COL["blue"])),
                     labels = c("0" = "Missing", "1" = "Observed")) +
  theme_project()
ggsave("plots/08_pattern_bili_ps.pdf", p_pattern, width = 10, height = 6)
cat("Saved: plots/08_pattern_bili_ps.pdf\n")

# ============================================
# Part 5: Method 4 - Multiple Imputation (MAR)
# ============================================

cat("\n=== Method 4: Multiple Imputation (MICE, Eq. 6.17) ===\n")

imputed_data <- mice(data[, c("chol", "age", "sex", "bili", "albumin")], 
                     m = 5, method = "pmm", seed = 123, printFlag = FALSE)

mi_model <- with(imputed_data, lm(chol ~ age + sex + bili + albumin))
mi_results <- pool(mi_model)
summary(mi_results)

pdf("plots/04_imputation_density.pdf", width = 10, height = 8)
set_base_par()
densityplot(imputed_data)
dev.off()
cat("Saved: plots/04_imputation_density.pdf\n")

# ============================================
# Part 6: Method 5 - Common Factor Model (factor-score proxy)
# ============================================

cat("\n=== Method 5: Common Factor Model (EFA factor-score proxy) ===\n")

outcomes <- c("chol", "platelet", "albumin", "bili")
imputed_for_fa <- mice(data[, c(outcomes, "age", "sex")], m = 1, method = "pmm", printFlag = FALSE)
data_imputed <- complete(imputed_for_fa)

fa_full <- fa(data_imputed[, outcomes], nfactors = 2, scores = "regression")
data$factor1 <- fa_full$scores[, 1]
data$factor2 <- fa_full$scores[, 2]

selection_eq_cfm <- chol_observed ~ age + sex + factor1 + factor2 + alk.phos
outcome_eq_cfm <- chol ~ age + sex + factor1 + factor2
heckman_cfm <- selection(selection_eq_cfm, outcome_eq_cfm, data = data)
summary(heckman_cfm)

# ============================================
# Part 7: Summary comparison (Age coefficient)
# ============================================

cat("\n=== Summary comparison across methods (Age coefficient) ===\n")

coef_naive <- coef(naive_model)["age"]
coef_heckman <- summary(heckman_model)$estimate["age", "Estimate"] # This picks the first "age", usually selection. Need to be specific.
# In sampleSelection summary, coefficients are named "outcome_age" or "age" depending on structure.
# Let's check the structure or access via specific indices if names are ambiguous, 
# but usually it's selection then outcome. 
# Better: extraction from the 2-step or ML object directly often separates them.
# For 'selection' object: coef(heckman_model) has names like "S:age", "O:age" or similar.
# Let's check names from the terminal output. 
# Terminal output shows:
# Probit selection equation: ... age ...
# Outcome equation: ... age ...
# The summary object usually has a matrix. 
# Let's try to be safer by checking names.

# Inspecting terminal output for 6_missing.R again.
# The summary prints "Probit selection equation" then "Outcome equation".
# If we look at `coef(heckman_model)`, it usually returns a vector with prefixes.
# Let's use specific indices or names if possible. 
# In 'sampleSelection', coefficients are typically S:var and O:var.
# Let's try to find "O:age" or "outcome_age".

coef_heckman <- coef(heckman_model)["O:age"]
if(is.na(coef_heckman)) coef_heckman <- coef(heckman_model)["age"] # Fallback if names differ (e.g. tobit2)
# Actually, looking at the terminal output for summary(heckman_model), it lists them separately.
# Let's use the coefficient vector directly which is safer.
coefs_heck <- coef(heckman_model)
coef_heckman <- coefs_heck[grep("^O:age|^outcome:age", names(coefs_heck))]
if(length(coef_heckman) == 0) coef_heckman <- coefs_heck["age"] # Fallback

coef_pmm <- coef(pmm_ipw_model)["age"]
coef_mi <- summary(mi_results)$estimate[2] # Assuming age is 2nd

# For CFM (Heckman selection model)
coefs_cfm <- coef(heckman_cfm)
coef_cfm <- coefs_cfm[grep("^O:age|^outcome:age", names(coefs_cfm))]
if(length(coef_cfm) == 0) coef_cfm <- tryCatch(summary(heckman_cfm)$estimate["age", "Estimate"], error = function(e) NA)


coef_comparison <- data.frame(
  Method = c("Naive", "Heckman", "PMM (IPW)", "MI", "CFM"),
  Coefficient = c(coef_naive, coef_heckman, coef_pmm, coef_mi, coef_cfm),
  Assumption = c("MCAR", "MNAR", "MNAR", "MAR", "MNAR")
)

print(coef_comparison)

p3 <- ggplot(coef_comparison, aes(x = reorder(Method, -Coefficient), y = Coefficient, fill = Assumption)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = P_COL["red"]) +
  labs(
    title = "Age Coefficient Across Five Missing-Data Methods",
    subtitle = "Different assumptions can lead to different conclusions",
    y = "Coefficient (Age)",
    x = NULL
  ) +
  scale_fill_manual(values = c("MCAR" = unname(P_COL["green"]), "MNAR" = unname(P_COL["orange"]), "MAR" = unname(P_COL["blue"]))) +
  theme_project()

ggsave("plots/09_age_coefficient_comparison.pdf", p3, width = 10, height = 6)
cat("Saved: plots/09_age_coefficient_comparison.pdf\n")

pdf("plots/05_residual_diagnostics.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
set_base_par()
plot(naive_model)
dev.off()
cat("Saved: plots/05_residual_diagnostics.pdf\n")

cat("\n========================================\n")
cat("Done.\n")
cat("========================================\n")
