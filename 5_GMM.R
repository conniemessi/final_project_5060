# ============================================
# Latent Subgroups Project (Chapter 5)
# Goal: Use Finite Mixture Models to identify hidden patient subgroups
# with different disease progression trajectories in PBC data
# ============================================
# 
# Theoretical Background (Chapter 5):
# - Finite Mixture Model (Eq 5.4): p(y) = π₁p₁(y) + ... + πₖpₖ(y)
# - Unobserved Heterogeneity (Eq 5.2-5.3): Captures latent subgroups
# - Latent Class Analysis: Identifies K distinct trajectory patterns

library(lcmm)
library(ggplot2)
library(dplyr)
library(gridExtra)

if(!dir.exists("plots_gmm")) dir.create("plots_gmm")
cat("Saving plots to: plots_gmm/\n\n")

data <- read.csv("pbc_longitudinal.csv")

# ============================================
# Part 1: Data overview and preprocessing
# ============================================

cat("=== Data overview ===\n")
cat("Total observations:", nrow(data), "\n")
cat("Number of patients:", length(unique(data$id)), "\n")
cat("Mean visits per patient:", round(nrow(data) / length(unique(data$id)), 2), "\n\n")

obs_per_patient <- data %>%
  group_by(id) %>%
  summarise(n_obs = n(), 
            max_day = max(day),
            mean_bili = mean(bili, na.rm=TRUE))

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

# Plot 1: Visit count distribution
pdf("plots_gmm/01_observation_distribution.pdf", width=10, height=5)
par(mfrow=c(1,2))
set_base_par()
hist(obs_per_patient$n_obs, breaks=20, 
     main="Number of visits per patient", 
     xlab="Visits", 
     ylab="Patients",
     col=P_COL["blue"], border="white")
hist(obs_per_patient$max_day, breaks=30,
     main="Maximum follow-up time",
     xlab="Days",
     ylab="Patients",
     col=P_COL["orange"], border="white")
dev.off()
cat("Saved: plots_gmm/01_observation_distribution.pdf\n")

# Plot 2: Spaghetti plot of bilirubin (sample of patients)
pdf("plots_gmm/02_spaghetti_plot.pdf", width=10, height=6)
set.seed(123)
sample_ids <- sample(unique(data$id), min(50, length(unique(data$id))))
data_sample <- data %>% filter(id %in% sample_ids)

ggplot(data_sample, aes(x=day, y=log(bili+0.1), group=id)) +
  geom_line(alpha=0.3, color=P_COL["blue"]) +
  geom_smooth(aes(group=1), method="loess", se=TRUE, 
              color=P_COL["red"], size=1.2) +
  labs(title="Bilirubin trajectories (50 randomly sampled patients)",
       subtitle="LOESS smooth (population trend)",
       x="Days since baseline",
       y="log(bilirubin + 0.1)") +
  theme_project()
dev.off()
cat("Saved: plots_gmm/02_spaghetti_plot.pdf\n")

# ============================================
# Part 2: Fit latent class mixed models
# ============================================

cat("\n=== Fitting latent class mixed models ===\n")
cat("This may take a few minutes...\n\n")

cat("Fitting model: K=1 (single class)\n")
m1 <- hlme(bili ~ day, 
           random = ~day, 
           subject = 'id', 
           data = data, 
           ng = 1,
           verbose = FALSE)

cat("Fitting model: K=2\n")
m2 <- hlme(bili ~ day, 
           mixture = ~day,
           random = ~day, 
           subject = 'id', 
           data = data, 
           ng = 2, 
           B = m1,
           verbose = FALSE)

cat("Fitting model: K=3\n")
m3 <- hlme(bili ~ day, 
           mixture = ~day, 
           random = ~day, 
           subject = 'id', 
           data = data, 
           ng = 3, 
           B = m1,
           verbose = FALSE)

cat("Fitting model: K=4\n")
m4 <- hlme(bili ~ day, 
           mixture = ~day, 
           random = ~day, 
           subject = 'id', 
           data = data, 
           ng = 4, 
           B = m1,
           verbose = FALSE)

cat("\nAll models fitted.\n\n")

# ============================================
# Part 3: Model comparison and selection
# ============================================

cat("=== Model comparison (Information Criteria) ===\n")
comparison <- summarytable(m1, m2, m3, m4)
print(comparison)

model_stats <- data.frame(
  K = 1:4,
  LogLik = c(m1$loglik, m2$loglik, m3$loglik, m4$loglik),
  AIC = c(m1$AIC, m2$AIC, m3$AIC, m4$AIC),
  BIC = c(m1$BIC, m2$BIC, m3$BIC, m4$BIC)
)

best_aic <- which.min(model_stats$AIC)
best_bic <- which.min(model_stats$BIC)

cat("\nBest model selection:\n")
cat("  By AIC: K =", best_aic, "(AIC =", round(model_stats$AIC[best_aic], 2), ")\n")
cat("  By BIC: K =", best_bic, "(BIC =", round(model_stats$BIC[best_bic], 2), ")\n\n")

# Plot 3: AIC/BIC curves
pdf("plots_gmm/03_model_comparison.pdf", width=10, height=5)
par(mfrow=c(1,2))
set_base_par()

plot(model_stats$K, model_stats$AIC, 
     type="b", pch=19, col=P_COL["blue"], lwd=2,
     xlab="Number of classes (K)", 
     ylab="AIC (lower is better)",
     main="AIC",
     xaxt="n")
axis(1, at=1:4)
points(best_aic, model_stats$AIC[best_aic], 
       pch=19, col=P_COL["red"], cex=1.6)
text(best_aic, model_stats$AIC[best_aic], 
     labels=paste0("Best: K=", best_aic), pos=3, col=P_COL["red"], font=2)

plot(model_stats$K, model_stats$BIC, 
     type="b", pch=19, col=P_COL["green"], lwd=2,
     xlab="Number of classes (K)", 
     ylab="BIC (lower is better)",
     main="BIC",
     xaxt="n")
axis(1, at=1:4)
points(best_bic, model_stats$BIC[best_bic], 
       pch=19, col=P_COL["red"], cex=1.6)
text(best_bic, model_stats$BIC[best_bic], 
     labels=paste0("Best: K=", best_bic), pos=3, col=P_COL["red"], font=2)

dev.off()
cat("Saved: plots_gmm/03_model_comparison.pdf\n")

# ============================================
# Part 4: 分析最佳模型 (假设选择m2或m3)
# ============================================

# 根据BIC选择最佳模型
best_model <- list(m1, m2, m3, m4)[[best_bic]]
K_best <- best_bic

cat("\n=== Best model analysis: K =", K_best, "===\n\n")

# 提取类别概率 (mixing proportions πk)
if(K_best > 1) {
  class_props <- best_model$pprob[1:K_best, 2]
  names(class_props) <- paste0("Class ", 1:K_best)
  
  cat("Mixing proportions (pi_k):\n")
  for(k in 1:K_best) {
    cat(sprintf("  Class %d: %.1f%% (n ≈ %d patients)\n",
                k, class_props[k]*100, 
                round(class_props[k] * length(unique(data$id)))))
  }
  cat("\n")
}

# Plot 4: Class proportions
if(K_best > 1) {
  pdf("plots_gmm/04_class_proportions.pdf", width=8, height=6)
  actual_n_classes <- length(class_props)

  set_base_par()
  barplot(class_props, 
          col=c(P_COL["blue"], P_COL["orange"], P_COL["green"], P_COL["purple"])[1:actual_n_classes],
          main=paste0("Finite Mixture Model: K = ", K_best),
          ylab="Proportion (pi_k)",
          xlab="Class",
          ylim=c(0, max(class_props)*1.2),
          border="white")
  text(x=seq(0.7, by=1.2, length.out=actual_n_classes), 
       y=class_props + 0.03, 
       labels=paste0(round(class_props*100, 1), "%"),
       font=2)
  dev.off()
  cat("Saved: plots_gmm/04_class_proportions.pdf\n")
}

# Plot 5: Posterior classification quality
if(K_best > 1) {
  pdf("plots_gmm/05_posterior_probabilities.pdf", width=10, height=6)
  plot(best_model, which="postprob", 
       main=paste0("Posterior classification probabilities (K=", K_best, ")"))
  dev.off()
  cat("Saved: plots_gmm/05_posterior_probabilities.pdf\n")
}

# ============================================
# Part 5: 可视化不同类别的轨迹
# ============================================

# Plot 6: Class-specific trajectories (lcmm)
if(K_best > 1) {
  pdf("plots_gmm/06_class_trajectories.pdf", width=12, height=8)
  plot(best_model, which="fit", var.time="day", 
       break.times=10,
       main=paste0("Class-specific bilirubin trajectories (K=", K_best, ")"))
  dev.off()
  cat("Saved: plots_gmm/06_class_trajectories.pdf\n")
}

# Plot 7: Class trajectories (ggplot)
if(K_best > 1) {
  posterior <- best_model$pprob
  patient_classes <- posterior[, c("id", "class")]
  
  data_with_class <- data %>%
    left_join(patient_classes, by="id")
  
  trajectory_summary <- data_with_class %>%
    filter(!is.na(bili)) %>%
    mutate(time_bin = cut(day, breaks=10)) %>%
    group_by(class, time_bin) %>%
    summarise(
      day_mid = mean(day),
      mean_bili = mean(bili),
      se_bili = sd(bili) / sqrt(n()),
      .groups = "drop"
    )
  
  p <- ggplot(trajectory_summary, aes(x=day_mid, y=log(mean_bili+0.1), 
                                       color=factor(class), 
                                       fill=factor(class))) +
    geom_line(size=1.2) +
    geom_ribbon(aes(ymin=log(mean_bili - se_bili + 0.1), 
                    ymax=log(mean_bili + se_bili + 0.1)),
                alpha=0.2, color=NA) +
    labs(
      title=paste0("Class-specific mean trajectories (K=", K_best, ")"),
      subtitle="Mean bilirubin over time (± SE)",
      x="Days since baseline",
      y="log(mean bilirubin + 0.1)",
      color="Class",
      fill="Class"
    ) +
    scale_color_manual(values = setNames(c(P_COL["blue"], P_COL["orange"], P_COL["green"], P_COL["purple"])[1:K_best], 1:K_best)) +
    scale_fill_manual(values = setNames(c(P_COL["blue"], P_COL["orange"], P_COL["green"], P_COL["purple"])[1:K_best], 1:K_best)) +
    theme_project()
  
  ggsave("plots_gmm/07_class_trajectories_ggplot.pdf", p, 
         width=10, height=6)
  cat("Saved: plots_gmm/07_class_trajectories_ggplot.pdf\n")
}

# ============================================
# Part 6: 类别特征分析
# ============================================

if(K_best > 1) {
  cat("\n=== Baseline characteristics by class ===\n")
  
  # 合并后验类别
  posterior <- best_model$pprob
  patient_classes <- posterior[, c("id", "class")]
  
  # 计算baseline特征
  baseline_data <- data %>%
    group_by(id) %>%
    filter(day == min(day)) %>%
    ungroup() %>%
    left_join(patient_classes, by="id")
  
  # 按类别汇总
  class_summary <- baseline_data %>%
    group_by(class) %>%
    summarise(
      N = n(),
      Age_mean = mean(age, na.rm=TRUE),
      Age_sd = sd(age, na.rm=TRUE),
      Female_pct = mean(sex=="f", na.rm=TRUE) * 100,
      Bili_mean = mean(bili, na.rm=TRUE),
      Bili_sd = sd(bili, na.rm=TRUE),
      Albumin_mean = mean(albumin, na.rm=TRUE),
      .groups = "drop"
    )
  
  print(class_summary)
  
# Plot 8: Baseline characteristics
  pdf("plots_gmm/08_class_characteristics.pdf", width=12, height=8)
  par(mfrow=c(2,2))
  set_base_par()
  
  # 获取实际的类别数量（修复bug）
  actual_n_classes <- length(unique(baseline_data$class))
  actual_class_labels <- sort(unique(baseline_data$class))
  
  boxplot(age ~ class, data=baseline_data,
          col=c(P_COL["blue"], P_COL["orange"], P_COL["green"], P_COL["purple"])[1:actual_n_classes],
          main="Baseline age",
          xlab="Class",
          ylab="Age (years)")
  
  boxplot(log(bili+0.1) ~ class, data=baseline_data,
          col=c(P_COL["blue"], P_COL["orange"], P_COL["green"], P_COL["purple"])[1:actual_n_classes],
          main="Baseline bilirubin",
          xlab="Class",
          ylab="log(bilirubin + 0.1)")
  
  boxplot(albumin ~ class, data=baseline_data,
          col=c(P_COL["blue"], P_COL["orange"], P_COL["green"], P_COL["purple"])[1:actual_n_classes],
          main="Baseline albumin",
          xlab="Class",
          ylab="Albumin (g/dL)")
  
  sex_table <- table(baseline_data$class, baseline_data$sex)
  barplot(sex_table, beside=TRUE,
          col=c(P_COL["blue"], P_COL["orange"], P_COL["green"], P_COL["purple"])[1:actual_n_classes],
          legend=paste("Class", actual_class_labels),
          main="Sex distribution",
          xlab="Sex",
          ylab="Patients")
  
  dev.off()
  cat("Saved: plots_gmm/08_class_characteristics.pdf\n")
}

# ============================================
# Part 7: 生存分析（如果有）
# ============================================

if(K_best > 1 && "status" %in% colnames(data)) {
  cat("\n=== Survival outcomes by class ===\n")
  
  survival_data <- data %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup() %>%
    left_join(patient_classes, by="id")
  
  survival_summary <- survival_data %>%
    group_by(class) %>%
    summarise(
      N = n(),
      Death_rate = mean(status==2, na.rm=TRUE) * 100,
      Median_futime = median(futime, na.rm=TRUE),
      .groups = "drop"
    )
  
  print(survival_summary)
  
# Plot 9: Survival outcomes
  pdf("plots_gmm/09_survival_outcomes.pdf", width=10, height=6)
  par(mfrow=c(1,2))
  set_base_par()
  
  # 获取实际的类别编号（修复bug：类别可能不连续）
  actual_classes <- survival_summary$class
  n_classes <- nrow(survival_summary)
  
  barplot(survival_summary$Death_rate,
          names.arg=paste("Class", actual_classes),
          col=c(P_COL["blue"], P_COL["orange"], P_COL["green"], P_COL["purple"])[1:n_classes],
          main="Mortality rate by class",
          ylab="Mortality (%)",
          ylim=c(0, max(survival_summary$Death_rate)*1.2))
  text(x=seq(0.7, by=1.2, length.out=n_classes),
       y=survival_summary$Death_rate + 2,
       labels=paste0(round(survival_summary$Death_rate, 1), "%"),
       font=2)
  
  barplot(survival_summary$Median_futime,
          names.arg=paste("Class", actual_classes),
          col=c(P_COL["blue"], P_COL["orange"], P_COL["green"], P_COL["purple"])[1:n_classes],
          main="Median follow-up time by class",
          ylab="Days")
  
  dev.off()
  cat("Saved: plots_gmm/09_survival_outcomes.pdf\n")
}

# ============================================
# Part 8: 诊断检查
# ============================================

# Plot 10: Residual diagnostics
pdf("plots_gmm/10_residual_diagnostics.pdf", width=12, height=8)
plot(best_model, which="residuals")
dev.off()
cat("Saved: plots_gmm/10_residual_diagnostics.pdf\n")

cat("Best model: K =", K_best, "\n")
cat("========================================\n")
