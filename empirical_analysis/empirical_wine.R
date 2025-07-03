required_packages <- c(
  "MASS", "sure", "DHARMa", "tidyverse", "ggplot2", "gridExtra", 
  "ordinal", "brant", "car", "performance", "stats", "generalhoslem", 
  "gofcat", "Matrix", "insight", "VGAM", "rms", "corrplot", "GGally"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("Installing missing packages:", paste(new_packages, collapse=", "), "\n")
  install.packages(new_packages, dependencies=TRUE)
}

for(pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only=TRUE))
}
cat("All required packages loaded successfully.\n")

set.seed(123)

cb_cols <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
             "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

project_dir <- getwd()
results_dir <- file.path(project_dir, "results_empirical_redwine_v3")
plots_dir <- file.path(project_dir, "plots_empirical_redwine_v3")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

start_time <- Sys.time()
cat("\n========== RED WINE QUALITY EMPIRICAL ANALYSIS ==========\n")
cat("Starting analysis at:", format(start_time), "\n")
cat("Results will be saved to:", results_dir, "\n")
cat("Plots will be saved to:", plots_dir, "\n\n")

cat("\n--- 1. Data Loading and Preparation ---\n")

wine_data <- read.csv("C:/Users/leoko/Documents/TU/bachelorarbeit/simulationStudy/winequality-red.csv", sep = ";", header = TRUE)

colnames(wine_data) <- gsub(" ", "_", colnames(wine_data))

cat("\nDataset dimensions:", nrow(wine_data), "rows,", ncol(wine_data), "columns\n")
cat("Column names:", paste(colnames(wine_data), collapse = ", "), "\n")

wine_data$quality <- as.ordered(wine_data$quality)
cat("\nQuality levels:", paste(levels(wine_data$quality), collapse = " < "), "\n")

missing_count <- sum(is.na(wine_data))
cat("Missing values:", missing_count, "\n")

summary_file <- file.path(results_dir, "data_summary_redwine.txt")
sink(summary_file)
cat("RED WINE QUALITY DATA SUMMARY\n")
cat("=============================\n\n")
print(summary(wine_data))
sink()
cat("Data summary saved to:", summary_file, "\n")

cat("\n--- 2. Exploratory Data Analysis ---\n")

cat("\nQuality distribution:\n")
quality_table <- table(wine_data$quality)
print(quality_table)
cat("\nProportions:\n")
print(prop.table(quality_table))

png(file.path(plots_dir, "wine_quality_distribution.png"), width = 800, height = 600, res = 100)
print(
  ggplot(wine_data, aes(x = quality)) +
    geom_bar(fill = cb_cols[3], alpha = 0.7) +
    labs(title = "Distribution of Red Wine Quality Ratings",
         x = "Quality Rating", y = "Count") +
    theme_minimal()
)
dev.off()

predictors <- wine_data[, -which(names(wine_data) == "quality")]
cor_matrix <- cor(predictors)
cb_corr_colors <- colorRampPalette(c(cb_cols[7], "white", cb_cols[3]))(100)

png(file.path(plots_dir, "wine_correlation_matrix.png"), width = 1000, height = 1000, res = 100)
corrplot(cor_matrix, method = "color", type = "upper", 
         order = "hclust", tl.cex = 0.8, tl.col = "black",
         addCoef.col = "black", number.cex = 0.7,
         col = cb_corr_colors)
dev.off()

png(file.path(plots_dir, "wine_predictors_by_quality.png"), width = 1600, height = 1200, res = 100)
par(mfrow = c(3, 4))
for(i in 1:(ncol(wine_data)-1)) {
  boxplot(wine_data[,i] ~ wine_data$quality, 
          main = names(wine_data)[i],
          xlab = "Quality", ylab = "",
          col = cb_cols[2:7])
}
dev.off()

skewness_file <- file.path(results_dir, "skewness_outliers_redwine.txt")
sink(skewness_file)
cat("SKEWNESS AND OUTLIER ANALYSIS\n")
cat("==============================\n\n")
for(var in names(predictors)) {
  cat(var, ":\n")
  cat("  Skewness:", round(e1071::skewness(wine_data[[var]]), 3), "\n")
  outliers <- boxplot.stats(wine_data[[var]])$out
  cat("  Number of outliers:", length(outliers), "\n")
  cat("  Outlier range:", 
      if(length(outliers) > 0) paste(round(range(outliers), 3), collapse = " - ") else "None", "\n\n")
}
sink()
cat("Skewness and outlier analysis saved to:", skewness_file, "\n")

cat("\n--- 3. Collinearity Diagnostics ---\n")

pred_matrix <- as.matrix(predictors)
n_vars <- ncol(pred_matrix)

pred_scaled <- scale(pred_matrix)

cat("\nCalculating condition indices...\n")

cor_eigen <- eigen(cor_matrix)
eigenvalues <- cor_eigen$values
condition_indices <- sqrt(max(eigenvalues) / eigenvalues)

cat("Calculating variance proportions...\n")

eigenvectors <- cor_eigen$vectors
variance_props <- eigenvectors^2

variance_props <- apply(variance_props, 2, function(x) x / sum(x))
rownames(variance_props) <- names(predictors)
colnames(variance_props) <- paste("PC", 1:n_vars, sep="")

cat("Calculating VIF values...\n")

vif_values <- numeric(n_vars)
names(vif_values) <- names(predictors)

for(i in 1:n_vars) {
  y <- pred_matrix[, i]
  X <- pred_matrix[, -i]
  
  lm_fit <- lm(y ~ X)
  r_squared <- summary(lm_fit)$r.squared
  
  vif_values[i] <- 1 / (1 - r_squared)
}

cat("Calculating condition number...\n")

condition_number <- sqrt(max(eigenvalues) / min(eigenvalues))

max_condition_index <- max(condition_indices)

collinearity_file <- file.path(results_dir, "collinearity_diagnostics_redwine.txt")
sink(collinearity_file)
cat("COLLINEARITY DIAGNOSTICS - RED WINE DATA\n")
cat("=========================================\n\n")

cat("CONDITION NUMBER:\n")
cat("-----------------\n")
cat("Condition Number (sqrt(λ_max/λ_min)):", round(condition_number, 2), "\n")
cat("Maximum Condition Index:", round(max_condition_index, 2), "\n")
cat("Interpretation: CN > 15 suggests weak dependencies, CN > 30 suggests strong collinearity\n\n")

cat("CONDITION INDICES:\n")
cat("------------------\n")
condition_df <- data.frame(
  Component = 1:n_vars,
  Eigenvalue = round(eigenvalues, 4),
  Condition_Index = round(condition_indices, 2)
)
print(condition_df)
cat("\nInterpretation: CI > 15 suggests weak dependencies, CI > 30 suggests strong collinearity\n\n")

cat("VARIANCE PROPORTIONS:\n")
cat("---------------------\n")
print(round(variance_props, 3))
cat("\nInterpretation: High proportions (>0.5) on high condition indices indicate problematic variables\n\n")

cat("VARIANCE INFLATION FACTORS:\n")
cat("---------------------------\n")
vif_df <- data.frame(
  Variable = names(vif_values),
  VIF = round(vif_values, 2),
  Tolerance = round(1/vif_values, 3)
)
print(vif_df)
cat("\nInterpretation: VIF > 5 suggests moderate collinearity, VIF > 10 suggests severe collinearity\n")
sink()

cat("Collinearity diagnostics saved to:", collinearity_file, "\n")

png(file.path(plots_dir, "wine_vif_plot.png"), width = 1000, height = 600, res = 100)
vif_plot_df <- data.frame(
  Variable = factor(names(vif_values), levels = names(vif_values)[order(vif_values)]),
  VIF = vif_values[order(vif_values)]
)

print(
  ggplot(vif_plot_df, aes(x = Variable, y = VIF)) +
    geom_col(fill = cb_cols[4], alpha = 0.7) +
    geom_hline(yintercept = 5, color = "red", linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = 10, color = "red", linetype = "solid", alpha = 0.7) +
    labs(title = "Variance Inflation Factors - Red Wine Predictors",
         x = "Predictor Variables", y = "VIF") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    annotate("text", x = 2, y = 5.5, label = "VIF = 5", color = "red", size = 3) +
    annotate("text", x = 2, y = 10.5, label = "VIF = 10", color = "red", size = 3)
)
dev.off()

png(file.path(plots_dir, "wine_condition_indices.png"), width = 800, height = 600, res = 100)
print(
  ggplot(condition_df, aes(x = Component, y = Condition_Index)) +
    geom_col(fill = cb_cols[5], alpha = 0.7) +
    geom_hline(yintercept = 15, color = "orange", linetype = "dashed") +
    geom_hline(yintercept = 30, color = "red", linetype = "dashed") +
    labs(title = "Condition Indices - Red Wine Data",
         x = "Principal Component", y = "Condition Index") +
    theme_minimal() +
    annotate("text", x = 2, y = 16, label = "CI = 15", color = "orange", size = 3) +
    annotate("text", x = 2, y = 31, label = "CI = 30", color = "red", size = 3)
)
dev.off()

cat("VIF and condition indices plots saved.\n")

model_formula <- quality ~ alcohol + volatile.acidity + sulphates + free.sulfur.dioxide + 
  chlorides

model_polr_wine <- polr(model_formula, data = wine_data, Hess = TRUE)
cat("POLR model fitted successfully\n")

model_clm_wine <- clm(model_formula, data = wine_data, link = "logit")
cat("CLM model fitted successfully\n")

cat("\n=== POLR Model Summary ===\n")
print(summary(model_polr_wine))

cat("\n=== CLM Model Summary ===\n")
print(summary(model_clm_wine))

cat("\n=== Model Comparison ===\n")
cat("POLR AIC:", AIC(model_polr_wine), "\n")
cat("CLM AIC:", AIC(model_clm_wine), "\n")
cat("POLR Log-Likelihood:", logLik(model_polr_wine), "\n")
cat("CLM Log-Likelihood:", logLik(model_clm_wine), "\n")

cat("\n--- Brant Test for Proportional Odds Assumption ---\n")

brant_results_file <- file.path(results_dir, "brant_test_wine.txt")
sink(brant_results_file)

cat("======================================\n")
cat("Brant Test for Wine Quality Model\n")
cat("======================================\n\n")

brant_standard <- capture.output({
  standard_result <- brant(model_polr_wine)
})
cat(paste(brant_standard, collapse = "\n"))

cat("\n\nBy Variable:\n")
brant_byvar <- capture.output({
  byvar_result <- brant(model_polr_wine, by.var = TRUE)
})
cat(paste(brant_byvar, collapse = "\n"))

sink()
cat("Brant test results saved to:", brant_results_file, "\n")

cat("\n--- Nominal Effects Test ---\n")

nominal_results_file <- file.path(results_dir, "nominal_test_wine.txt")
sink(nominal_results_file)

cat("======================================\n")
cat("Nominal Test for Wine Quality Model\n")
cat("======================================\n\n")

print(nominal_test(model_clm_wine))

sink()
cat("Nominal test results saved to:", nominal_results_file, "\n")

cat("\n=== BRANT TEST SUMMARY ===\n")
print(brant(model_polr_wine))

cat("\n=== NOMINAL TEST SUMMARY ===\n")
print(nominal_test(model_clm_wine))


model_formula_2 <- quality ~ alcohol + volatile.acidity + sulphates + free.sulfur.dioxide + 
  chlorides

model_clm_partial <- clm(model_formula_2, data = wine_data, link = "logit",
                         nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)
cat("Partial proportional odds CLM model fitted successfully\n")

model_clm_partial_conservative <- clm(model_formula_2, data = wine_data, link = "logit",
                                      nominal = ~ alcohol)

model_clm_partial_liberal <- clm(model_formula_2, data = wine_data, link = "logit",
                                 nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)

cat("Alternative partial proportional odds models fitted\n")

cat("\n=== MODEL COMPARISON ===\n")
cat("Standard CLM AIC:", AIC(model_clm_wine), "\n")
cat("Partial PO (conservative) AIC:", AIC(model_clm_partial_conservative), "\n") 
cat("Partial PO (liberal) AIC:", AIC(model_clm_partial_liberal), "\n")

lr_test_conservative <- anova(model_clm_wine, model_clm_partial_conservative)
lr_test_liberal <- anova(model_clm_wine, model_clm_partial_liberal)

print("Conservative partial model vs standard:")
print(lr_test_conservative)

print("Liberal partial model vs standard:")
print(lr_test_liberal)

set.seed(123)
nsim_value <- 100
residuals_wine_v3 <- sure::resids(model_clm_partial_liberal, method = "jitter", nsim = nsim_value)

qq_plot <- autoplot.resid(residuals_wine_v3, what = "qq",
                          distribution = qnorm,
                          qqpoint.color = cb_cols[3], 
                          qqline.color = cb_cols[2],  
                          alpha = 0.8) +
  labs(title = "Wine Quality - Surrogate Residuals QQ Plot")

p1 <- ggplot(data.frame(residual = residuals_wine_v3, 
                        alcohol = wine_data$alcohol),
             aes(x = alcohol, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) +  
  labs(title = "Residuals vs Alcohol", x = "Alcohol (%)", y = "Surrogate residual") +
  theme_minimal()

p2 <- ggplot(data.frame(residual = residuals_wine_v3, 
                        volatile_acidity = wine_data$volatile.acidity),
             aes(x = volatile_acidity, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) + 
  labs(title = "Residuals vs Volatile Acidity", x = "Volatile Acidity", y = "Surrogate residual") +
  theme_minimal()

p3 <- ggplot(data.frame(residual = residuals_wine_v3, 
                        sulphates = wine_data$sulphates),
             aes(x = sulphates, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) +  
  labs(title = "Residuals vs Sulphates", x = "Sulphates", y = "Surrogate residual") +
  theme_minimal()

p4 <- ggplot(data.frame(residual = residuals_wine_v3, 
                        free_sulfur = wine_data$free.sulfur.dioxide),
             aes(x = free_sulfur, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) +  
  labs(title = "Residuals vs Free Sulfur Dioxide", x = "Free Sulfur Dioxide", y = "Surrogate residual") +
  theme_minimal()

p5 <- ggplot(data.frame(residual = residuals_wine_v3, 
                        chlorides = wine_data$chlorides),
             aes(x = chlorides, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) +  
  labs(title = "Residuals vs Chlorides", x = "Chlorides", y = "Surrogate residual") +
  theme_minimal()

png(file.path(plots_dir, "surrogate_residuals_wine_qq.png"), width = 800, height = 600, res = 100)
print(qq_plot)
dev.off()

png(file.path(plots_dir, "surrogate_residuals_wine_sulphates_chlorides.png"), width = 1200, height = 600, res = 100)
grid.arrange(p3, p5, ncol = 2)
dev.off()

png(file.path(plots_dir, "surrogate_residuals_wine_alcohol_volatileacidity_freesulfur.png"), width = 1800, height = 600, res = 100)
grid.arrange(p1, p2, p4, ncol = 3)
dev.off()

cat("Surrogate residual plots saved to:", plots_dir, "\n")
cat("- QQ plot: surrogate_residuals_wine_qq.png\n")
cat("- Sulphates & Chlorides: surrogate_residuals_wine_sulphates_chlorides.png\n")
cat("- Alcohol, Volatile Acidity & Free Sulfur: surrogate_residuals_wine_alcohol_volatileacidity_freesulfur.png\n")

outlier_indices <- which(
  wine_data$sulphates > 1.5 | 
    wine_data$free.sulfur.dioxide > 60 |
    wine_data$chlorides > 0.4
)

wine_data_clean <- wine_data[-outlier_indices, ]
model_without_outliers <- clm(model_formula_2, data = wine_data_clean, link = "logit",
                              nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)

cat("\n--- 8. Influence and Leverage Diagnostics ---\n")

n <- nrow(wine_data)
p_reg <- length(coef(model_clm_partial_liberal))  
p_intercepts <- length(model_clm_partial_liberal$Theta)  
p_all <- p_reg + p_intercepts 

cutoff_cooksd <- 4/n
cutoff_dfbetas <- 2/sqrt(n)

cat("Sample size (n):", n, "\n")
cat("Number of parameters (p):", p_all, "\n")
cat("Cook's Distance cutoff:", cutoff_cooksd, "\n")
cat("DFBETAS cutoff:", cutoff_dfbetas, "\n\n")


sample_size <- 1599
sample_indices <- sample(1:n, sample_size)
cooksd <- numeric(n)
cooksd[] <- NA

cat("Calculating Cook's distance for sample of", sample_size, "observations...\n")

for(idx in seq_along(sample_indices)) {
  i <- sample_indices[idx]
  if(idx %% 10 == 0) cat("Progress:", round(idx/sample_size*100), "%\n")
  
  tryCatch({

    data_subset <- wine_data[-i, ]
    
    model_i <- clm(formula(model_clm_partial_liberal), data = data_subset, 
                   link = "logit", nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)
    
    orig_coef <- coef(model_clm_partial_liberal)
    coef_i <- coef(model_i)
    common_coefs <- intersect(names(orig_coef), names(coef_i))
    
    delta_coef <- orig_coef[common_coefs] - coef_i[common_coefs]
    vcov_subset <- vcov(model_clm_partial_liberal)[common_coefs, common_coefs]
    
    cooksd[i] <- t(delta_coef) %*% solve(vcov_subset) %*% delta_coef / length(common_coefs)
    
  }, error = function(e) {
    cooksd[i] <- NA
    cat("  Error for observation", i, ":", e$message, "\n")
  })
}

extreme_alcohol <- which(wine_data$alcohol > quantile(wine_data$alcohol, 0.95) | 
                           wine_data$alcohol < quantile(wine_data$alcohol, 0.05))
extreme_volatile_acidity <- which(wine_data$volatile.acidity > quantile(wine_data$volatile.acidity, 0.95) | 
                            wine_data$volatile.acidity < quantile(wine_data$volatile.acidity, 0.05))
extreme_sulphates <- which(wine_data$sulphates > quantile(wine_data$sulphates, 0.95) | 
                             wine_data$sulphates < quantile(wine_data$sulphates, 0.05))
extreme_free_sulfur_dioxide <- which(wine_data$free.sulfur.dioxide > quantile(wine_data$free.sulfur.dioxide, 0.95) | 
                               wine_data$free.sulfur.dioxide < quantile(wine_data$free.sulfur.dioxide, 0.05))
extreme_chlorides <- which(wine_data$chlorides > quantile(wine_data$chlorides, 0.95) | 
                             wine_data$chlorides < quantile(wine_data$chlorides, 0.05))

high_leverage_proxy <- unique(c(extreme_alcohol, extreme_volatile_acidity, extreme_sulphates, 
                                extreme_free_sulfur_dioxide, extreme_chlorides))

influence_file <- file.path(results_dir, "influence_diagnostics_wine.txt")
sink(influence_file)

cat("INFLUENCE AND LEVERAGE DIAGNOSTICS - WINE QUALITY\n")
cat("=================================================\n\n")
cat("Sample size (n):", n, "\n")
cat("Number of parameters (p):", p_all, "\n\n")

cat("Recommended cutoff values:\n")
cat("  Cook's Distance:", cutoff_cooksd, "\n")
cat("  DFBETAS:", cutoff_dfbetas, "\n\n")

cat("EXTREME VALUES (LEVERAGE PROXY)\n")
cat("-------------------------------\n")
cat("Observations in top/bottom 5% of any predictor:", length(high_leverage_proxy), "\n")
cat("Percentage of data:", round(length(high_leverage_proxy)/n*100, 1), "%\n\n")

if(!all(is.na(cooksd))) {
  cat("COOK'S DISTANCE (subset of", sample_size, "observations)\n")
  cat("---------------------\n")
  cooksd_valid <- cooksd[!is.na(cooksd)]
  if(length(cooksd_valid) > 0) {
    cat("  Min:", min(cooksd_valid), "\n")
    cat("  Median:", median(cooksd_valid), "\n")
    cat("  Max:", max(cooksd_valid), "\n")
    cat("  Influential points:", sum(cooksd_valid > cutoff_cooksd), "\n")
    
    influential_obs <- which(cooksd > cutoff_cooksd)
    if(length(influential_obs) > 0) {
      cat("\nMost influential observations:\n")
      for(obs in influential_obs[1:min(10, length(influential_obs))]) {
        cat("  Observation", obs, ": Cook's D =", round(cooksd[obs], 6), "\n")
      }
    }
  }
}

sink()
cat("Influence diagnostics saved to:", influence_file, "\n")

png(file.path(plots_dir, "influence_plots_wine_enhanced.png"), width = 1400, height = 1000)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))

if(!all(is.na(cooksd))) {
  cooksd_indices <- which(!is.na(cooksd))
  plot(cooksd_indices, cooksd[cooksd_indices], type = "h",
       main = "Cook's Distance - Wine Quality", xlab = "Observation Index", ylab = "Cook's Distance",
       col = ifelse(cooksd[cooksd_indices] > cutoff_cooksd, cb_cols[7], cb_cols[1]))
  abline(h = cutoff_cooksd, col = cb_cols[7], lty = 2, lwd = 2)
  text(x = max(cooksd_indices)*0.7, y = cutoff_cooksd*1.1, 
       labels = paste("Cutoff =", round(cutoff_cooksd, 6)), col = cb_cols[7])
}

predictors_info <- list(
  list(var = wine_data$alcohol, name = "alcohol", title = "Alcohol Distribution", ylabel = "Alcohol (%)"),
  list(var = wine_data$volatile.acidity, name = "volatile.acidity", title = "Volatile Acidity Distribution", ylabel = "Volatile Acidity"),
  list(var = wine_data$sulphates, name = "sulphates", title = "Sulphates Distribution", ylabel = "Sulphates"),
  list(var = wine_data$free.sulfur.dioxide, name = "free.sulfur.dioxide", title = "Free Sulfur Dioxide Distribution", ylabel = "Free Sulfur Dioxide"),
  list(var = wine_data$chlorides, name = "chlorides", title = "Chlorides Distribution", ylabel = "Chlorides")
)

for(i in 1:5) {
  pred_data <- predictors_info[[i]]
  extreme_indices <- get(paste0("extreme_", gsub("\\.", "_", pred_data$name)))
  
  p05 <- quantile(pred_data$var, 0.05)
  p95 <- quantile(pred_data$var, 0.95)
  
  plot(pred_data$var, type = "p", pch = 16, cex = 2, main = pred_data$title, ylab = pred_data$ylabel,
       xlab = "Observation Index", col = ifelse(1:n %in% extreme_indices, cb_cols[7], cb_cols[1]),
       cex.main = 3, cex.lab = 2.5, cex.axis = 2)
  
  abline(h = p05, col = cb_cols[2], lty = 2, lwd = 2)
  abline(h = p95, col = cb_cols[2], lty = 2, lwd = 2)
  
  text(x = n*0.02, y = p05, labels = paste("5th percentile:", round(p05, 3)), 
       col = cb_cols[2], pos = 3, cex = 2)
  text(x = n*0.02, y = p95, labels = paste("95th percentile:", round(p95, 3)), 
       col = cb_cols[2], pos = 1, cex = 2)

  text(x = n*0.7, y = max(pred_data$var)*0.9, 
       labels = paste("Extreme values:", length(extreme_indices)), 
       col = cb_cols[7], cex = 3, font = 2)
}

dev.off()

cat("Influence plots saved to:", plots_dir, "\n")


all_extreme_indices <- unique(c(extreme_alcohol, extreme_volatile_acidity, extreme_sulphates, 
                                extreme_free_sulfur_dioxide, extreme_chlorides))

wine_data_noextremes <- wine_data[-all_extreme_indices, ]

model_clm_partial_liberal_noextreme <- clm(model_formula_2, data = wine_data_noextremes, link = "logit",
                                 nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)

set.seed(123)
nsim_value <- 100
residuals_model_clm_partial_liberal_noextreme <- sure::resids(model_clm_partial_liberal_noextreme, method = "jitter", nsim = nsim_value)

qq_plot <- autoplot.resid(residuals_model_clm_partial_liberal_noextreme, what = "qq",
                          distribution = qnorm,
                          qqpoint.color = cb_cols[3],  
                          qqline.color = cb_cols[2],  
                          alpha = 0.8) +
  labs(title = "Wine Quality - Surrogate Residuals QQ Plot (Clean Data)")


p1 <- ggplot(data.frame(residual = residuals_model_clm_partial_liberal_noextreme, 
                        alcohol = wine_data_noextremes$alcohol),
             aes(x = alcohol, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) + 
  labs(title = "Residuals vs Alcohol", x = "Alcohol (%)", y = "Surrogate residual") +
  theme_minimal()

p2 <- ggplot(data.frame(residual = residuals_model_clm_partial_liberal_noextreme, 
                        volatile_acidity = wine_data_noextremes$volatile.acidity),
             aes(x = volatile_acidity, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) +  
  labs(title = "Residuals vs Volatile Acidity", x = "Volatile Acidity", y = "Surrogate residual") +
  theme_minimal()

p3 <- ggplot(data.frame(residual = residuals_model_clm_partial_liberal_noextreme, 
                        sulphates = wine_data_noextremes$sulphates),
             aes(x = sulphates, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) +  
  labs(title = "Residuals vs Sulphates", x = "Sulphates", y = "Surrogate residual") +
  theme_minimal()

p4 <- ggplot(data.frame(residual = residuals_model_clm_partial_liberal_noextreme, 
                        free_sulfur = wine_data_noextremes$free.sulfur.dioxide),
             aes(x = free_sulfur, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) +  
  labs(title = "Residuals vs Free Sulfur Dioxide", x = "Free Sulfur Dioxide", y = "Surrogate residual") +
  theme_minimal()

p5 <- ggplot(data.frame(residual = residuals_model_clm_partial_liberal_noextreme, 
                        chlorides = wine_data_noextremes$chlorides),
             aes(x = chlorides, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) + 
  labs(title = "Residuals vs Chlorides", x = "Chlorides", y = "Surrogate residual") +
  theme_minimal()

png(file.path(plots_dir, "surrogate_residuals_wine_cleandata_noextremes_qq.png"), width = 800, height = 600, res = 100)
print(qq_plot)
dev.off()

png(file.path(plots_dir, "surrogate_residuals_wine_cleandata_noextremes_sulphates_chlorides.png"), width = 1200, height = 600, res = 100)
grid.arrange(p3, p5, ncol = 2)
dev.off()

png(file.path(plots_dir, "surrogate_residuals_wine_cleandata_noextremes_alcohol_volatileacidity_freesulfur.png"), width = 1800, height = 600, res = 100)
grid.arrange(p1, p2, p4, ncol = 3)
dev.off()

cat("Surrogate residual plots saved to:", plots_dir, "\n")
cat("- QQ plot (clean data): surrogate_residuals_wine_cleandata_noextremes_qq.png\n")
cat("- Sulphates & Chlorides (clean data): surrogate_residuals_wine_cleandata_noextremes_sulphates_chlorides.png\n")
cat("- Alcohol, Volatile Acidity & Free Sulfur (clean data): surrogate_residuals_wine_cleandata_noextremes_alcohol_volatileacidity_freesulfur.png\n")


model_quadratic <- clm(quality ~ alcohol + volatile.acidity + free.sulfur.dioxide + 
                         sulphates + I(sulphates^2) + 
                         chlorides + I(chlorides^2), 
                       data = wine_data,
                       nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)

model_log <- clm(quality ~ alcohol + volatile.acidity + free.sulfur.dioxide + 
                   log(sulphates) + log(chlorides), 
                 data = wine_data,
                 nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)


model_hybrid <- clm(quality ~ alcohol + volatile.acidity + free.sulfur.dioxide + 
                      sulphates + I(sulphates^2) + log(chlorides), 
                    data = wine_data,
                    nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)

model_selective <- clm(quality ~ alcohol + volatile.acidity + free.sulfur.dioxide + 
                         sulphates + I(sulphates^2) + chlorides, 
                       data = wine_data,
                       nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)


model_selective_probit <- clm(quality ~ alcohol + volatile.acidity + free.sulfur.dioxide + 
                      sulphates + I(sulphates^2) + chlorides, 
                    data = wine_data,
                    link = "probit",
                    nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)

model_selective_cloglog <- clm(quality ~ alcohol + volatile.acidity + free.sulfur.dioxide + 
                       sulphates + I(sulphates^2) + chlorides, 
                     data = wine_data,
                     link = "cloglog", 
                     nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)

model_selective_cauchit <- clm(quality ~ alcohol + volatile.acidity + free.sulfur.dioxide + 
                       sulphates + I(sulphates^2) + chlorides, 
                     data = wine_data,
                     link = "cauchit",
                     nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)


model_probit <- clm(quality ~ alcohol + volatile.acidity + free.sulfur.dioxide + 
                                sulphates + chlorides, 
                              data = wine_data,
                              link = "probit",
                              nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)

model_cloglog <- clm(quality ~ alcohol + volatile.acidity + free.sulfur.dioxide + 
                                 sulphates + chlorides, 
                               data = wine_data,
                               link = "cloglog",
                               nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)

model_cauchit <- clm(quality ~ alcohol + volatile.acidity + free.sulfur.dioxide + 
                                 sulphates + chlorides, 
                               data = wine_data,
                               link = "cauchit",  
                               nominal = ~ alcohol + volatile.acidity + free.sulfur.dioxide)

summary_file <- file.path(results_dir, "model_summaries_wine.txt")
sink(summary_file)

cat("===============================================\n")
cat("MODEL SUMMARIES\n")
cat("===============================================\n\n")

cat("QUADRATIC MODEL:\n")
cat("===============================================\n")
print(summary(model_quadratic))

cat("\n\nLOG MODEL:\n")
cat("===============================================\n")
print(summary(model_log))

cat("\n\nHYBRID MODEL:\n")
cat("===============================================\n")
print(summary(model_hybrid))

cat("\n\nSELECTIVE MODEL:\n")
cat("===============================================\n")
print(summary(model_selective))

cat("\n\nSELECTIVE PROBIT MODEL:\n")
cat("===============================================\n")
print(summary(model_selective_probit))

cat("\n\nSELECTIVE CLOGLOG MODEL:\n")
cat("===============================================\n")
print(summary(model_selective_cloglog))

cat("\n\nSELECTIVE CAUCHIT MODEL:\n")
cat("===============================================\n")
print(summary(model_selective_cauchit))

cat("\n\nPROBIT MODEL:\n")
cat("===============================================\n")
print(summary(model_probit))

cat("\n\nCLOGLOG MODEL:\n")
cat("===============================================\n")
print(summary(model_cloglog))

cat("\n\nCAUCHIT MODEL:\n")
cat("===============================================\n")
print(summary(model_cauchit))


sink()

cat("Model summaries saved to 'model_summaries.txt'\n")


set.seed(123)
nsim_value <- 100
residuals_model_selective <- sure::resids(model_selective, method = "jitter", nsim = nsim_value)


p1 <- ggplot(data.frame(residual = residuals_model_selective, 
                        alcohol = wine_data$alcohol),
             aes(x = alcohol, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) + 
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) +
  labs(title = "Residuals vs Alcohol", x = "Alcohol (%)", y = "Surrogate residual") +
  theme_minimal()

p2 <- ggplot(data.frame(residual = residuals_model_selective, 
                        volatile_acidity = wine_data$volatile.acidity),
             aes(x = volatile_acidity, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) + 
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) +
  labs(title = "Residuals vs Volatile Acidity", x = "Volatile Acidity", y = "Surrogate residual") +
  theme_minimal()

p3 <- ggplot(data.frame(residual = residuals_model_selective, 
                        sulphates = wine_data$sulphates),
             aes(x = sulphates, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) + 
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) + 
  labs(title = "Residuals vs Sulphates", x = "Sulphates", y = "Surrogate residual") +
  theme_minimal()

p4 <- ggplot(data.frame(residual = residuals_model_selective, 
                        free_sulfur = wine_data$free.sulfur.dioxide),
             aes(x = free_sulfur, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) + 
  labs(title = "Residuals vs Free Sulfur Dioxide", x = "Free Sulfur Dioxide", y = "Surrogate residual") +
  theme_minimal()

p5 <- ggplot(data.frame(residual = residuals_model_selective, 
                        chlorides = wine_data$chlorides),
             aes(x = chlorides, y = residual)) +
  geom_point(color = cb_cols[3], alpha = 0.6) +  
  geom_smooth(method = "loess", color = cb_cols[2], se = TRUE) +  
  labs(title = "Residuals vs Chlorides", x = "Chlorides", y = "Surrogate residual") +
  theme_minimal()

png(file.path(plots_dir, "surrogate_residuals_wine_sel.png"), width = 800, height = 600, res = 100)
print(qq_plot)
dev.off()

png(file.path(plots_dir, "surrogate_residuals_wine_sel_sulphates_chlorides.png"), width = 1200, height = 600, res = 100)
grid.arrange(p3, p5, ncol = 2)
dev.off()

png(file.path(plots_dir, "surrogate_residuals_wine_sel_alcohol_volatileacidity_freesulfur.png"), width = 1800, height = 600, res = 100)
grid.arrange(p1, p2, p4, ncol = 3)
dev.off()