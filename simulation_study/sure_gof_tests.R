run_all_gof_tests <- function(model, model_name, model_type, results_dir, plots_dir, residuals_model, nsim = 100) {
  if (!requireNamespace("sure", quietly = TRUE)) {
    stop("Package 'sure' is needed for this function to work. Please install it.")
  }
  
  theme_set(
    theme_classic(base_size = 14) +
      theme(
        axis.text  = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = .5),
        legend.position = "bottom")
  )
  
  cb_cols <- c(
    black   = "#000000",
    orange  = "#E69F00",
    sky     = "#56B4E9",
    green   = "#009E73",
    yellow  = "#F0E442",
    blue    = "#0072B2",
    vermil  = "#D55E00",
    purple  = "#CC79A7"
  )
  
  par(
    fg  = cb_cols["black"],
    col = cb_cols["sky"],
    pch = 16
  )
  
  ks_test <- sure::gof(model, nsim = nsim, test = "ks")
  ad_test <- sure::gof(model, nsim = nsim, test = "ad")
  cvm_test <- sure::gof(model, nsim = nsim, test = "cvm")
  
  results <- list(
    ks_test = ks_test,
    ad_test = ad_test,
    cvm_test = cvm_test
  )
  
  ks_significant <- mean(ks_test < 0.05)
  ad_significant <- mean(ad_test < 0.05)
  cvm_significant <- mean(cvm_test < 0.05)
  
  summary_df <- data.frame(
    Test = c("Kolmogorov-Smirnov", "Anderson-Darling", "Cramer-Von Mises"),
    Mean_pvalue = c(mean(ks_test), mean(ad_test), mean(cvm_test)),
    Prop_significant = c(ks_significant, ad_significant, cvm_significant)
  )
  
  deeper_results_dir <- file.path(results_dir, "sure_gof_checks")
  if (!dir.exists(deeper_results_dir)) dir.create(deeper_results_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  deeper_plots_dir <- file.path(plots_dir, "sure_gof_plots")
  if (!dir.exists(deeper_plots_dir)) dir.create(deeper_plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  txt_filename <- paste0("gof_summary_", model_type, "_", model_name, ".txt")
  txt_path <- file.path(deeper_results_dir, txt_filename)
  
  summary_text <- paste0(
    "Goodness-of-Fit Test Results for ", model_type, " ", model_name, "\n",
    "Number of simulations: ", nsim, "\n",
    "Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
    "Test                  | Mean p-value | Proportion Significant\n",
    "--------------------- | ------------ | ---------------------\n",
    sprintf("%-21s | %12.4f | %12.4f\n", "Kolmogorov-Smirnov", mean(ks_test), ks_significant),
    sprintf("%-21s | %12.4f | %12.4f\n", "Anderson-Darling", mean(ad_test), ad_significant),
    sprintf("%-21s | %12.4f | %12.4f\n", "Cramer-Von Mises", mean(cvm_test), cvm_significant)
  )
  
  writeLines(summary_text, txt_path)
  
  png(file.path(deeper_plots_dir, paste0("gof_plots_", model_type, "_", model_name, ".png")), 
      width = 10, height = 8, units = "in", res = 300)
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  plot(seq(0, 1, length.out = 100), seq(0, 1, length.out = 100), 
       type = "l", lwd = 2, col = "black", 
       xlab = "Theoretical", ylab = "Empirical",
       main = "Combined P-value ECDFs")
  
  lines(sort(ks_test), (1:nsim)/nsim, col = cb_cols["blue"], lwd = 2)
  lines(sort(ad_test), (1:nsim)/nsim, col = cb_cols["vermil"], lwd = 2)
  lines(sort(cvm_test), (1:nsim)/nsim, col = cb_cols["green"], lwd = 2)
  
  legend("bottomright", 
         legend = c("Theoretical", "Kolmogorov-Smirnov", "Anderson-Darling", "Cramer-Von Mises"),
         col = c(cb_cols["black"], cb_cols["blue"], cb_cols["vermil"], cb_cols["green"]),
         lwd = 2, cex = 0.8)
  
  plot(ks_test, main = "Kolmogorov-Smirnov Test", col = cb_cols["blue"])
  plot(ad_test, main = "Anderson-Darling Test", col = cb_cols["vermil"])
  plot(cvm_test, main = "Cramer-Von Mises Test", col = cb_cols["green"])
  
  dev.off()
  
  return(results)
}