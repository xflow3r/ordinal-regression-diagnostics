run_col_tests <- function(model, data, model_name, model_type, 
                          results_dir, plots_dir) {
  
  deeper_results_dir <- file.path(results_dir, "collinearity_checks")
  dir.create(deeper_results_dir, showWarnings = FALSE, recursive = TRUE)
  
  txt_filename <- paste0(deeper_results_dir, "/collinearity_checks_", model_type, "_", model_name, ".txt")
  
  tm <- stats::terms(formula(model))
  tm_noresp <- stats::delete.response(tm)
  X <- stats::model.matrix(tm_noresp, data = data)
  
  col_diag <- tryCatch(
    performance::check_collinearity(model),
    error = function(e) NULL
  )
  
  cond_number <- tryCatch(
    base::kappa(X, exact = TRUE),
    error = function(e) NA_real_
  )
  

  eig_values <- tryCatch({
    mat <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    eigen(stats::cor(mat), symmetric = TRUE, only.values = TRUE)$values
  }, error = function(e) NA_real_)
  
  
  cor_matrix <- tryCatch({
    mat <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    stats::cor(mat)
  }, error = function(e) matrix(NA_real_, ncol = 0, nrow = 0))
  
  sink(txt_filename)
  
  cat("==============================================================\n")
  cat("COLLINEARITY DIAGNOSTICS FOR", toupper(model_name), "MODEL (", toupper(model_type), ")\n")
  cat("==============================================================\n\n")
  
  cat("1. VARIANCE INFLATION FACTORS (VIF)\n")
  cat("------------------------------------------\n")
  if (!is.null(col_diag)) {
    df_col <- as.data.frame(col_diag$VIF)
    print(df_col)
    
    cat("\nInterpretation guidelines:\n")
    cat("  - VIF < 5: No indication of multicollinearity\n")
    cat("  - VIF 5-10: Moderate multicollinearity\n")
    cat("  - VIF > 10: High multicollinearity\n\n")
    
    high_vif <- sum(df_col$VIF > 10, na.rm = TRUE)
    mod_vif <- sum(df_col$VIF > 5 & df_col$VIF <= 10, na.rm = TRUE)
    
    if (high_vif > 0) {
      cat("ALERT: ", high_vif, " variable(s) have high VIF (>10), indicating serious multicollinearity.\n")
    } else if (mod_vif > 0) {
      cat("NOTE: ", mod_vif, " variable(s) have moderate VIF (5-10), suggesting some multicollinearity.\n")
    } else {
      cat("All variables have acceptable VIF values (<5).\n")
    }
  } else {
    cat("VIF calculation failed or is not applicable for this model type.\n")
  }
  cat("\n")
  
  cat("2. GLOBAL CONDITION NUMBER\n")
  cat("------------------------------------------\n")
  cat("Condition number of design matrix:", cond_number, "\n\n")
  cat("Interpretation guidelines:\n")
  cat("  - < 10: No serious collinearity issues\n")
  cat("  - 10-30: Moderate collinearity\n")
  cat("  - > 30: Severe collinearity\n\n")
  
  if (!is.na(cond_number)) {
    if (cond_number > 30) {
      cat("ALERT: Condition number is >30, indicating severe multicollinearity.\n")
    } else if (cond_number > 10) {
      cat("NOTE: Condition number is between 10-30, suggesting moderate multicollinearity.\n")
    } else {
      cat("Condition number is <10, suggesting no serious multicollinearity.\n")
    }
  }
  cat("\n")
  
  cat("3. EIGENVALUES OF PREDICTOR CORRELATION MATRIX\n")
  cat("------------------------------------------\n")
  if (!all(is.na(eig_values))) {
    cat("Eigenvalues:", paste(round(eig_values, 4), collapse = ", "), "\n\n")
    min_eig <- min(eig_values)
    cat("Minimum eigenvalue:", min_eig, "\n\n")
    cat("Interpretation guidelines:\n")
    cat("  - Minimum eigenvalue near 0 indicates collinearity\n")
    cat("  - General rule: eigenvalues <0.01 may indicate problems\n\n")
    
    if (min_eig < 0.01) {
      cat("ALERT: Minimum eigenvalue is <0.01, suggesting collinearity.\n")
    } else if (min_eig < 0.05) {
      cat("NOTE: Minimum eigenvalue is low (<0.05), potential collinearity.\n")
    } else {
      cat("All eigenvalues appear reasonable, no strong indication of collinearity from eigenvalue analysis.\n")
    }
  } else {
    cat("Eigenvalue calculation failed.\n")
  }
  cat("\n")
  
  cat("4. CORRELATION MATRIX OF PREDICTORS\n")
  cat("------------------------------------------\n")
  if (nrow(cor_matrix) > 0 && ncol(cor_matrix) > 0) {

    pred_names <- colnames(X)[colnames(X) != "(Intercept)"]
    dimnames(cor_matrix) <- list(pred_names, pred_names)
    
    print(round(cor_matrix, 3))
    cat("\n")
    
    high_cors <- which(abs(cor_matrix) > 0.8 & abs(cor_matrix) < 1, arr.ind = TRUE)
    
    if (nrow(high_cors) > 0) {
      cat("High correlations (>0.8) found between predictors:\n")
      for (i in 1:nrow(high_cors)) {
        if (high_cors[i, 1] < high_cors[i, 2]) {
          var1 <- pred_names[high_cors[i, 1]]
          var2 <- pred_names[high_cors[i, 2]]
          corr <- cor_matrix[high_cors[i, 1], high_cors[i, 2]]
          cat(sprintf("  - %s and %s: %.3f\n", var1, var2, corr))
        }
      }
    } else {
      cat("No high correlations (>0.8) found between predictors.\n")
    }
  } else {
    cat("Correlation matrix calculation failed.\n")
  }
  cat("\n")
  
  cat("5. SUMMARY AND RECOMMENDATIONS\n")
  cat("------------------------------------------\n")
  
  has_vif_issue <- (!is.null(col_diag) && any(df_col$VIF > 5, na.rm = TRUE))
  has_cond_issue <- (!is.na(cond_number) && cond_number > 10)
  has_eig_issue <- (!all(is.na(eig_values)) && min(eig_values) < 0.05)
  has_cor_issue <- (nrow(cor_matrix) > 0 && ncol(cor_matrix) > 0 && 
                      any(abs(cor_matrix) > 0.8 & abs(cor_matrix) < 1))
  
  total_issues <- sum(has_vif_issue, has_cond_issue, has_eig_issue, has_cor_issue)
  
  if (total_issues == 0) {
    cat("No collinearity issues detected. The model appears to have well-behaved predictors.\n")
  } else {
    cat("Collinearity issues detected:\n")
    if (has_vif_issue) cat("  - High variance inflation factors\n")
    if (has_cond_issue) cat("  - High condition number\n")
    if (has_eig_issue) cat("  - Low eigenvalues\n")
    if (has_cor_issue) cat("  - High correlations between predictors\n")
  }
  
  cat("\n==============================================================\n")
  cat("END OF COLLINEARITY DIAGNOSTICS\n")
  cat("==============================================================\n")
  
  sink()
}