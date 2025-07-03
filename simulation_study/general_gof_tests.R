run_gof_tests <- function(model, data, model_name, model_type, 
                          results_dir, plots_dir) {
  
if (!dir.exists(results_dir)) dir.create(results_dir, recursive=TRUE)
if (!dir.exists(plots_dir))   dir.create(plots_dir,   recursive=TRUE)
  
  
deeper_results_dir <- file.path(results_dir, "general_gof_tests")
if (!dir.exists(deeper_results_dir))   dir.create(deeper_results_dir,   recursive=TRUE)

model_formula <- formula(model)
terms_list <- attr(terms(model_formula), "term.labels")
categorical_predictors <- character(0)

for (term in terms_list) {
  if (!grepl(":", term) && term %in% names(data)) {
    if (is.factor(data[[term]]) || is.ordered(data[[term]])) {
      categorical_predictors <- c(categorical_predictors, term)
    }
  }
}

txt_filename <- paste0(deeper_results_dir, "/general_gof_tests_", model_type, "_", model_name, ".txt")

res_lip1  <- lipsitz.test(model, g=10)
res_lip2  <- lipsitz(model, group=10)

if (length(categorical_predictors) > 0) {
  res_pr1   <- pulkrob.chisq(model, catvars=categorical_predictors)
  res_pr2   <- pulkrob.deviance(model, catvars=categorical_predictors)
  res_pr3   <- pulkroben(model, test="chisq")
  res_pr4   <- pulkroben(model, test="deviance")
}

obs <- data[[as.character(model_formula[[2]])]]
exp <- predict(model, newdata = data, type = "prob")
res_hlm1 <- generalhoslem::logitgof(obs, exp, g = 10, ord = TRUE)
res_hlm2 <- gofcat::hosmerlem(model, group = 10)

cat("Goodness-of-Fit Tests for Model:", model_name, "\n", file = txt_filename)
cat("Model Type:", model_type, "\n", file = txt_filename, append = TRUE)
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = txt_filename, append = TRUE)

cat("Categorical Predictor Check:\n", file = txt_filename, append = TRUE)
if (length(categorical_predictors) > 0) {
  cat("  Model has", length(categorical_predictors), "categorical predictor(s):\n", file = txt_filename, append = TRUE)
  for (var in categorical_predictors) {
    cat("  - ", var, "\n", file = txt_filename, append = TRUE)
  }
} else {
  cat("  Model has NO categorical predictors.\n", file = txt_filename, append = TRUE)
  cat("  Pulkstenis-Robinson tests could not be run.\n", file = txt_filename, append = TRUE)
}
cat("\n", file = txt_filename, append = TRUE)

cat("Lipsitz Tests:\n", file = txt_filename, append = TRUE)
if (!is.null(res_lip1)) {
  cat("  generalhoslem::lipsitz.test:\n", file = txt_filename, append = TRUE)
  cat("    Statistic:", res_lip1$statistic, "\n", file = txt_filename, append = TRUE)
  cat("    df:", res_lip1$parameter, "\n", file = txt_filename, append = TRUE)
  cat("    p-value:", res_lip1$p.value, "\n\n", file = txt_filename, append = TRUE)
} else {
  cat("  generalhoslem::lipsitz.test: Failed to run\n\n", file = txt_filename, append = TRUE)
}

if (!is.null(res_lip2)) {
  cat("  gofcat::lipsitz:\n", file = txt_filename, append = TRUE)
  cat("    LR:", res_lip2$LR, "\n", file = txt_filename, append = TRUE)
  cat("    df:", res_lip2$df, "\n", file = txt_filename, append = TRUE)
  cat("    p-value:", res_lip2$p.value, "\n\n", file = txt_filename, append = TRUE)
} else {
  cat("  gofcat::lipsitz: Failed to run\n\n", file = txt_filename, append = TRUE)
}

if (length(categorical_predictors) > 0) {
  cat("Pulkstenis-Robinson Tests:\n", file = txt_filename, append = TRUE)
  
  if (!is.null(res_pr1)) {
    cat("  generalhoslem::pulkrob.chisq:\n", file = txt_filename, append = TRUE)
    cat("    Statistic:", res_pr1$statistic, "\n", file = txt_filename, append = TRUE)
    cat("    df:", res_pr1$parameter, "\n", file = txt_filename, append = TRUE)
    cat("    p-value:", res_pr1$p.value, "\n\n", file = txt_filename, append = TRUE)
  } else {
    cat("  generalhoslem::pulkrob.chisq: Failed to run\n\n", file = txt_filename, append = TRUE)
  }
  
  if (!is.null(res_pr2)) {
    cat("  generalhoslem::pulkrob.deviance:\n", file = txt_filename, append = TRUE)
    cat("    Statistic:", res_pr2$statistic, "\n", file = txt_filename, append = TRUE)
    cat("    df:", res_pr2$parameter, "\n", file = txt_filename, append = TRUE)
    cat("    p-value:", res_pr2$p.value, "\n\n", file = txt_filename, append = TRUE)
  } else {
    cat("  generalhoslem::pulkrob.deviance: Failed to run\n\n", file = txt_filename, append = TRUE)
  }
  
  if (!is.null(res_pr3)) {
    cat("  gofcat::pulkroben (chisq):\n", file = txt_filename, append = TRUE)
    cat("    Statistic:", res_pr3$stat, "\n", file = txt_filename, append = TRUE)
    cat("    df:", res_pr3$df, "\n", file = txt_filename, append = TRUE)
    cat("    p-value:", res_pr3$p.value, "\n\n", file = txt_filename, append = TRUE)
  } else {
    cat("  gofcat::pulkroben (chisq): Failed to run\n\n", file = txt_filename, append = TRUE)
  }
  
  if (!is.null(res_pr4)) {
    cat("  gofcat::pulkroben (deviance):\n", file = txt_filename, append = TRUE)
    cat("    Statistic:", res_pr4$stat, "\n", file = txt_filename, append = TRUE)
    cat("    df:", res_pr4$df, "\n", file = txt_filename, append = TRUE)
    cat("    p-value:", res_pr4$p.value, "\n\n", file = txt_filename, append = TRUE)
  } else {
    cat("  gofcat::pulkroben (deviance): Failed to run\n\n", file = txt_filename, append = TRUE)
  }
}

cat("Hosmer-Lemeshow Tests:\n", file = txt_filename, append = TRUE)

if (!is.null(res_hlm1)) {
  cat("  generalhoslem::logitgof:\n", file = txt_filename, append = TRUE)
  cat("    Statistic:", as.numeric(res_hlm1$statistic), "\n",
      file = txt_filename, append = TRUE)
  cat("    df:",        as.numeric(res_hlm1$parameter), "\n",
      file = txt_filename, append = TRUE)
  cat("    p-value:",   res_hlm1$p.value, "\n\n",
      file = txt_filename, append = TRUE)
} else {
  cat("  generalhoslem::logitgof: Failed to run\n\n", file = txt_filename, append = TRUE)
}

if (!is.null(res_hlm2)) {
  cat("  gofcat::hosmerlem:\n", file = txt_filename, append = TRUE)
  cat("    Statistic:", res_hlm2$chi.sq, "\n", file = txt_filename, append = TRUE)
  cat("    df:", res_hlm2$df, "\n", file = txt_filename, append = TRUE)
  cat("    p-value:", res_hlm2$p.value, "\n", file = txt_filename, append = TRUE)
} else {
  cat("  gofcat::hosmerlem: Failed to run\n", file = txt_filename, append = TRUE)
}

message("Results saved to ", txt_filename)

}