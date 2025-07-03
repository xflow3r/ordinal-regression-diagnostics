run_nominal_test <- function(model, data, model_name, model_type, results_dir, plots_dir) {
 
  nominal_results_dir <- file.path(results_dir, "nominal_tests")
  
  if (!dir.exists(nominal_results_dir)) dir.create(nominal_results_dir, showWarnings = FALSE, recursive = TRUE)
  
  output_file <- file.path(nominal_results_dir, paste0("nominal_test_", model_name, ".txt"))
   
  sink(output_file)
  
  cat("\n\n======================================\n")
  cat("Nominal Test for", toupper(model_name), "Model\n")
  cat("======================================\n\n")
  
  print(nominal_test(model))
  
  sink()
}