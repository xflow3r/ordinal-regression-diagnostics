run_brant_test <- function(model, data, model_name, model_type, results_dir, plots_dir) {
  
  brant_results_dir <- file.path(results_dir, "brant_tests")
  
  if (!dir.exists(brant_results_dir)) dir.create(brant_results_dir, showWarnings = FALSE, recursive = TRUE)
  
  output_file <- file.path(brant_results_dir, paste0("brant_test_", model_name, ".txt"))
  
  header <- c(
    "\n\n======================================",
    paste("Brant Test for", toupper(model_name), "Model"),
    "======================================\n\n"
  )
  
  standard_output <- capture.output({
    standard_result <- brant(model)
  })
  
  byvar_header <- c("\n\nBy Variable:")
  byvar_output <- capture.output({
    byvar_result <- brant(model, by.var = TRUE)
  })
  
  writeLines(
    c(header, standard_output, byvar_header, byvar_output),
    output_file
  )
  
  invisible(list(
    standard = standard_result,
    by_variable = byvar_result
  ))
}