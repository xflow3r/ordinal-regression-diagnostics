
required_packages <- c(
  "MASS", "sure", "DHARMa", "tidyverse", "ggplot2", "gridExtra", 
  "simstudy", "ordinal", "brant", "magick", "VGAM", "rms", 
  "generalhoslem", "gofcat", "Matrix", "car", "performance", 
  "insight", "stats", "faraway", "broom", "data.table"
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

options(DHARMaSignalColor = "#E69F00")

setwd("C:/Users/leoko/Documents/TU/bachelorarbeit/simulationStudy")
project_dir <- getwd()

results_dir <- file.path(project_dir, "results")
plots_dir <- file.path(project_dir, "plots")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

start_time <- Sys.time()
cat("Starting simulation at:", format(start_time), "\n")
cat("Results will be saved to:", results_dir, "\n")
cat("Plots will be saved to:", plots_dir, "\n")

source("generate_data.R")

n <- 300
seed <- 123

data_correct <- generate_ideal_data(n = 500, seed = 123)
data_nonprop <- generate_PO_violation_data(n = 500, seed = 123)
data_outliers <- generate_outlier_data(n = 500, seed = 123)
data_wrong_link <- generate_link_misspecification_data(n = 500, seed = 123)

data_multicollinearity_weak <- generate_multicollinearity_data(n = 500, correlation = 0.5, seed = 123)
data_multicollinearity_strong <- generate_multicollinearity_data(n = 500, correlation = 0.9, seed = 123)

data_S5_outliers           <- generate_outlier_data(n = 500, seed = 123)
data_S6_nonlinear          <- generate_nonlinear_data(n = 500, seed = 123)
data_S7_omitted_variable   <- generate_omitted_variable_data(n = 500, seed = 123)
data_S8_heteroscedastic    <- generate_heteroscedastic_data(n = 500, seed = 123)
data_S9_link_misspec       <- generate_link_misspecification_data(n = 500, seed = 123)

end_time <- Sys.time()
cat("Finished data generation at:", format(end_time), "\n")

# Model Fitting Section
#-------------------------------------------------------------------------

model_correct <- polr(y ~ x1 + x2, data = data_correct, method = "logistic", Hess = TRUE)
model_wrong_link <- polr(y ~ x1 + x2, data = data_wrong_link, method = "cloglog", Hess = TRUE)
model_nonprop <- polr(y ~ x1 + x2, data = data_nonprop, method = "logistic", Hess = TRUE)
model_outliers <- polr(y ~ x1 + x2, data = data_outliers, method = "logistic", Hess = TRUE)

model_multicollinearity_weak <- polr(y ~ x1 + x2 + x3 + x4, data = data_multicollinearity_weak, method = "logistic", Hess = TRUE)
model_multicollinearity_strong <- polr(y ~ x1 + x2 + x3 + x4, data = data_multicollinearity_strong, method = "logistic", Hess = TRUE)

model_nonlinear <- polr(y ~ x1 + x2, data = data_S6_nonlinear, method = "logistic", Hess = TRUE)
model_omitted_variable <- polr(y ~ x1 + x2, data = data_S7_omitted_variable, method = "logistic", Hess = TRUE)
model_heteroscedastic <- polr(y ~ x1 + x2, data = data_S8_heteroscedastic, method = "logistic", Hess = TRUE)

clm_correct <- clm(y ~ x1 + x2, data = data_correct, link = "logit")
clm_wrong_link <- clm(y ~ x1 + x2, data = data_wrong_link, link = "cloglog")
clm_nonprop <- clm(y ~ x1 + x2, data = data_nonprop, link = "logit")
clm_outliers <- clm(y ~ x1 + x2, data = data_outliers, link = "logit")

clm_multicollinearity_weak <- clm(y ~ x1 + x2 + x3 + x4, data = data_multicollinearity_weak, link = "logit")
clm_multicollinearity_strong <- clm(y ~ x1 + x2 + x3 + x4, data = data_multicollinearity_strong, link = "logit")

clm_nonlinear <- clm(y ~ x1 + x2, data = data_S6_nonlinear, link = "logit")
clm_omitted_variable <- clm(y ~ x1 + x2, data = data_S7_omitted_variable, link = "logit")
clm_heteroscedastic <- clm(y ~ x1 + x2, data = data_S8_heteroscedastic, link = "logit")


vgam_correct <- vglm(y ~ x1 + x2, family = cumulative(link = "logit", parallel = TRUE, reverse  = FALSE), data = data_correct)
vgam_wrong_link <- vglm(y ~ x1 + x2, family = cumulative(link = "cloglog", parallel = TRUE, reverse  = FALSE), data = data_wrong_link)
vgam_nonprop <- vglm(y ~ x1 + x2, family = cumulative(link = "logit", parallel = TRUE, reverse  = FALSE), data = data_nonprop)
vgam_outliers <- vglm(y ~ x1 + x2, family = cumulative(link = "logit", parallel = TRUE, reverse  = FALSE), data = data_outliers)

vgam_nonlinear <- vglm(y ~ x1 + x2, family = cumulative(link = "logit", parallel = TRUE, reverse = FALSE), data = data_S6_nonlinear)
vgam_omitted_variable <- vglm(y ~ x1 + x2, family = cumulative(link = "logit", parallel = TRUE, reverse = FALSE), data = data_S7_omitted_variable)
vgam_heteroscedastic <- vglm(y ~ x1 + x2, family = cumulative(link = "logit", parallel = TRUE, reverse = FALSE), data = data_S8_heteroscedastic)

rms_correct <- orm(y ~ x1 + x2, data = data_correct, family = "logistic")
rms_wrong_link <- orm(y ~ x1 + x2, data = data_wrong_link, family = "cloglog")
rms_nonprop <- orm(y ~ x1 + x2, data = data_nonprop, family = "logistic")
rms_outliers <- orm(y ~ x1 + x2, data = data_outliers, family = "logistic")


rms_nonlinear <- orm(y ~ x1 + x2, data = data_S6_nonlinear, family = "logistic")
rms_omitted_variable <- orm(y ~ x1 + x2, data = data_S7_omitted_variable, family = "logistic")
rms_heteroscedastic <- orm(y ~ x1 + x2, data = data_S8_heteroscedastic, family = "logistic")



# ---- Collinearity Checks ----
cat("Starting Collinearity Tests at:", format(Sys.time()), "\n")
start_time <- Sys.time()


source("collinearity_checks.R")
run_col_tests(
  model = model_correct,
  data = data_correct,
  model_name = "correct",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_col_tests(
  model = model_wrong_link,
  data = data_wrong_link,
  model_name = "wrong_link",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_col_tests(
  model = model_nonprop,
  data = data_nonprop,
  model_name = "nonprop",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_col_tests(
  model = model_outliers,
  data = data_outliers,
  model_name = "outliers",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_col_tests(
  model = model_multicollinearity_weak,
  data = data_multicollinearity_weak,
  model_name = "multicollinearity_weak",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_col_tests(
  model = model_multicollinearity_strong,
  data = data_multicollinearity_strong,
  model_name = "multicollinearity_strong",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_col_tests(
  model = model_nonlinear,
  data = data_S6_nonlinear,
  model_name = "nonlinear",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_col_tests(
  model = model_omitted_variable,
  data = data_S7_omitted_variable,
  model_name = "omitted_variable",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_col_tests(
  model = model_heteroscedastic,
  data = data_S8_heteroscedastic,
  model_name = "heteroscedastic",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

end_time <- Sys.time()
elapsed <- end_time - start_time
cat("Completed Collinearity Tests at:", format(end_time), "\n")
cat("Time taken:", format(elapsed), "\n")


# Influence and leverage Diagnostics - only for polr
cat("Starting Influence and leverage Tests at:", format(Sys.time()), "\n")
start_time <- Sys.time()

source("influence_leverage_tests.R")
run_influence_leverage_tests(
  model = model_correct,
  data = data_correct,
  model_name = "correct",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_influence_leverage_tests(
  model = model_outliers,
  data = data_outliers,
  model_name = "outliers",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_influence_leverage_tests(
  model = model_wrong_link,
  data = data_wrong_link,
  model_name = "wrong_link",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_influence_leverage_tests(
  model = model_nonprop,
  data = data_nonprop,
  model_name = "nonprop",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_influence_leverage_tests(
  model = model_nonlinear,
  data = data_S6_nonlinear,
  model_name = "nonlinear",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_influence_leverage_tests(
  model = model_omitted_variable,
  data = data_S7_omitted_variable,
  model_name = "omitted_variable",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_influence_leverage_tests(
  model = model_heteroscedastic,
  data = data_S8_heteroscedastic,
  model_name = "heteroscedastic",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)


end_time <- Sys.time()
elapsed <- end_time - start_time
cat("Completed Influence and leverage Tests at:", format(end_time), "\n")
cat("Time taken:", format(elapsed), "\n")

cat("Starting Surrogate Residuals tests and plots at:", format(Sys.time()), "\n")
start_time <- Sys.time()



source("plot_surrogate_residuals.R")
results_dir <- "results"
plots_dir <- "plots"


# --- Generate Surrogate Residuals + Plots ---

residuals_clm_correct <- generate_residuals_and_plots(clm_correct, data_correct, "correct", "clm", results_dir, plots_dir)
residuals_clm_wrong_link <- generate_residuals_and_plots(clm_wrong_link, data_wrong_link, "wrong_link", "clm", results_dir, plots_dir)
residuals_clm_nonprop <- generate_residuals_and_plots(clm_nonprop, data_nonprop, "nonprop", "clm", results_dir, plots_dir)
residuals_clm_outliers <- generate_residuals_and_plots(clm_outliers, data_outliers, "outliers", "clm", results_dir, plots_dir)

residuals_clm_multicollinearity_weak <- generate_residuals_and_plots(clm_multicollinearity_weak, data_multicollinearity_weak, "multicollinearity_weak", "clm", results_dir, plots_dir)
residuals_clm_multicollinearity_strong <- generate_residuals_and_plots(clm_multicollinearity_strong, data_multicollinearity_strong, "multicollinearity_strong", "clm", results_dir, plots_dir)

residuals_clm_nonlinear <- generate_residuals_and_plots(clm_nonlinear, data_S6_nonlinear, "nonlinear", "clm", results_dir, plots_dir)

residuals_clm_omitted_variable <- generate_residuals_and_plots(clm_omitted_variable, data_S7_omitted_variable, "omitted_variable", "clm", results_dir, plots_dir)

residuals_clm_heteroscedastic <- generate_residuals_and_plots(clm_heteroscedastic, data_S8_heteroscedastic, "heteroscedastic", "clm", results_dir, plots_dir)

all_residuals <- rbind(
  data.frame(model_name = "correct", model_type = "clm", residuals = residuals_clm_correct, 
             x1 = data_correct$x1, x2 = data_correct$x2, y = data_correct$y),
  data.frame(model_name = "wrong_link", model_type = "clm", residuals = residuals_clm_wrong_link, 
             x1 = data_wrong_link$x1, x2 = data_wrong_link$x2, y = data_wrong_link$y),
  data.frame(model_name = "nonprop", model_type = "clm", residuals = residuals_clm_nonprop, 
             x1 = data_nonprop$x1, x2 = data_nonprop$x2, y = data_nonprop$y),
  data.frame(model_name = "outliers", model_type = "clm", residuals = residuals_clm_outliers, 
             x1 = data_outliers$x1, x2 = data_outliers$x2, y = data_outliers$y),
  data.frame(model_name = "multicollinearity_weak", model_type = "clm", residuals = residuals_clm_multicollinearity_weak, 
             x1 = data_multicollinearity_weak$x1, x2 = data_multicollinearity_weak$x2, y = data_multicollinearity_weak$y),
  data.frame(model_name = "multicollinearity_strong", model_type = "clm", residuals = residuals_clm_multicollinearity_strong, 
             x1 = data_multicollinearity_strong$x1, x2 = data_multicollinearity_strong$x2, y = data_multicollinearity_strong$y),
  data.frame(model_name = "nonlinear", model_type = "clm", residuals = residuals_clm_nonlinear, 
             x1 = data_S6_nonlinear$x1, x2 = data_S6_nonlinear$x2, y = data_S6_nonlinear$y),
  data.frame(model_name = "omitted_variable", model_type = "clm", residuals = residuals_clm_omitted_variable, 
             x1 = data_S7_omitted_variable$x1, x2 = data_S7_omitted_variable$x2, y = data_S7_omitted_variable$y),
  data.frame(model_name = "heteroscedastic", model_type = "clm", residuals = residuals_clm_heteroscedastic, 
             x1 = data_S8_heteroscedastic$x1, x2 = data_S8_heteroscedastic$x2, y = data_S8_heteroscedastic$y)
)

write.csv(all_residuals, file.path(results_dir, "all_residuals.csv"), row.names = FALSE)

end_time <- Sys.time()
elapsed <- end_time - start_time
cat("Completed Surrogate Residuals tests and plots at:", format(end_time), "\n")
cat("Time taken:", format(elapsed), "\n")

# --- DHARMA
cat("Starting Dharma Tests and plots at:", format(Sys.time()), "\n")
start_time <- Sys.time()
source("plot_dharma.R")

dharma_residuals_correct_polr <- generate_dharma_residuals_and_plots(
  model = model_correct,
  data = data_correct,
  model_name = "correct",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

dharma_residuals_wrong_link_polr <- generate_dharma_residuals_and_plots(
  model = model_wrong_link,
  data = data_wrong_link,
  model_name = "wrong_link",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

dharma_residuals_nonprop_polr <- generate_dharma_residuals_and_plots(
  model = model_nonprop,
  data = data_nonprop,
  model_name = "nonprop",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

dharma_residuals_outliers_polr <- generate_dharma_residuals_and_plots(
  model = model_outliers,
  data = data_outliers,
  model_name = "outliers",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

dharma_residuals_nonlinear_clm <- generate_dharma_residuals_and_plots(
  model = model_nonlinear,
  data = data_S6_nonlinear,
  model_name = "nonlinear",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

dharma_residuals_omitted_variable_clm <- generate_dharma_residuals_and_plots(
  model = model_omitted_variable,
  data = data_S7_omitted_variable,
  model_name = "omitted_variable",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

dharma_residuals_heteroscedastic_clm <- generate_dharma_residuals_and_plots(
  model = model_heteroscedastic,
  data = data_S8_heteroscedastic,
  model_name = "heteroscedastic",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

end_time <- Sys.time()
elapsed <- end_time - start_time
cat("Completed Dharma plots and tests at:", format(end_time), "\n")
cat("Time taken:", format(elapsed), "\n")


# ---- Brant  Test ----
cat("Starting Brant Tests at:", format(Sys.time()), "\n")
start_time <- Sys.time()

source("brant_tests.R")

run_brant_test(
  model = model_correct,
  data = data_correct,
  model_name = "correct",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_brant_test(
  model = model_nonprop,
  data = data_nonprop,
  model_name = "nonprop",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_brant_test(
  model = model_multicollinearity_weak,
  data = data_multicollinearity_weak,
  model_name = "multicollinearity_weak",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_brant_test(
  model = model_multicollinearity_strong,
  data = data_multicollinearity_strong,
  model_name = "multicollinearity_strong",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_brant_test(
  model = model_nonlinear,
  data = data_S6_nonlinear,
  model_name = "nonlinear",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_brant_test(
  model = model_omitted_variable,
  data = data_S7_omitted_variable,
  model_name = "omitted_variable",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_brant_test(
  model = model_heteroscedastic,
  data = data_S8_heteroscedastic,
  model_name = "heteroscedastic",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

end_time <- Sys.time()
elapsed <- end_time - start_time
cat("Completed Brant at:", format(end_time), "\n")
cat("Time taken:", format(elapsed), "\n")


# ---- NOMINAL
cat("Starting Nominal tests at:", format(Sys.time()), "\n")
start_time <- Sys.time()

source("nominal_tests.R")

run_nominal_test(
  model = clm_correct,
  data = data_correct,
  model_name = "correct",
  model_type = "clm",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_nominal_test(
  model = clm_multicollinearity_weak,
  data = data_multicollinearity_weak,
  model_name = "multicollinearity_weak",
  model_type = "clm",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_nominal_test(
  model = clm_multicollinearity_strong,
  data = data_multicollinearity_strong,
  model_name = "multicollinearity_strong",
  model_type = "clm",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_nominal_test(
  model = clm_nonlinear,
  data = data_S6_nonlinear,
  model_name = "nonlinear",
  model_type = "clm",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_nominal_test(
  model = clm_omitted_variable,
  data = data_S7_omitted_variable,
  model_name = "omitted_variable",
  model_type = "clm",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_nominal_test(
  model = clm_heteroscedastic,
  data = data_S8_heteroscedastic,
  model_name = "heteroscedastic",
  model_type = "clm",
  results_dir = results_dir,
  plots_dir = plots_dir
)

end_time <- Sys.time()
elapsed <- end_time - start_time
cat("Completed Nominal tests at:", format(end_time), "\n")
cat("Time taken:", format(elapsed), "\n")

# ---- Performing General GOF tests ----
cat("Starting general GOF tests at:", format(Sys.time()), "\n")
start_time <- Sys.time()

source("general_gof_tests.R")

run_gof_tests(
  model = model_correct,
  data = data_correct,
  model_name = "correct",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_gof_tests(
  model = model_multicollinearity_weak,
  data = data_multicollinearity_weak,
  model_name = "multicollinearity_weak",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_gof_tests(
  model = model_multicollinearity_strong,
  data = data_multicollinearity_strong,
  model_name = "multicollinearity_strong",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_gof_tests(
  model = model_nonlinear,
  data = data_S6_nonlinear,
  model_name = "nonlinear",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_gof_tests(
  model = model_omitted_variable,
  data = data_S7_omitted_variable,
  model_name = "omitted_variable",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

run_gof_tests(
  model = model_heteroscedastic,
  data = data_S8_heteroscedastic,
  model_name = "heteroscedastic",
  model_type = "polr",
  results_dir = results_dir,
  plots_dir = plots_dir
)

end_time <- Sys.time()
elapsed <- end_time - start_time
cat("Completed general GOF tests at:", format(end_time), "\n")
cat("Time taken:", format(elapsed), "\n")
