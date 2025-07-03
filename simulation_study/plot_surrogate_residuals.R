generate_residuals_and_plots <- function(model, data, model_name, model_type, results_dir, plots_dir) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    install.packages("gridExtra")
  }
  
  set.seed(123)
  
  library(ggplot2)
  library(gridExtra)
  library(sure)
  library(VGAM)
  library(ordinal)
  library(magick)
  
  cb_cols <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  theme_set(
    theme_classic(base_size = 14) +
      theme(axis.text  = element_text(size = 12),
            axis.title = element_text(face = "bold"),
            plot.title = element_text(face = "bold", hjust = .5),
            legend.position = "bottom")
  )
  
  nsim_value <- 2
  residuals_model <- resids(model, method = "latent", nsim = nsim_value)
  
  residuals_df <- data.frame(
    model_name = model_name,
    model_type = model_type,
    residuals = residuals_model,
    x1 = data$x1,
    x2 = data$x2,
    y = data$y
  )
  
  deeper_results_dir <- file.path(results_dir, "surrogate_residuals")
  dir.create(deeper_results_dir, showWarnings = FALSE, recursive = TRUE)
  
  csv_filename <- paste0("residuals_", model_type, "_", model_name, ".csv")
  write.csv(residuals_df, file.path(deeper_results_dir, csv_filename), row.names = FALSE)
  
  qq_plot <- autoplot.resid(residuals_model, what = "qq",
                            distribution = qnorm,
                            qqpoint.color = cb_cols[3],
                            qqline.color  = cb_cols[2],
                            alpha = 0.8) +
    labs(title = sprintf("%s %s – QQ plot", model_type, model_name))
  
  fitted_plot <- autoplot.resid(residuals_model, what = "fitted",
                                fit = model,
                                color = cb_cols[3],
                                smooth = TRUE,
                                smooth.color = cb_cols[2],
                                alpha = 0.6,
                                size = 1.3) +
    labs(title = sprintf("%s %s – fitted", model_type, model_name),
         y = "Surrogate residual")
  
  x1_plot <- autoplot.resid(
    residuals_model,
    what         = "covariate",
    x            = data$x1,
    color        = cb_cols[3],      
    alpha        = 0.7,             
    size         = 1.3,             
    smooth       = TRUE,
    smooth.color = cb_cols[2],      
    smooth.size  = 1.1              
  ) +
    labs(x = "x1", y = "Surrogate residual") +
    theme_classic(base_size = 14)
  
  x2_plot <- autoplot.resid(
    residuals_model,
    what          = "covariate",
    x             = data$x2,
    fill          = cb_cols[5],   
    color         = cb_cols[3],    
    smooth.color  = cb_cols[1],    
    smooth.size   = 1.0,          
    jitter.width  = 0.2,
    alpha         = 0.6,
    size          = 2
  ) +
    scale_x_continuous(breaks = c(0, 1),
                       labels = c("0", "1")) +
    labs(x = "x2 (0 / 1)", y = "Surrogate residual") +
    theme_classic(base_size = 14)
  
  plot_grid <- grid.arrange(
    qq_plot + ggtitle(paste(model_type, model_name, "- QQ Plot")),
    fitted_plot + ggtitle(paste(model_type, model_name, "- Fitted Plot")),
    x1_plot + ggtitle(paste(model_type, model_name, "- x1 Plot")),
    x2_plot + ggtitle(paste(model_type, model_name, "- x2 Plot")),
    ncol = 2
  )
  
  plot_filename <- paste0("plots_", model_type, "_", model_name, ".png")
  
  deeper_plots_dir <- file.path(plots_dir, "surrogate_residual_plots")
  dir.create(deeper_results_dir, showWarnings = FALSE, recursive = TRUE)
  
  ggsave(file.path(deeper_plots_dir, plot_filename), plot_grid, width = 12, height = 10)
  
  source("sure_gof_tests.R")
  gof_results <- run_all_gof_tests(
    model = model,
    model_name = model_name,
    model_type = model_type,
    results_dir = results_dir,
    plots_dir = plots_dir,
    residuals_model = residuals_model,
    nsim = 100  
  )
  
  return(residuals_model)
}