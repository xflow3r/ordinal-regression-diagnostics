generate_dharma_residuals_and_plots <- function(model, data, model_name, model_type, results_dir, plots_dir) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    install.packages("gridExtra")
  }
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    install.packages("DHARMa")
  }
  
  set.seed(123)
  
  library(ggplot2)
  library(gridExtra)
  library(DHARMa)
  stopifnot(packageVersion("DHARMa") >= "0.4.5")
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
  
  theme_set(
    theme_classic(base_size = 14) +
      theme(
        axis.text  = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = .5),
        legend.position = "bottom")
  )
  
  options(DHARMaSignalColor = "orange")

  if (model_type == "polr") {
   
    nsim <- 250  
    
    sim_list <- simulate(model, nsim = nsim)      
    sim_mat <- sapply(sim_list, as.integer)       
    sim_mat <- matrix(sim_mat, nrow = nrow(sim_mat)) 
    
    obs <- as.integer(data$y)                     
    
    fit_mostlikely <- apply(predict(model, type = "probs"),
                            1, which.max)
    
    dharma_res <- createDHARMa(
      simulatedResponse = sim_mat,
      observedResponse = obs,
      fittedPredictedResponse = fit_mostlikely,
      integerResponse = TRUE,
      seed = 123
    )
  }
  else if (model_type == "clm") {
    # didnt work
  }
  else if (model_type == "vgam" || model_type == "vglm") {
    nsim <- 250
    
    preds <- predictvglm(model, type = "response")
    
    sim_list <- list()
    for (i in 1:nsim) {
      sim_list[[i]] <- factor(apply(preds, 1, function(p) sample(1:length(p), 1, prob = p)),
                              levels = 1:length(levels(data$y)),
                              labels = levels(data$y),
                              ordered = TRUE)
    }
    
    sim_mat <- sapply(sim_list, as.integer)
    sim_mat <- matrix(sim_mat, nrow = nrow(sim_mat))
    
    obs <- as.integer(data$y)
    
    fit_mostlikely <- apply(preds, 1, which.max)
    
    dharma_res <- createDHARMa(
      simulatedResponse = sim_mat,
      observedResponse = obs,
      fittedPredictedResponse = fit_mostlikely,
      integerResponse = TRUE,
      seed = 123
    )
  }
  else if (model_type == "rms") {
    nsim <- 250
    
    preds <- predict(model, type = "fitted.ind")
    
    sim_list <- list()
    for (i in 1:nsim) {
      sim_list[[i]] <- factor(apply(preds, 1, function(p) sample(1:length(p), 1, prob = p)),
                              levels = 1:length(levels(data$y)),
                              labels = levels(data$y),
                              ordered = TRUE)
    }
    
    sim_mat <- sapply(sim_list, as.integer)
    sim_mat <- matrix(sim_mat, nrow = nrow(sim_mat))
    
    obs <- as.integer(data$y)
    
    fit_mostlikely <- apply(preds, 1, which.max)
    
    dharma_res <- createDHARMa(
      simulatedResponse = sim_mat,
      observedResponse = obs,
      fittedPredictedResponse = fit_mostlikely,
      integerResponse = TRUE,
      seed = 123
    )
  }
  else {
    stop("Unsupported model type: ", model_type)
  }
  
  residuals_model <- dharma_res$scaledResiduals
  
  residuals_df <- data.frame(
    model_name = model_name,
    model_type = paste0("dharma_", model_type),
    residuals = residuals_model,
    x1 = data$x1,
    x2 = data$x2,
    y = data$y
  )
  
  deeper_results_dir <- file.path(results_dir, "dharma_residuals")
  if (!dir.exists(deeper_results_dir)) dir.create(deeper_results_dir, showWarnings = FALSE, recursive = TRUE)
  
  csv_filename <- paste0("dharma_residuals_", model_type, "_", model_name, ".csv")
  write.csv(residuals_df, file.path(deeper_results_dir, csv_filename), row.names = FALSE)
  
  
  deeper_plots_dir <- file.path(plots_dir, "dharma_residual_plots")
  if (!dir.exists(deeper_plots_dir)) dir.create(deeper_plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  png(file.path(deeper_plots_dir,
                sprintf("dharma_standard_%s_%s.png", model_type, model_name)),
      width = 10, height = 5, units = "in", res = 300)
  plot(
    dharma_res,
    quantreg   = TRUE,
    main       = sprintf("%s %s", model_type, model_name),
    col        = cb_cols["sky"],    
    colSmooth  = cb_cols["orange"], 
    pch        = 16,                
    cex        = 0.7,               
    lwd        = 2                  
  )
  dev.off()
  
  png(file.path(deeper_plots_dir,
                sprintf("dharma_qq_%s_%s.png", model_type, model_name)),
      width = 5, height = 5, units = "in", res = 300)
  
  plotQQunif(
    dharma_res,
    testUniformity = TRUE,
    testOutliers   = TRUE,
    testDispersion = TRUE,
    main           = sprintf("%s %s – QQ plot", model_type, model_name),
    col            = cb_cols["sky"],
    pch            = 16,
    cex            = 0.8,
    colLine        = cb_cols["orange"],
    lwdLine        = 2
  )
  dev.off()
  
  png(file.path(deeper_plots_dir, paste0("dharma_fitted_", model_type, "_", model_name, ".png")), 
      width = 5, height = 5, units = "in", res = 300)
  DHARMa::plotResiduals(
    dharma_res,
    form        = fit_mostlikely,                 
    quantreg    = TRUE,
    main        = "Residuals vs fitted values",
    col         = cb_cols["sky"],                 
    pch         = 16,
    cex         = 0.7,
    colSmooth   = cb_cols["orange"],              
    lwdSmooth   = 2,
    ltySmooth   = 1
  )
  dev.off()
  
  png(file.path(deeper_plots_dir, paste0("dharma_x1_", model_type, "_", model_name, ".png")), 
      width = 5, height = 5, units = "in", res = 300)
  DHARMa::plotResiduals(
    dharma_res,
    form        = data$x1,                 
    quantreg    = TRUE,
    main        = "Residuals vs x1",
    col         = cb_cols["sky"],                 
    pch         = 16,
    cex         = 0.7,
    colSmooth   = cb_cols["orange"],              
    lwdSmooth   = 2,
    ltySmooth   = 1
  )
  dev.off()
  
  png(file.path(deeper_plots_dir, paste0("dharma_x2_", model_type, "_", model_name, ".png")), 
      width = 5, height = 5, units = "in", res = 300)
  
  DHARMa::plotResiduals(
    dharma_res,
    form        = data$x2,                    
    quantreg    = TRUE,
    main        = sprintf("%s %s – Residuals vs x2", 
                          model_type, model_name),
    xlab        = "x2",                        
    ylab        = "DHARMa scaled residual",    
    col         = cb_cols["sky"],              
    pch         = 16,
    cex         = 0.7,
    colSmooth   = cb_cols["orange"],           
    lwdSmooth   = 2,
    ltySmooth   = 1
  )
  dev.off()
  
  png(file.path(deeper_plots_dir,
                sprintf("dharma_hist_%s_%s.png", model_type, model_name)),
      width = 5, height = 5, units = "in", res = 300)
  

  hist(
    dharma_res,
    breaks = seq(-0.02, 1.02, len = 53),          
    col    = c(cb_cols["orange"],                
               rep("lightgrey", 50),             
               cb_cols["orange"]),               
    border = cb_cols["black"],
    main   = sprintf("%s %s – residuals histogram", model_type, model_name),
    xlab   = "DHARMa scaled residual"
  )
  dev.off()
  
  
  return(residuals_model)
}