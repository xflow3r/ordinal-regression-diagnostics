library(simstudy)
library(MASS)

default_beta1 <- 1.5
default_beta2 <- -1.0
default_beta3 <- 0.8
default_baseprobs <- c(0.15, 0.20, 0.30, 0.25, 0.10)

generate_ideal_data <- function(n = 500, beta1 = default_beta1, beta2 = default_beta2, 
                                baseprobs = default_baseprobs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  def <- defData(varname = "x1", formula = 0, variance = 1) 
  def <- defData(def, varname = "x2", dist = "binary", formula = 0.5)
  def <- defData(def, varname = "linpred", formula = paste(beta1, "* x1 +", beta2, "* x2"), variance = 0)
  data <- genData(n, def)
  data <- genOrdCat(data, adjVar = "linpred", baseprobs = baseprobs, 
                    catVar = "y", asFactor = TRUE)
  data$y <- factor(data$y, levels = 1:length(baseprobs), 
                   labels = as.character(1:length(baseprobs)), ordered = TRUE)
  data$linpred <- NULL
  return(data)
}


generate_multicollinearity_data <- function(n = 500, beta1 = default_beta1, beta2 = default_beta2, 
                                            beta3 = 0.6, beta4 = -0.4,
                                            correlation = 0.5, baseprobs = default_baseprobs, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  def <- defData(varname = "temp", formula = 0, variance = 1)
  data <- genData(n, def)
  
  cor_matrix <- matrix(c(
    1.0,         correlation, 0.0, 0.0,    
    correlation, 1.0,         0.0, 0.0,    
    0.0,         0.0,         1.0, 0.3,    
    0.0,         0.0,         0.3, 1.0     
  ), nrow = 4, byrow = TRUE)
  
  data <- addCorData(data, idname = "id", mu = c(0, 0, 0, 0), sigma = c(1, 1, 1, 1), 
                     corMatrix = cor_matrix, cnames = "x1,x2,x3,x4")
  
  data$linpred <- beta1 * data$x1 + beta2 * data$x2 + beta3 * data$x3 + beta4 * data$x4
  
  data <- genOrdCat(data, adjVar = "linpred", baseprobs = baseprobs, 
                    catVar = "y", asFactor = TRUE)
  
  data$y <- factor(data$y, levels = 1:length(baseprobs), 
                   labels = as.character(1:length(baseprobs)), ordered = TRUE)
  
  data$temp <- NULL
  data$linpred <- NULL
  return(data)
}

generate_PO_violation_data <- function(n = 500, beta1 = default_beta1, beta2 = default_beta2, 
                                       baseprobs = default_baseprobs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  def <- defData(varname = "x1", formula = 0, variance = 1)
  def <- defData(def, varname = "x2", dist = "binary", formula = 0.5)
  data <- genData(n, def)
  data$linpred <- beta1 * data$x1 + beta2 * data$x2

  np_adjust <- matrix(c(2, 1, 0, -1, -2), nrow = 1)

  data <- genOrdCat(data, adjVar = "linpred", baseprobs = baseprobs, catVar = "y", 
                    asFactor = TRUE, npVar = "x2", npAdj = np_adjust)
  data$y <- factor(data$y, levels = 1:length(baseprobs), 
                   labels = as.character(1:length(baseprobs)), ordered = TRUE)
  data$linpred <- NULL
  return(data)
}

generate_outlier_data <- function(n = 500, beta1 = default_beta1, beta2 = default_beta2, 
                                  outlier_pct = 0.05, outlier_shift = 5, 
                                  baseprobs = default_baseprobs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  def <- defData(varname = "x1", formula = 0, variance = 1)
  def <- defData(def, varname = "x2", dist = "binary", formula = 0.5)
  def <- defData(def, varname = "linpred", formula = paste(beta1, "* x1 +", beta2, "* x2"), variance = 0)
  data <- genData(n, def)
  
  data <- genOrdCat(data, adjVar = "linpred", baseprobs = baseprobs, 
                    catVar = "y", asFactor = TRUE)
  
  data$y <- factor(data$y, levels = 1:length(baseprobs), 
                   labels = as.character(1:length(baseprobs)), ordered = TRUE)
  
  n_outliers <- round(n * outlier_pct)
  if (n_outliers > 0) {
    outlier_indices <- sample(1:n, n_outliers)
    for (i in outlier_indices) {

      data$x1[i] <- data$x1[i] + rnorm(1, mean = outlier_shift, sd = 0.5)
      
      if (data$y[i] == as.character(1)) {
        data$y[i] <- as.character(length(baseprobs))
      } else if (data$y[i] == as.character(length(baseprobs))) {
        data$y[i] <- "1"
      } else {
        data$y[i] <- sample(c("1", as.character(length(baseprobs))), 1)
      }
    }
    data$y <- factor(data$y, levels = as.character(1:length(baseprobs)), ordered = TRUE)
  }
  
  data$linpred <- NULL
  return(data)
}

generate_nonlinear_data <- function(n = 500, beta1 = default_beta1, beta2 = default_beta2, 
                                    nonlinear_coeff = 0.7, baseprobs = default_baseprobs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  def <- defData(varname = "x1", formula = 0, variance = 1)
  def <- defData(def, varname = "x2", dist = "binary", formula = 0.5)
  data <- genData(n, def)
  
  data$linpred <- beta1 * data$x1 + beta2 * data$x2 + nonlinear_coeff * (data$x1^2)
  
  data <- genOrdCat(data, adjVar = "linpred", baseprobs = baseprobs, 
                    catVar = "y", asFactor = TRUE)
  
  data$y <- factor(data$y, levels = 1:length(baseprobs), 
                   labels = as.character(1:length(baseprobs)), ordered = TRUE)
  
  data$linpred <- NULL
  return(data)
}

generate_omitted_variable_data <- function(n = 500, beta1 = default_beta1, beta2 = default_beta2, beta3 = default_beta3, 
                                           x3_correlation = 0.5, baseprobs = default_baseprobs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  def <- defData(varname = "x2", dist = "binary", formula = 0.5)
  data <- genData(n, def)
  
  data <- addCorData(data, idname = "id", mu = c(0, 0), sigma = c(1, 1), 
                     rho = x3_correlation, corstr = "cs", cnames = "x1,x3")
  
  data$linpred <- beta1 * data$x1 + beta2 * data$x2 + beta3 * data$x3
  
  data <- genOrdCat(data, adjVar = "linpred", baseprobs = baseprobs, 
                    catVar = "y", asFactor = TRUE)
  
  data$y <- factor(data$y, levels = 1:length(baseprobs), 
                   labels = as.character(1:length(baseprobs)), ordered = TRUE)
  
  data$x3 <- NULL
  data$linpred <- NULL
  return(data)
}

generate_heteroscedastic_data <- function(n = 500, beta1 = default_beta1, beta2 = default_beta2, 
                                          base_scale = 1.0, scale_factor = 5.0, 
                                          baseprobs = default_baseprobs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  def <- defData(varname = "x1", formula = 0, variance = 1)
  def <- defData(def, varname = "x2", dist = "binary", formula = 0.5)
  data <- genData(n, def)
  
  scale_param <- base_scale * exp(scale_factor * abs(data$x1))
  
  linear_pred <- beta1 * data$x1 + beta2 * data$x2
  
  latent_y <- numeric(n)
  for(i in 1:n) {
    latent_y[i] <- linear_pred[i] + rlogis(1, location = 0, scale = scale_param[i])
  }
  
  cum_probs <- cumsum(baseprobs)
  thresholds <- qlogis(cum_probs[-length(cum_probs)])
  
  data$y <- cut(latent_y, 
                breaks = c(-Inf, thresholds, Inf), 
                labels = 1:length(baseprobs),
                ordered = TRUE)
  
  return(data)
}


generate_link_misspecification_data <- function(n = 500, beta1 = default_beta1, beta2 = default_beta2, 
                                                baseprobs = default_baseprobs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  def <- defData(varname = "x1", formula = 0, variance = 1)
  def <- defData(def, varname = "x2", dist = "binary", formula = 0.5)
  data <- genData(n, def)
  
  eta <- beta1 * data$x1 + beta2 * data$x2
  
  cum_probs <- cumsum(baseprobs)[-length(baseprobs)]
  thresholds <- qnorm(cum_probs)
  
  Z <- eta + rnorm(n, 0, 1)
  
  y_num <- findInterval(Z, vec = c(-Inf, thresholds, Inf), rightmost.closed = TRUE)
  
  data$y <- factor(y_num, levels = 1:length(baseprobs), 
                   labels = as.character(1:length(baseprobs)), ordered = TRUE)
  
  return(data)
}


generate_scenario_data <- function(scenario, n = 500, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  scenario <- toupper(scenario)
  if (scenario == "S1") {
    return(generate_ideal_data(n = n, seed = seed))
  } else if (scenario == "S2") {
    return(generate_multicollinearity_data(n = n, correlation = 0.5, seed = seed))
  } else if (scenario == "S3") {
    return(generate_multicollinearity_data(n = n, correlation = 0.9, seed = seed))
  } else if (scenario == "S4") {
    return(generate_PO_violation_data(n = n, seed = seed))
  } else if (scenario == "S5") {
    return(generate_outlier_data(n = n, seed = seed))
  } else if (scenario == "S6") {
    return(generate_nonlinear_data(n = n, seed = seed))
  } else if (scenario == "S7") {
    return(generate_omitted_variable_data(n = n, seed = seed))
  } else if (scenario == "S8") {
    return(generate_heteroscedastic_data(n = n, seed = seed))
  } else if (scenario == "S9") {
    return(generate_link_misspecification_data(n = n, seed = seed))
  } else {
    stop("Invalid scenario identifier. Use 'S1' through 'S9'.")
  }
}


generate_all_datasets <- function(scenarios = paste0("S", 1:9), n = 500, seed = 123) {
  datasets <- list()
  for (sc in scenarios) {
    cat("Generating data for scenario", sc, "...\n")
    datasets[[sc]] <- generate_scenario_data(sc, n = n, seed = seed)
  }
  return(datasets)
}
