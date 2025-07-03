run_influence_leverage_tests <- function(model, data, model_name, model_type,
                                         results_dir, plots_dir) {

  if (!inherits(model, "polr")) {
    stop("This function only supports models fitted with MASS::polr")
  }
  
  inf_plots_dir <- file.path(plots_dir, "influence_leverage")
  if (!dir.exists(inf_plots_dir)) dir.create(inf_plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  deeper_results_dir <- file.path(results_dir, "influence_leverage")
  dir.create(deeper_results_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  cb_cols <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  
  output_file <- file.path(deeper_results_dir, paste0(model_name, "_", model_type, "_influence_diagnostics.txt"))
  
  sink(output_file)
  
  cat("===================================================================\n")
  cat("INFLUENCE AND LEVERAGE DIAGNOSTICS FOR", toupper(model_name), "MODEL\n")
  cat("Model Type: POLR\n")
  cat("===================================================================\n\n")
  
  tryCatch({
    n <- nrow(data)
    p <- length(coef(model))
    
    p_all <- length(insight::get_parameters(model))
    
    cat("Sample size (n):", n, "\n")
    cat("Number of parameters (p):",  p_all, "\n\n")
    
    cutoff_cooksd <- 4/n
    cutoff_dfbetas <- 2/sqrt(n)
    cutoff_dffits <- 2*sqrt(( p_all+1)/(n- p_all-1))
    
    cat("Recommended cutoff values:\n")
    cat("  Cook's Distance:", cutoff_cooksd, "\n")
    cat("  DFBETAS:", cutoff_dfbetas, "\n")
    cat("  DFFITS:", cutoff_dffits, "\n\n")
    
    cat("LEVERAGE (HAT VALUES)\n")
    cat("---------------------\n")
    
    tryCatch({
      require(insight)
      
      X <- insight::get_modelmatrix(model)
      
      H <- X %*% solve(crossprod(X)) %*% t(X)
      hat_values <- diag(H)
      
      cat("Summary statistics:\n")
      cat("  Min:", min(hat_values), "\n")
      cat("  1st Qu:", quantile(hat_values, 0.25), "\n")
      cat("  Median:", median(hat_values), "\n")
      cat("  Mean:", mean(hat_values), "\n")
      cat("  3rd Qu:", quantile(hat_values, 0.75), "\n")
      cat("  Max:", max(hat_values), "\n\n")
      
      avg_leverage <- (p_all+1)/n
      high_leverage <- which(hat_values > 2*avg_leverage)
      cat("Average leverage:", avg_leverage, "\n")
      cat("High leverage points (> 2 * avg):", length(high_leverage), "observations\n")
      if (length(high_leverage) > 0 && length(high_leverage) <= 20) {
        cat("Observations with high leverage:\n")
        for (i in high_leverage) {
          cat("  Obs", i, ": Hat =", hat_values[i], "\n")
        }
      } else if (length(high_leverage) > 20) {
        cat("Top 20 observations with highest leverage:\n")
        top_leverage <- order(hat_values, decreasing = TRUE)[1:20]
        for (i in top_leverage) {
          cat("  Obs", i, ": Hat =", hat_values[i], "\n")
        }
      }
      cat("\n")
      
    }, error = function(e) {
      cat("Could not calculate hat values manually.\n")
      cat("Error message:", conditionMessage(e), "\n\n")
      high_leverage <- integer(0)
    })
    
    cat("COOK'S DISTANCE (manual calculation)\n")
    cat("----------------------------------\n")
    
    tryCatch({
      orig_coef <- coef(model)
      orig_vcov <- vcov(model)
      
      cooksd <- numeric(n)
      
      if (n > 100) {
        cat("Calculating Cook's distance for all", n, "observations (this may take a while)...\n")
        progress_step <- ceiling(n/10)
      } else {
        cat("Calculating Cook's distance for all", n, "observations.\n")
        progress_step <- n
      }
      
      for (i in 1:n) {
        if (n > 100 && i %% progress_step == 0) {
          cat("  Progress:", round(i/n*100), "% complete\n")
        }
        
        tryCatch({
          data_subset <- data[-i, ]
          
          model_i <- MASS::polr(formula(model), data = data_subset, Hess = TRUE, method = model$method)
          
          coef_i <- coef(model_i)
          
          common_coefs <- intersect(names(orig_coef), names(coef_i))
          
          delta_coef <- orig_coef[common_coefs] - coef_i[common_coefs]
          
          vcov_subset <- orig_vcov[common_coefs, common_coefs]
          
          cooksd[i] <- t(delta_coef) %*% solve(vcov_subset) %*% delta_coef / length(common_coefs)
          
        }, error = function(e) {
          cooksd[i] <- NA
          cat("  Warning: Failed to calculate Cook's distance for observation", i, "\n")
          cat("  Error message:", conditionMessage(e), "\n")
        })
      }
      
      cooksd_valid <- cooksd[!is.na(cooksd)]
      indices_valid <- which(!is.na(cooksd))
      
      if (length(cooksd_valid) > 0) {
        cat("\nSummary statistics:\n")
        cat("  Min:", min(cooksd_valid), "\n")
        cat("  1st Qu:", quantile(cooksd_valid, 0.25), "\n")
        cat("  Median:", median(cooksd_valid), "\n")
        cat("  Mean:", mean(cooksd_valid), "\n")
        cat("  3rd Qu:", quantile(cooksd_valid, 0.75), "\n")
        cat("  Max:", max(cooksd_valid), "\n\n")
        
        influential_cooksd <- indices_valid[cooksd_valid > cutoff_cooksd]
        cat("Influential points based on Cook's distance (>", cutoff_cooksd, "):", 
            length(influential_cooksd), "observations\n")
        
        if (length(influential_cooksd) > 0 && length(influential_cooksd) <= 20) {
          cat("Observations with high Cook's distance:\n")
          for (i in influential_cooksd) {
            cat("  Obs", i, ": Cook's D =", cooksd[i], "\n")
          }
        } else if (length(influential_cooksd) > 20) {
          cat("Top 20 observations with highest Cook's distance:\n")
          top_cooksd <- order(cooksd, decreasing = TRUE, na.last = TRUE)[1:20]
          for (i in top_cooksd) {
            if (!is.na(cooksd[i])) {
              cat("  Obs", i, ": Cook's D =", cooksd[i], "\n")
            }
          }
        }
      } else {
        cat("Could not calculate Cook's distance for any observation.\n")
      }
      cat("\n")
      
    }, error = function(e) {
      cat("Could not calculate Cook's distance.\n")
      cat("Error message:", conditionMessage(e), "\n\n")
      influential_cooksd <- integer(0)
    })
    
    cat("DFBETAS (leave-one-out approach)\n")
    cat("-----------------------------\n")
    
    tryCatch({
      orig_coef <- coef(model)
      orig_vcov <- vcov(model)
      
      dfbetas <- matrix(NA, nrow = n, ncol = length(orig_coef))
      colnames(dfbetas) <- names(orig_coef)
      max_dfbetas <- numeric(n)
      
      if (exists("cooksd") && length(cooksd) == n) {
        cat("Reusing model fits from Cook's distance calculations for DFBETAS...\n\n")
        
        for (i in 1:n) {
          if (!is.na(cooksd[i])) {
            data_subset <- data[-i, ]
            model_i <- MASS::polr(formula(model), data = data_subset, Hess = TRUE, method = model$method)
            
            coef_i <- coef(model_i)
            
            common_coefs <- intersect(names(orig_coef), names(coef_i))
            delta_coef <- orig_coef[common_coefs] - coef_i[common_coefs]
            
            vcov_subset <- orig_vcov[common_coefs, common_coefs]
            
            dfbetas_i <- delta_coef / sqrt(diag(vcov_subset))
            dfbetas[i, common_coefs] <- dfbetas_i
            max_dfbetas[i] <- max(abs(dfbetas_i))
          }
        }
      } else {
        cat("Calculating DFBETAS for all observations...\n")
        
        if (n > 100) {
          progress_step <- ceiling(n/10)
        } else {
          progress_step <- n
        }
        
        for (i in 1:n) {
          if (n > 100 && i %% progress_step == 0) {
            cat("  Progress:", round(i/n*100), "% complete\n")
          }
          
          tryCatch({
            data_subset <- data[-i, ]
            
            model_i <- MASS::polr(formula(model), data = data_subset, Hess = TRUE, method = model$method)

            coef_i <- coef(model_i)
            
            common_coefs <- intersect(names(orig_coef), names(coef_i))
            delta_coef <- orig_coef[common_coefs] - coef_i[common_coefs]
            
            vcov_subset <- orig_vcov[common_coefs, common_coefs]
            
            dfbetas_i <- delta_coef / sqrt(diag(vcov_subset))
            dfbetas[i, common_coefs] <- dfbetas_i
            max_dfbetas[i] <- max(abs(dfbetas_i))
            
          }, error = function(e) {
            cat("Error refitting model without observation", i, ":", conditionMessage(e), "\n")
          })
        }
      }
      
      max_dfbetas_valid <- max_dfbetas[!is.na(max_dfbetas)]
      dfbetas_indices <- which(!is.na(max_dfbetas))
      
      if (length(max_dfbetas_valid) > 0) {
        cat("\nSummary of maximum absolute DFBETAS values:\n")
        cat("  Min:", min(max_dfbetas_valid), "\n")
        cat("  1st Qu:", quantile(max_dfbetas_valid, 0.25), "\n")
        cat("  Median:", median(max_dfbetas_valid), "\n")
        cat("  Mean:", mean(max_dfbetas_valid), "\n")
        cat("  3rd Qu:", quantile(max_dfbetas_valid, 0.75), "\n")
        cat("  Max:", max(max_dfbetas_valid), "\n\n")
        
        influential_dfbetas <- dfbetas_indices[max_dfbetas[dfbetas_indices] > cutoff_dfbetas]
        
        cat("Influential points based on max |DFBETAS| (>", cutoff_dfbetas, "):", 
            length(influential_dfbetas), "observations\n")
        
        if (length(influential_dfbetas) > 0 && length(influential_dfbetas) <= 20) {
          cat("Observations with high DFBETAS values:\n")
          for (i in influential_dfbetas) {
            cat("  Obs", i, ": Max |DFBETAS| =", max_dfbetas[i], "\n")
            most_affected <- which.max(abs(dfbetas[i,]))
            coef_name <- colnames(dfbetas)[most_affected]
            cat("    Most affected coefficient:", coef_name, "with DFBETA =", 
                dfbetas[i, most_affected], "\n")
          }
        } else if (length(influential_dfbetas) > 20) {
          cat("Top 20 observations with highest max |DFBETAS|:\n")
          top_dfbetas <- order(max_dfbetas, decreasing = TRUE, na.last = TRUE)[1:20]
          for (i in top_dfbetas) {
            if (!is.na(max_dfbetas[i])) {
              cat("  Obs", i, ": Max |DFBETAS| =", max_dfbetas[i], "\n")
              most_affected <- which.max(abs(dfbetas[i,]))
              coef_name <- colnames(dfbetas)[most_affected]
              cat("    Most affected coefficient:", coef_name, "with DFBETA =", 
                  dfbetas[i, most_affected], "\n")
            }
          }
        }
        cat("\n")
      } else {
        cat("Could not calculate DFBETAS using leave-one-out approach.\n\n")
      }
      
    }, error = function(e) {
      cat("Could not calculate DFBETAS.\n")
      cat("Error message:", conditionMessage(e), "\n\n")
      influential_dfbetas <- integer(0)
    })
    
    cat("DFFITS (leave-one-out approach)\n")
    cat("-----------------------------\n")
    
    tryCatch({
      dffits <- numeric(n)
      
      if (exists("hat_values") && exists("dfbetas") && all(dim(dfbetas) == c(n, length(orig_coef)))) {
        cat("Calculating DFFITS using leverage and DFBETAS values...\n\n")
        
        for (i in 1:n) {
          if (!is.na(max_dfbetas[i]) && !is.na(hat_values[i])) {
            dffits[i] <- max_dfbetas[i] * sqrt(hat_values[i] / (1 - hat_values[i]))
          }
        }
        
        dffits_valid <- dffits[!is.na(dffits)]
        dffits_indices <- which(!is.na(dffits))
        
        if (length(dffits_valid) > 0) {
          cat("Summary of DFFITS values:\n")
          cat("  Min:", min(dffits_valid), "\n")
          cat("  1st Qu:", quantile(dffits_valid, 0.25), "\n")
          cat("  Median:", median(dffits_valid), "\n")
          cat("  Mean:", mean(dffits_valid), "\n")
          cat("  3rd Qu:", quantile(dffits_valid, 0.75), "\n")
          cat("  Max:", max(dffits_valid), "\n\n")
          
          influential_dffits <- dffits_indices[dffits[dffits_indices] > cutoff_dffits]
          
          cat("Influential points based on |DFFITS| (>", cutoff_dffits, "):", 
              length(influential_dffits), "observations\n")
          
          if (length(influential_dffits) > 0 && length(influential_dffits) <= 20) {
            cat("Observations with high DFFITS values:\n")
            for (i in influential_dffits) {
              cat("  Obs", i, ": |DFFITS| =", dffits[i], "\n")
            }
          } else if (length(influential_dffits) > 20) {
            cat("Top 20 observations with highest |DFFITS|:\n")
            top_dffits <- order(dffits, decreasing = TRUE, na.last = TRUE)[1:20]
            for (i in top_dffits) {
              if (!is.na(dffits[i])) {
                cat("  Obs", i, ": |DFFITS| =", dffits[i], "\n")
              }
            }
          }
          cat("\n")
        } else {
          cat("Could not calculate DFFITS values.\n\n")
        }
      } else {
        cat("Could not calculate DFFITS because leverage or DFBETAS calculations failed.\n\n")
      }
      
    }, error = function(e) {
      cat("Could not calculate DFFITS.\n")
      cat("Error message:", conditionMessage(e), "\n\n")
      influential_dffits <- integer(0)
    })
    
    cat("COVARIANCE RATIO (COVRATIO)\n")
    cat("-------------------------\n")
    
    tryCatch({
      covratio <- numeric(n)
      
      if (exists("cooksd") && length(cooksd) == n) {
        cat("Calculating COVRATIO using existing model fits...\n\n")
        
        vcov_full <- vcov(model)
        
        for (i in 1:n) {
          if (!is.na(cooksd[i])) {
            data_subset <- data[-i, ]
            
            model_i <- MASS::polr(formula(model), data = data_subset, Hess = TRUE, method = model$method)
            
            vcov_i <- vcov(model_i)
            
            common_coefs <- intersect(names(coef(model)), names(coef(model_i)))
            
            vcov_full_sub <- vcov_full[common_coefs, common_coefs]
            vcov_i_sub <- vcov_i[common_coefs, common_coefs]
            
            covratio[i] <- det(vcov_i_sub) / det(vcov_full_sub)
          }
        }
      } else {
        cat("Calculating COVRATIO for all", n, "observations...\n")
        if (n > 100) {
          progress_step <- ceiling(n/10)
        } else {
          progress_step <- n
        }
        
        vcov_full <- vcov(model)
        
        for (i in 1:n) {
          if (n > 100 && i %% progress_step == 0) {
            cat("  Progress:", round(i/n*100), "% complete\n")
          }
          
          tryCatch({
            data_subset <- data[-i, ]
            
            model_i <- MASS::polr(formula(model), data = data_subset, Hess = TRUE, method = model$method)
            
            vcov_i <- vcov(model_i)
            
            common_coefs <- intersect(names(coef(model)), names(coef(model_i)))
            
            vcov_full_sub <- vcov_full[common_coefs, common_coefs]
            vcov_i_sub <- vcov_i[common_coefs, common_coefs]
            
            covratio[i] <- det(vcov_i_sub) / det(vcov_full_sub)
            
          }, error = function(e) {
            covratio[i] <- NA
            cat("  Warning: Failed to calculate COVRATIO for observation", i, "\n")
            cat("  Error message:", conditionMessage(e), "\n")
          })
        }
      }
      
      covratio_valid <- covratio[!is.na(covratio)]
      covratio_indices <- which(!is.na(covratio))
      
      if (length(covratio_valid) > 0) {
        cat("\nSummary statistics for COVRATIO:\n")
        cat("  Min:", min(covratio_valid), "\n")
        cat("  1st Qu:", quantile(covratio_valid, 0.25), "\n")
        cat("  Median:", median(covratio_valid), "\n")
        cat("  Mean:", mean(covratio_valid), "\n")
        cat("  3rd Qu:", quantile(covratio_valid, 0.75), "\n")
        cat("  Max:", max(covratio_valid), "\n\n")
        
        cutoff_covratio_min <- 1 - 3 * p_all / n
        cutoff_covratio_max <- 1 + 3 * p_all / n
        
        cat("Recommended cutoff values for COVRATIO:\n")
        cat("  Lower bound:", cutoff_covratio_min, "\n")
        cat("  Upper bound:", cutoff_covratio_max, "\n\n")
        
        influential_covratio_low <- covratio_indices[covratio[covratio_indices] < cutoff_covratio_min]
        influential_covratio_high <- covratio_indices[covratio[covratio_indices] > cutoff_covratio_max]
        influential_covratio <- c(influential_covratio_low, influential_covratio_high)
        
        cat("Influential points based on COVRATIO (outside [", 
            cutoff_covratio_min, ",", cutoff_covratio_max, "]):", 
            length(influential_covratio), "observations\n")
        
        if (length(influential_covratio) > 0 && length(influential_covratio) <= 20) {
          cat("Observations with extreme COVRATIO values:\n")
          for (i in influential_covratio) {
            cat("  Obs", i, ": COVRATIO =", covratio[i], 
                ifelse(covratio[i] < cutoff_covratio_min, " (precision decreasing)", " (precision increasing)"), "\n")
          }
        } else if (length(influential_covratio) > 20) {
          cat("Top 20 observations with most extreme COVRATIO values:\n")
          sorted_by_extreme <- order(abs(covratio - 1), decreasing = TRUE, na.last = TRUE)[1:20]
          for (i in sorted_by_extreme) {
            if (!is.na(covratio[i])) {
              cat("  Obs", i, ": COVRATIO =", covratio[i], 
                  ifelse(covratio[i] < cutoff_covratio_min, " (precision decreasing)", 
                         ifelse(covratio[i] > cutoff_covratio_max, " (precision increasing)", "")), "\n")
            }
          }
        }
        cat("\n")
      } else {
        cat("Could not calculate COVRATIO for any observation.\n\n")
      }
      
    }, error = function(e) {
      cat("Could not calculate COVRATIO.\n")
      cat("Error message:", conditionMessage(e), "\n\n")
      influential_covratio <- integer(0)
    })

    cat("=============================================\n")
    cat("SUMMARY OF INFLUENCE DIAGNOSTICS ANALYSIS\n")
    cat("=============================================\n\n")
    
    influential_obs <- unique(c(
      if(exists("high_leverage")) high_leverage else integer(0),
      if(exists("influential_cooksd")) influential_cooksd else integer(0),
      if(exists("influential_dfbetas")) influential_dfbetas else integer(0),
      if(exists("influential_dffits")) influential_dffits else integer(0),
      if(exists("influential_covratio")) influential_covratio else integer(0)
    ))
    
    if (length(influential_obs) > 0) {
      cat("Total influential observations detected:", length(influential_obs), "\n")
      if (length(influential_obs) <= 20) {
        cat("List of all influential observations:", paste(sort(influential_obs), collapse = ", "), "\n\n")
      } else {
        cat("Large number of influential observations detected (", length(influential_obs), " in total).\n\n")
      }
      
      cat("Diagnostic summary for each influential observation:\n")
      cat("--------------------------------------------------\n")
      cat(sprintf("%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                  "Obs", "High Leverage", "High Cook's D", "High DFBETAS", "High DFFITS", 
                  "Extreme COVRATIO"))
      
      for (i in sort(influential_obs)) {
        is_high_leverage <- if(exists("high_leverage")) i %in% high_leverage else FALSE
        is_high_cooksd <- if(exists("influential_cooksd")) i %in% influential_cooksd else FALSE
        is_high_dfbetas <- if(exists("influential_dfbetas")) i %in% influential_dfbetas else FALSE
        is_high_dffits <- if(exists("influential_dffits")) i %in% influential_dffits else FALSE
        is_extreme_covratio <- if(exists("influential_covratio")) i %in% influential_covratio else FALSE
        
        cat(sprintf("%-10d %-15s %-15s %-15s %-15s %-15s\n", 
                    i, 
                    ifelse(is_high_leverage, "YES", "NO"),
                    ifelse(is_high_cooksd, "YES", "NO"),
                    ifelse(is_high_dfbetas, "YES", "NO"),
                    ifelse(is_high_dffits, "YES", "NO"),
                    ifelse(is_extreme_covratio, "YES", "NO")))
      }
      cat("\n")
    } else {
      cat("No influential observations detected by any method.\n\n")
    }
    
  }, error = function(e) {
    cat("ERROR: Failed to calculate influence diagnostics\n")
    cat("Error message:", conditionMessage(e), "\n")
  })
  
  sink()
  
  cat("Influence and leverage diagnostics saved to:", output_file, "\n")
  
  plot_data <- data.frame(
    index = 1:n,
    leverage = hat_values,
    cooks_d = cooksd,
    max_dfbetas = max_dfbetas,
    dffits = dffits,
    covratio = covratio
  )
  
  for (j in 1:ncol(dfbetas)) {
    plot_data[[paste0("dfbetas_", colnames(dfbetas)[j])]] <- dfbetas[,j]
  }
  
  plot_data$combined_influence <- scale(abs(plot_data$leverage)) + 
    scale(abs(plot_data$cooks_d)) + 
    scale(abs(plot_data$max_dfbetas)) + 
    scale(abs(plot_data$dffits)) + 
    scale(abs(plot_data$covratio - 1))
  
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = leverage, y = cooks_d)) +
    ggplot2::geom_point(color = cb_cols[3], alpha = 0.7) +
    ggplot2::geom_hline(yintercept = cutoff_cooksd, linetype = "dashed", color = cb_cols[7]) +
    ggplot2::geom_vline(xintercept = 2*(p_all+1)/n, linetype = "dashed", color = cb_cols[7]) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(cooks_d > cutoff_cooksd | leverage > 2*(p_all+1)/n, 
                                                   as.character(index), "")), 
                       vjust = -0.5, hjust = -0.2, size = 3) +
    ggplot2::labs(title = "Leverage vs. Cook's Distance",
                  subtitle = paste("Model:", model_name, "- Type:", model_type),
                  x = "Leverage (hat values)",
                  y = "Cook's Distance",
                  caption = paste("Horizontal line at Cook's D =", round(cutoff_cooksd, 4),
                                  "\nVertical line at Leverage =", round(2*(p_all+1)/n, 4))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
  
  ggplot2::ggsave(file.path(inf_plots_dir, paste0(model_name, "_", model_type, "_leverage_cooks.png")), 
                  p1, width = 8, height = 6, bg = "white")
  
  dfbetas_data <- reshape2::melt(plot_data, id.vars = "index", 
                                 measure.vars = grep("dfbetas_", names(plot_data), value = TRUE),
                                 variable.name = "coefficient", value.name = "dfbetas")
  
  dfbetas_data$coefficient <- gsub("dfbetas_", "", dfbetas_data$coefficient)
  
  p2 <- ggplot2::ggplot(dfbetas_data, ggplot2::aes(x = index, y = dfbetas, color = coefficient)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = c(-cutoff_dfbetas, cutoff_dfbetas), 
                        linetype = "dashed", color = cb_cols[7]) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(abs(dfbetas) > cutoff_dfbetas, 
                                                   as.character(index), "")), 
                       vjust = -0.5, hjust = -0.2, size = 3) +
    ggplot2::scale_color_manual(values = cb_cols) +
    ggplot2::labs(title = "DFBETAS for Each Coefficient",
                  subtitle = paste("Model:", model_name, "- Type:", model_type),
                  x = "Observation Index",
                  y = "DFBETAS",
                  caption = paste("Horizontal lines at ±", round(cutoff_dfbetas, 4))) +
    ggplot2::theme_minimal() +
    ggplot2::facet_wrap(~coefficient, scales = "free_y") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
  
  ggplot2::ggsave(file.path(inf_plots_dir, paste0(model_name, "_", model_type, "_dfbetas.png")), 
                  p2, width = 10, height = 8, bg = "white")
  
  p3 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = index, y = dffits)) +
    ggplot2::geom_point(color = cb_cols[4], alpha = 0.7) +
    ggplot2::geom_hline(yintercept = c(-cutoff_dffits, cutoff_dffits), 
                        linetype = "dashed", color = cb_cols[7]) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(abs(dffits) > cutoff_dffits, 
                                                   as.character(index), "")), 
                       vjust = -0.5, hjust = -0.2, size = 3) +
    ggplot2::labs(title = "DFFITS vs. Observation Index",
                  subtitle = paste("Model:", model_name, "- Type:", model_type),
                  x = "Observation Index",
                  y = "DFFITS",
                  caption = paste("Horizontal lines at ±", round(cutoff_dffits, 4))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
  
  ggplot2::ggsave(file.path(inf_plots_dir, paste0(model_name, "_", model_type, "_dffits.png")), 
                  p3, width = 8, height = 6, bg = "white")
  
  p4 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = index, y = covratio)) +
    ggplot2::geom_point(color = cb_cols[5], alpha = 0.7) +
    ggplot2::geom_hline(yintercept = c(cutoff_covratio_min, cutoff_covratio_max), 
                        linetype = "dashed", color = cb_cols[7]) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted", color = cb_cols[1]) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(covratio < cutoff_covratio_min | 
                                                     covratio > cutoff_covratio_max, 
                                                   as.character(index), "")), 
                       vjust = -0.5, hjust = -0.2, size = 3) +
    ggplot2::labs(title = "COVRATIO vs. Observation Index",
                  subtitle = paste("Model:", model_name, "- Type:", model_type),
                  x = "Observation Index",
                  y = "COVRATIO",
                  caption = paste("Horizontal lines at", round(cutoff_covratio_min, 4), 
                                  "and", round(cutoff_covratio_max, 4),
                                  "\nDotted line at COVRATIO = 1")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
  
  ggplot2::ggsave(file.path(inf_plots_dir, paste0(model_name, "_", model_type, "_covratio.png")), 
                  p4, width = 8, height = 6, bg = "white")
  
  p5 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = index, y = combined_influence)) +
    ggplot2::geom_point(color = cb_cols[6], alpha = 0.7) +
    ggplot2::geom_hline(yintercept = mean(plot_data$combined_influence, na.rm = TRUE) + 
                          2 * sd(plot_data$combined_influence, na.rm = TRUE), 
                        linetype = "dashed", color = cb_cols[7]) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(combined_influence > 
                                                     mean(combined_influence, na.rm = TRUE) + 
                                                     2 * sd(combined_influence, na.rm = TRUE), 
                                                   as.character(index), "")), 
                       vjust = -0.5, hjust = -0.2, size = 3) +
    ggplot2::labs(title = "Combined Influence Index",
                  subtitle = paste("Model:", model_name, "- Type:", model_type),
                  x = "Observation Index",
                  y = "Standardized Combined Influence",
                  caption = "Standardized sum of all influence measures\nDashed line at mean + 2*SD") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
  
  ggplot2::ggsave(file.path(inf_plots_dir, paste0(model_name, "_", model_type, "_combined_influence.png")), 
                  p5, width = 8, height = 6, bg = "white")
  
  
  return(output_file)
}