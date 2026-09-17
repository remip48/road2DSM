run_summary_gam <- function(run_models,
                            response,
                            calibdata,
                            log1p_trans,
                            best_stacking) {


  ### Modify and improve spline plot
  build_spline_plot <- function(model, model4plotting, which_soap, log1p_trans) {
    # Extract gratia plot to modify it
    plot <- gratia::draw(model4plotting, rug = T, select = -which_soap)

    # Extract edfs for each covariates
    gam_summary <- summary(model)
    edf_values <- gam_summary$s.table[, "edf"]

    # Modify each subplots
    for (n in 1:length(plot)){
      current_covariate <- gsub("s\\((.*)\\)", "\\1", plot[[n]]$labels$title)

      # Transform aspect and slope to degrees
      if (current_covariate == "aspect") {
        aspect_breaks <- seq(0, 2 * pi, by = pi / 2)
        aspect_labels <- aspect_breaks * (180 / pi) - 180

        plot[[n]] <- plot[[n]] +
          scale_x_continuous(name = "aspect (°)", breaks = aspect_breaks, labels = aspect_labels)
      }
      if (current_covariate == "slope") {
        slope_range <- range(model4plotting$model$slope, na.rm = TRUE)
        slope_degree_breaks <- seq(0, 30, by = 5)
        slope_breaks <- tan(slope_degree_breaks * (pi / 180))
        slope_labels <- slope_degree_breaks

        plot[[n]] <- plot[[n]] +
          scale_x_continuous(name = "slope (°)", breaks = slope_breaks, labels = slope_labels)
      }

      plot[[n]] <- plot[[n]] +
        ylab(sprintf("s(%s, %.2f)", current_covariate, edf_values[n+1])) +
        # labs(caption = NULL) +
        labs(title = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
        )

      # Log-scale for log-transformed covariates
      if(current_covariate %in% log1p_trans){
        plot[[n]] <- plot[[n]] +
          scale_x_continuous(trans = "log1p")
      }

      # Add units for each covariate
      if(current_covariate == "EKE"){
        plot[[n]]$labels$x = paste0(plot[[n]]$labels$x, " (m/s)")
      }
      if(current_covariate == "NPPV"){
        plot[[n]]$labels$x = bquote(.(plot[[n]]$labels$x) ~ "(mg/m"^3*"/day)")
      }
      # if(current_covariate %in% c("slope","aspect")){
      #   plot[[n]]$labels$x = paste0(plot[[n]]$labels$x, " (°)")
      # }
      if(current_covariate %in% c("mSST","gradSST")){
        plot[[n]]$labels$x = paste0(plot[[n]]$labels$x, " (°C)")
      }
      if(current_covariate %in% c("dist_coast","dist_50m","dist_200m","dist_2000m")){
        plot[[n]]$labels$x = paste0(plot[[n]]$labels$x, " (m)")
      }
      if(current_covariate == "bathy"){
        plot[[n]]$labels$x = "depth (m)"
      }
    }

    out <- plot
  }


  for (i in 1:length(run_models$best_models)) {

    cat("\n\n")
    cat(paste0("## Model ", i, "\n\n"))

    print(summary(run_models$best_models[[i]]))

    which_soap <- which(vapply(run_models[["best_models"]][[i]][["smooth"]], inherits, logical(1), what = "soap.film"))

    if (length(which_soap) >= 1) {
      for (j in which_soap) {
        (plot(run_models$best_models[[i]], select = j))
      }

      spline_plot <- build_spline_plot(run_models$best_models[[i]], run_models$best_models4plotting[[i]], which_soap, log1p_trans)
      print(spline_plot)

      if (i == best_stacking & save_plots) {
        if (!dir.exists(paste0(prediction_folder, "/Plots/", output_file))) {
          dir.create(paste0(prediction_folder, "/Plots/", output_file))
        }

        ggsave(paste0(prediction_folder, "/Plots/", output_file, "/Spline Plot.png"), plot = spline_plot,
               width = 9, height = 5, dpi = 300)
      }

    } else {
      if (grepl("t2\\(X,Y,", paste(stringr::str_remove_all(as.character(run_models[["best_models"]][[1]][["formula"]]),
                                                           " "),
                                   collapse = ""))) {
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = F, select = 1)))
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = F, select = 2)))
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = F, select = -c(1,2))))
      } else if (grepl("X,Y", paste(stringr::str_remove_all(as.character(run_models[["best_models"]][[1]][["formula"]]),
                                                            " "),
                                    collapse = ""))) {
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = F, select = 1)))
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = F, select = -1)))
      } else {
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = F)))
      }
    }

    mgcv::qq.gam(run_models$best_models[[i]], rep = 1000)

    cat("\n\n")
    cat(paste0("#### ASPE & Ratio of the number of observed ", ifelse(str_detect(response, "group"), "groups", "individuals"),
               " / number of predicted ", ifelse(str_detect(response, "group"), "groups", "individuals"),
               "\n"))

    p <- predict(run_models$best_models[[i]], newdata=calibdata, type='response')
    dens <-as.numeric(p)
    obs_n <- calibdata %>%
      pull(response)
    ASPE <- (sum((obs_n - dens)^2, na.rm=TRUE) / nrow(calibdata))
    cat("ASPE =", ASPE, "\n")

    ratio <- obs_n/dens

    print(summary(ratio))

    print(ggplot2::ggplot() +
            ggplot2::geom_histogram(data = data.frame(Ratio = ratio) %>%
                                      dplyr::filter(!is.na(Ratio)), ggplot2::aes(x = Ratio)) +
            ggplot2::scale_y_sqrt(name = "Count"))

    checks <- as.data.frame(calibdata)
    checks$new <- checks[, response]
    checks$Value <- "Observed value"
    checks <- checks %>%
      rbind(checks %>%
              mutate(new = unname(as.numeric(p)),
                     Value = "Predicted value")) %>%
      mutate(Value = as.factor(Value))

    color1 <- viridis::viridis(256)[1]
    color2 <- viridis::viridis(256)[256]

    rgb1 <- col2rgb(color1) / 255
    rgb2 <- col2rgb(color2) / 255
    blended_rgb <- (rgb1 + rgb2) / 2
    blended_color <- rgb(blended_rgb[1], blended_rgb[2], blended_rgb[3])

    rootg <- ggplot() +
      geom_histogram(data = checks %>%
                       dplyr::filter(as.character(Value) == "Observed value"), aes(x = new, fill = Value), alpha = 1, binwidth = 1) +
      scale_fill_manual(values = color1) +
      new_scale_fill() +
      geom_histogram(data = checks %>%
                       dplyr::filter(as.character(Value) == "Predicted value"), aes(x = new, fill = Value), alpha = .5, binwidth = 1) +
      scale_fill_manual(values = color2) +
      new_scale_fill() +
      geom_tile(aes(x = Inf, y = Inf, fill = "Observed + Predicted value")) +
      scale_fill_manual(name = "",
                        values = c("Observed + Predicted value" = blended_color)) +
      scale_y_sqrt() +
      ylab("Count") +
      xlab("Number of individuals") +
      theme(legend.title = element_blank())

    if (max(checks$new, na.rm = T) > 20) {
      rootg <- rootg +
        scale_x_sqrt()
    }

    print(rootg)

    if (i == best_stacking & save_plots) {
      ggsave(paste0(prediction_folder, "/Plots/", output_file, "/Rootogram Plot.png"), plot = rootg,
             width = 8, height = 5, dpi = 300)
    }
  }


  if (save_plots) {
    i <- best_stacking
    which_soap <- which(vapply(run_models$best_models[[i]]$smooth, inherits, logical(1), what = "soap.film"))

    for (j in which_soap) {
      if (!dir.exists(paste0(prediction_folder, "/Plots/", output_file))) {
        dir.create(paste0(prediction_folder, "/Plots/", output_file))
      }

      png(paste0(prediction_folder, "/Plots/", output_file, "/Soap Plot.png"), width = 2000, height = 2000, res = 300)
      plot(run_models$best_models[[i]], select = j)
      dev.off()
    }

    png(filename = paste0(prediction_folder, "/Plots/", output_file, "/QQ Plot.png"),
        width = 1500, height = 1500, res = 300)
    mgcv::qq.gam(run_models$best_models[[i]], rep = 1000)
    dev.off()
  }

  invisible()
}
