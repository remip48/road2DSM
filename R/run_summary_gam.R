run_summary_gam <- function(run_models,
                            response,
                            calibdata) {

  for (i in 1:length(run_models$best_models)) {

    cat("\n\n")
    cat(paste0("## Model ", i, "\n\n"))

    print(summary(run_models$best_models[[i]]))

    which_soap <- which(vapply(run_models[["best_models"]][[i]][["smooth"]], inherits, logical(1), what = "soap.film"))

    if (length(which_soap) >= 1) {
      for (j in which_soap) {
        (plot(run_models$best_models4plotting[[i]], select = j))
      }
      try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = T, select = -which_soap)))
    } else {
      if (grepl("t2\\(X,Y,", paste(stringr::str_remove_all(as.character(run_models[["best_models"]][[1]][["formula"]]),
                                                           " "),
                                   collapse = ""))) {
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = T, select = 1)))
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = T, select = 2)))
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = T, select = -c(1,2))))
      } else if (grepl("X,Y", paste(stringr::str_remove_all(as.character(run_models[["best_models"]][[1]][["formula"]]),
                                                            " "),
                                    collapse = ""))) {
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = T, select = 1)))
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = T, select = -1)))
      } else {
        try(print(gratia::draw(run_models$best_models4plotting[[i]], rug = T)))
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

    rootg <- ggplot2::ggplot() +
      ggplot2::geom_histogram(data = checks %>%
                                dplyr::filter(as.character(Value) == "Observed value"), ggplot2::aes(x = new, fill = Value), alpha = 1, binwidth = 1) +
      ggplot2::scale_fill_manual(values = viridis::viridis(256)[1]) +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_histogram(data = checks %>%
                                dplyr::filter(as.character(Value) == "Predicted value"), ggplot2::aes(x = new, fill = Value), alpha = .5, binwidth = 1) +
      ggplot2::scale_fill_manual(values = viridis::viridis(256)[256]) +
      ggplot2::theme(legend.title = ggplot2::element_blank())

    if (max(checks$new, na.rm = T) > 20) {
      rootg <- rootg +
        scale_x_sqrt()
    }

    print(rootg)

  }

  invisible()
}
