run_summary_nimble <- function(run_models,
                               response,
                               calibdata) {

  for (i in 1:length(run_models$best_models)) {

    cat("\n\n")
    cat(paste0("## Model ", i, "\n\n"))

    print(run_models$best_models[[i]]$WAIC)

    samples <- do.call(rbind, mcmc_glmm$samples)

    var <- colnames(samples)[stringr::str_detect(colnames(samples), "beta_X")]
    var <- do.call("c", lapply(var, function(v) {
      out <- stringr::str_remove_all(v, "beta_X")
      out <- stringr::str_split_1(out, fixed("["))[1]
    })) %>%
      unique() %>%
      as.numeric() %>%
      sort()

    terms <- do.call("c", lapply(1:length(var), function(v) {
      run_models$sm_X[[i]][[v]]$term
    }))

    splines <- purrr::map_dfr(var, function(x) {
      beta_cols <- grep(paste0("^beta_X", x, "\\["), colnames(samples))
      beta_samples <- samples[, beta_cols]

      X <- run_models$data_constant[[i]][[paste0("X_X", x)]]

      smooth_draws <- beta_samples %*% t(X)

      smooth_mean <- apply(smooth_draws, 2, mean)

      smooth_low <- apply(smooth_draws, 2, quantile, 0.025)

      smooth_high <- apply(smooth_draws, 2, quantile, 0.975)

      data.frame(variable = terms[match(x, var)],
                 mean = smooth_mean,
                 low = smooth_low,
                 high = smooth_high,
                 x = calibdata %>%
                   dplyr::pull(get(terms[match(x, var)])))
    })

    print(ggplot() +
            geom_ribbon(data = splines, aes(x = x, ymin = low, ymax = high), fill = "midnightblue", alpha = .5) +
            geom_point(data = splines, aes(x = x, y = mean)) +
            facet_wrap(~ variable, scales = "free_x"))

    cat("\n\n")
    cat(paste0("#### ASPE & Ratio of the number of observed ", ifelse(str_detect(response, "group"), "groups", "individuals"),
               " / number of predicted ", ifelse(str_detect(response, "group"), "groups", "individuals"),
               "\n"))

    samples_mat <- as.matrix(samples)

    intercept_s <- samples_mat[, "intercept"]

    beta_X   <- map(var, function(v) {
      K <- ncol(run_models$data_constant[[i]][[paste0("X_X", v)]])

      get_cols(samples_mat, paste0("beta_X", v), K)
    })

    n_iter <- nrow(samples_mat)
    n_sim  <- 1000   # match the GAM side's simulation budget / keep memory sane
    if (n_iter > n_sim) {
      keep <- sort(sample(seq_len(n_iter), n_sim))
      intercept_s <- intercept_s[keep]

      for (v in 1:length(beta_X)) {
        beta_X[[v]] <- beta_X[[v]][keep, , drop = FALSE]
      }
      n_iter <- n_sim
    }

    eta <- matrix(intercept_s, nrow = nrow(calibdata), ncol = n_iter, byrow = TRUE)

    nch <- run_models$best_models[[i]]$samples
    Rhat_var <- vector("list", length(beta_X))

    for (v in 1:length(beta_X)) { # for every covariate
      pred_mat <- mgcv::PredictMat(run_models$sm_X[[i]][[v]],
                                   data.frame(X1 = calibdata %>%
                                                dplyr::pull(get(terms[v]))) %>%
                                     dplyr::rename(!!terms[v] := X1))

      pred_mat <- sweep(pred_mat, 2, colMeans(run_models$sm_X[[i]][[v]]$X), "-")

      eta <- eta +  pred_mat %*% t(beta_X[[v]])

      ## estimate Rhat
      # smooth_chains <- vector("list", nch)
      #
      # for (ch in 1:nch) {
      #   beta <- run_models$best_models[[i]]$samples[[ch]][, grep("^", paste0("beta_X",
      #                                                                        var[v]),
      #                                                            "\\[", colnames(run_models$best_models[[i]]$samples[[ch]]))]
      #
      #   smooth_chains[[ch]] <- pred_mat %*% t(beta)
      # }
      #
      # Rhat_var[[v]] <- map_dfr(seq_len(nrow(pred_mat)), function(ix) {
      #   chains <- mcmc.list(lapply(smooth_chains, function(x) mcmc(x[ix, ])))
      #
      #   data.frame(ix = ix,
      #              Rhat = gelman.diag(chains)$psrf[1],
      #              ESS = effectiveSize(chains))
      # }) %>%
      #   dplyr::mutate(vi = v,
      #                 variable = run_models_noESWg0$sm_X[[i]][[v]]$term,
      #                 x = calibdata %>% dplyr::pull(run_models$sm_X[[i]][[v]]$term))
    }

    # print(ggplot() +
    #         geom_line(data = do.call("rbind", Rhat_var), aes(x = x, y = Rhat)) +
    #         geom_hline(yintercept = 1.00) +
    #         geom_hline(yintercept = 1.01) +
    #         geom_hline(yintercept = 1.05) +
    #         facet_grid(model ~ variable, scales = "free"))

    dens_pred  <- exp(eta)
    abund_pred <- dens_pred * calibdata$effort_km2
    p <- rowMeans(abund_pred, na.rm = TRUE)

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
