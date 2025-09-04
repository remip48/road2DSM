backward_selection <- function(variable,
                               seg_data,
                               response,
                               ...complexity...) { # use actually the code from CoastalFutues

  covariates <- variable

  seg_data_init <- seg_data %>%
    dplyr::select(all_of(covariates), effort_km2, ppho, X, Y) %>%
    drop_na() %>%
    as.data.frame()

  seg_data_scale <- seg_data_init

  rescale2 <- function (ynew, y = NULL)
  {
    assert_that(is.numeric(ynew) && is.vector(ynew) && length(ynew) >=
                  3)
    assert_that((is.numeric(y) && is.vector(y) && length(y) >=
                   3) || is.null(y))
    if (is.null(y)) {
      out <- (ynew - mean(ynew, na.rm = TRUE))/(sd(ynew, na.rm = TRUE))
    }
    else {
      out <- (ynew - mean(y, na.rm = TRUE))/(sd(y, na.rm = TRUE))
    }
    return(out)
  }

  seg_data_scale[, covariates] <- apply(seg_data_scale[, covariates], 2, rescale2)

  seg_data_scale <- as.data.frame(seg_data_scale)

  model <- paste0(response, " ~ 1 + ",
                  "te(X, Y, bs = c('cs', 'cs')) + ",
                  paste(paste0("s(", covariates, ", bs = 'cs', k = 10) + "), collapse = ""),
                  "offset(I(log(effort_km2)))")

  model_it <- mgcv::gam(as.formula(model),
                        method = "REML",
                        family = nb(),
                        data = seg_data_scale
  )

  print(gratia::draw(model_it, rug = F))

  csm <- summary(model_it)

  table_var <- data.frame(var = c("XY", covariates),
                          pvalues = abs(csm$s.pv),
                          chi = csm$chi.sq) %>%
    dplyr::filter(var != "XY") %>%
    arrange(pvalues, -chi)

  rho <- cor(seg_data_scale %>%
               dplyr::select(all_of(covariates)))

  diag(rho) <- 0

  all_thresholds <- c(.99, .9, .8, .7, .5)

  print(max(rho))

  ### correlations
  for (thresh_cor in all_thresholds) {
    while (any(rho > thresh_cor)) { #  | any(table_var$pvalues > .05)
      cat("Threshold:", thresh_cor, ". Max correlation:", print(round(max(rho), 2)), ". Max P-value:", max(round(table_var$pvalues, 3)), ".\n")

      rc <- which(map_dbl(1:nrow(rho), function(r) {return(max(abs(rho[r, ]), na.rm = T))}) > tcor)

      to_rem <- table_var %>%
        dplyr::filter(var %in% covariates[rc]) %>%
        dplyr::slice_max(pvalues, n = 1) %>%
        dplyr::slice_min(chi, n = 1)

      remove_v <- c(remove_v, to_rem %>%
                      pull(var))

      cat("__________________________________________________________________\ncorrelated    ---   ",
          to_rem %>%
            pull(var) %>%
            paste(., collapse = " / "), "(P:",
          to_rem  %>%
            pull(pvalues) %>%
            paste(., collapse = " / "), ", cor:", round(rc[match(to_rem %>%
                                                                   pull(var),
                                                                 covariates)], 2),
          ") removed   ---  ",
          length(covariates) - 1, "left.\nRemoved:",
          paste(remove_v, collapse = ", "), "\n")

      to_rem <- to_rem %>%
        pull(var)

      covariates <- covariates[!(covariates %in% to_rem)]

      model <- paste0("ppho ~ 1 + ",
                      "te(X, Y, bs = c('cs', 'cs')) + ",
                      paste(paste0("s(", covariates, ", bs = 'cs', k = 10) + "), collapse = ""),
                      "offset(I(log(effort_km2)))")

      model_it <- mgcv::gam(as.formula(model),
                            method = "REML",
                            family = nb(),
                            data = seg_data_scale
      )

      print(gratia::draw(model_it, rug = F))
      Sys.sleep(2)

      csm <- summary(model_it)

      table_var <- data.frame(var = c("XY", covariates),
                              pvalues = abs(csm$s.pv),
                              chi = csm$chi.sq) %>%
        dplyr::filter(var != "XY") %>%
        arrange(pvalues, -chi)

      rho <- cor(seg_data_scale %>%
                   dplyr::select(all_of(covariates)))

      diag(rho) <- 0

      print(table_var)

      gc()
    }
  }

  ### p values
  while (any(table_var$pvalues > .05)) {
    print(table_var[nrow(table_var), ])
    cat("Removed    ---   ", nrow(table_var) - 1, "left\n")

    table_var <- table_var[-nrow(table_var), ]

    covariates <- table_var$var

    model <- paste0("ppho ~ 1 + ",
                    "te(X, Y, bs = c('cs', 'cs')) + ",
                    paste(paste0("s(", covariates, ", bs = 'cs', k = 10) + "), collapse = ""),
                    "offset(I(log(effort_km2)))")

    model_it <- mgcv::gam(as.formula(model),
                          method = "REML",
                          family = nb(),
                          data = seg_data_scale
    )

    print(gratia::draw(model_it, rug = F))
    Sys.sleep(2)

    csm <- summary(model_it)

    table_var <- data.frame(var = c("XY", covariates),
                            pvalues = abs(csm$s.pv),
                            chi = csm$chi.sq) %>%
      dplyr::filter(var != "XY") %>%
      arrange(pvalues, -chi)

    print(table_var)

    gc()
  }

  # while (nrow(table_var) > 8) {
  #   print(table_var[nrow(table_var), ])
  #   cat("Removed    ---   ", nrow(table_var) - 1, "left\n")
  #
  #   table_var <- table_var[-nrow(table_var), ]
  #
  #   covariates <- table_var$var
  #
  #   model <- paste0("ppho ~ 1 + ",
  #                   paste(paste0("s(", covariates, ", bs = 'cs', k = 10) + "), collapse = ""),
  #                   "te(X, Y, bs = c('cs', 'cs')) + ",
  #                   "offset(I(log(effort_km2)))")
  #
  #   model_it <- mgcv::gam(as.formula(model),
  #                         method = "REML",
  #                         family = nb(),
  #                         data = seg_data_scale
  #   )
  #
  #   print(gratia::draw(model_it, rug = F))
  #   Sys.sleep(2)
  #
  #   csm <- summary(model_it)
  #
  #   table_var <- data.frame(var = c(covariates, "XY"),
  #                           pvalues = abs(csm$s.pv),
  #                           chi = csm$chi.sq) %>%
  #     dplyr::filter(var != "XY") %>%
  #     arrange(pvalues, -chi)
  #
  #   print(table_var)
  #
  #   gc()
  # }

  model_actual <- mgcv::gam(as.formula(model),
                            method = "REML",
                            family = nb(),
                            data = seg_data_init)

  final_file <- list(n_best = 1,

                     all_fits_binded = data.frame(model = model,
                                                  index = 1,
                                                  Convergence = ifelse(model_it$converged, 1, 0),
                                                  AIC = model_it$aic, ResDev = model_it$deviance, NulDev = model_it$null.deviance,
                                                  ExpDev = 100 * round(1 - model_it$deviance/model_it$null.deviance,
                                                                       3),
                                                  RMSE = qpcR::RMSE(model_it),
                                                  looic = 0,
                                                  se_looic = 0,
                                                  stacking_weights = 1),

                     best_models = list(model_it),

                     all_models_tried = "Backward selection",

                     best_models4plotting = list(model_actual),

                     splines_by = NULL,

                     random = NULL)
}
