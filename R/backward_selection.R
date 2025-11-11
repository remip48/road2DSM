#' Title
#'
#' @param variable
#' @param seg_data
#' @param response
#' @param spline_to_add
#' @param soap
#' @param offset_effort
#' @param treshold_discrete
#' @param max_correlation
#' @param Pvalue_max
#'
#' @return
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @export
#'
#' @examples
backward_selection <- function(variable,
                               smoother = function(variable,
                                                   bs,
                                                   complexity) {
                                 return(paste0("s(", variable, ", bs = '", bs, "', k = ", complexity, ")"))
                               }, # input can be a modified function to modify the smoothers
                               bs = "cs",
                               complexity = 10,
                               seg_data,
                               response,
                               spline_to_add = NULL, # can be static or anything else that wont be part of the selection and in all models
                               soap = list(bnd = NULL, knots = NULL, coordinates = c("X", "Y"), by = NULL),
                               offset_effort = "effort_km2",
                               treshold_discrete = 25,
                               max_correlation = .5,
                               nb_max_pred = Inf,
                               Pvalue_max = .05) { # use actually the code from CoastalFutues

  covariates <- variable

  bnd <- soap$bnd
  knots <- soap$knots

  assertthat::assert_that(is.factor(soap$by))

  if (!is.null(spline_to_add)) {
    cat(spline_to_add, "will be added to the model, be values wont be scaled contrary to the predictors used in 'variable'.\n")
  }

  seg_data <- seg_data %>%
    st_drop_geometry() %>%
    as.data.frame()

  seg_data_init <- seg_data %>%
    dplyr::select(all_of(c(covariates, offset_effort, response)), label) %>%
    tidyr::drop_na() %>%
    left_join(seg_data %>%
                dplyr::select(-all_of(c(covariates, offset_effort, response))),
              by = "label") %>%
    as.data.frame()

  seg_data_scale <- seg_data_init

  rescale2 <- function (ynew, y = NULL)
  {
    # assert_that(is.numeric(ynew) && is.vector(ynew) && length(ynew) >=
    #               3)
    # assert_that((is.numeric(y) && is.vector(y) && length(y) >=
    #                3) || is.null(y))
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
                  ifelse(all(!is.null(bnd)),
                         paste0("s(", soap$coordinates[1], ", ", soap$coordinates[2],
                                ifelse(all(!is.null(soap$by)),
                                       paste0(", by = ", soap$by),
                                       ""),
                                ", bs = 'so', xt = list(bnd = bnd)) + "),
                         ""),
                  # paste(paste0("s(", covariates, ", bs = 'cs', k = 10)"), collapse = " + "),
                  paste(smoother(variable = covariates, bs = bs, complexity = complexity), collapse = " + "),
                  ifelse(!is.null(spline_to_add),
                         paste0(" + ", spline_to_add),
                         ""),
                  ifelse(!is.null(offset_effort),
                         paste0("+ offset(I(log(", offset_effort, ")))"),
                         ""))

  # cat(model)

  if (length(covariates) > treshold_discrete) {
    model_it <- mgcv::bam(as.formula(model),
                          method = "fREML",
                          nthreads = parallel::detectCores() - 1,
                          discrete=T,
                          family = mgcv::nb(),
                          knots = knots,
                          data = seg_data_scale
    )
  } else {
    model_it <- mgcv::gam(as.formula(model),
                          method = "REML",
                          family = mgcv::nb(),
                          knots = knots,
                          data = seg_data_scale
    )
  }

  # try(print(gratia::draw(model_it, rug = F)))

  csm <- summary(model_it)

  var <- na.omit(c(ifelse(all(!is.null(bnd)),
                          paste0("soap ", soap$coordinates[1], ",", soap$coordinates[2]),
                          NA),
                   covariates))

  table_var <- data.frame(pvalues = abs(csm$s.pv),
                          chi = csm$chi.sq) %>%
    dplyr::mutate(var = na.omit(c(var,
                                  ifelse(n() > length(var),
                                         1:(n() - length(var)),
                                         NA)))) %>%
    dplyr::filter(var %in% covariates) %>%
    arrange(pvalues, -chi)

  rho <- cor(seg_data_scale %>%
               dplyr::select(all_of(covariates)))

  diag(rho) <- 0

  all_thresholds <- c(.99, .9, .8, .7, .5)
  all_thresholds <- unique(c(all_thresholds[all_thresholds > max_correlation], max_correlation))

  # print(max(rho))

  remove_v <- c()

  ### correlations
  for (thresh_cor in all_thresholds) {
    while (any(rho > thresh_cor)) { #  | any(table_var$pvalues > .05)
      cat("Threshold:", thresh_cor, ". Max correlation:", print(round(max(rho), 2)), ". Max P-value:", max(round(table_var$pvalues, 3)), ".\n")

      rc <- which(map_dbl(1:nrow(rho), function(r) {return(max(abs(rho[r, ]), na.rm = T))}) > thresh_cor)

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
          paste0(to_rem  %>%
                   pull(pvalues) %>%
                   round(., 3), ", cor:", round(max(rho[match(to_rem %>%
                                                                pull(var),
                                                              covariates), ]), 3),
                 ") removed   ---  "),
          length(covariates) - 1, "left.\nRemoved:",
          paste(remove_v, collapse = ", "), "\n")

      to_rem <- to_rem %>%
        pull(var)

      covariates <- covariates[!(covariates %in% to_rem)]

      model <- paste0(response, " ~ 1 + ",
                      ifelse(all(!is.null(bnd)),
                             paste0("s(", soap$coordinates[1], ", ", soap$coordinates[2], ", bs = 'so', xt = list(bnd = bnd)) + "),
                             ""),
                      # paste(paste0("s(", covariates, ", bs = 'cs', k = 10)"), collapse = " + "),
                      paste(smoother(variable = covariates, bs = bs, complexity = complexity), collapse = " + "),
                      ifelse(!is.null(offset_effort),
                             paste0("+ offset(I(log(", offset_effort, ")))"),
                             ""),
                      ifelse(!is.null(spline_to_add),
                             paste0(" + ", spline_to_add),
                             ""))

      if (length(covariates) > treshold_discrete) {
        model_it <- mgcv::bam(as.formula(model),
                              method = "fREML",
                              nthreads = parallel::detectCores() - 1,
                              discrete=T,
                              family = mgcv::nb(),
                              knots = knots,
                              data = seg_data_scale
        )
      } else {
        model_it <- mgcv::gam(as.formula(model),
                              method = "REML",
                              family = mgcv::nb(),
                              knots = knots,
                              data = seg_data_scale
        )
      }

      # try(print(gratia::draw(model_it, rug = F)))
      # Sys.sleep(2)

      csm <- summary(model_it)

      var <- na.omit(c(ifelse(all(!is.null(bnd)),
                              paste0("soap ", soap$coordinates[1], ",", soap$coordinates[2]),
                              NA),
                       covariates))

      table_var <- data.frame(pvalues = abs(csm$s.pv),
                              chi = csm$chi.sq) %>%
        dplyr::mutate(var = na.omit(c(var,
                                      ifelse(n() > length(var),
                                             1:(n() - length(var)),
                                             NA)))) %>%
        dplyr::filter(var %in% covariates) %>%
        arrange(pvalues, -chi)

      rho <- cor(seg_data_scale %>%
                   dplyr::select(all_of(covariates)))

      diag(rho) <- 0

      print(table_var %>%
              dplyr::select(-var))

      gc()
    }
  }

  ### p values
  while (any(table_var$pvalues > Pvalue_max) | nrow(table_var) > nb_max_pred) {
    print(table_var[nrow(table_var), ] %>%
            dplyr::select(-var))
    if (all(table_var$pvalues <= Pvalue_max)) {
      cat("All covariates are now significant: removing covariates to reach the maximum number of predictors allowed", paste0("(", nb_max_pred, ")\n"))
    }
    cat("Removed    ---   ", nrow(table_var) - 1, "left\n")

    table_var <- table_var[-nrow(table_var), ]

    covariates <- table_var$var

    model <- paste0(response, " ~ 1 + ",
                    ifelse(all(!is.null(bnd)),
                           paste0("s(", soap$coordinates[1], ", ", soap$coordinates[2], ", bs = 'so', xt = list(bnd = bnd)) + "),
                           ""),
                    # paste(paste0("s(", covariates, ", bs = 'cs', k = 10)"), collapse = " + "),
                    paste(smoother(variable = covariates, bs = bs, complexity = complexity), collapse = " + "),
                    ifelse(!is.null(offset_effort),
                           paste0("+ offset(I(log(", offset_effort, ")))"),
                           ""),
                    ifelse(!is.null(spline_to_add),
                           paste0(" + ", spline_to_add),
                           ""))

    if (length(covariates) > treshold_discrete) {
      model_it <- mgcv::bam(as.formula(model),
                            method = "fREML",
                            nthreads = parallel::detectCores() - 1,
                            discrete=T,
                            family = mgcv::nb(),
                            knots = knots,
                            data = seg_data_scale
      )
    } else {
      model_it <- mgcv::gam(as.formula(model),
                            method = "REML",
                            family = mgcv::nb(),
                            knots = knots,
                            data = seg_data_scale
      )
    }

    # try(print(gratia::draw(model_it, rug = F)))
    # Sys.sleep(2)

    csm <- summary(model_it)

    var <- na.omit(c(ifelse(all(!is.null(bnd)),
                            paste0("soap ", soap$coordinates[1], ",", soap$coordinates[2]),
                            NA),
                     covariates))

    table_var <- data.frame(pvalues = abs(csm$s.pv),
                            chi = csm$chi.sq) %>%
      dplyr::mutate(var = na.omit(c(var,
                                    ifelse(n() > length(var),
                                           1:(n() - length(var)),
                                           NA)))) %>%
      dplyr::filter(var %in% covariates) %>%
      arrange(pvalues, -chi)

    print(table_var %>%
            dplyr::select(-var))

    gc()
  }

  model_actual <- mgcv::gam(as.formula(model),
                            method = "REML",
                            family = mgcv::nb(),
                            knots = knots,
                            data = seg_data_init)

  cat("\nFinal model:\n",
      model,
      "\n")

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

  return(final_file)
}

