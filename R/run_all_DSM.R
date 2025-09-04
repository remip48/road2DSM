run_all_DSM <- function (segdata_obs,
                         response = "ppho",
                         predictors,
                         spatial_options = list(by = NULL, complexity = NA),
                         soap = list(xt = NULL, knots = NULL),
                         max_cor = 0.5,
                         splines_bs = "cs",
                         complexity = 4,
                         use_loo = FALSE,
                         random = NULL,
                         offset_effort = NULL,
                         method = "REML",
                         nb_min_pred = 1,
                         nb_max_pred = 3,
                         k = 5,
                         intermediate_model_save = NULL,
                         load_saved_models = F,
                         fit_all_once = T,
                         save_list_models = "models.rds",
                         list_models_to_do = NULL,
                         force_one_off = NULL,
                         not_together = NULL,
                         force_include = NULL,
                         fit_with_actual_data = T,
                         use_select = F,
                         spline_to_add = NULL,
                         dataset_4correlation = NULL,
                         data_valid_loo = NULL,
                         list_knots = NULL,
                         ncores = NULL,
                         outfile = "logs.txt",
                         first_try = F,
                         first_try_AIC = F,
                         by_complexity = NULL,
                         by_te = F,
                         use_ti = F,
                         all_in_te = F,
                         list_knots2 = NULL,
                         splines_by = NULL,
                         no_by = NULL,
                         no_by2 = NULL,
                         month_spline_bs  = F,
                         likelihood = "negbin",
                         fit_models = T,
                         verbose = FALSE)
{
  tab = T
  weighted = FALSE

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

  combn2 <- function (x, m, FUN = NULL, simplify = TRUE, force_include = NULL, force_one_off = NULL, not_together = NULL, ...)
  {
    fun_filter <- function(all_mods, force_include = NULL, force_one_off = NULL, not_together = NULL) {

      if (!is.null(force_include)) {
        check <- map_dbl(force_include, function(p) {
          return(ifelse(p %in% all_mods, p, NA))
        })
        if (any(is.na(check))) {
          return(rep(NA, length(all_mods)))
        }
        # for (p in force_include) {
        # if (!(p %in% all_mods)) {
        #   return(rep(NA, length(all_mods)))
        # }
        # }
      }

      if (!is.null(force_one_off)) {
        # check <- apply(map_dfc(1:length(force_one_off), function(l) {
        check <- do.call("c", map(1:length(force_one_off), function(l) {
          # out <- data.frame(new = map_dbl(force_one_off[[l]], function(p) {
          #   return(ifelse(p %in% all_mods, p, NA))
          # }))
          out <- map_dbl(force_one_off[[l]], function(p) {
            return(ifelse(p %in% all_mods, p, NA))
          })
          # colnames(out) <- paste0("var", l)
          return(length(which(!is.na(out))))
        }))
        # , 2, function(ii) {
        #   return(ifelse(length(which(!is.na(ii))) == 1,
        #                 1,
        #                 0))
        # })

        # if (any(check != 1)) {
        if (any(check < 1)) {
          return(rep(NA, length(all_mods)))
        }
        # for (l in 1:length(force_one_off)) {
        #
        #   rem <- map_dfc(force_one_off[[l]],
        #                  function(p) {
        #                    out <- data.frame(new_col = (p %in% all_mods))
        #                    colnames(out) <- paste0("var", p)
        #                    return(out)
        #                  }) %>%
        #     mutate(id = 1:n())
        #
        #   if(any(rem[-length(rem)]) & !all(rem[-length(rem)])) {
        #   } else {
        #     return(rep(NA, length(all_mods)))
        #   }
        # }
      }

      if (!is.null(not_together)) {
        # check <- apply(map_dfc(1:length(not_together), function(l) {
        check <- do.call("c", map(1:length(not_together), function(l) {
          # out <- data.frame(new = map_dbl(not_together[[l]], function(p) {
          #   return(ifelse(p %in% all_mods, p, NA))
          # }))
          out <- map_dbl(not_together[[l]], function(p) {
            return(ifelse(p %in% all_mods, p, NA))
          })
          # colnames(out) <- paste0(out, l)
          return(length(which(!is.na(out))))
        }))
        # , 2, function(ii) {
        #   return(ifelse(length(which(!is.na(ii))) %in% c(0, 1),
        #                 1,
        #                 0))
        # })

        if (any(check > 1)) {
          return(rep(NA, length(all_mods)))
        }
        # for (l in 1:length(not_together)) {
        #   rem <- map_dfc(not_together[[l]],
        #                  function(p) {
        #                    out <- data.frame(new_col = (p %in% all_mods))
        #                    colnames(out) <- paste0("var", p)
        #                    return(out)
        #                  }) %>%
        #     mutate(id = 1:n())
        #
        #   if(!all(rem[-length(rem)])) {
        #   } else {
        #     return(rep(NA, length(all_mods)))
        #   }
        # }
      }

      return(all_mods)
    }
    stopifnot(length(m) == 1L, is.numeric(m))
    if (m < 0)
      stop("m < 0", domain = NA)
    if (is.numeric(x) && length(x) == 1L && x > 0 && trunc(x) ==
        x)
      x <- seq_len(x)
    n <- length(x)
    if (n < m)
      stop("n < m", domain = NA)
    x0 <- x
    if (simplify) {
      if (is.factor(x))
        x <- as.integer(x)
    }
    m <- as.integer(m)
    e <- 0
    h <- m
    a <- seq_len(m)
    nofun <- is.null(FUN)
    if (!nofun && !is.function(FUN))
      stop("'FUN' must be a function or NULL")
    len.r <- length(r <- if (nofun) x[a] else FUN(x[a], ...))
    # count <- as.integer(round(choose(n, m)))
    count <- round(choose(n, m), 0)
    if (simplify) {
      dim.use <- if (nofun)
        c(m, count)
      else {
        d <- dim(r)
        if (length(d) > 1L)
          c(d, count)
        else if (len.r != 1L)
          c(len.r, count)
        else c(d, count)
      }
    }
    if (simplify) {
      out <- matrix(r, nrow = len.r, ncol = count)
    }  else {
      out <- vector("list", count)
      # out[[1L]] <- r
      out[[1L]] <- fun_filter(all_mods = r,
                              force_include = force_include, force_one_off = force_one_off, not_together = not_together)
    }
    if (m > 0) {
      i <- 2L
      nmmp1 <- n - m + 1L
      while (a[1L] != nmmp1) {
        if (e < n - h) {
          h <- 1L
          e <- a[m]
          j <- 1L
        }
        else {
          e <- a[m - h]
          h <- h + 1L
          j <- 1L:h
        }
        a[m - h + j] <- e + j
        r <- x[a]#if (nofun)

        # else FUN(x[a], ...)
        # if (simplify) {
        # if (is.null(filter_comb)) {
        #   out[, i] <- r
        # } else {
        new_val <- fun_filter(all_mods = r,
                              force_include = force_include, force_one_off = force_one_off, not_together = not_together)

        # cat(i, length(out[, i]), length(new_val), "\n")

        out[, i] <- new_val
        # }
        # }
        # else out[[i]] <- r
        i <- i + 1L
      }
    }
    # library(tidyr)
    out <- out[, !is.na(out[1, ])]

    if (simplify) {
      if (is.factor(x0)) {
        levels(out) <- levels(x0)
        class(out) <- class(x0)
      }
      # dim(out) <- dim.use
    }
    return(out)
  }

  assert_that(is.logical(use_loo))
  assert_that(is.logical(verbose))
  segdata_obs <- segdata_obs %>% units::drop_units()
  knots <- soap$knots
  # xt <- soap$xt
  bnd <- soap$xt

  if (is.null(by_complexity)) {
    by_complexity <- complexity
  }

  segdata_obs <- as.data.frame(segdata_obs)

  get_k_best_models <- function(tab_model, k = 5, use_loo) {

    if(use_loo) {

      writeLines("\t* Using leave-one out for model selection")
      tab_model <- tab_model[order(tab_model$looic, decreasing = FALSE), ]

    } else {

      writeLines("\t* Using AIC for model selection")
      tab_model <- tab_model[order(tab_model$AIC, decreasing = FALSE), ]

    }

    tab_model <- tab_model[1:min(c(k, nrow(tab_model))), "index"]
    return(tab_model)
  }

  for (p in c(predictors, unique(unlist(splines_by)))) {
    if (!is.null(p) & !is.na(p)) {
      NAs <- which(is.na(segdata_obs %>%
                           pull(p)))
      if (length(NAs) > 0) {
        cat(length(NAs), "removed for", p, "\n")
        segdata_obs <- segdata_obs[-NAs, ]
      }
    }
  }
  cat(nrow(segdata_obs), "remained\n")

  if (!is.null(splines_by) & (!by_te & !use_ti)) {
    assertthat::assert_that(is.character(splines_by))
    assertthat::assert_that(segdata_obs %has_name% splines_by)
    if (!by_te) {
      list_is_splines_by <- lapply(splines_by, function(x) {
        is.factor(segdata_obs %>%
                    pull(get(x)))
      })
      df_is_splines_by <- do.call(rbind, list_is_splines_by)
      assertthat::assert_that(all(df_is_splines_by[, 1] ==
                                    TRUE), msg = "element of splines_by is not a factor")
    }
    if (all(levels(segdata_obs[, splines_by]) %in% segdata_obs[,
                                                               splines_by]) == FALSE) {
      cli::cli_alert_warning("One or multiple levels of {.val splines_by} are empty, consider relevel factor of {.val splines_by} : \n \t {.code segdata_obs$session <- droplevels(segdata_obs$session)}")
    }
  }
  if (!is.null(spatial_options$by)) {
    assertthat::assert_that(is.character(spatial_options$by))
    assertthat::assert_that(segdata_obs %has_name% spatial_options$by)
    list_is_spatial_by <- lapply(spatial_options$by, function(x) {
      is.factor(segdata_obs[, x])
    })
    df_is_spatial_by <- do.call(rbind, list_is_spatial_by)
    assertthat::assert_that(all(df_is_spatial_by[, 1] ==
                                  TRUE), msg = "element of spatial_options$by is not a factor")
    if (all(levels(segdata_obs[, spatial_options$by]) %in%
            segdata_obs[, spatial_options$by]) == FALSE) {
      cli::cli_alert_warning("One or multiple levels of {.val spatial_options$by} are empty, consider relevel factor of {.val spatial_options$by} : \n \t {.code segdata_obs$session <- droplevels(segdata_obs$session)}")
    }
  }
  assertthat::assert_that(segdata_obs %has_name% predictors)
  X <- segdata_obs
  if (all(is.null(splines_by))) {
    smoothers <- paste("s(", predictors, ", k = ", complexity,
                       ifelse(!is.null(splines_bs) & !is.na(splines_bs),
                              paste0(", bs = '", splines_bs, "'"),
                              ""), ")", sep = "")
  } else {
    if ((length(splines_bs) - 1) != ifelse(is.list(splines_by), length(splines_by[[1]]), length(splines_by))) {
      cat("WARNING: Check that splines_bs and splines_by are appropriate.\n")
    }
    smoothers <- do.call("c",
                         map(predictors, function(p) {
                           if (by_te & !use_ti) {
                             if (p %in% no_by) {
                               return(paste("s(", p, ", k = ", complexity,
                                            ifelse(all(!is.null(splines_bs)) & all(!is.na(splines_bs)),
                                                   paste0(", bs = '", splines_bs[1], "'"),
                                                   ""), ")", sep = ""))
                             } else if (p %in% no_by2) {
                               return(paste(do.call("c", map(splines_by, function(s) {
                                 paste("te(", p, ", ", paste(s[1], collapse = ", "), ", k = c(", complexity, ", ",
                                       paste(rep(by_complexity, length(s[1])), collapse = ", "),
                                       ")",
                                       ifelse(all(!is.null(splines_bs)) & all(!is.na(splines_bs)),
                                              paste0(", bs = c(", paste(paste0("'", splines_bs[1:2], "'"), collapse = ", "), ")"),
                                              ""),
                                       ifelse(!is.null(list_knots2), paste0(", xt = ", list_knots2), ""),
                                       ")", sep = "")
                               })), collapse = " + "))
                             } else {
                               return(paste(do.call("c", map(splines_by, function(s) {
                                 paste("te(", p, ", ", paste(s, collapse = ", "), ", k = c(", complexity, ", ", paste(rep(by_complexity, length(s)), collapse = ", "),
                                       ")",
                                       ifelse(all(!is.null(splines_bs)) & all(!is.na(splines_bs)),
                                              paste0(", bs = c(", paste(paste0("'", splines_bs, "'"), collapse = ", "), ")"),
                                              ""),
                                       ifelse(!is.null(list_knots), paste0(", xt = ", list_knots), ""),
                                       ")", sep = "")
                               })), collapse = " + "))
                             }
                           } else if (use_ti) {
                             if (p %in% no_by | p %in% no_by2) {
                               return(paste(paste("s(", p, ", k = ", complexity,
                                                  ifelse(all(!is.null(splines_bs)) & all(!is.na(splines_bs)),
                                                         paste0(", bs = '", splines_bs[1], "'"),
                                                         ""), ")", sep = "")))

                             } else {
                               return(paste(paste("ti(", p, ", k = ", complexity,
                                                  ifelse(all(!is.null(splines_bs)) & all(!is.na(splines_bs)),
                                                         paste0(", bs = '", splines_bs[1], "'"),
                                                         ""), ")", sep = ""),
                                            " + ",
                                            paste(do.call("c", map(1:length(splines_by), function(s) {
                                              paste("ti(", p, ", ",
                                                    splines_by[s], ", k = c(", by_complexity, ", ", by_complexity, ")",
                                                    ", bs = c('", splines_bs[1], "',", " '", splines_bs[s + 1], "'))", sep = "")
                                            })), collapse = " + ")))
                             }
                           } else {
                             if (p %in% no_by | p %in% no_by2) {
                               return(paste(paste("s(", p, ", k = ", complexity,
                                                  ifelse(all(!is.null(splines_bs)) & all(!is.na(splines_bs)),
                                                         paste0(", bs = '", splines_bs[1], "'"),
                                                         ""), ")", sep = "")))

                             } else {
                               return(paste(paste("s(", p, ", k = ", complexity,
                                                  ifelse(all(!is.null(splines_bs)) & all(!is.na(splines_bs)),
                                                         paste0(", bs = '", splines_bs[1], "'"),
                                                         ""), ")", sep = ""),
                                            " + ",
                                            paste(do.call("c", map(1:length(splines_by), function(s) {
                                              paste("s(", p, ", ", ifelse(splines_bs[s + 1] != "fs", "by = ", ""),
                                                    splines_by[s], ", k = ", by_complexity,
                                                    ", bs = '", splines_bs[s + 1], "')", sep = "")
                                            })), collapse = " + ")))
                             }
                           }
                         }))
  }

  if (by_te | use_ti) {
    segdata_obs[, c(c(predictors, unique(na.omit(unlist(splines_by))))#, "lon", "lat"
    )] <- apply(segdata_obs[,
                            c(c(predictors, unique(na.omit(unlist(splines_by))))
                              # , "lon", "lat"
                            )], 2, rescale2)
  } else {
    segdata_obs[, c(c(predictors)#, "lon", "lat"
    )] <- apply(segdata_obs[,
                            c(c(predictors)
                              # , "lon", "lat"
                            )], 2, rescale2)
  }

  if (is.null(soap$xt) && is.null(soap$knots)) {
    if (is.null(spatial_options$by) && is.na(spatial_options$complexity)) {
      intercept <- "~ 1"
    }
    else if (!is.null(spatial_options$by)) {
      intercept <- paste0("~ te(X, Y, bs = c('cs', 'cs'), by = ",
                          spatial_options$by, ", k = ", spatial_options$complexity,
                          ")")
    }
    else {
      intercept <- paste0("~ te(X, Y, bs = c('cs', 'cs'), k = ",
                          spatial_options$complexity, ")")
    }

  } else {
    if (is.null(spatial_options$by) && is.na(spatial_options$complexity)) {
      # intercept <- "~ 1"
      intercept <- paste0("~ 1 + s(X, Y, bs = 'so', xt = list(bnd = bnd))")
    }
    else if (!is.null(spatial_options$by)) {
      intercept <- paste0("~ s(X, Y, bs = 'so', by = ",
                          spatial_options$by, ", k = ", spatial_options$complexity,
                          ", xt = list(bnd = bnd))")
    }
    else {
      intercept <- paste0("~ s(X, Y, bs = 'so', k = ",
                          spatial_options$complexity, ", xt = list(bnd = bnd))")
    }
  }

  if (!is.null(month_spline_bs) & !is.na(month_spline_bs) & month_spline_bs != F) {
    assert_that(is.numeric(segdata_obs$month))
    intercept <- paste0(intercept, " + s(month, bs = 'cc'", ifelse(is.numeric(month_spline_bs),
                                                                   paste0(", k = ", month_spline_bs),
                                                                   ""), ", xt = list(range = c(0.5, 12.5)))")
  }

  if (use_ti) {
    intercept <- paste0(intercept, " + ti(", splines_by, ", bs = '", splines_bs[2], "'", ", k = ", complexity, ")")
  }

  if (!is.null(random)) {
    if (!all(random %in% names(segdata_obs))) {
      stop("Check random effect: no matching column in table \"segdata_obs\"")
    }
    else {
      for (kk in 1:length(random)) {
        if (!is.factor(segdata_obs %>%
                       pull(random[kk]))) {
          cat("Transform", random[kk], "to factor\n")
          segdata_obs[, random[kk]] <- factor(c(segdata_obs[,
                                                            random[kk]])[[1]], levels = c(unique(segdata_obs[,
                                                                                                             random[kk]]), "new_level"))
          X[, random[kk]] <- factor(c(X[, random[kk]])[[1]],
                                    levels = c(unique(X[, random[kk]][[1]]), "new_level"))
        }
        intercept <- paste(intercept, " + s(", random[kk],
                           ", bs = 're')", sep = "")
      }
    }
  }

  if (all(!is.null(spline_to_add))) {
    intercept <- paste0(intercept, " + ", spline_to_add)
  }

  # all_x <- lapply(1:nb_max_pred, combn, x = length(predictors))
  # all_x <- list(all_x[[1]],
  #               all_x[[length(all_x)]])

  # all_x <- lapply(c(1, nb_max_pred), combn, x = length(predictors))
  if (!is.null(force_include)) {
    if (any(!(force_include %in% predictors))) {
      stop("1 of force_include is not in predictors")
    }
    force_includei <- do.call("c",
                              lapply(force_include, function(f) {
                                return(match(f, predictors))
                              }))
  } else {
    force_includei <- force_include
  }

  if (!is.null(force_one_off)) {
    force_one_offi <- lapply(force_one_off, function(f) {
      do.call("c", lapply(f, function(fo) {
        return(match(fo, predictors))
      }))
    })
  } else {
    force_one_offi <- force_one_off
  }

  if (!is.null(not_together)) {
    not_togetheri <- lapply(not_together, function(f) {
      unname(na.omit(do.call("c", lapply(f, function(fo) {
        return(ifelse(fo %in% predictors,
                      match(fo, predictors),
                      NA))}))
      ))
    })

    length0 <- map_dbl(not_togetheri, length)

    if (any(length0 == 0)) {
      not_togetheri <- not_togetheri[-which(length0 %in% 0:1)]
    }

  } else {
    not_togetheri <- not_together
  }

  r <- ifelse(response == "ind",
              "count",
              response)

  ### initial way
  # all_x <- map(unique(c(1, nb_min_pred:nb_max_pred)), function(ii) {
  #
  #   if (ii == 1) {
  #     return(combn(x = length(predictors), m = 1))
  #   } else {
  #     out <- combn2(x = length(predictors), m = ii,
  #                   force_include = force_includei,
  #                   force_one_off = force_one_offi,
  #                   not_together = not_togetheri)
  #     if (length(out) == 0 | nrow(out) == 0 | all(is.na(out)) | all(is.null(out)) | all(is.nan(out))) {
  #       return(NA)
  #     } else {
  #       return(out)
  #     }
  #   }
  # })

  # if (!is.null(save_list_models)) {
  # if (!file.exists(save_list_models)) {
  if (all(is.null(dataset_4correlation))) {
    rho <- X %>%
      dplyr::select(all_of(predictors)) %>%
      drop_na() %>%
      cor()
  } else {
    rho <- dataset_4correlation %>%
      dplyr::select(all_of(predictors)) %>%
      drop_na() %>%
      cor()
  }

  diag(rho) <- 0

  if (is.null(list_models_to_do)) {
    ## first version
    all_x <- map(unique(c(1, nb_min_pred:nb_max_pred)), function(ii) {

      if (ii == 1) {
        return(combn(x = length(predictors), m = 1))
      } else {
        out <- combn2(x = length(predictors), m = ii,
                      force_include = force_includei,
                      force_one_off = force_one_offi,
                      not_together = not_togetheri)
        if (length(out) == 0 | nrow(out) == 0 | all(is.na(out)) | all(is.null(out)) | all(is.nan(out))) {
          return(NA)
        } else {
          return(out)
        }
      }
    })

    rem <- c()
    for (i in 1:length(all_x)) {
      if (all(is.na(all_x[[i]]))) {
        rem <- c(rem, i)
      }
    }

    if (length(rem) > 0) {
      all_x <- all_x[-rem]
    }

    if (nb_max_pred == 1) {
      rm_combn <- c(rep(0, length(predictors) + 1))
    } else {

      n_cores <- ifelse(parallel, ifelse(is.null(ncores),
                                         detectCores() - 1,
                                         ncores), 1)
      clust <- makeCluster(n_cores, outfile = outfile)
      doParallel::registerDoParallel(clust)

      rm_combn <- list()
      # rm_combn <- lapply(all_x[-1], function(mat) {
      for (m in 2:length(all_x)) {
        mat <- all_x[[m]]

        if (ncol(mat) > 1000) {
          print("PARALLELISE")
          out <- do.call("c", foreach(i = 1:ncol(mat),
                                      # .noexport = ls()[!(ls() %in% c("mat", "list_NA", "predictors", "X"))],
                                      .noexport = ls()[!(ls() %in% c("mat", "rho"))],
                                      .packages = "dplyr"
          ) %dopar% {
            cat(i, "IN PARALLEL /", ncol(mat), "\n")

            # rem <- unique(as.numeric(do.call("c", lapply(predictors[mat[, i]], function(p) {
            #   return(as.character(which(is.na(c(as.data.frame(X[, p]))[[1]]))))
            # }))))
            # rem <- list_NA %>%
            #   dplyr::filter(pred %in% predictors[mat[, i]]) %>%
            #   pull(id) %>%
            #   unique()
            #
            # if (length(rem) > 0) {
            #   rho <- cor(X[-rem, predictors[mat[, i]]])
            # } else {
            #   rho <- cor(X[, predictors[mat[, i]]])
            # }
            #
            # diag(rho) <- 0
            return(max(abs(as.numeric(rho[mat[, i], mat[, i]]))))
          })
        } else {
          # print("SERIAL")
          out <- sapply(1:ncol(mat), function(i) {
            # cat(i, "/", ncol(mat), "\n")

            # rem <- unique(as.numeric(do.call("c", lapply(predictors[mat[, i]], function(p) {
            #   return(as.character(which(is.na(c(as.data.frame(X[, p]))[[1]]))))
            # }))))
            # rem <- list_NA %>%
            #   dplyr::filter(pred %in% predictors[mat[, i]]) %>%
            #   pull(id) %>%
            #   unique()
            #
            # if (length(rem) > 0) {
            #   rho <- cor(X[-rem, predictors[mat[, i]]])
            # } else {
            #   rho <- cor(X[, predictors[mat[, i]]])
            # }
            #
            # diag(rho) <- 0
            return(max(abs(as.numeric(rho[mat[, i], mat[, i]]))))
          })
        }
        # return(out)
        rm_combn[[m-1]] <- out
      }
      # )
      # }
      # rm_combn <- c(c(rep(0, length(predictors) + 1)), unlist(rm_combn))
      if (nb_min_pred == 1) {
        rm_combn <- c(c(rep(0, length(predictors))), unlist(rm_combn))
      } else {
        rm_combn <- unlist(rm_combn)
      }
    }
    # mlist <- function(n, y, predictors) {
    #   paste(y, apply(X = combn(predictors, n), MARGIN = 2,
    #                  paste, collapse = " + "), sep = paste(intercept,
    #                                                        "+", sep = " "))
    # }
    r <- ifelse(response == "ind",
                "count",
                response)

    if (nb_min_pred > 1) {
      all_x <- all_x[-1]
    }
  } else {
    remx <- do.call("c",
                    map(list_models_to_do, function(m) {
                      if (length(m) > 1) {
                        if (all(is.null(dataset_4correlation))) {
                          rho <- cor(X %>%
                                       dplyr::select(all_of(m)))
                        } else {
                          rho <- cor(dataset_4correlation %>%
                                       dplyr::select(all_of(m)) %>%
                                       drop_na())
                        }

                        diag(rho) <- 0

                        return(max(abs(as.numeric(rho))))
                      } else {
                        return(0)
                      }
                    }))

    cat(length(which(remx >= max_cor)), "models removed due to correlated variables\n")

    list_models_to_do <- list_models_to_do[remx <= max_cor]

    all_x <- list(do.call("cbind",
                          map(list_models_to_do, function(m) {

                            out <- match(m, predictors)

                            if ((all(force_includei %in% out) | all(is.null(force_includei))) &
                                (all(map_lgl(force_one_offi, function(t) {
                                  any(t %in% out)
                                })) | all(is.null(force_one_offi))) &
                                (all(map_dbl(not_togetheri, function(t) {
                                  length(which(t %in% out))
                                }) < 2) | all(is.null(not_togetheri)))
                            ) {
                              return(out)
                            } else {
                              return(NULL)
                            }
                          })))
  }

  n_cores <- ifelse(parallel, ifelse(is.null(ncores),
                                     detectCores() - 1,
                                     ncores), 1)
  clust <- makeCluster(n_cores, outfile = outfile)
  doParallel::registerDoParallel(clust)

  all_mods <- do.call("c", map(all_x, function(mat) {
    r <- r
    intercept <- intercept
    smoothers <- smoothers
    all_in_te <- all_in_te
    predictors <- predictors
    splines_bs <- splines_bs
    complexity <- complexity
    return(do.call("c",
                   (foreach(i = 1:ncol(mat),
                            .noexport = ls()[!(ls() %in% c("intercept", "r", "smoothers", "mat", "all_in_te", "predictors", "splines_bs", "complexity"))],
                            .packages = "dplyr"
                   ) %dopar% {
                     return(paste(paste(r, intercept, sep = ""), paste(ifelse(all_in_te,
                                                                              paste0(ifelse(length(predictors[mat[, i]]) > 1, "te(", "s("), paste(predictors[mat[, i]], collapse = ", "), ", bs = c(",
                                                                                     paste(paste0("'", rep(splines_bs[1], length(predictors[mat[, i]])), "'"), collapse = ", "),
                                                                                     "), k = c(",
                                                                                     paste(rep(complexity, length(predictors[mat[, i]])), collapse = ", "), "))"),
                                                                              paste(smoothers[mat[, i]], collapse = " + ")), collapse = " + "), sep = " + "))
                   })
    ))
  }))

  if (is.null(list_models_to_do)) {
    cat(length(all_mods), "models for", length(rm_combn), "correlations\n")
    cat(length(which(rm_combn > max_cor)), "models removed due to correlated variables\n")
    all_mods <- all_mods[rm_combn <= max_cor]

  }

  stopCluster(clust)
  gc()

  # all_mods <- all_mods[which(rm_combn < max_cor)] ## used previously

  all_mods <- paste0(all_mods, "+ offset(I(log(", offset_effort, ")))")

  if (weighted) {
    covariable <- predictors
    w <- lapply(all_x, function(tab) {
      sapply(1:ncol(tab), function(j) {
        make_cfact_2(calibration_data = segdata_obs,
                     test_data = segdata_obs, var_name = covariable[tab[,
                                                                        j]], percent = FALSE, near_by = TRUE)
      })
    })
    w <- cbind(rep(1, nrow(segdata_obs)), do.call("cbind",
                                                  w))
    w <- w[, which(rm_combn < max_cor)]
  } else {
    w <- matrix(1, nrow = nrow(segdata_obs), ncol = #length(which(rm_combn <
                  # max_cor))
                  length(all_mods)
    )
  }
  my_dsm_fct <- function(x, tab = TRUE, segdata_obs, loo = FALSE, method,
                         bnd = soap$xt, knots = soap$knots, verbose = F) {
    if (verbose) {
      glue("\t\t* Fitting model with formula {x}\n")
    }
    test <- try({if (!is.null(knots)) {
      model <- mgcv::gam(as.formula(all_mods[x]),
                         data = segdata_obs,
                         method = method,
                         select = use_select,
                         drop.unused.levels = F,
                         knots = knots,
                         family = "nb")
    } else {
      model <- mgcv::gam(as.formula(all_mods[x]),
                         data = segdata_obs,
                         method = method,
                         select = use_select,
                         drop.unused.levels = F,
                         family = "nb")
    }})

    if (all(class(test) != "try-error")) {
      if (loo) {
        tab <- FALSE
        beta <- mvtnorm::rmvnorm(1000, mean = model$coefficients,
                                 sigma = model$Vp)
        Z <- predict(model, newdata = model$model, off.set = model$offset,
                     type = "lpmatrix")
        mu <- exp(as.matrix(beta %*% t(Z)))
        y <- model$model[[1]]
        lppd = switch(likelihood, negbin = {
          apply(mu, 1, function(iter) {
            w[, x] * dnbinom(y, size = model$family$getTheta(trans = TRUE),
                             mu = exp(model$offset) * iter, log = TRUE)
          })
        }, poisson = {
          apply(mu, 1, function(iter) {
            w[, x] * dpois(y, lambda = exp(model$offset) *
                             iter, log = TRUE)
          })
        }, tweedie = {
          apply(mu, 1, function(iter) {
            w[, x] * log(tweedie::dtweedie(y, xi = model$family$getTheta(trans = TRUE),
                                           mu = exp(model$offset) * iter, phi = model$sig2))
          })
        })
        out <- loo::loo.matrix(t(lppd), save_psis = TRUE)
      }
      else {
        out <- model
      }
      if (tab) {
        return(data.frame(model = all_mods[x], index = x,
                          Convergence = ifelse(model$converged, 1, 0),
                          AIC = model$aic, ResDev = model$deviance, NulDev = model$null.deviance,
                          ExpDev = 100 * round(1 - model$deviance/model$null.deviance,
                                               3), RMSE = qpcR::RMSE(model)))
      }
      else {
        return(out)
      }
    } else {
      return(NULL)
    }

  }

  if (!fit_models) {
    return(all_mods)
  }

  if (!load_saved_models) {
    writeLines("\t* Fitting all possible models, please wait")
    print(paste0("Example 1: ", all_mods[1]))

    cat("\n", length(all_mods), "models to fit\n")
    n_cores <- ifelse(parallel, ifelse(is.null(ncores),
                                       detectCores() - 1,
                                       ncores), 1)
    clust <- makeCluster(n_cores, outfile = outfile)
    doParallel::registerDoParallel(clust)

    if (first_try != F | first_try_AIC != F) {
      BAM_try <- foreach(x = 1:length(all_mods),
                         .combine = rbind,
                         .noexport = ls()[!(ls() %in% c("segdata_obs", "all_mods", "knots", "bnd", "method"))], # my_dsm_fct
                         .packages = c("qpcR", "mgcv", "dplyr")
      ) %dopar% {
        # out <- my_dsm_fct(x, segdata_obs = segdata_obs, all_mods = all_mods)
        cat(x, "MODELS IN PARALLEL /", length(all_mods), ":", all_mods[x], "\n")

        test <- try({if (!is.null(knots)) {
          bnd <- bnd
          model <- mgcv::bam(as.formula(all_mods[x]),
                             data = segdata_obs,
                             method = method,
                             # drop.unused.levels = F,
                             knots = knots,
                             discrete = T,
                             family = "nb")
        } else {
          model <- mgcv::bam(as.formula(all_mods[x]),
                             data = segdata_obs,
                             method = method,
                             # drop.unused.levels = F,
                             discrete = T,
                             family = "nb")
        }})
        if (all(class(test) != "try-error")) {
          if (FALSE) {
            tab <- FALSE
            beta <- mvtnorm::rmvnorm(1000, mean = model$coefficients,
                                     sigma = model$Vp)
            Z <- predict(model, newdata = model$model, off.set = model$offset,
                         type = "lpmatrix")
            mu <- exp(as.matrix(beta %*% t(Z)))
            y <- model$model[[1]]
            lppd = switch(likelihood, negbin = {
              apply(mu, 1, function(iter) {
                w[, x] * dnbinom(y, size = model$family$getTheta(trans = TRUE),
                                 mu = exp(model$offset) * iter, log = TRUE)
              })
            }, poisson = {
              apply(mu, 1, function(iter) {
                w[, x] * dpois(y, lambda = exp(model$offset) *
                                 iter, log = TRUE)
              })
            }, tweedie = {
              apply(mu, 1, function(iter) {
                w[, x] * log(tweedie::dtweedie(y, xi = model$family$getTheta(trans = TRUE),
                                               mu = exp(model$offset) * iter, phi = model$sig2))
              })
            })
            out <- loo::loo.matrix(t(lppd), save_psis = TRUE)
          }
          else {
            out <- model
          }

          if (T) {
            return(data.frame(model = all_mods[x], index = x,
                              Convergence = ifelse(model$converged, 1, 0),
                              AIC = model$aic, ResDev = model$deviance, NulDev = model$null.deviance,
                              ExpDev = 100 * round(1 - model$deviance/model$null.deviance,
                                                   3), RMSE = qpcR::RMSE(model)))
          }
          else {
            return(out)
          }
        } else {
          return(data.frame(model = all_mods[x], index = x,
                            Convergence = NA,
                            AIC = NA, ResDev = NA, NulDev = NA,
                            ExpDev = NA, RMSE = NA))
        }
      } %>%
        arrange(AIC)

      if (!is.na(first_try) & !is.null(first_try) & first_try != F) {
        BAM_try <- BAM_try[1:min(first_try, nrow(BAM_try)), ] %>%
          dplyr::group_by(index) %>%
          dplyr::mutate(new_variable = map_chr(model, function(m) {
            model <- str_split_1(m, fixed("+"))
            model <- model[length(model) - 1]
            model <- str_split_1(str_split_1(model, fixed("("))[2], fixed(","))
            model <- model[!str_detect(model, "k = ") & !str_detect(model, "bs = ") &
                             !str_detect(model, "xt = ")]
            # model <- paste(model, collapse = ", ")
            model <- str_remove_all(model, " ")
            model <- paste(model[model %in% predictors], collapse = ", ")
            return(model)
          })) %>%
          ungroup()
      }

      if (!is.na(first_try_AIC) & !is.null(first_try_AIC) & first_try_AIC != F) {
        BAM_try <- BAM_try %>%
          dplyr::filter(AIC <= (min(AIC, na.rm = T) + first_try_AIC)) %>%
          dplyr::group_by(index) %>%
          dplyr::mutate(new_variable = map_chr(model, function(m) {
            model <- str_split_1(m, fixed("+"))
            model <- model[length(model) - 1]
            model <- str_split_1(str_split_1(model, fixed("("))[2], fixed(","))
            model <- model[!str_detect(model, "k = ") & !str_detect(model, "bs = ") &
                             !str_detect(model, "xt = ")]
            # model <- paste(model, collapse = ", ")
            model <- str_remove_all(model, " ")
            model <- paste(model[model %in% predictors], collapse = ", ")
            return(model)
          })) %>%
          ungroup()
      }

      cat("Selected:\n")

      print(paste0(nrow(BAM_try), ": ",
                   BAM_try %>%
                     as.data.frame() %>%
                     dplyr::select(-c(model, Convergence, ResDev, NulDev, RMSE)) %>%
                     dplyr::mutate(dAIC = round(AIC - min(AIC, na.rm = T), 1),
                                   toprint = paste0(new_variable, " (", dAIC, ")")) %>%
                     pull(toprint) %>%
                     paste(., collapse = ",   ")))
    } else {
      BAM_try <- data.frame(model = all_mods)
    }

    stopCluster(clust)
    gc()

    clust <- makeCluster(min(n_cores, nrow(BAM_try)), outfile = outfile)
    doParallel::registerDoParallel(clust)

    if (fit_all_once) {
      all_models_fitted <- foreach(x = 1:length(all_mods),
                                   .noexport = ls()[!(ls() %in% c("segdata_obs", "all_mods", "knots", "bnd", "method", "BAM_try",
                                                                  "use_select"))], # my_dsm_fct
                                   .packages = c("qpcR", "mgcv", "dplyr")
      ) %dopar% {
        if (all_mods[x] %in% BAM_try$model) {
          cat(x, "FROM ALL MODELS IN PARALLEL /", length(all_mods), ":", all_mods[x], "\n")

          test <- try({if (!is.null(knots)) {
            bnd <- bnd
            model <- mgcv::gam(as.formula(all_mods[x]),
                               data = segdata_obs,
                               method = method,
                               select = use_select,
                               # drop.unused.levels = F,
                               knots = knots,
                               family = "nb")
          } else {
            model <- mgcv::gam(as.formula(all_mods[x]),
                               data = segdata_obs,
                               method = method,
                               select = use_select,
                               # drop.unused.levels = F,
                               family = "nb")
          }})
          if (all(class(test) != "try-error")) {
            return(model)
          } else {
            return(NULL)
          }
        } else {
          return(NULL)
        }
      }

      all_fits <- foreach(x = 1:length(all_mods),
                          .noexport = ls()[!(ls() %in% c("all_models_fitted", "all_mods"))], # my_dsm_fct
                          .packages = c("qpcR", "mgcv")
      ) %dopar% {
        cat(x, "MODELS IN PARALLEL /", length(all_mods), "\n")
        model <- all_models_fitted[[x]]

        if (all(!is.null(model))) {
          if (T) {
            return(data.frame(model = all_mods[x], index = x,
                              Convergence = ifelse(model$converged, 1, 0),
                              AIC = model$aic, ResDev = model$deviance, NulDev = model$null.deviance,
                              ExpDev = 100 * round(1 - model$deviance/model$null.deviance,
                                                   3), RMSE = qpcR::RMSE(model)))
          }
          else {
            return(out)
          }
        } else {
          return(data.frame(model = all_mods[x], index = x,
                            Convergence = NA,
                            AIC = NA, ResDev = NA, NulDev = NA,
                            ExpDev = NA, RMSE = NA))
        }
      }
      stopCluster(clust)
      gc()
    } else {
      all_fits <- foreach(x = 1:length(all_mods),
                          .noexport = ls()[!(ls() %in% c("segdata_obs", "all_mods", "knots", "bnd", "method", "BAM_try",
                                                         "use_select"))], # my_dsm_fct
                          .packages = c("qpcR", "mgcv", "dplyr")
      ) %dopar% {
        if (all_mods[x] %in% BAM_try$model) {
          cat(x, "MODELS IN PARALLEL /", length(all_mods), ":", all_mods[x], "\n")

          test <- try({if (!is.null(knots)) {
            bnd <- bnd
            model <- mgcv::gam(as.formula(all_mods[x]),
                               data = segdata_obs,
                               method = method,
                               select = use_select,
                               # drop.unused.levels = F,
                               knots = knots,
                               family = "nb")
          } else {
            model <- mgcv::gam(as.formula(all_mods[x]),
                               data = segdata_obs,
                               method = method,
                               select = use_select,
                               # drop.unused.levels = F,
                               family = "nb")
          }})
          if (all(class(test) != "try-error")) {
            if (FALSE) {
              tab <- FALSE
              beta <- mvtnorm::rmvnorm(1000, mean = model$coefficients,
                                       sigma = model$Vp)
              Z <- predict(model, newdata = model$model, off.set = model$offset,
                           type = "lpmatrix")
              mu <- exp(as.matrix(beta %*% t(Z)))
              y <- model$model[[1]]
              lppd = switch(likelihood, negbin = {
                apply(mu, 1, function(iter) {
                  w[, x] * dnbinom(y, size = model$family$getTheta(trans = TRUE),
                                   mu = exp(model$offset) * iter, log = TRUE)
                })
              }, poisson = {
                apply(mu, 1, function(iter) {
                  w[, x] * dpois(y, lambda = exp(model$offset) *
                                   iter, log = TRUE)
                })
              }, tweedie = {
                apply(mu, 1, function(iter) {
                  w[, x] * log(tweedie::dtweedie(y, xi = model$family$getTheta(trans = TRUE),
                                                 mu = exp(model$offset) * iter, phi = model$sig2))
                })
              })
              out <- loo::loo.matrix(t(lppd), save_psis = TRUE)
            }
            else {
              out <- model
            }

            if (T) {
              return(data.frame(model = all_mods[x], index = x,
                                Convergence = ifelse(model$converged, 1, 0),
                                AIC = model$aic, ResDev = model$deviance, NulDev = model$null.deviance,
                                ExpDev = 100 * round(1 - model$deviance/model$null.deviance,
                                                     3), RMSE = qpcR::RMSE(model)))
            }
            else {
              return(out)
            }
          } else {
            return(data.frame(model = all_mods[x], index = x,
                              Convergence = NA,
                              AIC = NA, ResDev = NA, NulDev = NA,
                              ExpDev = NA, RMSE = NA))
          }
        } else {
          return(data.frame(model = all_mods[x], index = x,
                            Convergence = NA,
                            AIC = NA, ResDev = NA, NulDev = NA,
                            ExpDev = NA, RMSE = NA))
        }
      }
      stopCluster(clust)
      gc()
    }

    all_fits <- do.call("rbind", all_fits)

    if (!is.null(intermediate_model_save)) {
      if (fit_all_once) {
        saveRDS(object = list(all_models_fitted, all_fits), file = intermediate_model_save)
      } else {
        saveRDS(object = all_fits, file = intermediate_model_save)
      }
    }
  } else {
    saved_models <- readRDS(intermediate_model_save)
    if(fit_all_once) {
      all_models_fitted <- saved_models[[1]]
      all_fits <- saved_models[[2]]
    } else {
      all_fits <- saved_models
    }
  }

  if (use_loo) {
    writeLines("\t* Estimating loocv on all models: please wait")
    gc()

    n_cores <- ifelse(parallel, ifelse(is.null(ncores),
                                       detectCores() - 1,
                                       ncores), 1)
    clust <- makeCluster(n_cores, outfile = outfile)
    doParallel::registerDoParallel(clust)

    if (fit_all_once) {
      all_psis <- foreach(x = 1:length(all_mods),
                          .noexport = ls()[!(ls() %in% c("all_models_fitted", "all_mods", "response", "offset_effort", "data_valid_loo",
                                                         # "tab",
                                                         # "knots", "bnd",
                                                         "w", "likelihood"))], # my_dsm_fct
                          .packages = c("qpcR", "mgcv", "mvtnorm", "loo", "dplyr")
      ) %dopar% {
        test <- try({{
          cat(x, "LOO IN PARALLEL /", length(all_mods), "\n")

          model <- all_models_fitted[[x]]

          # rm(all_models_fitted)
          # gc()

          beta <- mvtnorm::rmvnorm(1000, mean = model$coefficients,
                                   sigma = model$Vp)
          if (all(is.null(data_valid_loo))) {
            offset <- model$offset
            Z <- predict(model, newdata = model$model, off.set = offset,
                         type = "lpmatrix")
            mu <- exp(as.matrix(beta %*% t(Z)))
            y <- model$model[[1]]
          } else {
            offset <- data_valid_loo %>% pull(get(offset_effort))
            Z <- predict(model, newdata = data_valid_loo,
                         type = "lpmatrix")
            mu <- exp(as.matrix(beta %*% t(Z)))
            y <- data_valid_loo %>% pull(get(response))
          }
          lppd = switch(likelihood, negbin = {
            apply(mu, 1, function(iter) {
              w[, x] * dnbinom(y, size = model$family$getTheta(trans = TRUE),
                               mu = exp(offset) * iter, log = TRUE)
            })
          }, poisson = {
            apply(mu, 1, function(iter) {
              w[, x] * dpois(y, lambda = exp(offset) *
                               iter, log = TRUE)
            })
          }, tweedie = {
            apply(mu, 1, function(iter) {
              # w[, x] * log(tweedie::dtweedie(y, xi = model$family$getTheta(trans = TRUE),
              #                                mu = exp(offset) * iter, phi = model$sig2))
              stop("HERE")
            })
          })
          out <- loo::loo.matrix(t(lppd), save_psis = TRUE)

        } })

        if (all(class(test) == "try-error")) {
          return(NULL)
        } else {
          return(out)
        }
      }
      stopCluster(clust)
      gc()
    } else {
      all_psis <- foreach(x = 1:length(all_mods),
                          .noexport = ls()[!(ls() %in% c("segdata_obs", "all_mods", "tab", "knots", "bnd", "w", "method", "response", "offset_effort",
                                                         "data_valid_loo",
                                                         "likelihood", "use_select"))], # my_dsm_fct
                          .packages = c("qpcR", "mgcv", "mvtnorm", "loo", "dplyr")
      ) %dopar% {
        cat(x, "LOO IN PARALLEL /", length(all_mods), "\n")

        test <- try({if (!is.null(knots)) {
          bnd <- bnd
          model <- mgcv::gam(as.formula(all_mods[x]),
                             data = segdata_obs,
                             method = method,
                             select = use_select,
                             # drop.unused.levels = F,
                             knots = knots,
                             family = "nb")
        } else {
          model <- mgcv::gam(as.formula(all_mods[x]),
                             data = segdata_obs,
                             method = method,
                             select = use_select,
                             # drop.unused.levels = F,
                             family = "nb")
        }})

        if (all(class(test) != "try-error")) {
          if (TRUE) {
            tab <- FALSE
            beta <- mvtnorm::rmvnorm(1000, mean = model$coefficients,
                                     sigma = model$Vp)
            if (all(is.null(data_valid_loo))) {
              offset <- model$offset
              Z <- predict(model, newdata = model$model, off.set = offset,
                           type = "lpmatrix")
              mu <- exp(as.matrix(beta %*% t(Z)))
              y <- model$model[[1]]
            } else {
              offset <- data_valid_loo %>% pull(get(offset_effort))
              Z <- predict(model, newdata = data_valid_loo,
                           type = "lpmatrix")
              mu <- exp(as.matrix(beta %*% t(Z)))
              y <- data_valid_loo %>% pull(get(response))
            }
            lppd = switch(likelihood, negbin = {
              apply(mu, 1, function(iter) {
                w[, x] * dnbinom(y, size = model$family$getTheta(trans = TRUE),
                                 mu = exp(offset) * iter, log = TRUE)
              })
            }, poisson = {
              apply(mu, 1, function(iter) {
                w[, x] * dpois(y, lambda = exp(offset) *
                                 iter, log = TRUE)
              })
            }, tweedie = {
              apply(mu, 1, function(iter) {
                # w[, x] * log(tweedie::dtweedie(y, xi = model$family$getTheta(trans = TRUE),
                #                                mu = exp(offset) * iter, phi = model$sig2))
                stop("HERE")
              })
            })
            out <- loo::loo.matrix(t(lppd), save_psis = TRUE)
          }
          else {
            out <- model
          }
          if (tab) {
            return(data.frame(model = all_mods[x], index = x,
                              Convergence = ifelse(model$converged, 1, 0),
                              AIC = model$aic, ResDev = model$deviance, NulDev = model$null.deviance,
                              ExpDev = 100 * round(1 - model$deviance/model$null.deviance,
                                                   3), RMSE = qpcR::RMSE(model)))
          }
          else {
            return(out)
          }
        } else {
          return(NULL)
        }
      }
      stopCluster(clust)
      gc()
    }
    # }
    loo_ic_coefs <- do.call(rbind, lapply(all_psis, function(x) {
      if (is.null(x)) {
        return(NA)
      } else {
        return(x$estimates["looic", ])
      }
    }))
    colnames(loo_ic_coefs) <- c("looic", "se_looic")
    all_fits <- cbind(all_fits, loo_ic_coefs)
    all_fits$stacking_weights <- NA
    all_fits_best <- all_fits %>% slice_min(looic, n = k)
    index_order_best <- get_k_best_models(tab_model = all_fits_best,
                                          k = k, use_loo = TRUE)
    get_elpd_loo <- do.call("cbind", lapply(index_order_best,
                                            function(l) {
                                              all_psis[[l]]$pointwise[, "elpd_loo"]
                                            }))
    loow <- as.numeric(loo::stacking_weights(get_elpd_loo))
    all_fits$stacking_weights[index_order_best] <- loow
    all_fits_best_sw <- all_fits %>% slice_min(looic, n = k) %>% arrange(looic)
    index_order_sw <- all_fits_best_sw$index
    if (fit_all_once) {
      n_cores <- ifelse(parallel, ifelse(is.null(ncores),
                                         detectCores() - 1,
                                         ncores), 1)
      clust <- makeCluster(n_cores, outfile = outfile)
      doParallel::registerDoParallel(clust)

      best <- foreach(x = index_order_sw,
                      .noexport = ls()[!(ls() %in% c("X", "all_mods", "knots", "bnd", "index_order_sw", "method",
                                                     "use_select"))], # my_dsm_fct
                      .packages = c("qpcR", "mgcv", "dplyr")
      ) %dopar% {
        # out <- my_dsm_fct(x, segdata_obs = segdata_obs, all_mods = all_mods)
        cat(match(x, index_order_sw), "BEST MODELS IN PARALLEL /", length(index_order_sw), "\n")

        if (!is.null(knots)) {
          bnd <- bnd
          model <- mgcv::gam(as.formula(all_mods[x]),
                             data = X,
                             method = method,
                             select = use_select,
                             # drop.unused.levels = F,
                             knots = knots,
                             family = "nb")
        } else {
          model <- mgcv::gam(as.formula(all_mods[x]),
                             data = X,
                             method = method,
                             select = use_select,
                             # drop.unused.levels = F,
                             family = "nb")
        }

        # summary(model)
        # gratia::draw(model)

        if (F) {
          tab <- FALSE
          beta <- mvtnorm::rmvnorm(1000, mean = model$coefficients,
                                   sigma = model$Vp)
          Z <- predict(model, newdata = model$model, off.set = model$offset,
                       type = "lpmatrix")
          mu <- exp(as.matrix(beta %*% t(Z)))
          y <- model$model[[1]]
          lppd = switch(likelihood, negbin = {
            apply(mu, 1, function(iter) {
              w[, x] * dnbinom(y, size = model$family$getTheta(trans = TRUE),
                               mu = exp(model$offset) * iter, log = TRUE)
            })
          }, poisson = {
            apply(mu, 1, function(iter) {
              w[, x] * dpois(y, lambda = exp(model$offset) *
                               iter, log = TRUE)
            })
          }, tweedie = {
            apply(mu, 1, function(iter) {
              w[, x] * log(tweedie::dtweedie(y, xi = model$family$getTheta(trans = TRUE),
                                             mu = exp(model$offset) * iter, phi = model$sig2))
            })
          })
          out <- loo::loo.matrix(t(lppd), save_psis = TRUE)
        }
        else {
          out <- model
        }
        if (F) {
          return(data.frame(model = all_mods[x], index = x,
                            Convergence = ifelse(model$converged, 1, 0),
                            AIC = model$aic, ResDev = model$deviance, NulDev = model$null.deviance,
                            ExpDev = 100 * round(1 - model$deviance/model$null.deviance,
                                                 3), RMSE = qpcR::RMSE(model)))
        }
        else {
          return(out)
        }
      }
      best_std <- lapply(index_order_sw, function(x) {
        return(all_models_fitted[[x]])
      })

      stopCluster(clust)
      gc()
    }

    all_fits <- all_fits %>% arrange(looic)
  } else {
    all_fits_best <- all_fits %>% slice_min(AIC, n = k)
    all_fits$stacking_weights <- NA
    index_order_best <- as.numeric(get_k_best_models(tab_model = all_fits_best,
                                                     k = k, use_loo = FALSE))
    writeLines("\t* Estimating loocv on k best models: please wait")

    all_fits$stacking_weights[index_order_best] <- 1
    all_fits_best_sw <- all_fits %>% slice_min(AIC, n = k) %>% arrange(AIC)
    index_order_sw <- as.numeric(all_fits_best_sw$index)

    n_cores <- ifelse(parallel, ifelse(is.null(ncores),
                                       detectCores() - 1,
                                       ncores), 1)
    clust <- makeCluster(n_cores, outfile = outfile)
    doParallel::registerDoParallel(clust)

    if (fit_all_once) {
      if (fit_with_actual_data) {
        best <- foreach(x = index_order_sw,
                        .noexport = ls()[!(ls() %in% c("X", "all_mods", "knots", "bnd", "index_order_sw", "method",
                                                       "use_select"))], # my_dsm_fct
                        .packages = c("qpcR", "mgcv", "dplyr")
        ) %dopar% {
          # out <- my_dsm_fct(x, segdata_obs = segdata_obs, all_mods = all_mods)
          cat(match(x, index_order_sw), "BEST MODELS IN PARALLEL /", length(index_order_sw), "\n")

          if (!is.null(knots)) {
            bnd <- bnd
            model <- mgcv::gam(as.formula(all_mods[x]),
                               data = X,
                               method = method,
                               select = use_select,
                               # drop.unused.levels = F,
                               knots = knots,
                               family = "nb")
          } else {
            model <- mgcv::gam(as.formula(all_mods[x]),
                               data = X,
                               method = method,
                               select = use_select,
                               # drop.unused.levels = F,
                               family = "nb")
          }

          # summary(model)
          # gratia::draw(model)

          if (F) {
            tab <- FALSE
            beta <- mvtnorm::rmvnorm(1000, mean = model$coefficients,
                                     sigma = model$Vp)
            Z <- predict(model, newdata = model$model, off.set = model$offset,
                         type = "lpmatrix")
            mu <- exp(as.matrix(beta %*% t(Z)))
            y <- model$model[[1]]
            lppd = switch(likelihood, negbin = {
              apply(mu, 1, function(iter) {
                w[, x] * dnbinom(y, size = model$family$getTheta(trans = TRUE),
                                 mu = exp(model$offset) * iter, log = TRUE)
              })
            }, poisson = {
              apply(mu, 1, function(iter) {
                w[, x] * dpois(y, lambda = exp(model$offset) *
                                 iter, log = TRUE)
              })
            }, tweedie = {
              apply(mu, 1, function(iter) {
                w[, x] * log(tweedie::dtweedie(y, xi = model$family$getTheta(trans = TRUE),
                                               mu = exp(model$offset) * iter, phi = model$sig2))
              })
            })
            out <- loo::loo.matrix(t(lppd), save_psis = TRUE)
          }
          else {
            out <- model
          }
          if (F) {
            return(data.frame(model = all_mods[x], index = x,
                              Convergence = ifelse(model$converged, 1, 0),
                              AIC = model$aic, ResDev = model$deviance, NulDev = model$null.deviance,
                              ExpDev = 100 * round(1 - model$deviance/model$null.deviance,
                                                   3), RMSE = qpcR::RMSE(model)))
          }
          else {
            return(out)
          }
        }
      } else {
        best <- NULL
      }
      best_std <- lapply(index_order_sw, function(x) {
        return(all_models_fitted[[x]])
      })
    } else {
      best <- foreach(x = index_order_sw,
                      .noexport = ls()[!(ls() %in% c("X", "all_mods", "knots", "bnd", "index_order_sw", "method",
                                                     "use_select"))], # my_dsm_fct
                      .packages = c("qpcR", "mgcv", "dplyr")
      ) %dopar% {
        # out <- my_dsm_fct(x, segdata_obs = segdata_obs, all_mods = all_mods)
        cat(match(x, index_order_sw), "BEST MODELS IN PARALLEL /", length(index_order_sw), "\n")

        if (!is.null(knots)) {
          bnd <- bnd
          model <- mgcv::gam(as.formula(all_mods[x]),
                             data = X,
                             method = method,
                             select = use_select,
                             # drop.unused.levels = F,
                             knots = knots,
                             family = "nb")
        } else {
          model <- mgcv::gam(as.formula(all_mods[x]),
                             data = X,
                             method = method,
                             select = use_select,
                             # drop.unused.levels = F,
                             family = "nb")
        }

        # summary(model)
        # gratia::draw(model)

        if (F) {
          tab <- FALSE
          beta <- mvtnorm::rmvnorm(1000, mean = model$coefficients,
                                   sigma = model$Vp)
          Z <- predict(model, newdata = model$model, off.set = model$offset,
                       type = "lpmatrix")
          mu <- exp(as.matrix(beta %*% t(Z)))
          y <- model$model[[1]]
          lppd = switch(likelihood, negbin = {
            apply(mu, 1, function(iter) {
              w[, x] * dnbinom(y, size = model$family$getTheta(trans = TRUE),
                               mu = exp(model$offset) * iter, log = TRUE)
            })
          }, poisson = {
            apply(mu, 1, function(iter) {
              w[, x] * dpois(y, lambda = exp(model$offset) *
                               iter, log = TRUE)
            })
          }, tweedie = {
            apply(mu, 1, function(iter) {
              w[, x] * log(tweedie::dtweedie(y, xi = model$family$getTheta(trans = TRUE),
                                             mu = exp(model$offset) * iter, phi = model$sig2))
            })
          })
          out <- loo::loo.matrix(t(lppd), save_psis = TRUE)
        }
        else {
          out <- model
        }
        if (F) {
          return(data.frame(model = all_mods[x], index = x,
                            Convergence = ifelse(model$converged, 1, 0),
                            AIC = model$aic, ResDev = model$deviance, NulDev = model$null.deviance,
                            ExpDev = 100 * round(1 - model$deviance/model$null.deviance,
                                                 3), RMSE = qpcR::RMSE(model)))
        }
        else {
          return(out)
        }
      }
      best_std <- foreach(x = index_order_sw,
                          .noexport = ls()[!(ls() %in% c("segdata_obs", "all_mods", "knots", "bnd", "index_order_sw", "method",
                                                         "use_select"))], # my_dsm_fct
                          .packages = c("qpcR", "mgcv", "dplyr")
      ) %dopar% {
        # out <- my_dsm_fct(x, segdata_obs = segdata_obs, all_mods = all_mods)
        cat(match(x, index_order_sw), "BEST MODELS_STD IN PARALLEL /", length(index_order_sw), "\n")

        if (!is.null(knots)) {
          bnd <- bnd
          model <- mgcv::gam(as.formula(all_mods[x]),
                             data = segdata_obs,
                             method = method,
                             select = use_select,
                             # drop.unused.levels = F,
                             knots = knots,
                             family = "nb")
        } else {
          model <- mgcv::gam(as.formula(all_mods[x]),
                             data = segdata_obs,
                             method = method,
                             select = use_select,
                             # drop.unused.levels = F,
                             family = "nb")
        }

        # summary(model)
        # gratia::draw(model)

        if (F) {
          tab <- FALSE
          beta <- mvtnorm::rmvnorm(1000, mean = model$coefficients,
                                   sigma = model$Vp)
          Z <- predict(model, newdata = model$model, off.set = model$offset,
                       type = "lpmatrix")
          mu <- exp(as.matrix(beta %*% t(Z)))
          y <- model$model[[1]]
          lppd = switch(likelihood, negbin = {
            apply(mu, 1, function(iter) {
              w[, x] * dnbinom(y, size = model$family$getTheta(trans = TRUE),
                               mu = exp(model$offset) * iter, log = TRUE)
            })
          }, poisson = {
            apply(mu, 1, function(iter) {
              w[, x] * dpois(y, lambda = exp(model$offset) *
                               iter, log = TRUE)
            })
          }, tweedie = {
            apply(mu, 1, function(iter) {
              w[, x] * log(tweedie::dtweedie(y, xi = model$family$getTheta(trans = TRUE),
                                             mu = exp(model$offset) * iter, phi = model$sig2))
            })
          })
          out <- loo::loo.matrix(t(lppd), save_psis = TRUE)
        }
        else {
          out <- model
        }
        if (FALSE) {
          return(data.frame(model = all_mods[x], index = x,
                            Convergence = ifelse(model$converged, 1, 0),
                            AIC = model$aic, ResDev = model$deviance, NulDev = model$null.deviance,
                            ExpDev = 100 * round(1 - model$deviance/model$null.deviance,
                                                 3), RMSE = qpcR::RMSE(model)))
        }
        else {
          return(out)
        }
      }
    }

    all_fits <- all_fits %>% arrange(AIC)
  }
  if (parallel) {
    try(stopCluster(clust))
  }
  return(list(n_best = k, all_fits_binded = all_fits, best_models = best_std, all_models_tried = all_mods,
              best_models4plotting = best, splines_by = splines_by,
              random = random))
}
