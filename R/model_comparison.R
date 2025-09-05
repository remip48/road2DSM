#' Title
#'
#' @param run_models
#' @param seg_data
#' @param variable
#' @param effort_column
#' @param version_preds
#' @param output_file
#' @param log1p_trans
#' @param grid_folder
#' @param static_grid
#' @param prediction_folder
#' @param block_file
#' @param data_file
#' @param sub_area_analysis_file
#' @param study_area
#' @param correct_bias
#' @param save_results_bias_corrected
#' @param use_threshold
#' @param quantile_mgcv_fixed
#' @param threshold
#' @param breaks_plot
#' @param labels_plot
#' @param corr_groupsize
#' @param response
#' @param subspecies
#' @param filter_year_month_not_in
#' @param run_all
#' @param outfile
#' @param save_posterior_distribution
#' @param n_cores
#'
#' @return
#' @export
#'
#' @examples
model_comparison <- function(run_models, # output from run_all_DSM
                                      seg_data, # segments used for run_all_DSM. Should not be an sf object.
                                      variable,
                                      effort_column,
                                      version_preds = as.character(lubridate::today()),
                                      output_file = paste0(version_preds, "_model_comparison"), # without extension, will be the name of html file
                                      log1p_trans = NULL,
                                      grid_folder,
                                      static_grid, # grid with variables (static) that were not included in the extract_grids. Must have geometry, and use the exact same grid that the one used for extract_grid.
                                      prediction_folder,
                                      block_file = NULL, # add here the sf file containing your sub-blocks for your area if you want to use groupsizes as response (containing the column Name for the sub-names), otherwise let NULL
                                      data_file = NULL, # add here the intial GPS points to calculate group sizes if needed and print observations. Must contain the column AU for which groupsizes are averaged per AU, and the column Platform (aerial or ship).
                                      sub_area_analysis_file = NULL, # path for sf file containing sub-areas if you want to investigate abundance per sub-area in addition to globally. Must contain the column Name for each sub-area.
                                      study_area = NULL, # in case you want to predict only on a part of your prediction grid
                                      correct_bias = T,
                                      save_results_bias_corrected = F,
                                      use_threshold = T,
                                      quantile_mgcv_fixed = "mgcv", # "mgcv", "fixed" or quantile
                                      threshold = 1.5, ## either a multiplication factor to use with max(mgcv_density) if quantile_mgcv_fixed = "mgcv",
                                      # or a (fixed) maximum overall allowed density if quantile_mgcv_fixed = "fixed",
                                      # or a quantile (0 to 1) to use on the density values if quantile_mgcv_fixed = "quantile"
                                      breaks_plot = c(-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.25, 1.5, 1000),
                                      labels_plot = c("0.00 - 0.10", "0.11 - 0.20", "0.21 - 0.30", "0.31 - 0.40", "0.41 - 0.50",
                                                       "0.51 - 0.75", "0.76 - 1.00", "1.01 - 1.25", "1.26 - 1.5", "> 1.50"),
                                      corr_groupsize, # correction factor for groupsize, if groupsize is used, to multiply the predicted densities.
                                      response, # should contain "group" if it is modelling the number of groups rather than of individuals. But please
                                      # use n_SpeciesCode for the number of individuals, or n_group_SpeciesCode for the number of groups!
                                      subspecies, # in case response contains several species. The groupsize estimate will then account for all the subspecies
                                      filter_year_month_not_in = "0000-00", # year and month that should not be used for prediction
                                      run_all = F,
                                      outfile = "log.txt",
                                      save_posterior_distribution = F,
                                      n_cores = NULL) {

  ## this function is adapted to the type of data of ITAW.

  cat("Running version", version_preds, "-", ifelse(run_all,
                                                    "running EVERYTHING.",
                                                    "loading EXISTING files if any."), "\n")
  if (correct_bias) {
    cat("Bias will be corrected.\n")
    if (save_posterior_distribution) {
      cat("Posterior distribution of cell predictions will be saved under",
          paste0(version_preds, "_for_uncertainties.RData."),
          "They may take large space on the hard disk!")
    }
    if (use_threshold) {
      cat("Threshold applied to the pseudo-posterior distribution:", case_when(quantile_mgcv_fixed == "mgcv" ~ paste0(threshold, " * max(mgcv_predictions) inds/km²"),
                                                                               quantile_mgcv_fixed == "fixed" ~ paste0(threshold, " inds/km²"),
                                                                               quantile_mgcv_fixed == "quantile" ~ paste0(threshold*100, "% quantile of pseudo posterior distribution")
      ), "\n")
    }
  }

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

  if (all(!is.null(study_area))) {
    study_area <- study_area %>%
      st_transform(crs = 3035) %>%
      group_by() %>%
      dplyr::summarise(do_union = F) %>%
      st_cast("MULTIPOLYGON")
  }
  # library(rmarkdown)
  # library(mgcv)
  # library(gratia)
  # library(dplyr)
  # library(ggplot2)
  # library(ggnewscale)
  # library(tidyr)
  # library(viridis)
  # library(knitr)
  # library(units)
  # library(doParallel)
  # library(stringr)
  # library(sf)
  # library(purrr)
  # library(rnaturalearth)

  if (is.null(n_cores)) {
    n_cores <- detectCores() - 1
    cat("Parallel processing using", n_cores, "cores.\n")
  }

  ##############
  setup_chunk <- quote({

    knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE
    )
    # library(mgcv)
    # library(gratia)
    # library(dplyr)
    # library(ggplot2)
    # library(ggnewscale)
    # library(tidyr)
    # library(viridis)
    # library(knitr)
    # library(units)
    # library(doParallel)
    # library(stringr)
    # library(sf)
    # library(purrr)
    # library(rnaturalearth)
  })
  ##############
  obsn_chunk <- quote({
    cat("Variable used:", paste(sort(variable), collapse = ", "), ".<br><br>")
    log1p_trans <- log1p_trans[log1p_trans %in% variable]

    calibdata <- seg_data %>%
      as.data.frame()

    if (length(log1p_trans) > 0 & all(!is.na(log1p_trans))) {
      for (k in log1p_trans) {
        newcol <- calibdata %>%
          pull(k) %>%
          log1p()

        calibdata <- calibdata %>%
          dplyr::select(-all_of(k)) %>%
          mutate(new = newcol)

        colnames(calibdata)[colnames(calibdata) == "new"] <- k
      }
    }

    calibdata <- calibdata %>%
      dplyr::select(X, Y, year, all_of(c(variable, response, effort_column))) %>%
      drop_na()

    calibdata[, variable] <- apply(calibdata[, variable], 2, rescale2)

    do_plot <- F
    cat("There were", length(which(seg_data %>% pull(response) > 0)),
        "segments with sightings of", paste0(response, "."),
        "In total, this dataset contains", sum(seg_data %>% dplyr::pull(response)),
        "sightings of", response, "and", length(which(seg_data %>% pull(response) == 0)), "segments of 10 km with 0 sightings", "\n")

    cat(paste0("\n\n", ifelse(length(log1p_trans) > 0 & all(!is.na(log1p_trans)),
                              paste(log1p_trans, collapse = ", "),
                              "No covariates")), "were log-transformed.")
  })
  ##############
  chunk_best_models <- quote({
    kable(run_models$all_fits_binded[1:length(run_models$best_models), ] %>%
            dplyr::select(-c(index, ResDev, NulDev, Convergence)) %>%
            dplyr::mutate(model = 1:n(),
                          stacking_weights = round(stacking_weights, 3),
                          AIC = round(AIC, 2),
                          RMSE = round(RMSE, 3),
                          looic = round(looic, 3),
                          se_looic = round(se_looic, 3)) %>%
            as.data.frame())
  })
  ##############
  # Construct R Markdown text
  rmd_text <- c(
    "---",
    paste0("title: \"", "Results from the density surface models", "\""),
    paste0("date: '`r Sys.Date()`'"),
    "output: html_document",
    "---",
    "",
    "```{r setup, include=FALSE}",
    paste(deparse(setup_chunk), collapse = "\n"),
    "```",
    "``` {r obsn, echo = F, results='asis'}",
    paste(deparse(obsn_chunk), collapse = "\n"),
    "```",
    "```{r best_models, echo=FALSE}",
    paste(deparse(chunk_best_models), collapse = "\n"),
    "```"
  )

  ###################
  for (i in 1:length(run_models$best_models)) {
    chunk_modeli <- quote({
      summary(run_models$best_models[[i]])
      gratia::draw(run_models$best_models4plotting[[i]], rug = F)
      qq.gam(run_models$best_models4plotting[[i]], rep = 1000)

      cat("\n\n#### ASPE & Ratio of the number of observed", ifelse(str_detect(response, "group"), "groups", "individuals"),
          "/ number of predicted", ifelse(str_detect(response, "group"), "groups", "individuals"), "\n")

      p <- predict(run_models$best_models[[i]], newdata=calibdata, type='response')
      dens <-as.numeric(p)
      obs_n <- calibdata %>%
        pull(response)
      ASPE <- (sum((obs_n - dens)^2, na.rm=TRUE) / nrow(calibdata))
      cat("ASPE =", ASPE, "\n")

      ratio <- obs_n/dens

      summary(ratio)

      ggplot() +
        geom_histogram(data = data.frame(Ratio = ratio) %>%
                         dplyr::filter(!is.na(Ratio)), aes(x = Ratio)) +
        scale_y_sqrt(name = "Count")

      checks <- as.data.frame(calibdata)
      checks$new <- checks[, response]
      checks$Value <- "Observed value"
      checks <- checks %>%
        rbind(checks %>%
                mutate(new = unname(as.numeric(p)),
                       Value = "Predicted value")) %>%
        mutate(Value = as.factor(Value))

      rootg <- ggplot() +
        geom_histogram(data = checks %>%
                         dplyr::filter(as.character(Value) == "Observed value"), aes(x = new, fill = Value), alpha = 1, binwidth = 1) +
        scale_fill_manual(values = viridis::viridis(256)[1]) +
        new_scale_fill() +
        geom_histogram(data = checks %>%
                         dplyr::filter(as.character(Value) == "Predicted value"), aes(x = new, fill = Value), alpha = .5, binwidth = 1) +
        scale_fill_manual(values = viridis::viridis(256)[256]) +
        theme(legend.title = element_blank())

      if (max(checks$new, na.rm = T) > 20) {
        rootg <- rootg +
          scale_x_sqrt()
      }

      print(rootg)
    })
##################
    rmd_text <- c(
      rmd_text,
      paste0("```{r best_models", i, ", echo=FALSE}"),
      paste(deparse(chunk_modeli), collapse = "\n"),
      "```"
    )
  }

  ##############
  run_preds_chunk <- quote({
    to_runm <- run_models$all_fits_binded[1:length(run_models$best_models), ] %>%
      dplyr::mutate(id = 1:n()) %>%
      dplyr::filter(stacking_weights >= .1) %>%
      pull(id)

    fit_models <- run_models$all_fits_binded[to_runm, ]
    fit_models$stacking_weights <- fit_models$stacking_weights / sum(fit_models$stacking_weights)

    models <- run_models$best_models[to_runm]

    ls <- list.files(grid_folder)
    ls <- ls[!str_detect(ls, static_grid)]
    ls <- do.call("c", lapply(ls, function(l) {
      if (any(str_detect(l, fixed(paste0(as.character(unique(seg_data$year)),
                                         "-")))) &
          !any(filter_year_month_not_in, function(ym) {
            str_detect(l, fixed(ym))
          })) {
        return(l)
      } else {
        return(NULL)
      }
    }))

    static <- read_sf(paste0(grid_folder, "/", static_grid)) %>%
      st_transform(crs = 3035)

    if (all(!is.null(study_area))) {
      static <- st_intersection(static, study_area)
    }

    static <- static %>%
      dplyr::mutate(areakm2 = units::drop_units(st_area(.)) / 10^6)

    cc <- st_coordinates(st_centroid(static))

    static <- static %>%
      dplyr::mutate(X = cc[, 1],
                    Y = cc[, 2])

    rm(cc)

    calibdata <- seg_data %>%
      as.data.frame()

    if (length(log1p_trans) > 0 & all(!is.na(log1p_trans))) {

      for (k in log1p_trans) {
        newcol <- calibdata %>%
          pull(k) %>%
          log1p()

        calibdata <- calibdata %>%
          dplyr::select(-all_of(k)) %>%
          mutate(new = newcol)

        colnames(calibdata)[colnames(calibdata) == "new"] <- k
      }
    }

    calibdata <- calibdata %>%
      dplyr::select(X, Y, year, all_of(c(variable, response, effort_column))) %>%
      drop_na()

    if (!dir.exists(paste0(folder_predictions, "/", response))) {
      dir.create(paste0(folder_predictions, "/", response))
    }

    if (!file.exists(paste0(folder_predictions, "/", response, "/", version_preds,
                            "_average_predictions.gpkg")) | run_all) {

      cl <- makeCluster(n_cores, outfile = outfile)
      registerDoParallel(cl)

      run <- foreach::foreach(f = ls,
                              .packages = c("sf", "dplyr", "purrr", "stringr"),
                              .noexport = ls()[!(ls() %in% c("static", "ls", "variable", "year", "calibdata", "run_all", "rescale2",
                                                             "models", "to_runm", "version_preds", "response", "log1p_trans", "effort_column", "grid_folder", "folder_predictions"))]
      ) %dopar% {
        cat(match(f, ls), "/", length(ls), "\n")

        current_predgrid <- readRDS(paste0(grid_folder, "/", f)) %>%
          st_drop_geometry() %>%
          mutate(year = as.numeric(sub(".*_(\\d{4})-\\d{2}-\\d{2}\\.rds$", "\\1", f)))

        gridi <- static %>%
          dplyr::select(id, areakm2, all_of(variable[variable %in% colnames(static)]), X, Y) %>%
          left_join(current_predgrid,
                    by = "id")

        date <- str_sub(str_split_1(f, "_")[3], 1, 10)

        grid <- gridi %>%
          st_drop_geometry()

        grid$new <- grid$areakm2
        colnames(grid)[which(colnames(grid) == "new")] <- effort_column

        grid <- grid %>%
          as.data.frame()

        if (length(log1p_trans) > 0 & all(!is.na(log1p_trans))) {
          for (k in log1p_trans) {
            newcol <- grid %>%
              pull(k) %>%
              log1p()

            grid <- grid %>%
              dplyr::select(-all_of(k)) %>%
              mutate(new = newcol)

            colnames(grid)[colnames(grid) == "new"] <- k
          }
        }

        new_var <- map_dfc(variable, function(v) {

          gridv <- grid %>%
            pull(v)

          ref <- calibdata %>%
            pull(v)

          out <- data.frame(new = rescale2(ynew = gridv, y = ref))

          colnames(out) <- v

          return(out)

        })

        grid[, variable] <- new_var

        run <- map(1:length(models), function(i) {

          if (!dir.exists(paste0(folder_predictions, "/", response, "/model",
                                 to_runm[i]))) {
            dir.create(paste0(folder_predictions, "/", response, "/model",
                              to_runm[i]))
          }

          if (!file.exists(paste0(folder_predictions, "/", response, "/model",
                                  to_runm[i], "/", version_preds, "_prediction_", date, ".gpkg")) | run_all) {
            pred <- predict(models[[i]], grid, type = "response")

            new_pred <- gridi %>%
              mutate(n_pred = unname(as.numeric(pred)),
                     density_pred = n_pred / areakm2)

            saveRDS(new_pred %>%
                      st_drop_geometry(), paste0(folder_predictions, "/", response, "/model",
                                                 to_runm[i], "/", version_preds, "_prediction_", date, ".rds"))
          }

          return(NULL)

        })

        return(NULL)

      }
    }

    to_runmi <- NA

    ls <- list.files(paste0(folder_predictions, "/", response, "/model",
                            to_runm[1]))

    ls <- ls[str_detect(ls, "_prediction_") & str_detect(ls, version_preds)]

    dates <- unique(do.call("c", map(ls, function(l) {
      return(str_sub(str_split_1(l, "_")[3], 1, 10))
    })))

    kable(data.frame(date = sort(dates)) %>%
            dplyr::mutate(Year = as.character(lubridate::year(date))) %>%
            group_by(Year) %>%
            dplyr::summarise(From = min(date),
                             To = max(date)) %>%
            ungroup(),
          caption = "Dates used for prediction:")

    to_runm_files <- paste0("model", na.omit(c(to_runmi, to_runm)))

    sp <- last(str_split_1(response, "_"))

    if (!is.null(block_file)) {
      blocks <- read_sf(block_file) %>%
        st_transform(crs = 3035)

      if (all(!is.null(study_area))) {
        blocks <- st_intersection(blocks, study_area)
      }

    } else {
      blocks <- static %>%
        st_transform(crs = 3035) %>%
        st_buffer(units::set_units(100, m)) %>%
        group_by() %>%
        dplyr::summarise(do_union = T) %>%
        st_cast("MULTIPOLYGON")

      if (!("Name" %in% colnames(blocks))) {
        blocks$Name <- "Full_area"
      }
    }

    blocksti <- blocks

    if (!file.exists(paste0(folder_predictions, "/", response, "/", version_preds,
                            "_average_predictions.gpkg")) | run_all) {

      mod_pred <- foreach::foreach(d = dates,
                                   # .combine = rbind,
                                   .packages = c("sf", "dplyr", "purrr", "stringr", "ggplot2", "png"),
                                   .noexport = ls()[!(ls() %in% c("to_runm_files", "group_sizes_all", "do_plot", "response", "WorkDir",
                                                                  "version_preds", "folder_predictions"))]) %dopar% {

                                                                    run <- map_dfr(to_runm_files, function(i) {

                                                                      f <- list.files(paste0(folder_predictions, "/", response, "/",
                                                                                             i))

                                                                      f <- f[str_detect(f, d) & str_detect(f, version_preds)]

                                                                      out <- readRDS(paste0(folder_predictions, "/", response, "/",
                                                                                            i, "/", f)) %>%
                                                                        # st_drop_geometry() %>%
                                                                        dplyr::select(id, X, Y, areakm2, n_pred, density_pred) %>%
                                                                        dplyr::mutate(model = i)

                                                                      return(out)
                                                                    })

                                                                    return(run %>%
                                                                             st_drop_geometry())

                                                                  }

      mod_predictions <- do.call("rbind", mod_pred) %>%
        group_by(id, X, Y, areakm2, model) %>%
        dplyr::summarise(Avg_density = mean(density_pred, na.rm = T),
                         SD_density = sd(density_pred, na.rm = T),
                         CV_density = SD_density / Avg_density,
                         SE_density = SD_density / sqrt(n()))

      pred_grid <- static %>%
        dplyr::select(id) %>%
        left_join(mod_predictions,
                  by = "id")

      pred_grid_cent <- pred_grid %>%
        dplyr::filter(model == unique(model)[1])

      ######

      ######

      final <- map_dfr(unique(blocksti$Name), function(b) {

        t <- blocksti %>%
          dplyr::filter(Name == b) %>%
          group_by(Name) %>%
          dplyr::summarise(do_union = F) %>%
          st_cast("MULTIPOLYGON") %>%
          rename(AU = Name) %>%
          st_make_valid()

        out <- pred_grid_cent %>%
          st_make_valid() %>%
          arrange(id) %>%
          st_intersection(t) %>%
          dplyr::mutate(block = b,
                        area_km2_cropped = units::drop_units(st_area(.)) / 10^6) %>%
          dplyr::select(id, X, Y, block, area_km2_cropped, areakm2) %>%
          dplyr::mutate(id_cropped = paste(id, block, 1:n(), sep = "_")) %>%
          dplyr::rename(areakm2_precropped = areakm2)

        return(out)

      })

      final <- final %>%
        arrange(id_cropped) %>%
        dplyr::mutate(AU = str_sub(block, 1, 2)) %>%
        left_join(pred_grid %>%
                    st_drop_geometry() %>%
                    dplyr::select(id, model, Avg_density),
                  by = "id") %>%
        rename(areakm2 = area_km2_cropped)

      write_sf(final,
               append = F,
               paste0(folder_predictions, "/", response, "/", version_preds,
                      "_average_predictions.gpkg"))
    } else {
      final <- read_sf(paste0(folder_predictions, "/", response, "/", version_preds,
                              "_average_predictions.gpkg"))
    }

    if (!is.null(data_file)) {
      data <- read_sf(data_file)

      if (all(!is.na(subspecies))) {

        data <- data %>%
          dplyr::mutate(species = ifelse(species %in% subspecies, sp, species))
      }

      if (sp != "ppho") {
        data$podsize <- data$SIG_nmb
      }

      strats <- unique(data$AU)

      if (str_detect(response, "group")) {
        groupsizes <- T

        group_sizes_all <- data %>%
          st_drop_geometry() %>%
          dplyr::filter(species == sp &
                          !(subj %in% c("ll", "lx", "xl", "xx", "pp", "lp", "pl", "xp", "px","")) & sub_sig %in% c("m", "g")) %>%
          dplyr::select(species, podsize) %>%
          dplyr::mutate(podsize = podsize * ifelse(!is.na(corr_groupsize),
                                                   corr_groupsize,
                                                   1))

        group_sizes <- data %>%
          st_drop_geometry() %>%
          dplyr::mutate(AU = AU,
                        podsize = podsize * ifelse(!is.na(corr_groupsize),
                                                   corr_groupsize,
                                                   1)) %>%
          dplyr::filter(species == sp &
                          !(subj %in% c("ll", "lx", "xl", "xx", "pp", "lp", "pl", "xp", "px","")) & sub_sig %in% c("m", "g")) %>%
          group_by(AU) %>%
          dplyr::summarise(n_groups = n(),
                           mean = mean(podsize, na.rm = T),
                           median = median(podsize, na.rm = T),
                           SD = sd(podsize, na.rm = T),
                           max = max(podsize, na.rm = T))

        group_sizes <- group_sizes %>%
          as.data.frame()

        if (any(!(strats %in% group_sizes$AU))) {
          group_sizes <- rbind(group_sizes,
                               map_dfr(strats[which(!(strats %in% group_sizes$AU))],
                                       function(st) {
                                         return(data.frame(AU = st,
                                                           n_groups = 0,
                                                           mean = 0,
                                                           median = 0,
                                                           SD = 0,
                                                           max = 0))
                                       }))
        }
      } else {
        groupsizes <- F
        # for_plotting_gp_size <- NA

        group_sizes_all <- data.frame(mean = ifelse(!is.na(corr_groupsize),
                                                    corr_groupsize,
                                                    1), SD = ifelse(!is.na(corr_groupsize),
                                                                    corr_groupsize,
                                                                    1), max = ifelse(!is.na(corr_groupsize),
                                                                                     corr_groupsize,
                                                                                     1))
        group_sizes <- data.frame(mean = ifelse(!is.na(corr_groupsize),
                                                corr_groupsize,
                                                1), SD = ifelse(!is.na(corr_groupsize),
                                                                corr_groupsize,
                                                                1), max = ifelse(!is.na(corr_groupsize),
                                                                                 corr_groupsize,
                                                                                 1)) %>%
          group_by(mean, SD, max) %>%
          dplyr::reframe(AU = strats)
      }

    } else {
      groupsizes <- F
      # for_plotting_gp_size <- NA

      strats <- unique(blocksti$Name)

      group_sizes_all <- data.frame(mean = ifelse(!is.na(corr_groupsize),
                                                  corr_groupsize,
                                                  1), SD = ifelse(!is.na(corr_groupsize),
                                                                  corr_groupsize,
                                                                  1), max = ifelse(!is.na(corr_groupsize),
                                                                                   corr_groupsize,
                                                                                   1))
      group_sizes <- data.frame(mean = ifelse(!is.na(corr_groupsize),
                                              corr_groupsize,
                                              1), SD = ifelse(!is.na(corr_groupsize),
                                                              corr_groupsize,
                                                              1), max = ifelse(!is.na(corr_groupsize),
                                                                               corr_groupsize,
                                                                               1)) %>%
        group_by(mean, SD, max) %>%
        dplyr::reframe(AU = strats)
    }

    final <- final %>%
      dplyr::mutate(AU = block)
  })
  ##############
  rplot_chunk <- quote({
    ggplot() +
      geom_sf(data = final %>%
                dplyr::filter(areakm2 > 0) %>%
                dplyr::filter(!is.na(Avg_density)) %>%
                left_join(group_sizes %>%
                            dplyr::select(AU, mean),
                          by = "AU") %>%
                mutate(Avg_density = Avg_density * mean), aes(fill = cut(Avg_density,
                                                                         breaks = breaks_plot,
                                                                         labels = labels_plot
                )), color = NA) +
      facet_wrap(~ model) +
      scale_fill_viridis_d(drop = F, name = "Average density\n(ind/km2)") +
      labs(title = "Predictions")
  })
  ##############
  rallplot_chunk <- quote({
    if (!is.null(data_file)) {
      for_plotting_gp_size <- data %>%
        st_drop_geometry() %>%
        dplyr::mutate(AU = AU) %>%
        dplyr::filter(species == sp & !(subj %in% c("ll", "lx", "xl", "xx", "pp", "lp", "pl", "xp", "px","")) & sub_sig %in% c("m", "g")) %>%
        dplyr::select(species, podsize, lon, lat, AU) %>%
        dplyr::filter(podsize > 0) %>%
        dplyr::mutate(#Platform = "Plane",
          podsize = podsize * ifelse(!is.na(corr_groupsize),
                                     corr_groupsize,
                                     1))
    } else {
      for_plotting_gp_size <- seg_data %>%
        dplyr::mutate(podsize = get(response)) %>%
        dplyr::filter(podsize > 0) %>%
        st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
        st_transform(crs = 3035) %>%
        dplyr::mutate(X = st_coordinates(.)[,1],
                      Y = st_coordinates(.)[,2]) %>%
        st_drop_geometry()
    }

    if (!("Platform" %in% colnames(for_plotting_gp_size))) {
      for_plotting_gp_size$Platform <- "All"
    }

    blocksg <- blocks %>%
      st_transform(crs = 3035)

    print(ggplot() +
            geom_sf(data = blocksg, fill = NA, color = "red") +
            # facet_wrap(~ Platform) +
            geom_point(data = for_plotting_gp_size %>%
                         arrange(podsize),
                       aes(x = X, y = Y, color = cut(podsize,
                                                     breaks = c(-1, 1, 5, 10, 20, 50, 1000),
                                                     labels = c("1", "2 - 5", "6 - 10",
                                                                "11 - 20", "21 - 50", "> 50"))), size = 1) +
            scale_color_viridis_d(name = "Nb of sightings") +
            labs(title = ifelse(is.null(data_file) & groupsizes, "Number of groups observed per platform",
                                paste0("Number of sightings per segment",
                                       ifelse(!is.na(corr_groupsize) & corr_groupsize != 1,
                                              paste0("\nmultiplied by correction factor = ", corr_groupsize),
                                              "")))))
    # }

    if (groupsizes) {

      gp_print <- group_sizes

      if (!(all(group_sizes$AU == "Full_area"))) {
        gp_print <- rbind(group_sizes_all %>%
                            dplyr::mutate(AU = "All",
                                          n_groups = sum(group_sizes$n_groups)))
      }

      gp_print <- gp_print %>%
        dplyr::mutate(mean = round(mean, 2), SD = round(SD, 2))
      colnames(gp_print) <- c("Block", "Number of groups", "Mean groupsize", "Median groupsize", "SD groupsize", "Maximum groupsize")
      print(kable(gp_print, caption = "Groupsizes"))
    }
  })

  rmd_text <- c(
    rmd_text,
    "```{r run, echo=FALSE}",
    paste(deparse(run_preds_chunk), collapse = "\n"),
    "```",
    "```{r plots_pred, echo=FALSE, out.width = '1500px', out.height = '700px'}",
    paste(deparse(rplot_chunk), collapse = "\n"),
    "```",
    "```{r abund_pred, echo=FALSE, results='asis'}",
    paste(deparse(rallplot_chunk), collapse = "\n"),
    "```"
  )

  if (correct_bias) {
    n <- 1000

    rmd_text <- c(
      rmd_text,
      "### CV estimate",
      paste0("We have the predictions of abundance per grid cell (", nrow(static),
             "), per draw (", n, ") and per day (", length(unique(dates)),
             "). For the estimates of CV and abundances presented in this section, we follow the method:"),
      "* For the estimates of total abundance:",
      "1.	We sum the cells per day and per draw",
      "2.	We average the daily abundances per simulation",
      "3.	We calculate the mean & SD in this pseudo-posterior distribution",
      "* For the estimates of CV per cell:",
      "1.	We average the daily abundance per cell and draw",
      "2.	We calculate the mean & SD in this pseudo-posterior distribution of cell abundances"
    )

    chunk_bias <- quote({
      library(MASS)

      calibdata <- seg_data %>%
        as.data.frame()

      if (length(log1p_trans) > 0 & all(!is.na(log1p_trans))) {
        for (k in log1p_trans) {
          newcol <- calibdata %>%
            pull(k) %>%
            log1p()

          calibdata <- calibdata %>%
            dplyr::select(-all_of(k)) %>%
            mutate(new = newcol)

          colnames(calibdata)[colnames(calibdata) == "new"] <- k
        }
      }

      calibdata <- calibdata %>%
        dplyr::select(X, Y, year, all_of(c(variable, response, effort_column))) %>%
        drop_na()

      final_gp <- final %>%
        dplyr::filter(model == unique(model)[1]) %>%
        dplyr::mutate(area_no_BS = areakm2) %>%
        st_drop_geometry() %>%
        left_join(group_sizes %>%
                    dplyr::select(AU, mean),
                  by = "AU") %>%
        group_by(id) %>%
        dplyr::summarise(groupsize = sum(mean * areakm2) / unique(areakm2_precropped),
                         groupsize_no_BS = sum(mean * area_no_BS) / unique(areakm2_precropped),
                         ratio = groupsize_no_BS / groupsize,
                         areakm2_precropped = unique(areakm2_precropped)#,
        ) %>%
        ungroup() %>%
        arrange(id) %>%
        as.data.frame()

      final_corrected <- final %>%
        dplyr::mutate(area_no_BS = areakm2) %>%
        st_drop_geometry() %>%
        arrange(model, id) %>%
        left_join(group_sizes %>%
                    dplyr::select(AU, mean),
                  by = "AU") %>%
        group_by(model, id_cropped) %>%
        dplyr::summarise(id = id,
                         groupsize = sum(mean * areakm2) / unique(areakm2_precropped),
                         groupsize_no_BS = sum(mean * area_no_BS) / unique(areakm2_precropped),
                         ratio = groupsize_no_BS / groupsize,
                         areakm2_precropped = unique(areakm2_precropped),
                         Avg_density = mean(Avg_density, na.rm = T)) %>%
        ungroup() %>%
        as.data.frame()

      if (quantile_mgcv_fixed == "mgcv") {
        final_threshold <- max(final_corrected$Avg_density*threshold, na.rm = T)
      } else if (quantile_mgcv_fixed == "fixed") {
        final_threshold <- threshold
      }

      if (!is.null(sub_area_analysis_file)) {
        if (!file.exists(paste0(folder_predictions, "/", response, "/", version_preds,
                                "_AU_grid_predictions.rds")) | run_all) {
          hp <- read_sf(sub_area_analysis_file) %>%
            st_transform(crs = 3035) %>%
            st_cast("MULTIPOLYGON")

          cl <- makeCluster(min(c(length(unique(hp$Name)), detectCores() - 1)), outfile = "log.txt")
          registerDoParallel(cl)

          perHP <- do.call("rbind", foreach::foreach(n = unique(hp$Name),
                                                     .packages = c("dplyr", "sf"),
                                                     .noexport = ls()[!(ls() %in% c("final", "hp"))]
          ) %dopar% {

            bb <- st_bbox(hp %>%
                            dplyr::filter(Name == n))

            print(n)

            out <- final %>%
              dplyr::filter(model == unique(model)[1]) %>%
              dplyr::filter(as.numeric(X) >= (bb[1] - 2 * mean(sqrt(final$areakm2))) & as.numeric(X) <= (bb[3] + 2 * mean(sqrt(final$areakm2))) &
                              as.numeric(Y) >= (bb[2] - 2 * mean(sqrt(final$areakm2))) & as.numeric(Y) <= (bb[4] + 2 * mean(sqrt(final$areakm2)))) %>%
              st_intersection(hp %>%
                                dplyr::filter(Name == n) %>%
                                group_by() %>%
                                dplyr::summarise(do_union = F) %>%
                                st_cast("MULTIPOLYGON")) %>%
              mutate(in_hp_areakm2 = units::drop_units(st_area(.)) / 10^6,
                     Name = n)

            return(out)

          })

          list_id_perHP <- map_dfr(unique(hp$Name), function(n) {
            temp <- perHP %>%
              st_drop_geometry() %>%
              dplyr::filter(Name == n) %>%
              group_by(id) %>%
              dplyr::summarise(in_hp_areakm2 = sum(in_hp_areakm2, na.rm = T))

            out <- final_gp %>%
              left_join(temp %>%
                          dplyr::select(id, in_hp_areakm2),
                        by = "id") %>%
              dplyr::mutate(ratio_area_in = in_hp_areakm2 / areakm2_precropped,
                            Name = n)

            return(out)
          })

          list_id_perHP <- as.data.frame(list_id_perHP)

          saveRDS(list_id_perHP, paste0(folder_predictions, "/", response, "/", version_preds,
                                        "_AU_grid_predictions.rds"))
        } else {
          list_id_perHP <- readRDS(paste0(folder_predictions, "/", response, "/", version_preds,
                                          "_AU_grid_predictions.rds"))
        }
      }

      all_models <- run_models

      for (i in to_runm) {
        for_model <- i
        model_folder <- i
        GridDir <- paste0(folder_predictions, "/",
                          response, "/model", model_folder)
        set.seed(11)
        gam_mod <- all_models$best_models[[for_model]]
        model <- gam_mod
        covariates <- do.call("c", lapply(all_models$best_models[[i]]$smooth, function(s) {return(s$term)}))
        covariates <- covariates[!(covariates %in% c("X", "Y"))]
        log1p_trans <- log1p_trans[log1p_trans %in% covariates]
        ls <- list.files(GridDir)
        ls <- ls[str_detect(ls, "_prediction_") & str_detect(ls, version_preds)]
        dd <- do.call("c", lapply(ls, function(l) {
          return(str_sub(str_split_1(l, "_")[3], 1, 10))
        }))

        final_correctedm <- final_corrected %>%
          dplyr::filter(model == paste0("model", i))

        if (!file.exists(paste0(GridDir, "/", version_preds, "_Abundance_CV_nonbias.RData")) | run_all) {

          # source(paste0(WorkDir_CV, "/20211105_FctPredictions.R"))
          # source(paste0(WorkDir_CV, "SCANS-IV_model", "/08_OSPAR_Matthieu_CV_estimate/20211105_FctPredictions.R"))

          if (!file.exists(paste0(GridDir, "/", version_preds, "_for_uncertainties.RData")) | run_all) {
            cl <- makeCluster(detectCores() - 1)
            registerDoParallel(cl)

            to_pred <- (foreach(d = dd,
                                .packages = c("sf", "dplyr"),
                                .noexport = ls()[!(ls() %in% c("dd", "GridDir", "effort_column", "log1p_trans", "covariates", "ls",
                                                               "calibdata", "response", "rescale2"))]
            ) %dopar% {

              f_d <- ls[which(d == dd)]

              grid <- readRDS(file = paste(GridDir, f_d, sep = "/")) %>%
                mutate(date = d
                )

              colnames(grid)[colnames(grid) == "areakm2"] <- effort_column

              if (length(log1p_trans) > 0 & all(!is.na(log1p_trans))) {
                for (k in log1p_trans) {
                  newcol <- grid %>%
                    pull(k) %>%
                    log1p()

                  grid <- grid %>%
                    dplyr::select(-all_of(k)) %>%
                    mutate(new = newcol)

                  colnames(grid)[colnames(grid) == "new"] <- k
                }
              }

              for (k in covariates) {

                gridv <- grid %>%
                  pull(k)

                ref <- calibdata %>%
                  pull(k)

                newcol = rescale2(ynew = gridv, y = ref)

                grid <- grid %>%
                  dplyr::select(-all_of(k)) %>%
                  mutate(new = newcol)

                colnames(grid)[colnames(grid) == "new"] <- k
              }

              return(grid %>%
                       as.data.frame())

            })

            stopCluster(cl)
            gc()

            parallel <- T

            ##########
            beta <- function(fitted_model, n_sim = 200, seed_id = 19790428) {
              if(!is.null(seed_id)) {
                set.seed(seed_id)
              } else {
                set.seed(sample.int(1e6, size = 1))
              }
              theta <- mvtnorm::rmvnorm(n_sim, mean = fitted_model$coefficients, sigma = fitted_model$Vp)
              return(theta)
            }

            n_cores <- floor(nc_cores / 4)
            cat("Parallel processing with", n_cores, "for pseudo-posterior distribution estimate.\n")

            cl <- makeCluster(n_cores)
            registerDoParallel(cl)

            theta <- beta(fitted_model = gam_mod, n_sim = 1e3)

            allX <- foreach(X = to_pred,
                            .packages = c("mgcv"),
                            .noexport = ls()[!(ls() %in% c("gam_mod"))]) %dopar% {
                              return(predict.gam(gam_mod, newdata = X, type = "lpmatrix"))
                            }

            stopCluster(cl)
            gc()

            cl <- makeCluster(floor(detectCores() / 4))
            registerDoParallel(cl)

            allpred <- foreach(l = allX,
                               .noexport = ls()[!(ls() %in% c("theta"))]) %dopar% { return(exp(l %*% t(theta))) }

            stopCluster(cl)
            rm(theta)
            rm(allX)
            rm(to_pred)
            gc()

            ##########

            if (save_posterior_distribution){save(allpred,
                 file = paste0(GridDir, "/", version_preds, "_for_uncertainties.RData"))}

          } else {
            load(paste0(GridDir, "/", version_preds, "_for_uncertainties.RData"))
          }

          gc()

          if (use_threshold) {

            if (quantile_mgcv_fixed == "quantile") {
              qq_with_cv_no_BSn <- do.call("c", map(1:length(allpred), function(x) {
                return(quantile(c(allpred[[x]]), threshold, na.rm = T))
              }))

              q999 <- median(qq_with_cv_no_BSn)

            } else {
              q999 <- final_threshold
            }

            allpred <- map(1:length(allpred), function(x) {
              infv <- which(allpred[[x]] >= q999)

              if (length(infv) > 0) {
                # allpred[[x]][infv] <- NA
                allpred[[x]][infv] <- q999
              }

              return(allpred[[x]])
            })
          } else {
            q999 <- NULL
          }

          with_cv_no_BSn <- map_dfr(1:length(allpred), function(x) {
            out <- data.frame(Abundance = colSums(allpred[[x]] * final_gp[, "areakm2_precropped"] * final_gp[, "groupsize"] * final_gp[, "ratio"], na.rm = T),
                              day = x) %>%
              dplyr::mutate(simulation = 1:n())

            return(out)
          })

          cat("Parallel processing with", round(n_cores / 2, 0), "for BIAS estimate and correction.\n")
          clust <- makeCluster(round(n_cores / 2, 0))
          registerDoParallel(clust)

          per_cell <- (foreach(sim = 1:n,
                               .packages = c("dplyr", "purrr"),
                               .noexport = ls()[!(ls() %in% c("allpred", "final_gp"))]
          ) %dopar% {
            out <- map_dfr(allpred, function(x) {
              return(data.frame(Cell_abundance = x[ ,sim] *
                                  final_gp[, "areakm2_precropped"]
              ) %>%
                dplyr::mutate(
                  cell = final_gp$id))
            }) %>%
              as.data.frame()

            out <- out %>%
              group_by(cell) %>%
              dplyr::summarise(
                mean_Cell_abundance = mean(Cell_abundance, na.rm = T),
                Median_Cell_abundance = median(Cell_abundance, na.rm = T)#,
              ) %>%
              ungroup() %>%
              dplyr::mutate(groupsize = final_gp[, "groupsize"],
                            MaxInf_Cell_abundance = 0)

            return(out)
          }
          )
          per_cell <- do.call("rbind", per_cell)

          stopCluster(clust)
          gc()

          if (!is.null(sub_area_analysis_file)) {

            with_cp_perHP <- map_dfr(unique(list_id_perHP$Name[!is.na(list_id_perHP$Name)]), function(n) {
              n_list_id_perHP <- list_id_perHP %>%
                dplyr::filter(Name == n) %>%
                dplyr::mutate(ratio_area_in = ifelse(is.na(ratio_area_in), 0, ratio_area_in)) #%>%

              cell_in <- which(n_list_id_perHP$ratio_area_in > 0)

              out <- map_dfr(1:length(allpred), function(x) {

                out <- data.frame(Abundance = colSums(allpred[[x]][cell_in, ] * n_list_id_perHP[cell_in, "areakm2_precropped"] * n_list_id_perHP[cell_in, "groupsize"] * n_list_id_perHP[cell_in, "ratio_area_in"], na.rm = T),
                                  day = x) %>%
                  dplyr::mutate(simulation = 1:n(),
                                AU = n)

                return(out)
              })

              return(out)
            })

            save(with_cv_no_BSn, per_cell, with_cp_perHP, q999,
                 file = paste0(GridDir, "/", version_preds, "_Abundance_CV_nonbias.RData"))
          } else {
            save(with_cv_no_BSn, per_cell, q999,
                 file = paste0(GridDir, "/", version_preds, "_Abundance_CV_nonbias.RData"))
          }

        } else {
          load(file = paste0(GridDir, "/", version_preds, "_Abundance_CV_nonbias.RData"))

        }

        with_cv <- with_cv_no_BSn

        cat("\n\n#### Model", i, "\n")

        summary_abundi <- with_cv %>% ## with_cv: total abundance per day and per simulation. Here gives the average abundance per simulation.
          group_by(simulation) %>%
          dplyr::summarise(median = median(Abundance, na.rm = T),
                           Abundance = mean(Abundance, na.rm = T)) %>%
          ungroup()

        ## Calculation of total abundance based on the mgcv response output.
        mle <- final_correctedm %>%
          dplyr::mutate(#ratio = ifelse(species == "ppho", 1, ratio),
            abundance = Avg_density * areakm2_precropped * groupsize * ratio) %>%
          group_by() %>%
          dplyr::summarise(total = sum(abundance, na.rm = T)) %>%
          pull(total)

        cat("MGCV abundance:", round(mle, 0), "inds.\n<br>")

        bias <- 2 * (mean(summary_abundi$Abundance, na.rm = T) - mle)
        # bias <- 2 * (median(summary_abundi$median, na.rm = T) - mle)

        if (use_threshold) {
          cat("Bias on total abundance estimate:", round(bias, 0), "inds.",
              "\n<br>Threshold used:", round(q999, 4), ifelse(groupsizes,
                                                              "groups/km².",
                                                              "inds/km².\n"))
        } else {
          cat("Bias on total abundance estimate:", round(bias, 0), "inds.",
              "\n<br>No threshold used for bias estimate.\n")
        }

        ### here substract on mean or median
        summary_abundi <- summary_abundi %>%
          dplyr::mutate(median = median - bias,
                        Abundance = Abundance - bias)

        summary_abund <- summary_abundi %>%
          # group_by(simulation) %>%
          # dplyr::summarise(Abundance = mean(Abundance, na.rm = T)) %>%
          group_by() %>%
          dplyr::summarise(Low95 = quantile(Abundance, 0.025, na.rm = T),
                           Up95 = quantile(Abundance, 0.975, na.rm = T),
                           SD = sd(Abundance, na.rm = T),
                           Median = median(Abundance, na.rm = T),
                           Abundance = mean(Abundance, na.rm = T),
                           CV = round(SD / Abundance, 3),
                           Model = i) %>%
          dplyr::mutate(Formula = "Mean") %>%
          rbind(summary_abundi %>%
                  # group_by(simulation) %>%
                  # dplyr::summarise(Abundance = mean(Abundance, na.rm = T)) %>%
                  group_by() %>%
                  dplyr::summarise(Low95 = quantile(median, 0.025, na.rm = T),
                                   Up95 = quantile(median, 0.975, na.rm = T),
                                   SD = mad(median, na.rm = T),
                                   Median = median(median, na.rm = T),
                                   Abundance = mean(median, na.rm = T),
                                   CV = round(SD / Median, 3),
                                   Model = i) %>%
                  dplyr::mutate(Formula = "Median")) %>%
          dplyr::mutate(Abundance = round(Abundance, 0),
                        Median = round(Median, 0),
                        Low95 = round(Low95, 0),
                        Up95 = round(Up95, 0),
                        SD = round(SD, 1),
                        n_simulation = n) %>%
          dplyr::rename(mean_abundance = Abundance,
                        median_abundance = Median) %>%
          dplyr::select(Model, mean_abundance, median_abundance, Low95, Up95, CV, SD, Formula, n_simulation)

        summary_abund <- summary_abund %>%
          # dplyr::filter(Formula == "Median") %>%
          dplyr::filter(Formula == "Mean") %>%
          # dplyr::select(-c("Formula", "mean_abundance"))
          dplyr::select(-c("Formula", "median_abundance"))

        # colnames(summary_abund)[colnames(summary_abund) == "median_abundance"] <- "Abundance"
        colnames(summary_abund)[colnames(summary_abund) == "mean_abundance"] <- "Abundance"

        print(kable(summary_abund, caption = "Non-biased final abundance estimates:"))

        ######### STOPPED HERE: need to save the results if save_results_bias_corrected in a proper file

        print(ggplot() +
                geom_histogram(data = #summary_abundi %>%
                                 # dplyr::mutate(Formula = "Mean") %>%
                                 # rbind(
                                 summary_abundi %>%
                                 # dplyr::mutate(Formula = "Median",
                                 dplyr::mutate(Formula = "Mean",
                                               # Abundance = median)#)
                                               Abundance = Abundance)#)
                               , aes(x = Abundance), show.legend = F, fill = "midnightblue") +
                # scale_fill_viridis_d() +
                # scale_y_sqrt() +
                scale_x_continuous(name = "Abundance") +
                # facet_wrap(~ Formula, scales = "free") +
                labs(title = paste0("Histogram of the abundances from the ", n, " simulations")))

        ## Calculation of total abundance based on the mgcv response output.
        final_abund_mgcvi <- final_correctedm %>%
          dplyr::mutate(abundance = Avg_density * areakm2_precropped * groupsize * ratio) %>%
          group_by(id) %>%
          dplyr::summarise(abundance = sum(abundance, na.rm = T)) %>%
          ungroup() %>%
          arrange(id)

        if (!is.null(sub_area_analysis_file)) {
          final_abund_mgcv <- final_abund_mgcvi %>%
            pull(abundance)

          summary_abund_HPi <- with_cp_perHP %>%
            group_by(simulation, AU) %>%
            dplyr::summarise(median = median(Abundance, na.rm = T),
                             Abundance = mean(Abundance, na.rm = T)) %>%
            ungroup()

          mle_perHP <- map_dfr(unique(list_id_perHP$Name[!is.na(list_id_perHP$Name)]), function(n) {
            n_list_id_perHP <- list_id_perHP %>%
              dplyr::filter(Name == n) %>%
              dplyr::mutate(ratio_area_in = ifelse(is.na(ratio_area_in), 0, ratio_area_in)) %>%
              pull(ratio_area_in)

            cell_in <- which(n_list_id_perHP > 0)

            out <- data.frame(Abundance_mgcv = sum(final_abund_mgcv[cell_in] * n_list_id_perHP[cell_in], na.rm = T)) %>%
              dplyr::mutate(AU = n)

            return(out)
          })

          summary_abund_HPi <- map_dfr(unique(summary_abund_HPi$AU), function(au) {

            out <- summary_abund_HPi %>%
              dplyr::filter(AU == au)

            mle_AU <- mle_perHP %>%
              dplyr::filter(AU == au) %>%
              pull(Abundance_mgcv)

            bias <- 2 * (mean(out$Abundance, na.rm = T) - mle_AU)
            # bias <- 2 * (median(out$median, na.rm = T) - mle_AU)

            out <- out %>%
              dplyr::mutate(median = median - bias,
                            Abundance = Abundance - bias)

            return(out)
          })

          summary_abund_HP <- summary_abund_HPi %>%
            group_by(AU) %>%
            dplyr::summarise(Low95 = quantile(Abundance, 0.025, na.rm = T),
                             Up95 = quantile(Abundance, 0.975, na.rm = T),
                             SD = sd(Abundance, na.rm = T),
                             Median = median(Abundance, na.rm = T),
                             Abundance = mean(Abundance, na.rm = T),
                             CV = SD / Abundance,
                             Model = i) %>%
            dplyr::mutate(Formula = "Mean") %>%
            rbind(summary_abund_HPi %>%
                    group_by(AU) %>%
                    dplyr::summarise(Low95 = quantile(median, 0.025, na.rm = T),
                                     Up95 = quantile(median, 0.975, na.rm = T),
                                     SD = mad(median, na.rm = T),
                                     Median = median(median, na.rm = T),
                                     Abundance = mean(median, na.rm = T),
                                     CV = SD / Median,
                                     Model = i) %>%
                    dplyr::mutate(Formula = "Median")) %>%
            dplyr::select(AU, Model, Abundance, Median, Low95, Up95, CV, SD, Formula) %>%
            ungroup() %>%
            dplyr::mutate(Abundance = round(Abundance, 0),
                          Median = round(Median, 0),
                          Low95 = round(Low95, 0),
                          Up95 = round(Up95, 0),
                          SD = round(SD, 1),
                          CV = round(CV, 3),
                          n_simulation = n) %>%
            ungroup() %>%
            dplyr::rename(mean_abundance = Abundance,
                          median_abundance = Median) %>%
            dplyr::select(Model, mean_abundance, median_abundance, Low95, Up95, CV, SD, Formula, n_simulation)

          summary_abund_HP <- summary_abund_HP %>%
            dplyr::filter(Formula == "Mean") %>%
            dplyr::select(-c("Formula", "median_abundance"))

          colnames(summary_abund_HP)[colnames(summary_abund_HP) == "mean_abundance"] <- "Abundance"

          print(kable(summary_abund_HP))

          print(ggplot() +
                  geom_histogram(data = #summary_abund_HPi %>%
                                   # dplyr::mutate(Formula = "Mean") %>%
                                   # rbind(
                                   summary_abund_HPi %>%
                                   dplyr::mutate(Formula = "Mean",
                                                 Abundance = Abundance)#)
                                 , aes(x = Abundance), show.legend = F, fill = "midnightblue") +
                  # scale_fill_viridis_d() +
                  # scale_y_sqrt() +
                  facet_wrap(~ AU, scales = "free") +
                  # facet_grid(AU ~ Formula, scales = "free") +
                  scale_x_continuous(name = "Abundance") +
                  labs(title = paste0("Histogram of the abundances from the ", n, " simulations per Assessment Unit")))
        }

        if (!("groupsize" %in% colnames(per_cell))) { # for the modelled groupsize models where groupsize was directly integrated in the creation of per_cell object
          per_cell$groupsize <- 1
        }

        per_cell <- per_cell %>% ## abundance per cell, per simulation (already averaged over time)
          left_join(final_gp %>%
                      dplyr::rename(cell = id) %>%
                      dplyr::select(cell, ratio),
                    by = "cell") %>%
          dplyr::filter(ratio != 0) %>%
          dplyr::mutate(Median_Cell_abundance = Median_Cell_abundance * ratio * groupsize,
                        mean_Cell_abundance = mean_Cell_abundance * ratio * groupsize)

        bias_percell <- per_cell %>%
          group_by(cell) %>%
          dplyr::summarise(mean_ab = mean(mean_Cell_abundance, na.rm = T)) %>%
          ungroup() %>%
          left_join(final_abund_mgcvi %>%
                      dplyr::rename(cell = id),
                    by = "cell") %>%
          dplyr::mutate(bias = 2 * (mean_ab - abundance))

        per_cell <- per_cell %>%
          left_join(bias_percell, by = "cell") %>%
          dplyr::mutate(Median_Cell_abundance = Median_Cell_abundance - bias,
                        mean_Cell_abundance = mean_Cell_abundance - bias,
                        Median_Cell_abundance = ifelse(Median_Cell_abundance < 0, 1e-16, Median_Cell_abundance),
                        mean_Cell_abundance = ifelse(mean_Cell_abundance < 0, 1e-16, mean_Cell_abundance))

        summary_abund_percell <- per_cell %>%
          group_by(cell) %>%
          dplyr::summarise(SD_med = sd(Median_Cell_abundance, na.rm = T),
                           SD = sd(mean_Cell_abundance, na.rm = T),
                           mad_med = mad(Median_Cell_abundance, na.rm = T),
                           mad = mad(mean_Cell_abundance, na.rm = T),
                           var = SD^2,
                           var_med = SD_med^2,
                           groupsize = unique(groupsize),
                           Abundance_med = mean(Median_Cell_abundance, na.rm = T),
                           Abundance = mean(mean_Cell_abundance, na.rm = T),
                           Median_med = median(Median_Cell_abundance, na.rm = T),
                           Median = median(mean_Cell_abundance, na.rm = T),
                           Max = max(MaxInf_Cell_abundance, na.rm = T),
                           Low95 = quantile(mean_Cell_abundance, 0.025, na.rm = T),
                           Up95 = quantile(mean_Cell_abundance, 0.975, na.rm = T)
          ) %>%
          ungroup() %>%
          dplyr::mutate(CVaa = (SD / Abundance),
                        CVam = (SD_med / Abundance_med),
                        CVma = (mad / Median),
                        CVmm = (mad_med / Median_med),
                        Model = i) %>%
          dplyr::select(cell, SD, SD_med, var, var_med, Model, Abundance, Abundance_med, Median, Median_med, Low95, Up95,
                        groupsize,
                        CVaa, CVam, CVma, CVmm, Max, mad, mad_med) %>%
          dplyr::rename(id = cell)

        cat("\n<br>Sum of the cell abundances:", sum(summary_abund_percell$Abundance, na.rm = T), "inds\n<br>")

        pred_grid_cent <- static %>%
          dplyr::select(id, areakm2)

        bias_percell <- bias_percell %>%
          mutate(id = cell)

        if (save_results_bias_corrected) {
          if (file.exists(paste0(GridDir, "/", version_preds, "_biascorrected_results_Model_", i, ".shp")) & !run_all) {
            cat("Model", i, ": File already exists! Results won't be saved.", paste0("(", GridDir, "/",
                                                                                     version_preds, "_biascorrected_results_Model_", i, ".shp)"),
                "\n")
          } else {
            cat("Saving bias-corrected results from Model", i, "under", paste0(GridDir, "/", version_preds, "_biascorrected_results_Model_", i, ".shp"),
                "\n")

            write_sf(pred_grid_cent %>%
                       left_join(summary_abund_percell %>%
                                   dplyr::select(id, CVma, CVaa),
                                 by = "id") %>%
                       dplyr::filter(!is.na(CVaa)) %>%
                       dplyr::select(id, areakm2, Abundance, CVaa, Low95, Up95) %>%
                       left_join(bias_percell %>%
                                   dplyr::select(id, bias), by = "id") %>%
                       dplyr::rename(CV = CVaa,
                                     Corrected_bias = bias),
                     paste0(GridDir, "/", version_preds, "_biascorrected_results_Model_", i, ".shp"),
                     append = F
                       #ifelse(run_all, FALSE, NA)
                       )
          }
        }

        p1 <- ggplot() +
          geom_sf(data = pred_grid_cent %>%
                    # dplyr::filter(!is.na(Avg_density)) %>%
                    left_join(summary_abund_percell %>%
                                dplyr::mutate(Median = Median_med) %>%
                                # dplyr::mutate(Median = Median_med * groupsize,
                                #               Abundance = Abundance * groupsize) %>%
                                dplyr::select(id, Abundance, Median),
                              by = "id") %>%
                    dplyr::filter(!is.na(Abundance)), aes(fill = cut(Abundance / areakm2,
                                                                     # breaks = c(-1,.4,.8,1.2,1.5,2,2.5,3,3.5,100),
                                                                     breaks = breaks_plot,
                                                                     # labels = c("0.00 - 0.40", "0.41 - 0.80", "0.81 - 1.20",
                                                                     #            "1.21 - 1.50", "1.51 - 2.00", "2.01 - 2.50",
                                                                     #            "2.51 - 3.00", "3.01 - 3.50", "> 3.50")
                                                                     labels = labels_plot
                    )), color = NA) +
          scale_fill_viridis_d(drop = F, name = "Density\n(ind/km2)", #\naveraged between days,\nmedian between simulation)"
          ) +
          labs(title = ifelse(groupsizes,
                              # paste0("Model ", i, ": predictions multiplied by\naverage groupsize per SCANS-IV block"),
                              paste0("Model ", i, ": predictions"),
                              paste0("Model ", i, ": predictions")))

        print(p1)

        p2 <- ggplot() +
          geom_sf(data = pred_grid_cent %>%
                    # dplyr::filter(!is.na(Avg_density)) %>%
                    left_join(summary_abund_percell %>%
                                dplyr::mutate(CVma = CVmm) %>%
                                dplyr::select(id, CVma, CVaa),
                              by = "id") %>%
                    dplyr::filter(!is.na(CVaa)), aes(fill = cut(CVaa,
                                                                # breaks = c(-1, seq(.1, 1, .1), 2, Inf),
                                                                breaks = c(-1, 0.1, 0.25, .5, 1, 1.5, 2, Inf),
                                                                # labels = c("0.00 - 0.10", "0.11 - 0.20", "0.21 - 0.30",
                                                                #            "0.31 - 0.40", "0.41 - 0.50", "0.51 - 0.60",
                                                                #            "0.61 - 0.70", "0.71 - 0.80", "0.81 - 0.90",
                                                                #            "0.91 - 1.00", "1.01 - 2.00", "> 2.00")
                                                                labels = c("0.00 - 0.10", "0.11 - 0.25", "0.26 - 0.50", "0.51 - 1.00",
                                                                           "1.01 - 1.50",
                                                                           "1.51 - 2.00", "> 2.00"))), color = NA) +
          scale_fill_viridis_d(name = "CV", drop = F) +
          labs(title = paste0("Model ", i, ": CV"))

        if (any(summary_abund_percell$Max == Inf)) {
          library(ggnewscale)

          infvalue <- summary_abund_percell %>%
            dplyr::filter(Max == Inf)

          p2 <- p2 +
            new_scale_fill() +
            geom_sf(data = pred_grid_cent %>%
                      dplyr::filter(id %in% infvalue$id) %>%
                      left_join(infvalue %>%
                                  dplyr::select(id, Max),
                                by = "id") %>%
                      dplyr::mutate(Infinite_value = ""), aes(fill = factor(Infinite_value)), color = NA) +
            scale_fill_manual(values = "red", name = paste0(nrow(infvalue), " cells with an\ninfinite value\nin one of\nthe simulations~days"))
        }

        print(p2)

        test <- pred_grid_cent %>%
          left_join(bias_percell, by = "id")

        bias_plot <- ggplot() +
          geom_sf(data = test , aes(fill = bias), color = NA) +
          scale_fill_viridis_c(name = "Bias") +
          labs(title = paste0("Model ", i, ": Bias"))

        print(bias_plot)
      }
    })
  }

  # Write Rmd file
  tmp_md <- tempfile(fileext = ".Rmd")   # temporary file
  writeLines(rmd_text, tmp_md)
  # writeLines(rmd_text, file_name)

  # Render to HTML
  rmarkdown::render(
    # input = file_name,
    input = tmp_md,
    output_file = paste0(output_file, ".html"),
    quiet = TRUE
  )

  message("HTML report generated: ", paste0(output_file, ".html"))
  browseURL(output_file)
  # invisible(output_file)
}
