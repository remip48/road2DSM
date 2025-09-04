gap_analysis <- function(run_models, # output from run_all_DSM
                         seg_data, # segments used for run_all_DSM. Should not be an sf object.
                         variable,
                         crs = crs,
                         grid_resolution = NULL,
                         version_preds = as.character(lubridate::today()),
                         output_file = paste0(version_preds, "_gap_analysis"), # without extension, will be the name of html file
                         grid_folder,
                         static_grid, # grid with variables (static) that were not included in the extract_grids. Must have geometry, and use the exact same grid that the one used for extract_grid.
                         save_results_dsmextra = paste0(getwd(), "/", version_preds, "_gap_analysis_results"), # put NULL if you dont want to save results
                         study_area = NULL, # in case you want to predict only on a part of your prediction grid
                         filter_year_month_not_in = "0000-00", # year and month that should not be used for prediction
                         run_all = F,
                         outfile = "log.txt",
                         n_cores = NULL) {

  if (is.null(grid_resolution)) {
    stop("Please provide a value for grid_resolution used")
  }

  if (!is.null(save_results_dsmextra)) {
    cat("Results will be saved under", save_results_dsmextra, "during analysis.",
        ifelse(run_all,
               "EXISTING files will be loaded! run_all = FALSE.\n",
               "run EVERYTHING and save it! run_all = TRUE.\n"))
  }

  ## this function is adapted to the type of data of ITAW.
  # library(dsmextra)     # Extrapolation toolkit for ecological models
  # library(raster)       # Geographic data analysis and modeling
  # library(tidyverse)    # Packages for data science
  # library(magrittr)
  # library(colorspace)

  cat("Create ETRScrs regular grid. If other crs is targetted by the grid, please adapt the function script.\n")

  static <- read_sf(paste0(grid_folder, "/", static_grid)) %>%
    st_transform(crs = crs) %>%
    dplyr::mutate(X = st_coordinates(st_centroid(.))[,1],
                  Y = st_coordinates(st_centroid(.))[,2]) %>%
    dplyr::mutate(area = units::drop_units(st_area(.)))

  sfcrs <- sf::st_crs(static)$proj4string

  ####
  grid_ref <- (static %>%
    dplyr::slice_max(area))[1,]

  listX <- seq(grid_ref$X - grid_resolution * ceiling(abs(grid_ref$X - min(static$X)) / grid_resolution),
               grid_ref$X + grid_resolution * ceiling(abs(grid_ref$X - max(static$X)) / grid_resolution),
               grid_resolution)
  listY <- seq(grid_ref$Y - grid_resolution * ceiling(abs(grid_ref$Y - min(static$Y)) / grid_resolution),
               grid_ref$Y + grid_resolution * ceiling(abs(grid_ref$Y - max(static$Y)) / grid_resolution),
               grid_resolution)

  static <- static %>%
    group_by(id) %>%
    dplyr::mutate(X = listX[which.min(abs(unique(X) - listX))],
                  Y = listY[which.min(abs(unique(Y) - listY))],
                  x = X,
                  y = Y) %>%
    ungroup()

  ####
  if (all(!is.null(study_area))) {
    study_area <- study_area %>%
      st_transform(crs = crs) %>%
      group_by() %>%
      dplyr::summarise(do_union = F) %>%
      st_cast("MULTIPOLYGON")

    static <- st_intersection(static, study_area) %>%
      # dplyr::mutate(area = units::drop_units(st_area(.))) %>%
      st_cast()
  }

  study_area <- static %>%
    group_by() %>%
    dplyr::summarise(do_union = T) %>%
    st_cast("MULTIPOLYGON")

  static_sf <- static

  static <- static %>%
    st_drop_geometry()

  segs <- seg_data %>%
    dplyr::select(year, date, all_of(c(variable))) %>%
    drop_na() %>%
    as.data.frame()

  ls <- list.files(grid_folder)
  ls <- ls[ls != static_grid]

  covariates <- variable

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
  run_dsmextra_chunk <- quote({
    cat("\n\nCovariates used for gap analysis: all (", paste(covariates, collapse = ", "), ")<br><br>\n\n")

    if (!file.exists(paste0(save_results_dsmextra, "/", version_preds, "_dsmextra.RData"))) {
      # list_extrapolation <- list()
      # list_nearby <- list()

      # for (l in ls) {
      output <- foreach(l = ls,
                        .packages = c("dplyr", "sf", "raster", "dsmextra"),
                        .noexport = ls()[!(ls() %in% c("grid_folder", "static", "variable", "covariates", "sfcrs"))]
                        ) %dopar% {

                          predgrid <- readRDS(paste0(grid_folder, "/", l)) %>%
                            st_drop_geometry()

                          predgrid <- static %>%
                            dplyr::select(id, all_of(variable[variable %in% colnames(static)]), X, Y) %>%
                            left_join(predgrid,
                                      by = "id") %>%
                            dplyr::select(id, all_of(variable), X, Y) %>%
                            as.data.frame() %>%
                            drop_na()

                          extrapolation <- compute_extrapolation(samples = segs,
                                                                 covariate.names = covariates,
                                                                 prediction.grid = predgrid %>%
                                                                   dplyr::select(x, y, all_of(covariates)),
                                                                 coordinate.system = sfcrs
                          )

                          nearby <- compute_nearby(samples = segs,
                                                   prediction.grid = predgrid %>%
                                                     dplyr::select(x, y, all_of(covariates)),
                                                   coordinate.system = aftt_crs,
                                                   covariate.names = covariates,
                                                   nearby = 1)

                          values_extra <- predgrid %>%
                            dplyr::select(id, y, x) %>%
                            left_join(extrapolation$data$all, by = c("x", "y")) %>%
                            dplyr::mutate(date = str_remove_all(str_split_1(l, "_")[3], fixed(".rds"))) %>%
                            st_drop_geometry()

                           values_nearby <- predgrid %>%
                             dplyr::select(id, y, x) %>%
                             left_join(nearby$raster$perc_nearby %>%
                                         st_as_stars() %>%
                                         st_as_sf() %>%
                                         st_centroid() %>%
                                         dplyr::mutate(x = round(st_coordinates(.)[,1], 0),
                                                       y = round(st_coordinates(.)[,2], 0)), by = c("x", "y")) %>%
                            dplyr::mutate(date = str_remove_all(str_split_1(l, "_")[3], fixed(".rds"))) %>%
                            st_drop_geometry()

                           return(list(extrapolations = values_extra,
                                       nearby = values_nearby))
                        }

      values_extra <- do.call("rbind", map(output, function(x) {x$extrapolations})) %>%
        dplyr::mutate(year = as.numeric(lubridate::year(date)))

      values_nearby <- do.call("rbind", map(output, function(x) {x$nearby})) %>%
        dplyr::mutate(year = as.numeric(lubridate::year(date)))

      save(values_extra, values_nearby,
           file = paste0(save_results_dsmextra, "/", version_preds, "_dsmextra.RData"))
    } else {
      load(file = paste0(save_results_dsmextra, "/", version_preds, "_dsmextra.RData"))
    }

    values_extra_ref <- static_sf %>%
      dplyr::select(id) %>%
      left_join(values_extra, by = "id")

    values_nearby_ref <- static_sf %>%
      dplyr::select(id) %>%
      left_join(values_nearby, by = "id")

    for (yr in sort(unique(values_extra_ref$year))) {

      values_extra <- values_extra_ref %>%
        dplyr::filter(year == yr)

      ExDet <- values_extra_ref %>%
        pull(ExDet)

      print(ggplot() +
              geom_sf(data = countries, fill = "white", alpha = 0.1) +
              geom_sf(data = values_extra %>%
                        dplyr::filter(ExDet < 0),
                      aes(fill = -ExDet, color = -ExDet)) +
              scale_fill_continuous_sequential(palette = "Reds 3", begin = .15, name = "Univariate",
                                               limits = c(min(-1*ExDet[ExDet < 0], na.rm = T),
                                                          max(-1*ExDet[ExDet < 0], na.rm = T))) +
              scale_color_continuous_sequential(palette = "Reds 3", begin = .15, name = "Univariate",
                                                limits = c(min(-1*ExDet[ExDet < 0], na.rm = T),
                                                           max(-1*ExDet[ExDet < 0], na.rm = T))) +
              new_scale_fill() +
              new_scale_color() +
              geom_sf(data = values_extra %>%
                        dplyr::filter(ExDet >= 0 & ExDet <= 1),
                      aes(fill = ExDet, color = ExDet)) +
              scale_fill_continuous_sequential(palette = "Greens 3", begin = .15, name = "Analog",
                                               limits = c(min(ExDet[ExDet >= 0 & ExDet <= 1], na.rm = T),
                                                          max(ExDet[ExDet >= 0 & ExDet <= 1], na.rm = T))) +
              scale_color_continuous_sequential(palette = "Greens 3", begin = .15, name = "Analog",
                                                limits = c(min(ExDet[ExDet >= 0 & ExDet <= 1], na.rm = T),
                                                           max(ExDet[ExDet >= 0 & ExDet <= 1], na.rm = T))) +
              new_scale_fill() +
              new_scale_color() +
              geom_sf(data = values_extra %>%
                        dplyr::filter(ExDet > 1),
                      aes(fill = ExDet, color = ExDet)) +
              scale_fill_continuous_sequential(palette = "Purples 3", begin = .15, name = "Combinatorial",
                                               limits = c(min(ExDet[ExDet > 1], na.rm = T),
                                                          max(ExDet[ExDet > 1], na.rm = T))) +
              scale_color_continuous_sequential(palette = "Purples 3", begin = .15, name = "Combinatorial",
                                                limits = c(min(ExDet[ExDet > 1], na.rm = T),
                                                           max(ExDet[ExDet > 1], na.rm = T))) +
              facet_wrap(~ date, ncol = 10) +
              theme_bw() +
              theme(panel.background = element_rect(fill = "white"),
                    legend.position = "top",
                    legend.direction = "horizontal") +
              labs(title = paste(yr, ":", "Extrapolations and analog predictions")) +
              coord_sf(xlim = c(min(static$X), max(static$X)),
                       ylim = c(min(static$Y), max(static$Y)),
                       expand = F))

    }

    # cat("## Most influential covariates\n<br>")

    groups <- list(c("eke", "mlotst", "nppv"),
                   c("thetao", "gradthetao", "grad2thetao"),
                   c("dist_to_sandeel", "diffT", "gradSST"),
                   c("bathy", "slope", "distance_to_coast"))

    values_extra_ref <- values_extra_ref %>%
      dplyr::filter(mic != 0) %>%
      dplyr::mutate(mic = covariates[mic])

    list_mic <- unique(values_extra_ref$mic)

    for (yr in sort(unique(values_extra_ref$year))) {

      values_mic <- values_extra_ref %>%
        dplyr::filter(year == yr)

      print(ggplot() +
              # geom_sf(data = countries, fill = "white", alpha = 0.1) +
              #     geom_sf(data = values_extra %>%
              #               dplyr::filter(mic != 0),
              #         aes(fill = factor(covariates[mic])), color = NA) +
              #   scale_fill_discrete_qualitative(name = NULL) +
              geom_sf(data = countries, fill = "white", alpha = 0.1) +
              geom_sf(data = values_mic %>%
                        dplyr::filter(mic %in% groups[[1]]) %>%
                        dplyr::mutate(mic = factor(mic, levels = (groups[[1]])[groups[[1]] %in% list_mic])),
                      aes(fill = mic, color = mic)) +
              scale_fill_manual(values = sequential_hcl(
                n        = 4,       # number of discrete categories you need
                palette  = "Greens 3", # base palette name
                c        = 100,     # chroma (colorfulness)
                l        = c(35, 90) # lightness range (lower start = darker start)
              )[-4], name = NULL, drop = F) +
              scale_color_manual(values = sequential_hcl(
                n        = 4,       # number of discrete categories you need
                palette  = "Greens 3", # base palette name
                c        = 100,     # chroma (colorfulness)
                l        = c(35, 90) # lightness range (lower start = darker start)
              )[-4], name = NULL, drop = F) +

              new_scale_fill() +
              new_scale_color() +
              geom_sf(data = values_mic %>%
                        dplyr::filter(mic %in% groups[[2]]) %>%
                        dplyr::mutate(mic = factor(mic, levels = (groups[[2]])[groups[[2]] %in% list_mic])),
                      aes(fill = mic, color = mic)) +
              scale_fill_manual(values = sequential_hcl(
                n        = 4,       # number of discrete categories you need
                palette  = "Reds 3", # base palette name
                c        = 100,     # chroma (colorfulness)
                l        = c(40, 90) # lightness range (lower start = darker start)
              )[-4], name = NULL, drop = F) +
              scale_color_manual(values = sequential_hcl(
                n        = 4,       # number of discrete categories you need
                palette  = "Reds 3", # base palette name
                c        = 100,     # chroma (colorfulness)
                l        = c(40, 90) # lightness range (lower start = darker start)
              )[-4], name = NULL, drop = F) +

              new_scale_fill() +
              new_scale_color() +
              geom_sf(data = values_mic %>%
                        dplyr::filter(mic %in% groups[[3]]) %>%
                        dplyr::mutate(mic = factor(mic, levels = (groups[[3]])[groups[[3]] %in% list_mic])),
                      aes(fill = mic, color = mic)) +
              scale_fill_manual(values = sequential_hcl(
                n        = 4,       # number of discrete categories you need
                palette  = "YlOrBr", # base palette name
                c        = 100,     # chroma (colorfulness)
                l        = c(40, 95) # lightness range (lower start = darker start)
              )[-4], name = NULL, drop = F) +
              scale_color_manual(values = sequential_hcl(
                n        = 4,       # number of discrete categories you need
                palette  = "YlOrBr", # base palette name
                c        = 100,     # chroma (colorfulness)
                l        = c(40, 95) # lightness range (lower start = darker start)
              )[-4], name = NULL, drop = F) +

              new_scale_fill() +
              new_scale_color() +
              geom_sf(data = values_mic %>%
                        dplyr::filter(mic %in% groups[[4]]) %>%
                        dplyr::mutate(mic = factor(mic, levels = (groups[[4]])[groups[[4]] %in% list_mic])),
                      aes(fill = mic, color = mic)) +
              scale_fill_manual(values = sequential_hcl(
                n        = 4,       # number of discrete categories you need
                palette  = "Blues 3", # base palette name
                c        = 100,     # chroma (colorfulness)
                l        = c(16, 95) # lightness range (lower start = darker start)
              )[-4], name = NULL, drop = F) +
              scale_color_manual(values = sequential_hcl(
                n        = 4,       # number of discrete categories you need
                palette  = "Blues 3", # base palette name
                c        = 100,     # chroma (colorfulness)
                l        = c(16, 95) # lightness range (lower start = darker start)
              )[-4], name = NULL, drop = F) +

              facet_wrap(~ date, ncol = 10) +
              theme_bw() +
              theme(panel.background = element_rect(fill = "white"),
                    legend.position = "top",
                    legend.direction = "horizontal") +
              labs(title = paste(yr, ":", "Most influential covariates on extrapolations")) +
              coord_sf(xlim = c(min(static$X), max(static$X)),
                       ylim = c(min(static$Y), max(static$Y)),
                       expand = F))

    }

    # cat("## Nearby data\n<br>")

    for (yr in sort(unique(values_nearby_ref$year))) {

      values_nearby <- values_nearby_ref %>%
        dplyr::filter(year == yr)

      print(ggplot() +
              geom_sf(data = countries, fill = "white", alpha = 0.1) +
              geom_sf(data = values_nearby,
                      aes(fill = perc_nearby, color = perc_nearby)) +
              scale_fill_viridis_c(name = "% nearby",
                                   breaks = seq(0, floor(max(values_nearby_ref$perc_nearby, na.rm = T)), length = 4),
                                   labels = function(breaks) {round(breaks, 0)},
                                   limits = c(0, max(values_nearby_ref$perc_nearby, na.rm = T))) +
              scale_color_viridis_c(name = "% nearby",
                                    breaks = seq(0, floor(max(values_nearby_ref$perc_nearby, na.rm = T)), length = 4),
                                    limits = c(0, max(values_nearby_ref$perc_nearby, na.rm = T)),
                                    guide = "none") +
              facet_wrap(~ date, ncol = 10) +
              theme_bw() +
              theme(panel.background = element_rect(fill = "white"),
                    legend.position = "top",
                    legend.direction = "horizontal") +
              labs(title = paste(yr, ":", "Nearby data")) +
              coord_sf(xlim = c(min(static$X), max(static$X)),
                       ylim = c(min(static$Y), max(static$Y)),
                       expand = F))
    }
  })
  ##############

  rmd_text <- c(
    rmd_text,
    "```{r run, echo=FALSE}",
    paste(deparse(run_preds_chunk), collapse = "\n"),
    "```",
    "``` {r obsn, echo = F, eval=TRUE, out.width = '150%', results='asis', fig.width = 10, fig.height = 15, dpi = 300, fig.align = 'center'}",
    paste(deparse(run_dsmextra_chunk), collapse = "\n"),
    "```"
  )

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
