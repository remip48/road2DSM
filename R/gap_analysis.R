#' Run DSMextra R-package
#'
#' for gap analysis in the environmental space.
#'
#' @param seg_data
#' @param variable
#' @param crs
#' @param grid_resolution
#' @param version_preds
#' @param output_file
#' @param grid_folder
#' @param static_grid
#' @param save_results_dsmextra
#' @param study_area
#' @param filter_year_month_not_in
#' @param run_all
#' @param outfile
#' @param n_cores
#'
#' @return
#' @export
#'
#' @examples
gap_analysis <- function(seg_data, # segments used for run_all_DSM. Should not be an sf object.
                         variable,
                         crs = 3035,
                         grid_resolution = NULL,
                         version_preds = as.character(lubridate::today()),
                         output_file = paste0(version_preds, "_gap_analysis"), # without extension, will be the name of html file
                         grid_folder,
                         static_grid, # grid with variables (static) that were not included in the extract_grids. Must have geometry, and use the exact same grid that the one used for extract_grid.
                         save_results_dsmextra = getwd(), # put NULL if you dont want to save results. Path where to save results.
                         study_area = NULL, # in case you want to predict only on a part of your prediction grid
                         filter_year_month_not_in = "0000-00", # year and month that should not be used for prediction
                         run_all = F,
                         outfile = "log.txt",
                         n_cores = NULL) {

  if (is.null(grid_resolution)) {
    stop("Please provide a value c(res_x, res_y) for grid_resolution used")
  }

  if (!is.null(save_results_dsmextra)) {
    save_results_dsmextra <- gsub("\\\\", "/", save_results_dsmextra)

    if (!dir.exists(save_results_dsmextra)) {
      cat(save_results_dsmextra, "directory does not exists.", getwd(), "is used.\n")
      save_results_dsmextra <- getwd()
    }
    cat("Results will be saved under", paste0(save_results_dsmextra, "/", version_preds, "_dsmextra.RData"), "during analysis.\n",
        ifelse(run_all,
               "The function is running EVERYTHING! run_all = TRUE.\n\n",
               "EXISTING files will be loaded! run_all = FALSE.\n\n"))
  }

  seg_data <- seg_data %>%
    st_drop_geometry()

  ## this function is adapted to the type of data of ITAW.
  # pkg::fun(dsmextra)     # Extrapolation toolkit for ecological models
  # pkg::fun(raster)       # Geographic data analysis and modeling
  # pkg::fun(tidyverse)    # Packages for data science
  # pkg::fun(magrittr)
  # pkg::fun(colorspace)
  # # pkg::fun(patchwork)

  cat("Create regular grid on crs =", crs, ".If other crs is targetted, please modify in function arguments.\n")

  static <- read_sf(paste0(grid_folder, "/", static_grid)) %>%
    st_transform(crs = crs) %>%
    dplyr::mutate(X = st_coordinates(st_centroid(.))[,1],
                  Y = st_coordinates(st_centroid(.))[,2]) %>%
    dplyr::mutate(area = units::drop_units(st_area(.)))

  sfcrs <- sf::st_crs(static)$proj4string

  ####
  grid_ref <- (static %>%
                 dplyr::slice_max(area) #%>%
               # dplyr::mutate(X = round(X, 5),
               #               Y = round(Y, 5))
  )[1,]

  listX <- seq(grid_ref$X - grid_resolution[1] * ceiling(abs(grid_ref$X - min(static$X)) / grid_resolution[1]),
               grid_ref$X + grid_resolution[1] * ceiling(abs(grid_ref$X - max(static$X)) / grid_resolution[1]),
               grid_resolution[1])
  listY <- seq(grid_ref$Y - grid_resolution[2] * ceiling(abs(grid_ref$Y - min(static$Y)) / grid_resolution[2]),
               grid_ref$Y + grid_resolution[2] * ceiling(abs(grid_ref$Y - max(static$Y)) / grid_resolution[2]),
               grid_resolution[2])

  static_sf <- static %>%
    dplyr::select(id)

  static <- static_sf %>%
    left_join(static %>%
                st_drop_geometry() %>%
                group_by(id) %>%
                dplyr::mutate(X = listX[which.min(abs(unique(X) - listX))],
                              Y = listY[which.min(abs(unique(Y) - listY))],
                              x = X,
                              y = Y) %>%
                ungroup(), by = "id")

  ####
  if (all(!is.null(study_area))) {
    study_area <- study_area %>%
      st_as_sf() %>%
      st_transform(crs = crs) %>%
      group_by() %>%
      dplyr::summarise(do_union = F) %>%
      st_cast("MULTIPOLYGON")

    static <- static %>%
      dplyr::filter(c(st_intersects(.,
                                    study_area,
                                    sparse = F))) %>%
      # dplyr::mutate(area = units::drop_units(st_area(.))) %>%
      st_cast()
  }

  study_area <- static %>%
    group_by() %>%
    dplyr::summarise(do_union = T) %>%
    st_cast("MULTIPOLYGON")

  static_sf <- static %>%
    dplyr::select(id)

  static <- static %>%
    st_drop_geometry()

  seg_data <- seg_data %>%
    dplyr::select(year, date, all_of(c(variable))) %>%
    drop_na() %>%
    as.data.frame()

  ls <- list.files(grid_folder)
  ls <- ls[ls != static_grid]
  ls <- do.call("c", lapply(ls, function(l) {
    date <- str_remove_all(dplyr::last(str_split_1(l, "_")), fixed(".rds"))
    if (any(str_detect(date, fixed(filter_year_month_not_in)))) {
      return(NULL)
    } else {
      return(l)
    }
  }))

  if (is.null(n_cores)) {
    n_cores <- detectCores() - 1
    cat("Parallel processing using", n_cores, "cores.\n")
  }

  ##############
  setup_chunk <- quote({

    knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE
    )
    # # # pkg::fun(mgcv)
    # pkg::fun(gratia)
    # pkg::fun(dplyr)
    # pkg::fun(ggplot2)
    # pkg::fun(ggnewscale)
    # pkg::fun(tidyr)
    # pkg::fun(viridis)
    # pkg::fun(knitr)
    # pkg::fun(units)
    # pkg::fun(doParallel)
    # pkg::fun(stringr)
    # pkg::fun(sf)
    # pkg::fun(purrr)
    # # pkg::fun(rnaturalearth)
  })

  ##############
  run_dsmextra_chunk <- quote({
    cat("\n\nvariable used for gap analysis: all", paste0("(", paste(variable, collapse = ", "), ")"), "<br><br>\n\n")

    if (!file.exists(paste0(save_results_dsmextra, "/", version_preds, "_dsmextra.RData")) | run_all) {
      # list_extrapolation <- list()
      # list_nearby <- list()

      cl <- makeCluster(n_cores, outfile = outfile)
      registerDoParallel(cl)

      # for (l in ls) {
      output <- foreach(l = ls,
                        .packages = c("dplyr", "sf", "raster", "dsmextra", "stars", "tidyr", "stringr"),
                        .noexport = ls()[!(ls() %in% c("grid_folder", "static", "variable", "sfcrs", "listX", "listY", "seg_data"))]
      ) %dopar% {
        predgrid <- static %>%
          dplyr::select(id, all_of(variable[variable %in% colnames(static)]), x, y) %>%
          left_join(readRDS(paste0(grid_folder, "/", l)) %>%
                      st_drop_geometry() %>%
                      dplyr::select(id, all_of(colnames(.)[!(colnames(.) %in% c("id", variable[variable %in% colnames(static)], "X", "Y"))])),
                    by = "id") %>%
          dplyr::select(id, x, y, all_of(variable)) %>%
          as.data.frame() %>%
          drop_na()

        extrapolation <- compute_extrapolation(samples = seg_data,
                                               covariate.names = variable,
                                               prediction.grid = predgrid %>%
                                                 dplyr::select(x, y, all_of(variable)),
                                               coordinate.system = sfcrs
        )

        nearby <- compute_nearby(samples = seg_data,
                                 prediction.grid = predgrid %>%
                                   dplyr::select(x, y, all_of(variable)),
                                 coordinate.system = sfcrs,
                                 covariate.names = variable,
                                 nearby = 1)

        values_extra <- predgrid %>%
          dplyr::select(id, y, x) %>%
          left_join(extrapolation$data$all, by = c("x", "y")) %>%
          dplyr::mutate(date = str_remove_all(dplyr::last(str_split_1(l, "_")), fixed(".rds")),
                        year = as.numeric(lubridate::year(unique(date)))) %>%
          st_drop_geometry()

        values_nearby <- predgrid %>%
          dplyr::select(id, y, x) %>%
          left_join(nearby$raster$perc_nearby %>%
                      st_as_stars() %>%
                      st_as_sf() %>%
                      st_centroid() %>%
                      dplyr::mutate(x = st_coordinates(.)[,1],
                                    y = st_coordinates(.)[,2]) %>%
                      group_by(x, y) %>%
                      dplyr::mutate(x = listX[which.min(abs(unique(x) - listX))],
                                    y = listY[which.min(abs(unique(y) - listY))]) %>%
                      ungroup(), by = c("x", "y")) %>%
          dplyr::mutate(date = str_remove_all(dplyr::last(str_split_1(l, "_")), fixed(".rds")),
                        year = as.numeric(lubridate::year(unique(date)))) %>%
          st_drop_geometry()

        return(list(extrapolations = values_extra,
                    nearby = values_nearby))
      }

      stopCluster(cl)
      gc()

      values_extra <- do.call("rbind", lapply(output, function(x) {x$extrapolations}))

      values_nearby <- do.call("rbind", lapply(output, function(x) {x$nearby}))

      if (!is.null(save_results_dsmextra)) {
        save(values_extra, values_nearby,
             file = paste0(save_results_dsmextra, "/", version_preds, "_dsmextra.RData"))
      }
    } else {
      load(file = paste0(save_results_dsmextra, "/", version_preds, "_dsmextra.RData"))
    }

    values_extra_ref <- static_sf %>%
      left_join(values_extra, by = "id")

    values_nearby_ref <- static_sf %>%
      left_join(values_nearby, by = "id")

    ## entire period
    print(ggplot() +
            geom_sf(data = study_area, fill = "white", alpha = 0.1) +
            geom_sf(data = static_sf %>%
                      left_join(values_extra_ref %>%
                                  st_drop_geometry() %>%
                                  group_by(id) %>%
                                  dplyr::summarise(mean_ext = 100*length(which(ExDet < 0 | ExDet > 1)) / n()) %>%
                                  ungroup() %>%
                                  dplyr::filter(mean_ext > 0), by = "id"),
                    aes(fill = mean_ext), color = NA) +
            scale_fill_continuous_sequential(palette = "Reds 3", begin = .15, name = "Percentage of days\nwhere the cell is extrapolated.\n0 are excluded:") +
            theme_bw() +
            theme(panel.background = element_rect(fill = "white"),
                  legend.position = "top",
                  legend.direction = "horizontal") +
            labs(title = paste("Overall percentages of extrapolation")) +
            coord_sf(xlim = c(min(static$X), max(static$X)),
                     ylim = c(min(static$Y), max(static$Y)),
                     expand = F))

    print(ggplot() +
            geom_sf(data = study_area, fill = "white", alpha = 0.1) +
            geom_sf(data = static_sf %>%
                      left_join(values_extra_ref %>%
                                  st_drop_geometry() %>%
                                  dplyr::filter(mic != 0) %>%
                                  group_by(id) %>%
                                  dplyr::summarise(n_mic = length(unique(mic))) %>%
                                  ungroup(),
                                by = "id"),
                    aes(fill = n_mic), color = NA) +
            scale_fill_viridis_c(name = NULL) +
            theme_bw() +
            theme(panel.background = element_rect(fill = "white"),
                  legend.position = "top",
                  legend.direction = "horizontal") +
            labs(title = paste("Number of unique daily 'Most influential variable on extrapolations' over the entire period per cell:")) +
            coord_sf(xlim = c(min(static$X), max(static$X)),
                     ylim = c(min(static$Y), max(static$Y)),
                     expand = F))

    print(ggplot() +
            geom_sf(data = study_area, fill = "white", alpha = 0.1) +
            geom_sf(data = static_sf %>%
                      left_join(values_nearby_ref %>%
                                  st_drop_geometry() %>%
                                  group_by(id) %>%
                                  dplyr::reframe(value = c(mean(perc_nearby, na.rm = T),
                                                           median(perc_nearby, na.rm = T),
                                                           min(perc_nearby, na.rm = T),
                                                           max(perc_nearby, na.rm = T)),
                                                 type = c("Mean", "Median", "Min", "Max")) %>%
                                  dplyr::mutate(type = factor(type, levels = c("Mean", "Median", "Min", "Max"))),
                                by = "id"),
                    aes(fill = value), color = NA) +
            scale_fill_viridis_c(name = "% nearby"
            ) +
            theme_bw() +
            theme(panel.background = element_rect(fill = "white"),
                  legend.position = "top",
                  legend.direction = "horizontal") +
            labs(title = paste("Summary of the Nearby data per cell over the entire period")) +
            coord_sf(xlim = c(min(static$X), max(static$X)),
                     ylim = c(min(static$Y), max(static$Y)),
                     expand = F))
    ## per year/date

    ExDet <- values_extra_ref %>%
      pull(ExDet)

    for (yr in sort(unique(values_extra_ref$year))) {

      values_extra <- values_extra_ref %>%
        dplyr::filter(year == yr)

      print(ggplot() +
              geom_sf(data = study_area, fill = "white", alpha = 0.1) +
              geom_sf(data = values_extra %>%
                        dplyr::filter(ExDet < 0),
                      aes(fill = -ExDet), color = NA) +
              scale_fill_continuous_sequential(palette = "Reds 3", begin = .15, name = "Univariate",
                                               limits = c(min(-1*ExDet[ExDet < 0], na.rm = T),
                                                          max(-1*ExDet[ExDet < 0], na.rm = T))) +
              # scale_color_continuous_sequential(palette = "Reds 3", begin = .15, name = "Univariate",
              #                                   limits = c(min(-1*ExDet[ExDet < 0], na.rm = T),
              #                                              max(-1*ExDet[ExDet < 0], na.rm = T))) +
              new_scale_fill() +
              # new_scale_color() +
              geom_sf(data = values_extra %>%
                        dplyr::filter(ExDet >= 0 & ExDet <= 1),
                      aes(fill = ExDet), color = NA) +
              scale_fill_continuous_sequential(palette = "Greens 3", begin = .15, name = "Analog",
                                               limits = c(min(ExDet[ExDet >= 0 & ExDet <= 1], na.rm = T),
                                                          max(ExDet[ExDet >= 0 & ExDet <= 1], na.rm = T))) +
              # scale_color_continuous_sequential(palette = "Greens 3", begin = .15, name = "Analog",
              #                                   limits = c(min(ExDet[ExDet >= 0 & ExDet <= 1], na.rm = T),
              #                                              max(ExDet[ExDet >= 0 & ExDet <= 1], na.rm = T))) +
              new_scale_fill() +
              # new_scale_color() +
              geom_sf(data = values_extra %>%
                        dplyr::filter(ExDet > 1),
                      aes(fill = ExDet), color = NA) +
              scale_fill_continuous_sequential(palette = "Purples 3", begin = .15, name = "Combinatorial",
                                               limits = c(min(ExDet[ExDet > 1], na.rm = T),
                                                          max(ExDet[ExDet > 1], na.rm = T))) +
              # scale_color_continuous_sequential(palette = "Purples 3", begin = .15, name = "Combinatorial",
              #                                   limits = c(min(ExDet[ExDet > 1], na.rm = T),
              #                                              max(ExDet[ExDet > 1], na.rm = T))) +
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

    # cat("## Most influential variable\n<br>")

    # groups <- list(c("eke", "mlotst", "nppv"),
    #                c("thetao", "gradthetao", "grad2thetao"),
    #                c("dist_to_sandeel", "diffT", "gradSST"),
    #                c("bathy", "slope", "distance_to_coast"))

    values_extra_ref <- values_extra_ref %>%
      dplyr::filter(mic != 0) %>%
      dplyr::mutate(mic = variable[mic])

    list_mic <- unique(values_extra_ref$mic)

    for (yr in sort(unique(values_extra_ref$year))) {

      values_mic <- values_extra_ref %>%
        dplyr::filter(year == yr)

      print(ggplot() +
              geom_sf(data = study_area, fill = "white", alpha = 0.1) +
              geom_sf(data = values_mic %>%
                        dplyr::mutate(mic = factor(mic, levels = variable)),
                      aes(fill = mic, color = mic)) +
              scale_fill_viridis_d() +
              scale_color_viridis_d(direction = -1) +
              facet_wrap(~ date, ncol = 10) +
              theme_bw() +
              theme(panel.background = element_rect(fill = "white"),
                    legend.position = "top",
                    legend.direction = "horizontal") +
              labs(title = paste(yr, ":", "Most influential variable on extrapolations")) +
              coord_sf(xlim = c(min(static$X), max(static$X)),
                       ylim = c(min(static$Y), max(static$Y)),
                       expand = F))

      # print(ggplot() +
      #         geom_sf(data = study_area, fill = "white", alpha = 0.1) +
      #         geom_sf(data = values_mic %>%
      #                   dplyr::filter(mic %in% groups[[1]]) %>%
      #                   dplyr::mutate(mic = factor(mic, levels = (groups[[1]])[groups[[1]] %in% list_mic])),
      #                 aes(fill = mic, color = mic)) +
      #         scale_fill_manual(values = sequential_hcl(
      #           n        = 4,       # number of discrete categories you need
      #           palette  = "Greens 3", # base palette name
      #           c        = 100,     # chroma (colorfulness)
      #           l        = c(35, 90) # lightness range (lower start = darker start)
      #         )[-4], name = NULL, drop = F) +
      #         scale_color_manual(values = sequential_hcl(
      #           n        = 4,       # number of discrete categories you need
      #           palette  = "Greens 3", # base palette name
      #           c        = 100,     # chroma (colorfulness)
      #           l        = c(35, 90) # lightness range (lower start = darker start)
      #         )[-4], name = NULL, drop = F) +
      #
      #         new_scale_fill() +
      #         new_scale_color() +
      #         geom_sf(data = values_mic %>%
      #                   dplyr::filter(mic %in% groups[[2]]) %>%
      #                   dplyr::mutate(mic = factor(mic, levels = (groups[[2]])[groups[[2]] %in% list_mic])),
      #                 aes(fill = mic, color = mic)) +
      #         scale_fill_manual(values = sequential_hcl(
      #           n        = 4,       # number of discrete categories you need
      #           palette  = "Reds 3", # base palette name
      #           c        = 100,     # chroma (colorfulness)
      #           l        = c(40, 90) # lightness range (lower start = darker start)
      #         )[-4], name = NULL, drop = F) +
      #         scale_color_manual(values = sequential_hcl(
      #           n        = 4,       # number of discrete categories you need
      #           palette  = "Reds 3", # base palette name
      #           c        = 100,     # chroma (colorfulness)
      #           l        = c(40, 90) # lightness range (lower start = darker start)
      #         )[-4], name = NULL, drop = F) +
      #
      #         new_scale_fill() +
      #         new_scale_color() +
      #         geom_sf(data = values_mic %>%
      #                   dplyr::filter(mic %in% groups[[3]]) %>%
      #                   dplyr::mutate(mic = factor(mic, levels = (groups[[3]])[groups[[3]] %in% list_mic])),
      #                 aes(fill = mic, color = mic)) +
      #         scale_fill_manual(values = sequential_hcl(
      #           n        = 4,       # number of discrete categories you need
      #           palette  = "YlOrBr", # base palette name
      #           c        = 100,     # chroma (colorfulness)
      #           l        = c(40, 95) # lightness range (lower start = darker start)
      #         )[-4], name = NULL, drop = F) +
      #         scale_color_manual(values = sequential_hcl(
      #           n        = 4,       # number of discrete categories you need
      #           palette  = "YlOrBr", # base palette name
      #           c        = 100,     # chroma (colorfulness)
      #           l        = c(40, 95) # lightness range (lower start = darker start)
      #         )[-4], name = NULL, drop = F) +
      #
      #         new_scale_fill() +
      #         new_scale_color() +
      #         geom_sf(data = values_mic %>%
      #                   dplyr::filter(mic %in% groups[[4]]) %>%
      #                   dplyr::mutate(mic = factor(mic, levels = (groups[[4]])[groups[[4]] %in% list_mic])),
      #                 aes(fill = mic, color = mic)) +
      #         scale_fill_manual(values = sequential_hcl(
      #           n        = 4,       # number of discrete categories you need
      #           palette  = "Blues 3", # base palette name
      #           c        = 100,     # chroma (colorfulness)
      #           l        = c(16, 95) # lightness range (lower start = darker start)
      #         )[-4], name = NULL, drop = F) +
      #         scale_color_manual(values = sequential_hcl(
      #           n        = 4,       # number of discrete categories you need
      #           palette  = "Blues 3", # base palette name
      #           c        = 100,     # chroma (colorfulness)
      #           l        = c(16, 95) # lightness range (lower start = darker start)
      #         )[-4], name = NULL, drop = F) +
      #
      #         facet_wrap(~ date, ncol = 10) +
      #         theme_bw() +
      #         theme(panel.background = element_rect(fill = "white"),
      #               legend.position = "top",
      #               legend.direction = "horizontal") +
      #         labs(title = paste(yr, ":", "Most influential variable on extrapolations")) +
      #         coord_sf(xlim = c(min(static$X), max(static$X)),
      #                  ylim = c(min(static$Y), max(static$Y)),
      #                  expand = F))

    }

    # cat("## Nearby data\n<br>")

    for (yr in sort(unique(values_nearby_ref$year))) {

      values_nearby <- values_nearby_ref %>%
        dplyr::filter(year == yr)

      print(ggplot() +
              geom_sf(data = study_area, fill = "white", alpha = 0.1) +
              geom_sf(data = values_nearby,
                      aes(fill = perc_nearby), color = NA) +
              scale_fill_viridis_c(name = "% nearby",
                                   breaks = seq(0, floor(max(values_nearby_ref$perc_nearby, na.rm = T)*10)/10, length = 5),
                                   labels = function(breaks) {round(breaks, 0)},
                                   limits = c(0, max(values_nearby_ref$perc_nearby, na.rm = T))) +
              # scale_color_viridis_c(name = "% nearby",
              #                       breaks = seq(0, floor(max(values_nearby_ref$perc_nearby, na.rm = T)*10)/10, length = 5),
              #                       limits = c(0, max(values_nearby_ref$perc_nearby, na.rm = T)),
              #                       guide = "none") +
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
    "---",
    paste0("title: \"", "Gap analysis in the environmental space, based on dsmextra R-package:", "\""),
    paste0("date: '`r Sys.Date()`'"),
    "output: html_document",
    "---",
    ""
  )

  rmd_text <- c(
    rmd_text,
    "```{r setup, include=FALSE}",
    paste(deparse(setup_chunk), collapse = "\n"),
    "```",
    # "aaaa"
    "``` {r obsn, echo = F, eval=TRUE, out.width = '200%', results='asis', fig.width = 10, fig.height = 15, dpi = 100, fig.align = 'center'}",
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
    output_file = paste0(save_results_dsmextra, "/", output_file, ".html"),
    quiet = TRUE
  )

  message("HTML report generated: ", paste0(save_results_dsmextra, "/", output_file, ".html"))
  browseURL(paste0(save_results_dsmextra, "/", output_file, ".html"))
  invisible(paste0(save_results_dsmextra, "/", output_file, ".html"))
}
