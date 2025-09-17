#' Extract covariates on effort
#'
#' Based on grids extracted from nc files using extract_nc. You can also use this function to extract covariates from grids created by SDspace_on_sf by reproducing the output structure of extract_nc: the sf grids must be called as "file_set_X.shp" (where X is a number) in file_set_directory. The daily .rds grids with actual covariates to extract must be saved under folder called "file_set_X" where X is matching with the corresponding sf grid.
#'
#' @param effort
#' @param variable
#' @param file_set_directory
#' @param n_cores
#' @param outfile
#'
#' @return
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @export
#'
#' @examples
extract_effort <- function(effort,
                           variable,
                           rasters = list(bathy = NULL), # should be terra object
                           distance_to = list(distance_to_coast = NULL), # should be sf object
                           file_set_directory,
                           try_to_combine_filetsets = T,
                           n_cores = NULL,
                           outfile = "log.txt") {
  # pkg::fun(sf)
  # pkg::fun(dplyr)
  # pkg::fun(raster)
  # pkg::fun(stars)
  # pkg::fun(terra)
  # pkg::fun(lubridate)
  # pkg::fun(purrr)
  # pkg::fun(doParallel)
  # pkg::fun(stringr)

  if (try_to_combine_filetsets) {
    cat("Will try to combine the file_set together.\n")
  } else {
    cat("File_set won't be combined together before extraction. Might be longer if grids are actually similar.\n")
  }

  list_set <- list.files(file_set_directory)
  list_set <- list_set[str_detect(list_set, "file_set_")]
  list_set <- list_set[!is.na(as.numeric(str_remove_all(list_set, fixed("file_set_"))))]

  effort$date <- as.character(effort$date)

  if (any(is.na(effort$date))) {
    stop("There is NA in the 'dates' object")
  }

  list_date <- sort(unique(effort %>%
                             pull(date)))

  if (is.null(n_cores)) {
    n_cores <- detectCores() * 3 / 4
    cat("Parallel processing with", n_cores, "cores will be used.\n")
  }

  cl <- makeCluster(n_cores, outfile = outfile
  )
  registerDoParallel(cl)

  extract <- foreach(d = sort(list_date),
                     .noexport = ls()[!(ls() %in% c("effort", "variable", "list_set",
                                                    "file_set_directory", "try_to_combine_filetsets"))],
                     .packages = c("sf", "dplyr", "stringr")
  ) %dopar% {
    cat(d, "\n")

    out <- effort %>%
      filter(date == d) %>%
      st_simplify()

    out_grid <- lapply(list_set, function(f) {
      # cat(d, f, "\n")

      phy <- read_sf(paste0(file_set_directory, "/", f, ".shp")) %>%
        st_transform(crs = 3035) %>%
        st_simplify() %>%
        dplyr::mutate(lon_cent = round(lon_cent, 6),
                      lat_cent = round(lat_cent, 6))

      list_files <- list.files(paste(file_set_directory, f, sep = "/"))[str_detect(list.files(paste(file_set_directory, f, sep = "/")), fixed(d))]

      if (length(list_files) > 0) {
        for (lf in list_files) {
          temp <- readRDS(paste(file_set_directory, f, lf, sep = "/")) %>%
            dplyr::select(id, all_of(colnames(.)[!(colnames(.) %in% colnames(phy))]))

          phy <- phy %>%
            left_join(temp, by = "id") %>%
            st_cast()
        }

        all_cols <- do.call("c", lapply(colnames(phy), function(c) {
          if (any(str_detect(c, variable))) {c} else {NULL}
        }))

        for (c in all_cols[str_detect(all_cols, fixed(".center"))]) {
          if (str_replace_all(c, "center", "mean") %in% all_cols & any(is.na(phy %>% pull(get(c))))) {
            NAl <- which(is.na(phy %>% pull(get(c))))

            phy <- phy %>%
              dplyr::filter(is.na(get(c))) %>%
              dplyr::mutate(!!c := (phy[NAl, ] %>% pull(get(str_replace_all(c, "center", "mean"))))) %>%
              rbind(phy %>%
                      dplyr::filter(!is.na(get(c)))) %>%
              arrange(id)
          } else if (str_replace_all(c, "center", "mean_1p.mean") %in% all_cols & any(is.na(phy %>% pull(get(c))))) {
            NAl <- which(is.na(phy %>% pull(get(c))))

            phy <- phy %>%
              dplyr::filter(is.na(get(c))) %>%
              dplyr::mutate(!!c := (phy[NAl, ] %>% pull(get(str_replace_all(c, "center", "mean_1p.mean"))))) %>%
              rbind(phy %>%
                      dplyr::filter(!is.na(get(c)))) %>%
              arrange(id)
          }
        }

        return(phy %>%
                 # dplyr::mutate(file_set = f) %>%
                 # st_simplify() %>%
                 st_cast())
      } else {
        return(NULL)
      }
    })

    i <- 1
    while (i <= length(out_grid)) {
      if (is.null(nrow(out_grid[[i]]))) {
        out_grid <- out_grid[-i]
      } else {
        i <- i + 1
      }
    }

    if (length(out_grid) > 0) {

      if (try_to_combine_filetsets) {
        i <- 1
        while (i <= (length(out_grid) - 1)) {
          j <- i + 1
          while (j <= length(out_grid)) {
            if (all(unique(paste(floor(out_grid[[j]]$lat_cent*10^5)/10^5,
                                 floor(out_grid[[j]]$lon_cent*10^5)/10^5)) %in% unique(paste(floor(out_grid[[i]]$lat_cent*10^5)/10^5,
                                                                                             floor(out_grid[[i]]$lon_cent*10^5)/10^5)))) {
              # if (all(out_grid[[i]]$lon_cent == out_grid[[j]]$lon_cent) &
              #     all(out_grid[[i]]$lat_cent == out_grid[[j]]$lat_cent)) {

              out_grid[[i]] <- out_grid[[i]] %>%
                dplyr::mutate(flon = floor(lon_cent*10^5), flat = floor(lat_cent*10^5)) %>%
                left_join(out_grid[[j]] %>%
                            st_drop_geometry() %>%
                            dplyr::mutate(flon = floor(lon_cent*10^5), flat = floor(lat_cent*10^5)) %>%
                            dplyr::select(flon, flat, all_of(colnames(.)[!(colnames(.) %in% colnames(out_grid[[i]]))])),
                          by = c("flon", "flat")) %>%
                dplyr::select(-c(flon, flat)) %>%
                st_cast()

              out_grid <- out_grid[-j]
              # cat("Joined ")
              # } else if (any(unique(paste(floor(out_grid[[i]]$lat_cent*10^5)/10^5,
              #                             floor(out_grid[[i]]$lon_cent*10^5)/10^5)) %in% unique(paste(floor(out_grid[[i]]$lat_cent*10^5)/10^5,
              #                                                                                         floor(out_grid[[i]]$lon_cent*10^5)/10^5)))) {
              #   out_grid[[i]] <- out_grid[[i]] %>%
              #     dplyr::mutate(flon = floor(lon_cent*10^5), flat = floor(lat_cent*10^5)) %>%
              #     left_join(out_grid[[j]] %>%
              #                 st_drop_geometry() %>%
              #                 dplyr::mutate(flon = floor(lon_cent*10^5), flat = floor(lat_cent*10^5)) %>%
              #                 dplyr::select(flon, flat, all_of(colnames(.)[!(colnames(.) %in% colnames(out_grid[[i]]))])),
              #               by = c("flon", "flat")) %>%
              #     dplyr::select(-c(flon, flat)) %>%
              #     st_cast()
              #
              #   out_grid[[j]] <- out_grid[[j]] %>%
              #     dplyr::filter(!(paste(floor(lat_cent*10^5)/10^5,
              #                           floor(lon_cent*10^5)/10^5) %in% paste(floor(out_grid[[i]]$lat_cent*10^5)/10^5,
              #                                                                 floor(out_grid[[i]]$lon_cent*10^5)/10^5)))
              #
              #   j <- j + 1
            } else {
              j <- j + 1
            }
            # cat(i, j, "\n")
          }
          i <- i + 1
        }
      }

      bb <- st_bbox(out)

      out_grid <- lapply(out_grid, function(grid) {
        gridcc <- grid %>%
          st_coordinates(.) %>%
          as.data.frame() %>%
          group_by(L2) %>%
          dplyr::summarise(Xmin = min(X), Xmax = max(X),
                           Ymin = min(Y), Ymax = max(Y)) %>%
          ungroup()

        return(grid %>%
                 dplyr::mutate(Xmin = gridcc$Xmin,
                               Xmax = gridcc$Xmax,
                               Ymin = gridcc$Ymin,
                               Ymax = gridcc$Ymax) %>%
                 dplyr::filter(Xmax >= bb[1] & Xmin <= bb[3] &
                                 Ymax >= bb[2] & Ymin <= bb[4]))
      })

      all_cols <- do.call("c", lapply(out_grid, function(x) {
        do.call("c", lapply(colnames(x), function(c) {
          if (any(str_detect(c, variable))) {c} else {NULL}
        }))
      }))

      final <- out %>%
        group_by(label) %>%
        group_map(~ {

          c(do.call("c", lapply(out_grid, function(gridX) {
            ints <- gridX %>%
              st_intersection(., .x) %>%
              dplyr::mutate(area_km2 = as.numeric(st_area(.))) %>%
              st_drop_geometry()

            vals <- ints[, all_cols[all_cols %in% colnames(gridX)]]

            c(colSums(vals * ints$area_km2, na.rm = TRUE) / colSums(ints$area_km2 * !is.na(vals)))
          })),
          label = .y$label)

        }) %>%
        bind_rows() %>%
        dplyr::mutate(across(all_of(all_cols), ~ as.numeric(.x)))

      out <- out %>%
        left_join(final, by = "label") %>%
        st_cast() %>%
        st_drop_geometry()

      return(out)
    } else {
      return(NULL)
    }
  }

  stopCluster(cl)
  gc()

  extract <- do.call("rbind", extract)

  effort <- effort %>%
    left_join(extract %>%
                dplyr::select(label, all_of(colnames(extract)[!(colnames(extract) %in% colnames(effort))])), by = "label") %>%
    st_cast()

  crseffort <- st_crs(effort)

  for (i in 1:length(rasters)) {
    if (any(!is.null(rasters[[i]]))) {
      effort <- effort %>%
        st_transform(crs = terra::crs(rasters[[i]])) %>%
        dplyr::mutate(!!names(rasters)[i] := terra::extract(rasters[[i]], ., fun = mean, na.rm = TRUE)[, 2])
    }
  }

  for (i in 1:length(distance_to)) {
    if (any(!is.null(distance_to[[i]]))) {
      effort <- effort %>%
        st_transform(crs = crseffort) %>%
        dplyr::mutate(!!names(distance_to)[i] := c(units::drop_units(st_distance(st_centroid(.),
                                                                                 distance_to[[i]] %>%
                                                                                   st_transform(crs = crseffort) %>%
                                                                                   group_by() %>%
                                                                                   dplyr::summarise(do_union = F) %>%
                                                                                   st_cast()))))
    }
  }

  return(effort %>%
           st_cast())
}
