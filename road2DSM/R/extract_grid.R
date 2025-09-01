#' Title
#'
#' @param grid
#' @param variable
#' @param file_set_directory
#' @param writing_directory
#' @param dates
#' @param n_cores
#' @param outfile
#'
#' @return
#' @export
#'
#' @examples
extract_grid <- function(grid,
                         variable,
                         file_set_directory,
                         writing_directory,
                         dates,
                         n_cores = NULL,
                         outfile = "log.txt") {
  # library(sf)
  # library(dplyr)
  # library(raster)
  # library(stars)
  # library(terra)
  # library(lubridate)
  # library(purrr)
  # library(doParallel)
  # library(stringr)
  # library(timeDate)

  list_set <- list.files(file_set_directory)
  list_set <- list_set[str_detect(list_set, "file_set_")]

  if (any(is.na(dates))) {
    stop("There is NA in the object 'dates'.")
  }

  list_date <- sort(unique(dates))

  if (is.null(n_cores) | is.na(n_cores)) {
    n_cores <- detectCores() * 3 / 4
    cat("Parallel processing with", n_cores, "cores will be used.\n")
  }

  cl <- makeCluster(n_cores, outfile = outfile
  )
  registerDoParallel(cl)

  extract <- foreach(d = sort(list_date),
                     .noexport = ls()[!(ls() %in% c("grid", "variable", "writing_directory",
                                                    "file_set_directory"))],
                     .packages = c("sf", "dplyr", "stringr")
  ) %dopar% {

    out <- grid %>%
      st_simplify()

    out_grid <- lapply(unique(list_variable$file_set), function(f) {
      print(f)

      phy <- read_sf(paste0(file_set_directory, "/", f, ".shp")) %>%
        st_transform(crs = 3035) %>%
        st_simplify() %>%
        dplyr::mutate(lon_cent = round(lon_cent, 6),
                      lat_cent = round(lat_cent, 6))

      list_files <- list.files(paste(file_set_directory, f, sep = "/"))[str_detect(list.files(paste(file_set_directory, f, sep = "/")), fixed(d))]

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
        }
      }

      phycc <- phy %>%
        st_coordinates(.) %>%
        as.data.frame() %>%
        group_by(L2) %>%
        dplyr::summarise(Xmin = min(X), Xmax = max(X),
                         Ymin = min(Y), Ymax = max(Y)) %>%
        ungroup()

      grid <- phy %>%
        dplyr::mutate(file_set = f,
                      Xmin = phycc$Xmin,
                      Xmax = phycc$Xmax,
                      Ymin = phycc$Ymin,
                      Ymax = phycc$Ymax) %>%
        dplyr::filter(Xmax >= st_bbox(out)[1] & Xmin <= st_bbox(out)[3] &
                        Ymax >= st_bbox(out)[2] & Ymin <= st_bbox(out)[4]) %>%
        # st_simplify() %>%
        st_cast()

      return(grid)

    })

    i <- 1
    while (i <= (length(out_grid) - 1)) {
      j <- i + 1
      while (j <= length(out_grid)) {
        if (all(out_grid[[i]]$lon_cent == out_grid[[j]]$lon_cent) &
            all(out_grid[[i]]$lat_cent == out_grid[[j]]$lat_cent)) {

          out_grid[[i]] <- out_grid[[i]] %>%
            left_join(out_grid[[j]] %>%
                        st_drop_geometry() %>%
                        dplyr::select(lon_cent, lat_cent, all_of(colnames(.)[!(colnames(.) %in% colnames(out_grid[[i]]))])),
                      by = c("lon_cent", "lat_cent")) %>%
            st_cast()

          out_grid <- out_grid[-j]
          # cat("Joined ")
        } else {
          j <- j + 1
        }
        # cat(i, j, "\n")
      }
      i <- i + 1
    }

    all_cols <- do.call("c", lapply(out_grid, function(x) {
      do.call("c", lapply(colnames(x), function(c) {
        if (any(str_detect(c, variable))) {c} else {NULL}
      }))
    }))

    bb <- st_bbox(out)

    final <- lapply(out_grid, function(grid) {

      cols_grid <- all_cols[all_cols %in% colnames(grid)]

      return(grid[grid$Xmax >= bb[1] & grid$Xmin <= bb[3] &
                    grid$Ymax >= bb[2] & grid$Ymin <= bb[4], ] %>%
               dplyr::select(all_of(cols_grid)) %>%
               st_intersection(., out) %>%
               dplyr::mutate(area_km2 = as.numeric(st_area(.))) %>%
               st_drop_geometry() %>%
               dplyr::select(id, area_km2, all_of(cols_grid)) %>%
               group_by(id) %>%
               dplyr::summarise(across(all_of(cols_grid),
                                       ~ (sum(.x * area_km2, na.rm = T) /
                                            sum(area_km2[!is.na(.x)], na.rm = T)))
               ) %>%
               ungroup()
      )
    })

    for (i in 1:length(final)) {
      out <- out %>%
        left_join(final[[i]], by = "id") %>%
        st_cast()
    }

    saveRDS(out %>%
              st_drop_geometry() %>%
              dplyr::mutate(date = d),
            paste0(writing_directory, "/02_prediction_grids/20250826_predgrid_", d, ".rds"))

    return(NULL)
  }

  stopCluster(cl)
  gc()

  return(NULL)

}
