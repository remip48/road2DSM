#' Extract covariates from nc file
#'
#' On the grid saved in this nc file. SDspace, SDtime and mean can be extracted as well.
#' If SDtime is extracted: it will extract the SDtime of the center, SDspace and mean that are extracted (e.g. if only center and SDspace are extracted, only SDtime for those will be extracted and not for the mean)
#'
#' @param nc.path
#' @param list_variable
#' @param nc_files
#' @param all_pixel.radius
#' @param all_time.period
#' @param dates
#' @param lonmin
#' @param latmin
#' @param lonmax
#' @param latmax
#' @param n_cores
#' @param outfile
#'
#' @return
#' @export
#'
#' @examples
extract_nc <- function (nc.path, list_variable, nc_files, all_pixel.radius,
                        all_time.period, dates, lonmin = -Inf, latmin = -Inf, lonmax = Inf,
                        latmax = Inf, n_cores = NULL, outfile = "log.txt")
{
  if (any(is.na(dates))) {
    stop("There is NA in the 'dates' object")
  }
  order_dim <- function(infos_dimensions, lon, lat, depth,
                        time, use_depth = T) {
    out <- c()
    for (i in infos_dimensions$dim) {
      if (tolower(i) %in% c("lon", "longitude")) {
        out <- c(out, lon)
      }
      else if (tolower(i) %in% c("lat", "latitude")) {
        out <- c(out, lat)
      }
      else if (tolower(i) %in% c("time", "t")) {
        out <- c(out, time)
      }
      else if (tolower(i) %in% c("depth", "d") & use_depth) {
        out <- c(out, depth)
      }
    }
    return(out)
  }
  create_dim <- function(infos_dimensions, lon, lat, depth,
                         time) {
    out <- data.frame(init = NA)[-1]
    for (i in rev(infos_dimensions$dim)) {
      if (tolower(i) %in% c("lon", "longitude")) {
        out <- out %>% group_by(across(all_of(colnames(out)))) %>%
          dplyr::reframe(lon = lon)
      }
      else if (tolower(i) %in% c("lat", "latitude")) {
        out <- out %>% group_by(across(all_of(colnames(out)))) %>%
          dplyr::reframe(lat = lat)
      }
      else if (tolower(i) %in% c("depth", "d")) {
        out <- out %>% group_by(across(all_of(colnames(out)))) %>%
          dplyr::reframe(d = depth)
      }
      else if (tolower(i) %in% c("time", "t")) {
        out <- out %>% group_by(across(all_of(colnames(out)))) %>%
          dplyr::reframe(t = time)
      }
    }
    return(out)
  }
  cl <- makeCluster(ifelse(is.null(n_cores), detectCores() *
                             2/4, n_cores))
  registerDoParallel(cl)
  nc_filesi <- do.call("rbind", foreach(i = 1:nrow(nc_files),
                                        .packages = c("ncdf4", "stringr", "timeDate", "dplyr"),
                                        .noexport = ls()[!(ls() %in% c("nc.path", "nc_files"))]) %dopar%
                         {
                           out <- nc_open(file.path(nc.path, nc_files$file[i]))
                           names <- names(out$var)
                           dim <- names(out$dim)
                           time_var <- dim[str_detect(tolower(dim), "time")]
                           lt <- length(ncvar_get(out, time_var))
                           dt <- length(as.character(timeSequence(from = nc_files$date_start[i],
                                                                  to = nc_files$date_end[i], by = "day")))
                           if (lt != dt) {
                             cat("Number of days in", nc_files$file[i], "entered (",
                                 dt, ") is different from the NC file (", dt,
                                 ").\nFunction continues but stop it if you want to check.\n\n")
                           }
                           return(nc_files[i, ] %>% dplyr::group_by(file, date_start,
                                                                    date_end) %>% dplyr::reframe(variable = names))
                         })
  stopCluster(cl)
  gc()
  list_variable <- nc_filesi %>% dplyr::filter(variable %in%
                                                 list_variable) %>% dplyr::select(file, variable) %>%
    distinct() %>% arrange(file) %>% group_by(variable) %>%
    dplyr::summarise(allfile = paste(file, collapse = "///")) %>%
    ungroup() %>% dplyr::mutate(file_set = paste0("file_set_",
                                                  match(allfile, unique(allfile)))) %>% left_join(nc_filesi %>%
                                                                                                    dplyr::select(file, variable) %>% distinct(), by = "variable",
                                                                                                  relationship = "one-to-many") %>% dplyr::rename(file_id = file) %>%
    dplyr::select(-allfile)
  if (!("expr" %in% colnames(nc_files))) {
    nc_files$expr <- NA
  }
  predtype_ref <- nc_files %>% dplyr::rename(file_id = file) %>%
    left_join(list_variable %>% dplyr::select(file_set, file_id),
              by = "file_id") %>% group_by(file_set) %>% dplyr::reframe(predtype = str_remove_all(str_split_1(unique(type),
                                                                                                              fixed(",")), " "), expr = unique(expr)) %>% as.data.frame()
  ncinfoi_ref <- map_dfr(unique(list_variable$variable), function(i) {
    return(list_variable %>% dplyr::filter(variable == i) %>%
             left_join(nc_files, by = c(file_id = "file")) %>%
             dplyr::rename(nc.name = file_id) %>% dplyr::select(nc.name,
                                                                variable, file_set, date_start, date_end) %>% arrange(date_start) %>%
             dplyr::mutate(period = 1:n() - 1))
  }) %>% distinct()
  datesi <- data.frame(dates = dates)
  for (f in unique(list_variable$file_set)) {
    pred.type <- predtype_ref %>% dplyr::filter(file_set ==
                                                  f) %>% pull(predtype)
    variable <- list_variable %>% dplyr::filter(file_set ==
                                                  f) %>% pull(variable) %>% unique()
    cat("Run NC", f, paste0("(", ifelse(length(variable) ==
                                          1, "variable: ", "variables: "), paste(variable,
                                                                                 collapse = ", "), ") as"), paste(pred.type, collapse = ", "),
        "\n")
    ncinfoi <- ncinfoi_ref %>% dplyr::filter(file_set ==
                                               f)
    dates <- datesi %>% group_by(dates) %>% dplyr::mutate(period = if (any(ncinfoi$date_start <=
                                                                           unique(dates) & ncinfoi$date_end >= unique(dates))) {
      ncinfoi %>% dplyr::filter(date_start <= unique(dates) & date_end >=
                                  unique(dates)) %>% dplyr::pull(period) %>% unique()
    }
    else {
      NA
    }) %>% ungroup()
    infos_dim <- map_dfr(1:nrow(ncinfoi), function(i) {
      v <- nc_open(paste(nc.path, ncinfoi$nc.name[i], sep = "/"))$var
      out <- ncinfoi[i, ] %>% left_join(map_dfr(v, function(x) {
        map_dfr(x$dim, function(d) {
          return(data.frame(variable = x$name, dim = d$name,
                            n_values = length(d$vals)))
        })
      }) %>% dplyr::mutate(nc.name = ncinfoi$nc.name[i]),
      by = c("nc.name", "variable")) %>% ungroup()
      return(out)
    })
    if (!file.exists(paste0(nc.path, "/", f, ".shp"))) {
      datafile <- paste(nc.path, first(ncinfoi$nc.name[ncinfoi$file_set ==
                                                         f]), sep = "/")
      nc.data <- nc_open(datafile)
      namedims <- names(nc.data$dim)
      lat <- ncvar_get(nc.data, namedims[str_detect(namedims,
                                                    "lat")])
      lon <- ncvar_get(nc.data, namedims[str_detect(namedims,
                                                    "lon")])
      ref_c <- first(match(infos_dim %>% dplyr::mutate(id = 1:n()) %>%
                             group_by(nc.name, variable) %>% dplyr::mutate(n_dim = length(unique(dim))) %>%
                             ungroup() %>% dplyr::slice_max(n_dim) %>% dplyr::slice_min(id) %>%
                             pull(variable), ncinfoi$variable))
      data <- map_dfr(match(infos_dim %>% pull(variable) %>%
                              unique(), ncinfoi$variable), function(i) {
                                data.var <- try(ncvar_get(nc.data, ncinfoi$variable[i],
                                                          start = order_dim(infos_dim %>% dplyr::filter(variable ==
                                                                                                          ncinfoi$variable[ref_c] & nc.name == ncinfoi$nc.name[ref_c]),
                                                                            lon = 1, lat = 1, depth = 1, time = 1), count = order_dim(infos_dim %>%
                                                                                                                                        dplyr::filter(variable == ncinfoi$variable[ref_c] &
                                                                                                                                                        nc.name == ncinfoi$nc.name[ref_c]), lon = length(lon),
                                                                                                                                      lat = length(lat), depth = 1, time = 1),
                                                          verbose = FALSE))
                                if (all(class(data.var) == "try-error")) {
                                  cat("Trying to remove depth. If Time should be removed instead of depth, please adapt the function script.\n")
                                  data.var <- ncvar_get(nc.data, ncinfoi$variable[i],
                                                        start = order_dim(infos_dim %>% dplyr::filter(variable ==
                                                                                                        ncinfoi$variable[ref_c] & nc.name == ncinfoi$nc.name[ref_c]),
                                                                          lon = 1, lat = 1, depth = 1, time = 1, use_depth = F), count = order_dim(infos_dim %>%
                                                                                                                                                     dplyr::filter(variable == ncinfoi$variable[ref_c] &
                                                                                                                                                                     nc.name == ncinfoi$nc.name[ref_c]), lon = length(lon),
                                                                                                                                                   lat = length(lat), depth = 1, time = 1, use_depth = F),
                                                        verbose = FALSE)
                                }
                                if (dim(data.var)[1] == length(lon) & dim(data.var)[2] ==
                                    length(lat)) {
                                  if (lat[1] < lat[2]) {
                                    data.var <- data.var[, length(lat):1]
                                  }
                                  data.var <- t(data.var)
                                }
                                newgrid <- raster(data.var, xmn = min(lon) -
                                                    (sort(lon)[2] - sort(lon)[1])/2, xmx = max(lon) +
                                                    (sort(lon)[length(lon)] - sort(lon)[length(lon) -
                                                                                          1])/2, ymn = min(lat) - (sort(lat)[2] - sort(lat)[1])/2,
                                                  ymx = max(lat) + (sort(lat)[length(lat)] -
                                                                      sort(lat)[length(lat) - 1])/2, crs = "+proj=longlat +datum=WGS84") %>%
                                  rasterToPolygons(., dissolve = F) %>% st_as_sf()
                                data <- newgrid %>% st_transform(crs = 4326) %>%
                                  dplyr::mutate(lon_cent = st_coordinates(st_centroid(.))[,
                                                                                          1], lat_cent = st_coordinates(st_centroid(.))[,
                                                                                                                                        2]) %>% group_by(lon_cent) %>% dplyr::mutate(loni = lon[which.min(abs(unique(lon_cent) -
                                                                                                                                                                                                                lon))]) %>% ungroup() %>% group_by(lat_cent) %>%
                                  dplyr::mutate(lati = lat[which.min(abs(unique(lat_cent) -
                                                                           lat))]) %>% ungroup() %>% dplyr::mutate(lon_cent = loni,
                                                                                                                   lat_cent = lati) %>% dplyr::select(-c(layer,
                                                                                                                                                         loni, lati))
                                data <- data %>% dplyr::filter(lon_cent >= lonmin,
                                                               lon_cent <= lonmax, lat_cent >= latmin, lat_cent <=
                                                                 latmax)
                                return(data)
                              })
      data <- data %>% dplyr::filter(!duplicated(paste(lon_cent,
                                                       lat_cent)))
      data_cent <- data.frame(x = data$lon_cent, y = data$lat_cent) %>%
        st_as_sf(coords = c("x", "y"), crs = 4326) %>%
        st_transform(crs = 3035) %>% st_coordinates()
      data <- data %>% dplyr::mutate(X = data_cent[, 1],
                                     Y = data_cent[, 2])
      rm(data_cent)
      cat("Check grid plotted: are lon/lat ok ?\n")
      print(ggplot() + geom_point(data = data, aes(x = lon_cent,
                                                   y = lat_cent)))
      new_grid <- data %>% dplyr::mutate(id = 1:n()) %>%
        st_cast()
      new_grid <- new_grid %>% arrange(id)
      write_sf(new_grid, paste0(nc.path, "/", f, ".shp"),
               append = F)
      cat("Grid from the NC file_set saved under", paste0(nc.path,
                                                          "/", f, ".shp"), "\n")
    }
    else {
      new_grid <- read_sf(paste0(nc.path, "/", f, ".shp"))
    }
    data <- new_grid %>% st_drop_geometry() %>% dplyr::select(lon_cent,
                                                              lat_cent, id)
    Predictor.name_ref <- variable
    gc()
    if (!dir.exists(paste0(nc.path, "/", f))) {
      dir.create(paste0(nc.path, "/", f))
      cat("Extraction will be saved under", paste0(nc.path,
                                                   "/", f), "\n")
    }
    origin_all <- timeDate(ncinfoi %>% dplyr::slice_min(period) %>%
                             pull(date_start) %>% unique())
    origin_all <- strptime(origin_all, "%Y-%m-%d", tz = "GMT")
    all_days_period_ref <- data.frame(date = as.character(timeSequence(from = as.character(origin_all),
                                                                       to = as.character(last(dates$dates)), by = "day"))) %>%
      group_by(date) %>% dplyr::mutate(period = if (any(ncinfoi$date_start <=
                                                        date & ncinfoi$date_end >= date)) {
        ncinfoi %>% dplyr::filter(date_start <= date & date_end >=
                                    date) %>% pull(period) %>% unique()
      }
      else {
        NA
      }) %>% ungroup() %>% dplyr::filter(!is.na(period))
    expr <- predtype_ref %>% dplyr::filter(file_set == f) %>%
      pull(expr) %>% unique()
    if (is.null(n_cores)) {
      n_cores <- detectCores() * 2/4
    }
    cat("Parallel processing with", n_cores, "cores will be used. Check that computed can handle it.\n\n",
        "In case of large files, few iterations might return an error due to memory limits. The function can then be re-run to re-do only these few iterations\n")

    cl <- makeCluster(n_cores, outfile = outfile)
    registerDoParallel(cl)
    try({
      all_run <- foreach(dt = 1:nrow(dates), .packages = c("purrr", "rlang",
                                                           "dplyr", "timeDate", "ncdf4", "stringr", "data.table",
                                                           "collapse"), .noexport = ls()[!(ls() %in% c("data",
                                                                                                       "ncinfoi", "nc.path", "expr", "order_dim", "create_dim",
                                                                                                       "origin_all", "all_days_period_ref", "infos_dim",
                                                                                                       "f", "all_pixel.radius", "pixel.radius", "pred.type",
                                                                                                       "Predictor.name_ref", "all_time.period", "dates"))]) %dopar%
        {
          file_to_save <- paste0(nc.path, "/", f, "/",
                                 dates$dates[dt], "_variables.rds")
          if (file.exists(file_to_save)) {
            return(NULL)
          }
          else {
            print(dt)
            prd <- dates$period[dt]
            origin <- timeDate(ncinfoi %>% dplyr::filter(period ==
                                                           prd) %>% pull(date_start) %>% unique())
            origin <- strptime(origin, "%Y-%m-%d", tz = "GMT")
            SegDay <- length(as.character(timeSequence(from = as.character(origin),
                                                       to = as.character(dates$dates[dt]), by = "day")))
            all_days_period <- all_days_period_ref %>%
              dplyr::filter(date <= dates$dates[dt]) %>%
              dplyr::mutate(id = 1:n()) %>% dplyr::filter(id %in%
                                                            (n() - max(all_time.period) + 1):n())
            ncfile <- unique(ncinfoi$nc.name[ncinfoi$period ==
                                               prd])
            datafile <- paste(nc.path, ncfile, sep = "/")
            if (!file.exists(datafile)) {
              return("file did not exist")
            }
            nc.data <- nc_open(datafile)
            namedims <- names(nc.data$dim)
            lat <- ncvar_get(nc.data, namedims[str_detect(tolower(namedims),
                                                          "lat")])
            lon <- ncvar_get(nc.data, namedims[str_detect(tolower(namedims),
                                                          "lon")])
            if (any(str_detect(tolower(namedims), "depth"))) {
              depth <- ncvar_get(nc.data, namedims[str_detect(tolower(namedims),
                                                              "depth")])
            }
            else {
              depth = 0
            }
            time <- ncvar_get(nc.data, namedims[str_detect(tolower(namedims),
                                                           "time")])
            time <- 0:(length(time) - 1)
            day.index <- which((time + 1) == SegDay)
            if (length(day.index) == 0) {
              print(paste("Not OK! Missing dates in nc file for: ",
                          SegDay, sep = "  "))
            }
            nrows <- length(lon)
            ncols <- length(lat)
            ntimes <- length(time)
            time1 <- max(day.index - max(all_time.period) +
                           1, 1)
            numtimes <- min(max(all_time.period), day.index)
            Predictor.name <- Predictor.name_ref
            data.var_refi <- lapply(Predictor.name, function(pred) {
              data.var_ref <- try(ncvar_get(nc.data, pred,
                                            start = order_dim(infos_dim %>% dplyr::filter(variable ==
                                                                                            pred & nc.name == ncfile), lon = 1,
                                                              lat = 1, depth = 1, time = time1),
                                            count = order_dim(infos_dim %>% dplyr::filter(variable ==
                                                                                            pred & nc.name == ncfile), lon = nrows,
                                                              lat = ncols, depth = length(depth),
                                                              time = numtimes), verbose = FALSE))

              if (all(class(data.var_ref) == "try-error")) {
                # cat("Trying to remove depth. If Time should be removed instead of depth, please adapt the function script.\n")
                data.var_ref <- ncvar_get(nc.data, pred,
                                          start = order_dim(infos_dim %>% dplyr::filter(variable ==
                                                                                          pred & nc.name == ncfile), lon = 1,
                                                            lat = 1, depth = 1, time = time1, use_depth = F),
                                          count = order_dim(infos_dim %>% dplyr::filter(variable ==
                                                                                          pred & nc.name == ncfile), lon = nrows,
                                                            lat = ncols, depth = length(depth),
                                                            time = numtimes, use_depth = F), verbose = FALSE)
              }

              if (any(all_days_period$period < prd)) {
                last_period <- unique(all_days_period$period[all_days_period$period <
                                                               prd])
                ncfile1 <- ncinfoi$nc.name[ncinfoi$period ==
                                             last_period]
                datafile1 <- paste(nc.path, ncfile1,
                                   sep = "/")
                nc.data1 <- nc_open(datafile1)
                time11 <- ncvar_get(nc.data1, "time")
                time11 <- 0:(length(time11) - 1)
                SegDay1 <- length(as.character(timeSequence(from = ncinfoi$date_start[ncinfoi$period ==
                                                                                        last_period], to = last(all_days_period$date[all_days_period$period ==
                                                                                                                                       last_period]), by = "day")))
                day.index1 <- which((time11 + 1) == SegDay1)
                if (length(day.index1) == 0) {
                  print(paste("Not OK! Missing dates in nc file for: ",
                              SegDay1, sep = "  "))
                }
                time11 <- max(day.index1 - (max(all_time.period) -
                                              numtimes) + 1, 1)
                numtimes11 <- min(max(all_time.period) -
                                    numtimes, day.index1)
                data.var1 <- try(ncvar_get(nc.data1, pred,
                                           start = order_dim(infos_dim %>% dplyr::filter(variable ==
                                                                                           pred & nc.name == ncfile1), lon = 1,
                                                             lat = 1, depth = 1, time = time11),
                                           count = order_dim(infos_dim %>% dplyr::filter(variable ==
                                                                                           pred & nc.name == ncfile1), lon = nrows,
                                                             lat = ncols, depth = length(depth),
                                                             time = numtimes11), verbose = FALSE))

                if (all(class(data.var1) == "try-error")) {
                  data.var1 <- ncvar_get(nc.data1, pred,
                                         start = order_dim(infos_dim %>% dplyr::filter(variable ==
                                                                                         pred & nc.name == ncfile1), lon = 1,
                                                           lat = 1, depth = 1, time = time11, use_depth = F),
                                         count = order_dim(infos_dim %>% dplyr::filter(variable ==
                                                                                         pred & nc.name == ncfile1), lon = nrows,
                                                           lat = ncols, depth = length(depth),
                                                           time = numtimes11, use_depth = F), verbose = FALSE)
                }

                data.var_ref <- abind::abind(data.var1,
                                             data.var_ref, along = infos_dim %>%
                                               dplyr::filter(variable == pred &
                                                               nc.name == ncfile & n_values >
                                                               1) %>% dplyr::mutate(id = 1:n()) %>%
                                               dplyr::filter(dim == "time") %>%
                                               pull(id))
                nc_close(nc.data1)
              }
              out <- cbind(data.var_ref) %>% as.data.frame() %>%
                cbind(create_dim(infos_dim %>% dplyr::filter(variable ==
                                                               pred & nc.name == ncfile), lon = lon,
                                 lat = lat, depth = depth, time = 1:max(all_time.period)) %>%
                        as.data.frame()) %>% dplyr::rename(`:=`(!!paste0(pred),
                                                                data.var_ref))
              return(out)
            })
            nc_close(nc.data)
            data.var_ref <- data.var_refi[[1]]
            if (!("d" %in% colnames(data.var_ref))) {
              data.var_ref <- data.var_ref %>% dplyr::mutate(d = depth[1])
            }
            if (length(Predictor.name) > 1) {
              for (i in 2:length(data.var_refi)) {
                if (!("d" %in% colnames(data.var_refi[[i]]))) {
                  data.var_refi[[i]] <- data.var_refi[[i]] %>%
                    dplyr::mutate(d = depth[1])
                }
                data.var_ref <- data.var_ref %>% left_join(data.var_refi[[i]],
                                                           by = c("lon", "lat", "t", "d"))
              }
            }
            rm(data.var_refi)
            data.var_ref <- data.var_ref %>% as.data.frame()
            for (c in colnames(data.var_ref)) {
              data.var_ref[, c] <- c(data.var_ref[, c])
            }
            id_NA <- data.var_ref %>% dplyr::mutate(id = 1:n()) %>%
              dplyr::select(-c(lon, lat, t, d))
            id_NA <- id_NA %>% dplyr::filter(rowSums(is.na(across(colnames(.)))) >=
                                               (ncol(id_NA) - 1)) %>% pull(id)
            if (length(id_NA) > 0) {
              data.var_ref <- data.var_ref[-id_NA, ]
            }
            if (all(!is.na(expr)) & all(!is.null(expr))) {
              for (e in expr) {
                Predictor.name <- c(Predictor.name, str_remove_all(str_split_1(e,
                                                                               fixed("="))[1], " "))
                data.var_ref <- data.var_ref %>% mutate(!!parse_expr(e))
                colnames(data.var_ref)[ncol(data.var_ref)] <- last(Predictor.name)
              }
            }
            if (any(!(unique(paste(floor(data$lat_cent*10^5)/10^5,
                                   floor(data$lon_cent*10^5)/10^5)) %in% unique(paste(floor(data.var_ref$lat*10^5)/10^5,
                                                                                      floor(data.var_ref$lon*10^5)/10^5))))) {
              stop("lon or lat not in the data")
            }
            data.var_ref <- data %>% dplyr::mutate(lon_cent = floor(lon_cent*10^5)/10^5, lat_cent = floor(lat_cent*10^5)/10^5) %>%
              left_join(data.var_ref %>% dplyr::mutate(lon_cent = floor(lon*10^5)/10^5, lat_cent = floor(lat*10^5)/10^5), by = c("lon_cent",
                                                                                                                                 "lat_cent")) %>% dplyr::mutate(id_nc = id) %>%
              dplyr::select(t, d, lon, lat, id_nc, all_of(Predictor.name))
            numtimes <- max(all_time.period)
            data.var_ref_t1 <- data.var_ref %>% dplyr::filter(t ==
                                                                dplyr::first(na.omit(t))) %>% dplyr::select(lon,
                                                                                                            lat, id_nc) %>% as.data.frame()
            run_mean_SDspace <- any(c("mean", "SDspace") %in%
                                      pred.type)
            final <- lapply(all_pixel.radius, function(pixel.radius) {
              if (pixel.radius != all_pixel.radius[1] &
                  !run_mean_SDspace) {
                return(NULL)
              }
              res_lat <- mean(sort(lat)[-1] - sort(lat)[-length(lat)],
                              na.rm = T) * (pixel.radius + 0.5)
              res_lon <- mean(sort(lon)[-1] - sort(lon)[-length(lon)],
                              na.rm = T) * (pixel.radius + 0.5)
              outM <- map_dfr(unique(data$lon_cent),
                              function(l) {
                                data.var_ref_t1l <- data.var_ref_t1[abs(l -
                                                                          data.var_ref_t1$lon) <= res_lon,
                                ]
                                return(data[data$lon_cent == l, ] %>%
                                         group_by(id) %>% group_map(~{
                                           c(lat_cent = .x$lat_cent, lon_cent = .x$lon_cent,
                                             id = .y$id, id_nc = paste(data.var_ref_t1l$id_nc[abs(.x$lat_cent -
                                                                                                    data.var_ref_t1l$lat) <= res_lat],
                                                                       collapse = ","))
                                         }) %>% bind_rows())
                              }) %>% group_by(id) %>% dplyr::reframe(id_nc = str_split_1(id_nc,
                                                                                         ","), lon_cent = lon_cent, lat_cent = lat_cent) %>%
                dplyr::mutate(across(colnames(.), ~as.numeric(.x))) %>%
                left_join(data.var_ref, by = "id_nc",
                          relationship = "many-to-many") %>%
                dplyr::mutate(dist = abs(lat_cent - lat) +
                                abs(lon_cent - lon)) %>% dplyr::select(id,
                                                                       t, id_nc, dist, all_of(Predictor.name))
              if (any(c("SDspace", "mean") %in% pred.type)) {
                setDT(outM)
                if (all(c("SDspace", "mean") %in% pred.type)) {
                  outM <- outM[, c(list(id_nc = id_nc[which.min(dist)]),
                                   setNames(fmean(.SD, na.rm = TRUE),
                                            paste0(Predictor.name, ".mean")),
                                   setNames(fsd(.SD, na.rm = TRUE), paste0(Predictor.name,
                                                                           ".SDspace"))), by = list(id, t), .SDcols = Predictor.name]
                }
                else if ("mean" %in% pred.type) {
                  outM <- outM[, c(list(id_nc = id_nc[which.min(dist)]),
                                   setNames(fmean(.SD, na.rm = TRUE),
                                            paste0(Predictor.name, ".mean"))),
                               by = list(id, t), .SDcols = Predictor.name]
                }
                else if ("SDspace" %in% pred.type) {
                  outM <- outM[, c(list(id_nc = id_nc[which.min(dist)]),
                                   setNames(fsd(.SD, na.rm = TRUE), paste0(Predictor.name,
                                                                           ".SDspace"))), by = list(id, t), .SDcols = Predictor.name]
                }
                outM <- outM %>% as_tibble()
                if (run_mean_SDspace) {
                  colnames(outM)[str_detect(colnames(outM),
                                            fixed(".mean")) | str_detect(colnames(outM),
                                                                         fixed(".SDspace"))] <- paste0(colnames(outM)[str_detect(colnames(outM),
                                                                                                                                 fixed(".mean")) | str_detect(colnames(outM),
                                                                                                                                                              fixed(".SDspace"))], "_", pixel.radius,
                                                                                                       "p")
                }
                outM <- outM %>% left_join(data.var_ref %>%
                                             dplyr::select(id_nc, t, all_of(Predictor.name)),
                                           by = c("t", "id_nc")) %>% dplyr::select(-id_nc)
              }
              setDT(outM)
              cols1 <- colnames(outM)[stringr::str_detect(colnames(outM),
                                                          fixed(".mean")) | stringr::str_detect(colnames(outM),
                                                                                                fixed(".SDspace"))]
              cols2 <- Predictor.name
              outfinal <- map(all_time.period, function(time.period) {
                # print(time.period)
                if (time.period == 1) {
                  outf <- outM[t %in% (max(all_time.period) -
                                         time.period + 1):numtimes, -c("t")]
                  colnames(outf)[colnames(outf) %in%
                                   cols1] <- paste0(cols1, ".mean")
                  colnames(outf)[colnames(outf) %in%
                                   cols2] <- paste0(cols2, ".center")
                }
                else {
                  if (any(c("SDspace", "mean") %in% pred.type) &
                      any(c("center", "SDtime") %in% pred.type)) {
                    # outf <- outM[t %in% (max(all_time.period) - time.period + 1):numtimes,
                    #              c(setNames(lapply(.SD[,..cols1], mean, na.rm = TRUE), paste0(cols1, ".mean")),
                    #                setNames(lapply(.SD[,..cols1], sd, na.rm = TRUE), paste0(cols1, ".SDtime")),
                    #                setNames(lapply(.SD[,..cols2], mean, na.rm = TRUE), paste0(cols2, ".center")),
                    #                setNames(lapply(.SD[,..cols2], sd, na.rm = TRUE), paste0(cols2, ".SDtime"))), by = id]

                    res1 <- outM[
                      t %in% (max(all_time.period) - time.period + 1):numtimes,
                      c(
                        setNames(lapply(.SD, mean, na.rm = TRUE), paste0(cols1, ".mean")),
                        setNames(lapply(.SD, sd,   na.rm = TRUE), paste0(cols1, ".SDtime"))
                      ),
                      by = id,
                      .SDcols = cols1
                    ]

                    res2 <- outM[
                      t %in% (max(all_time.period) - time.period + 1):numtimes,
                      c(
                        setNames(lapply(.SD, mean, na.rm = TRUE), paste0(cols2, ".center")),
                        setNames(lapply(.SD, sd,   na.rm = TRUE), paste0(cols2, ".SDtime"))
                      ),
                      by = id,
                      .SDcols = cols2
                    ]

                    outf <- merge(res1, res2, by = "id")
                    rm(res1)
                    rm(res2)
                  }
                  else if (any(c("SDspace", "mean") %in%
                               pred.type)) {
                    outf <- outM[t %in% (max(all_time.period) -
                                           time.period + 1):numtimes, c(setNames(lapply(.SD, mean, na.rm = TRUE),
                                                                                 paste0(cols1, ".mean")), setNames(lapply(.SD, sd, na.rm = TRUE), paste0(cols1,
                                                                                                                                                         ".SDtime"))), by = id,
                                 .SDcols = cols1]
                  }
                  else if (any(c("center", "SDtime") %in%
                               pred.type)) {
                    outf <- outM[t %in% (max(all_time.period) -
                                           time.period + 1):numtimes, c(setNames(lapply(.SD, mean, na.rm = TRUE),
                                                                                 paste0(cols2, ".center")), setNames(lapply(.SD, sd, na.rm = TRUE), paste0(cols2,
                                                                                                                                                           ".SDtime"))), by = id,
                                 .SDcols = cols2]
                  }
                }
                outf <- outf %>% as_tibble()
                colnames(outf)[2:ncol(outf)] <- paste0(colnames(outf)[2:ncol(outf)],
                                                       "_", time.period, "d")
                return(outf)
              })
              out_final <- data
              for (i in 1:length(outfinal)) {
                if (any(!is.null(outfinal[[i]]))) {
                  out_final <- out_final %>% left_join(outfinal[[i]],
                                                       by = "id")
                }
              }
              return(out_final)
            })
            out <- data
            for (i in 1:length(final)) {
              if (any(!is.null(final[[i]]))) {
                out <- out %>% left_join(final[[i]] %>%
                                           dplyr::select(id, all_of(colnames(final[[i]])[!(colnames(final[[i]]) %in%
                                                                                             colnames(out))])), by = "id")
              }
            }
            saveRDS(out %>% mutate(date = dates$dates[dt]),
                    file_to_save)
            return(NULL)
          }
        }
    })
    try(stopCluster(cl))
    gc()
  }
  # return(NULL)
  invisible("Finished")
}
