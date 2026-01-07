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
#' @param Number_starting_name_file_set
#' @param resolution
#' @param name_dimension
#' @param vertical_variables
#' @param max_depth
#' @param rename_dimensions
#' @param crs_meter
#' @param SDspace_radius km or pixel, depending on the unit of values you put in all_pixel.radius (if SDspace is used)
#'
#' @return
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom data.table .SD .N :=
#' @import data.table
#'
#' @export
#'
#' @examples


extract_nc <- function (nc.path, list_variable, nc_files, all_pixel.radius,
                        resolution = "day", all_time.period, dates, lonmin = -Inf,
                        latmin = -Inf, lonmax = Inf, latmax = Inf, name_dimension = list(lon = "lon",
                                                                                         lat = "lat", time = NULL, depth = NULL), max_depth = Inf,
                        vertical_variables = NULL, Number_starting_name_file_set = 1,
                        crs_meter = 3035, SDspace_radius = "pixel",
                        rename_dimensions = list(lon = NULL, lat = NULL, time = NULL, depth = NULL),
                        n_cores = NULL, outfile = "log.txt")
{
  if (any(is.na(dates))) {
    stop("There is NA in the 'dates' object")
  }
  order_dim <- function(infos_dimensions, dimensions, use_1value = T) {
    out <- c()
    for (i in infos_dimensions$dim) {
      if (use_1value | infos_dimensions$n_values[infos_dimensions$dim ==
                                                 i] > 1) {
        out <- c(out, dimensions[[i]])
      }
    }
    return(out)
  }
  create_dim <- function(infos_dimensions, dimensions, name_dim) {
    out <- data.frame(init = NA)[-1]
    for (i in rev(infos_dimensions$dim)) {
      out <- out %>% group_by(across(all_of(colnames(out)))) %>%
        dplyr::reframe(`:=`(!!i, dimensions[[i]]))
      if (i == name_dim[["lon"]]) {
        colnames(out)[colnames(out) == i] <- "lon"
      }
      else if (i == name_dim[["lat"]]) {
        colnames(out)[colnames(out) == i] <- "lat"
      }
      else if (i == name_dim[["time"]]) {
        colnames(out)[colnames(out) == i] <- "t"
      }
      else if (i == name_dim[["depth"]]) {
        colnames(out)[colnames(out) == i] <- "d"
      }
    }
    return(out)
  }
  if (any(as.character(dates) < min(as.character(nc_files$date_start),
                                    na.rm = T))) {
    stop("At least one date to extract is prior to one NC file's starting date. The function will fail.")
  }
  if (any(as.character(lubridate::as_date(min(dates)) - ifelse(resolution ==
                                                               "day", lubridate::days(max(all_time.period - 1)), ifelse(resolution ==
                                                                                                                        "month", months(max(all_time.period - 1)), NA))) < min(as.character(nc_files$date_start),
                                                                                                                                                                               na.rm = T))) {
    warning("Are you sure the periods calculated for SDtime are inside NCDF file timeframes for every dates to extract ?")
  }
  if (resolution == "month" & any(as.numeric(lubridate::day(dates)) !=
                                  15)) {
    stop("Please label all dates as the 15th of the month for resolution = month")
  }
  if (resolution == "month" & any(as.numeric(lubridate::day(nc_files$date_start)) !=
                                  15)) {
    stop("Please label all nc_files$date_start as the 15th of the month for resolution = month")
  }
  if (resolution == "month" & any(as.numeric(lubridate::day(nc_files$date_end)) !=
                                  15)) {
    stop("Please label all nc_files$date_end as the 15th of the month for resolution = month")
  }
  cl <- makeCluster(ifelse(is.null(n_cores), detectCores() *
                             2/4, n_cores))
  registerDoParallel(cl)
  nc_filesi <- do.call("rbind", foreach(i = 1:nrow(nc_files),
                                        .packages = c("ncdf4", "stringr", "timeDate", "dplyr"),
                                        .noexport = ls()[!(ls() %in% c("nc.path", "nc_files",
                                                                       "resolution", "name_dimension"))]) %dopar% {
                                                                         out <- nc_open(file.path(nc.path, nc_files$file[i]))
                                                                         names <- names(out$var)
                                                                         dim <- names(out$dim)
                                                                         time_var <- dim[dim == name_dimension[["time"]]]
                                                                         if (length(time_var) == 1) {
                                                                           lt <- length(ncvar_get(out, time_var))
                                                                         }
                                                                         else {
                                                                           time_var <- names[names == name_dimension[["time"]]]
                                                                           if (length(time_var) == 1) {
                                                                             lt <- length(ncvar_get(out, time_var))
                                                                           }
                                                                           else {
                                                                             lt <- 1
                                                                           }
                                                                         }
                                                                         dt <- length(as.character(timeSequence(from = nc_files$date_start[i],
                                                                                                                to = nc_files$date_end[i], by = resolution)))
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
                                                  match(allfile, unique(allfile)) + Number_starting_name_file_set -
                                                    1)) %>% left_join(nc_filesi %>% dplyr::select(file,
                                                                                                  variable) %>% distinct(), by = "variable", relationship = "one-to-many") %>%
    dplyr::rename(file_id = file) %>% dplyr::select(-allfile)
  if (file.exists(paste0(nc.path, "/list_file_set.xlsx"))) {
    list_sets <- readxl::read_xlsx(paste0(nc.path, "/list_file_set.xlsx"))
    list_variablei <- list_variable
    list_variable <- list_variable %>% dplyr::select(-file_set) %>%
      dplyr::left_join(list_sets, by = "file_id")
    if (any(!(list_variable$file_id %in% list_sets$file_id))) {
      for (i in 1:nrow(list_variable)) {
        if (is.na(list_variable$file_set[i])) {
          file_ref <- list_variablei %>% dplyr::filter(file_set ==
                                                         list_variablei$file_set[i] & file_id !=
                                                         list_variable$file_id[i]) %>% pull(file_id) %>%
            unique()
          rep <- list_variable %>% dplyr::filter(file_id %in%
                                                   file_ref & !is.na(file_set)) %>% pull(file_set) %>%
            unique()
          if (length(rep) == 0) {
            list_variable$file_set[i] <- paste0("file_set_",
                                                (list_variable %>% pull(file_set) %>%
                                                   str_remove_all(., fixed("file_set_")) %>%
                                                   as.numeric() %>% max(., na.rm = T)) +
                                                  1)
          }
          else {
            list_variable$file_set[i] <- rep
          }
        }
      }
      writexl::write_xlsx(list_variable %>% dplyr::select(file_set,
                                                          file_id) %>% rbind(list_sets %>% dplyr::select(file_set,
                                                                                                         file_id)) %>% dplyr::arrange(file_set, file_id) %>%
                            dplyr::distinct(), paste0(nc.path, "/list_file_set.xlsx"))
    }
  }
  else {
    writexl::write_xlsx(list_variable %>% dplyr::select(file_set,
                                                        file_id) %>% dplyr::arrange(file_set, file_id) %>%
                          dplyr::distinct(), paste0(nc.path, "/list_file_set.xlsx"))
  }
  if (!("expr" %in% colnames(nc_files))) {
    nc_files$expr <- NA
  }
  predtype_ref <- nc_files %>% dplyr::rename(file_id = file) %>%
    left_join(list_variable %>% dplyr::select(file_set,
                                              file_id), by = "file_id") %>% group_by(file_set) %>%
    dplyr::reframe(predtype = str_remove_all(str_split_1(unique(type),
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
      ncinfoi %>% dplyr::filter(date_start <= unique(dates) &
                                  date_end >= unique(dates)) %>% dplyr::pull(period) %>%
        unique()
    }
    else {
      NA
    }) %>% ungroup()
    infos_dim <- map_dfr(1:nrow(ncinfoi), function(i) {
      v <- nc_open(paste(nc.path, ncinfoi$nc.name[i],
                         sep = "/"))$var
      out <- ncinfoi[i, ] %>% left_join(map_dfr(v, function(x) {
        map_dfr(x$dim, function(d) {
          if (d$name == name_dimension[["depth"]] & max_depth <
              Inf) {
            return(data.frame(variable = x$name, dim = d$name,
                              n_values = length(d$vals[d$vals <= max_depth])))
          }
          else {
            return(data.frame(variable = x$name, dim = d$name,
                              n_values = length(d$vals)))
          }
        })
      }) %>% dplyr::mutate(nc.name = ncinfoi$nc.name[i]),
      by = c("nc.name", "variable")) %>% ungroup()
      return(out)
    })

    for (i in 1:length(rename_dimensions)) {
      if (!is.null(rename_dimensions[[i]])) {
        infos_dim$dim[infos_dim$dim == rename_dimensions[[i]]] <- name_dimension[[names(rename_dimensions)[i]]]
      }
    }

    if (!file.exists(paste0(nc.path, "/", f, ".shp"))) {
      datafile <- paste(nc.path, first(ncinfoi$nc.name[ncinfoi$file_set ==
                                                         f]), sep = "/")
      nc_data <- nc_open(datafile)
      namedims <- names(nc_data$dim)
      test <- try({
        lat <- c(unique(as.vector(ncvar_get(nc_data,
                                            namedims[namedims == name_dimension[["lat"]]]))))
        lon <- c(unique(as.vector(ncvar_get(nc_data,
                                            namedims[namedims == name_dimension[["lon"]]]))))
      })
      if (all(class(test) == "try-error")) {
        namedims <- names(nc_data$var)
        test <- try({
          lat <- c(unique(as.vector(ncvar_get(nc_data,
                                              namedims[namedims == name_dimension[["lat"]]]))))
          lon <- c(unique(as.vector(ncvar_get(nc_data,
                                              namedims[namedims == name_dimension[["lon"]]]))))
        })
      }
      ref_c <- first(match(infos_dim %>% dplyr::mutate(id = 1:n()) %>%
                             group_by(nc.name, variable) %>% dplyr::mutate(n_dim = length(unique(dim))) %>%
                             ungroup() %>% dplyr::slice_max(n_dim) %>% dplyr::slice_min(id) %>%
                             pull(variable), ncinfoi$variable))
      data <- map_dfr(match(infos_dim %>% pull(variable) %>%
                              unique(), ncinfoi$variable), function(i) {
                                infos_dimi <- infos_dim %>% dplyr::filter(variable ==
                                                                            ncinfoi$variable[ref_c] & nc.name == ncinfoi$nc.name[ref_c])
                                data_var <- try(ncvar_get(nc_data, ncinfoi$variable[i],
                                                          start = order_dim(infos_dimi, dimensions = set_names(map(infos_dimi$dim,
                                                                                                                   function(x) {
                                                                                                                     1
                                                                                                                   }), infos_dimi$dim)), count = order_dim(infos_dimi,
                                                                                                                                                           dimensions = set_names(map(infos_dimi$dim,
                                                                                                                                                                                      function(x) {
                                                                                                                                                                                        ifelse(x %in% c(name_dimension[["lat"]],
                                                                                                                                                                                                        name_dimension[["lon"]]), infos_dimi %>%
                                                                                                                                                                                                 dplyr::filter(dim == x) %>% pull(n_values),
                                                                                                                                                                                               1)
                                                                                                                                                                                      }), infos_dimi$dim)), verbose = FALSE))
                                if (all(class(data_var) == "try-error")) {
                                  cat("Trying to remove depth. If Time should be removed instead of depth, please adapt the function script.\n")
                                  data_var <- ncvar_get(nc_data, ncinfoi$variable[i],
                                                        start = order_dim(infos_dimi, dimensions = set_names(map(infos_dimi$dim,
                                                                                                                 function(x) {
                                                                                                                   1
                                                                                                                 }), infos_dimi$dim)), count = order_dim(infos_dimi,
                                                                                                                                                         dimensions = set_names(map(infos_dimi$dim,
                                                                                                                                                                                    function(x) {
                                                                                                                                                                                      ifelse(x %in% c(name_dimension[["lat"]],
                                                                                                                                                                                                      name_dimension[["lon"]]), infos_dimi %>%
                                                                                                                                                                                               dplyr::filter(dim == x) %>% pull(n_values),
                                                                                                                                                                                             1)
                                                                                                                                                                                    }), infos_dimi$dim), use_1value = F),
                                                        verbose = FALSE)
                                }
                                if (dim(data_var)[1] == length(lon) & dim(data_var)[2] ==
                                    length(lat)) {
                                  if (lat[1] < lat[2]) {
                                    data_var <- data_var[, length(lat):1]
                                  }
                                  data_var <- t(data_var)
                                }
                                newgrid <- raster(data_var, xmn = min(lon) -
                                                    (sort(lon)[2] - sort(lon)[1])/2, xmx = max(lon) +
                                                    (sort(lon)[length(lon)] - sort(lon)[length(lon) -
                                                                                          1])/2, ymn = min(lat) - (sort(lat)[2] -
                                                                                                                     sort(lat)[1])/2, ymx = max(lat) + (sort(lat)[length(lat)] -
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
        st_transform(crs = crs_meter) %>% st_coordinates()
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
                                                              lat_cent, X, Y, id)
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
    llon <- unique(data$lon_cent)
    llat <- unique(data$lat_cent)
    dlon <- median(sort(llon)[-1] - sort(llon)[-length(llon)])
    dlat <- median(sort(llat)[-1] - sort(llat)[-length(llat)])
    cat("Parallel processing with", n_cores, "cores will be used. Check that computed can handle it.\n\n",
        "In case of large files, few iterations might return an error due to memory limits. The function can then be re-run to re-do only these few iterations\n")
    cl <- makeCluster(n_cores, outfile = outfile)
    registerDoParallel(cl)
    try({
      all_run <- foreach(dt = 1:nrow(dates), .packages = c("purrr",
                                                           "rlang", "dplyr", "timeDate", "ncdf4", "stringr",
                                                           "data.table", "collapse"), .noexport = ls()[!(ls() %in%
                                                                                                           c("data", "llon", "llat", "dlon", "dlat", "max_depth",
                                                                                                             "ncinfoi", "nc.path", "expr", "order_dim",
                                                                                                             "crs_meter", "SDspace_radius",
                                                                                                             "create_dim", "origin_all", "all_days_period_ref",
                                                                                                             "infos_dim", "name_dimension", "f", "all_pixel.radius",
                                                                                                             "pixel.radius", "pred.type", "vertical_variables",
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
            nc_data <- nc_open(datafile)
            namedims <- names(nc_data$dim)
            test <- try({
              lat <- c(unique(as.vector(ncvar_get(nc_data,
                                                  namedims[namedims == name_dimension[["lat"]]]))))
              lon <- c(unique(as.vector(ncvar_get(nc_data,
                                                  namedims[namedims == name_dimension[["lon"]]]))))
            })
            if (all(class(test) == "try-error")) {
              namedims <- names(nc_data$var)
              test <- try({
                lat <- c(unique(as.vector(ncvar_get(nc_data,
                                                    namedims[namedims == name_dimension[["lat"]]]))))
                lon <- c(unique(as.vector(ncvar_get(nc_data,
                                                    namedims[namedims == name_dimension[["lon"]]]))))
              })
            }

            # lat <- sort(lat)
            # lon <- sort(lon)

            namedims <- names(nc_data$dim)
            if (!is.null(name_dimension[["depth"]])) {
              depth <- try(ncvar_get(nc_data, namedims[namedims ==
                                                         name_dimension[["depth"]]]))
              if (all(class(depth) == "try-error")) {
                namedims <- names(nc_data$var)
                depth <- ncvar_get(nc_data, namedims[namedims ==
                                                       name_dimension[["depth"]]])
              }
            } else {
              depth = 0
            }
            depth <- depth[depth <= max_depth]
            time_var <- name_dimension[["time"]]
            if (length(time_var) > 0) {
              time <- ncvar_get(nc_data, time_var)
              time <- 0:(length(time) - 1)
            }
            else {
              time <- 0
            }
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
            data_var_refi <- lapply(Predictor.name,
                                    function(pred) {
                                      # print(pred)
                                      infos_dimi <- infos_dim %>% dplyr::filter(variable ==
                                                                                  pred & nc.name == ncfile)
                                      data_var_ref <- try(ncvar_get(nc_data, pred,
                                                                    start = order_dim(infos_dimi, dimensions = set_names(map(infos_dimi$dim,
                                                                                                                             function(x) {
                                                                                                                               ifelse(x %in% name_dimension[["time"]],
                                                                                                                                      time1, 1)
                                                                                                                             }), infos_dimi$dim)), count = order_dim(infos_dimi,
                                                                                                                                                                     dimensions = set_names(map(infos_dimi$dim,
                                                                                                                                                                                                function(x) {
                                                                                                                                                                                                  ifelse(x %in% name_dimension[["time"]],
                                                                                                                                                                                                         numtimes, ifelse(x %in% name_dimension[["depth"]],
                                                                                                                                                                                                                          length(depth), infos_dimi %>%
                                                                                                                                                                                                                            dplyr::filter(dim == x) %>%
                                                                                                                                                                                                                            pull(n_values)))
                                                                                                                                                                                                }), infos_dimi$dim)), verbose = FALSE))
                                      if (all(class(data_var_ref) == "try-error")) {
                                        data_var_ref <- try(ncvar_get(nc_data,
                                                                      pred, start = order_dim(infos_dimi,
                                                                                              dimensions = set_names(map(infos_dimi$dim,
                                                                                                                         function(x) {
                                                                                                                           ifelse(x %in% name_dimension[["time"]],
                                                                                                                                  time1, 1)
                                                                                                                         }), infos_dimi$dim), use_1value = F),
                                                                      count = order_dim(infos_dimi, dimensions = set_names(map(infos_dimi$dim,
                                                                                                                               function(x) {
                                                                                                                                 ifelse(x %in% name_dimension[["time"]],
                                                                                                                                        numtimes, ifelse(x %in% name_dimension[["depth"]],
                                                                                                                                                         length(depth), infos_dimi %>%
                                                                                                                                                           dplyr::filter(dim == x) %>%
                                                                                                                                                           pull(n_values)))
                                                                                                                               }), infos_dimi$dim), use_1value = F),
                                                                      verbose = FALSE))
                                      }
                                      if (any(all_days_period$period < prd)) {
                                        last_period <- unique(all_days_period$period[all_days_period$period <
                                                                                       prd])
                                        ncfile1 <- unique(ncinfoi$nc.name[ncinfoi$period ==
                                                                            last_period & ncinfoi$variable == pred])
                                        datafile1 <- paste(nc.path, ncfile1,
                                                           sep = "/")
                                        nc_data1 <- nc_open(datafile1)
                                        if (length(time_var) > 0) {
                                          time11 <- ncvar_get(nc_data1, time_var)
                                          time11 <- 0:(length(time11) - 1)
                                        }
                                        else {
                                          time11 <- 0
                                        }
                                        SegDay1 <- length(as.character(timeSequence(from = unique(ncinfoi$date_start[ncinfoi$period ==
                                                                                                                       last_period & ncinfoi$variable == pred]),
                                                                                    to = last(all_days_period$date[all_days_period$period == last_period]), by = "day")))
                                        day.index1 <- which((time11 + 1) ==
                                                              SegDay1)
                                        if (length(day.index1) == 0) {
                                          print(paste("Not OK! Missing dates in nc file for: ",
                                                      SegDay1, sep = "  "))
                                        }
                                        time11 <- max(day.index1 - (max(all_time.period) -
                                                                      numtimes) + 1, 1)
                                        numtimes11 <- min(max(all_time.period) -
                                                            numtimes, day.index1)
                                        infos_dimi1 <- infos_dim %>% dplyr::filter(variable ==
                                                                                     pred & nc.name == ncfile1)
                                        data_var1 <- try(ncvar_get(nc_data1,
                                                                   pred, start = order_dim(infos_dimi1,
                                                                                           dimensions = set_names(map(infos_dimi1$dim,
                                                                                                                      function(x) {
                                                                                                                        ifelse(x %in% name_dimension[["time"]],
                                                                                                                               time11, 1)
                                                                                                                      }), infos_dimi1$dim)), count = order_dim(infos_dimi1,
                                                                                                                                                               dimensions = set_names(map(infos_dimi1$dim,
                                                                                                                                                                                          function(x) {
                                                                                                                                                                                            ifelse(x %in% name_dimension[["time"]],
                                                                                                                                                                                                   numtimes11, ifelse(x %in%
                                                                                                                                                                                                                        name_dimension[["depth"]],
                                                                                                                                                                                                                      length(depth), infos_dimi1 %>%
                                                                                                                                                                                                                        dplyr::filter(dim ==
                                                                                                                                                                                                                                        x) %>% pull(n_values)))
                                                                                                                                                                                          }), infos_dimi1$dim)), verbose = FALSE))
                                        if (all(class(data_var1) == "try-error")) {
                                          data_var1 <- try(ncvar_get(nc_data1,
                                                                     pred, start = order_dim(infos_dimi1,
                                                                                             dimensions = set_names(map(infos_dimi1$dim,
                                                                                                                        function(x) {
                                                                                                                          ifelse(x %in% name_dimension[["time"]],
                                                                                                                                 time11, 1)
                                                                                                                        }), infos_dimi1$dim), use_1value = F),
                                                                     count = order_dim(infos_dimi1,
                                                                                       dimensions = set_names(map(infos_dimi1$dim,
                                                                                                                  function(x) {
                                                                                                                    ifelse(x %in% name_dimension[["time"]],
                                                                                                                           numtimes11, ifelse(x %in%
                                                                                                                                                name_dimension[["depth"]],
                                                                                                                                              length(depth), infos_dimi1 %>%
                                                                                                                                                dplyr::filter(dim ==
                                                                                                                                                                x) %>% pull(n_values)))
                                                                                                                  }), infos_dimi1$dim), use_1value = F),
                                                                     verbose = FALSE))
                                        }
                                        dimtime <- infos_dim %>% dplyr::filter(variable ==
                                                                                 pred & nc.name == ncfile & n_values >
                                                                                 1)
                                        if (any(dimtime$dim == name_dimension[["time"]])) {
                                          dimtime <- dimtime %>% dplyr::mutate(id = 1:n()) %>%
                                            dplyr::filter(dim == name_dimension[["time"]]) %>%
                                            pull(id)
                                        }
                                        else {
                                          dimtime <- length(dim(data_var_ref)) +
                                            1
                                        }
                                        data_var_ref <- abind::abind(data_var1,
                                                                     data_var_ref, along = dimtime)
                                        nc_close(nc_data1)
                                      }
                                      infos_dim2 <- infos_dim %>% dplyr::filter(variable ==
                                                                                  pred & nc.name == ncfile)

                                      new_name_dimension <- name_dimension

                                      if (length(time_var) == 0) {
                                        infos_dim2 <- infos_dim2 %>% group_by(nc.name,
                                                                              variable, file_set, date_start,
                                                                              date_end, period) %>% dplyr::reframe(dim = c(dim,
                                                                                                                           "time"), n_values = c(n_values,
                                                                                                                                                 1)) %>% ungroup

                                        new_name_dimension[["time"]] <- "time"
                                      }

                                      for (i in 1:length(new_name_dimension)) {
                                        if (is.null(new_name_dimension[i])) {
                                          new_name_dimension[[i]] <- "NAAAA"
                                        }
                                      }

                                      dimensions_list <- set_names(map(infos_dim2$dim,
                                                                       function(x) {
                                                                         if (x == new_name_dimension[["time"]]) {
                                                                           1:max(all_time.period)
                                                                         } else if (x == new_name_dimension[["depth"]]) {
                                                                           depth
                                                                         }
                                                                         else if (x == new_name_dimension[["lat"]]) {
                                                                           lat
                                                                         }
                                                                         else if (x == new_name_dimension[["lon"]]) {
                                                                           lon
                                                                         }
                                                                       }), infos_dim2$dim)

                                      dimensions <- create_dim(infos_dim2 %>%
                                                                 dplyr::filter(variable == pred & nc.name ==
                                                                                 ncfile), dimensions = dimensions_list,
                                                               name_dim = new_name_dimension) %>% as.data.frame()

                                      for (i in seq_along(dimensions_list)) {
                                        # print(i)
                                        col_name <- names(dimensions_list)[i]

                                        if (col_name == new_name_dimension[["lon"]]) {
                                          col_name <- "lon"
                                        }
                                        else if (col_name == new_name_dimension[["lat"]]) {
                                          col_name <- "lat"
                                        }
                                        else if (col_name == new_name_dimension[["time"]]) {
                                          col_name <- "t"
                                        }
                                        else if (col_name == new_name_dimension[["depth"]]) {
                                          col_name <- "d"
                                        }

                                        if (length(unique(dimensions_list[[i]])) > 1) {
                                          if (dimensions_list[[i]][2] > dimensions_list[[i]][1]) {
                                            dimensions <- dimensions %>%
                                              arrange(across(all_of(col_name)))
                                          } else {
                                            dimensions <- dimensions %>%
                                              arrange(desc(across(all_of(col_name))))
                                          }
                                        }
                                      }

                                      if (pred %in% vertical_variables) {
                                        out <- cbind(data_var_ref) %>% as.data.frame() %>%
                                          cbind(dimensions) %>% dplyr::filter(lon >=
                                                                                (min(llon) - dlon/2) & lon <= (max(llon) +
                                                                                                                 dlon/2) & lat >= (min(llat) - dlat/2) &
                                                                                lat <= (max(llat) + dlat/2)) %>%
                                          dplyr::arrange(d) %>% dplyr::group_by(lon,
                                                                                lat, t) %>% dplyr::summarise(va_data_var_ref = mean(data_var_ref,
                                                                                                                                    na.rm = TRUE), bot_data_var_ref = dplyr::last(na.omit(data_var_ref)),
                                                                                                             data_var_ref = dplyr::first(na.omit(data_var_ref)),
                                                                                                             d = 1, .groups = "drop") %>% dplyr::mutate(d_data_var_ref = data_var_ref -
                                                                                                                                                          bot_data_var_ref) %>% dplyr::rename_with(.fn = ~c(pred,
                                                                                                                                                                                                            paste0("VertAv_", pred), paste0("Bottom_",
                                                                                                                                                                                                                                            pred), paste0("VertDiff_", pred),
                                                                                                                                                                                                            "d"), .cols = c("data_var_ref",
                                                                                                                                                                                                                            "va_data_var_ref", "bot_data_var_ref",
                                                                                                                                                                                                                            "d_data_var_ref", "d"))
                                      }
                                      else {
                                        out <- cbind(data_var_ref) %>% as.data.frame() %>%
                                          cbind(dimensions) %>% dplyr::filter(lon >=
                                                                                (min(llon) - dlon/2) & lon <= (max(llon) +
                                                                                                                 dlon/2) & lat >= (min(llat) - dlat/2) &
                                                                                lat <= (max(llat) + dlat/2)) %>%
                                          dplyr::group_by(lon, lat, t) %>%
                                          dplyr::summarise(data_var_ref = dplyr::first(na.omit(data_var_ref)),
                                                           d = 1) %>% dplyr::ungroup() %>%
                                          dplyr::rename(`:=`(!!paste0(pred),
                                                             data_var_ref))
                                      }
                                      return(out)
                                    })
            nc_close(nc_data)
            data_var_ref <- data_var_refi[[1]]
            if (!("d" %in% colnames(data_var_ref))) {
              data_var_ref <- data_var_ref %>% dplyr::mutate(d = depth[1])
            }
            if (length(Predictor.name) > 1) {
              for (i in 2:length(data_var_refi)) {
                if (!("d" %in% colnames(data_var_refi[[i]]))) {
                  data_var_refi[[i]] <- data_var_refi[[i]] %>%
                    dplyr::mutate(d = depth[1])
                }
                data_var_ref <- data_var_ref %>% left_join(data_var_refi[[i]],
                                                           by = c("lon", "lat", "t", "d"))
              }
            }
            rm(data_var_refi)
            data_var_ref <- data_var_ref %>% as.data.frame()
            for (c in colnames(data_var_ref)) {
              data_var_ref[, c] <- c(data_var_ref[,
                                                  c])
            }
            id_NA <- data_var_ref %>% dplyr::mutate(id = 1:n()) %>%
              dplyr::select(-c(lon, lat, t, d))
            id_NA <- id_NA %>% dplyr::filter(rowSums(is.na(across(colnames(.)))) >=
                                               (ncol(id_NA) - 1)) %>% pull(id)
            if (length(id_NA) > 0) {
              data_var_ref <- data_var_ref[-id_NA, ]
            }
            if (all(!is.na(expr))) {
              for (e in expr) {
                Predictor.name <- c(Predictor.name,
                                    str_remove_all(str_split_1(e, fixed("="))[1],
                                                   " "))
                data_var_ref <- data_var_ref %>% mutate(!!parse_expr(e))
                colnames(data_var_ref)[ncol(data_var_ref)] <- last(Predictor.name)
              }
            }
            data_var_ref <- data_var_ref %>% group_by(lon) %>%
              dplyr::mutate(lon = llon[which.min(abs(unique(lon) -
                                                       llon))]) %>% ungroup() %>% group_by(lat) %>%
              dplyr::mutate(lat = llat[which.min(abs(unique(lat) -
                                                       llat))]) %>% ungroup()
            if (!(all(unique(paste(floor(data$lat_cent *
                                         10^5)/10^5, floor(data$lon_cent * 10^5)/10^5)) %in%
                      unique(paste(floor(data_var_ref$lat *
                                         10^5)/10^5, floor(data_var_ref$lon *
                                                           10^5)/10^5))) | all(unique(paste(floor(data_var_ref$lat *
                                                                                                  10^5)/10^5, floor(data_var_ref$lon * 10^5)/10^5)) %in%
                                                                               unique(paste(floor(data$lat_cent * 10^5)/10^5,
                                                                                            floor(data$lon_cent * 10^5)/10^5))))) {
              stop("Error in the function and from the data: lon or lat not in the data")
            }

            Predictor.name <- c(Predictor.name, na.omit(do.call("c",
                                                                lapply(vertical_variables, function(c) {
                                                                  if (c %in% Predictor.name) {
                                                                    paste0(c("VertAv_", "Bottom_", "VertDiff_"),
                                                                           c)
                                                                  } else {
                                                                    NA
                                                                  }
                                                                }))))
            data_var_ref <- data %>% dplyr::mutate(lon_cent = floor(lon_cent *
                                                                      10^5)/10^5, lat_cent = floor(lat_cent *
                                                                                                     10^5)/10^5) %>% left_join(data_var_ref %>%
                                                                                                                                 dplyr::mutate(lon_cent = floor(lon * 10^5)/10^5,
                                                                                                                                               lat_cent = floor(lat * 10^5)/10^5),
                                                                                                                               by = c("lon_cent", "lat_cent")) %>%
              dplyr::mutate(id_nc = id,
                            lon = lon_cent, lat = lat_cent) %>%
              dplyr::select(t, d, lon, lat, X, Y, id_nc, all_of(Predictor.name))

            numtimes <- max(all_time.period)

            data_var_ref_t1 <- data_var_ref %>% dplyr::filter(t ==
                                                                dplyr::first(na.omit(t))) %>% dplyr::select(lon,
                                                                                                            lat, X, Y, id_nc) %>% as.data.frame()

            data_var_ref <- data_var_ref %>%
              dplyr::select(-c(X, Y))

            XX <- sort(unique(data$X))
            YY <- sort(unique(data$Y))

            run_mean_SDspace <- any(c("mean", "SDspace") %in%
                                      pred.type)
            if (!run_mean_SDspace) {
              all_pixel.radius <- 0
            }
            final <- lapply(all_pixel.radius, function(pixel.radius) {
              if (pixel.radius != all_pixel.radius[1] &
                  !run_mean_SDspace) {
                return(NULL)
              }

              if (SDspace_radius == "km") {
                res_lat <- pixel.radius
                res_lon <- pixel.radius

                outM <- map_dfr(unique(data$X),
                                function(l) {
                                  data_var_ref_t1l <- data_var_ref_t1[abs(l -
                                                                            data_var_ref_t1$X) <= res_lon,
                                  ]
                                  return(data[data$X == l, ] %>%
                                           group_by(id) %>% group_map(~{
                                             c(lat_cent = .x$lat_cent, lon_cent = .x$lon_cent,
                                               X = .x$X, Y = .x$Y,
                                               id = .y$id,
                                               id_nc = paste(data_var_ref_t1l$id_nc[sqrt((.x$X - data_var_ref_t1l$X)^2 +
                                                                                           (.x$Y - data_var_ref_t1l$Y)^2) <= pixel.radius],
                                                             collapse = ","))
                                           }) %>% bind_rows())
                                })
              } else if (SDspace_radius == "pixel") {
                res_lat <- mean(sort(lat)[-1] - sort(lat)[-length(lat)],
                                na.rm = T) * (pixel.radius + 0.5)
                res_lon <- mean(sort(lon)[-1] - sort(lon)[-length(lon)],
                                na.rm = T) * (pixel.radius + 0.5)

                outM <- map_dfr(unique(data$lon_cent),
                                function(l) {
                                  data_var_ref_t1l <- data_var_ref_t1[abs(l -
                                                                            data_var_ref_t1$lon) <= res_lon,
                                  ]
                                  return(data[data$lon_cent == l, ] %>%
                                           group_by(id) %>% group_map(~{
                                             c(lat_cent = .x$lat_cent, lon_cent = .x$lon_cent,
                                               id = .y$id, id_nc = paste(data_var_ref_t1l$id_nc[abs(.x$lat_cent -
                                                                                                      data_var_ref_t1l$lat) <= res_lat],
                                                                         collapse = ","))
                                           }) %>% bind_rows())
                                })
              } else {
                stop("SDspace_radius is different from 'km' and 'pixel': please choose one of them even if SDspace is not calculated.")
              }

              outM <- outM %>%
                group_by(id) %>% dplyr::reframe(id_nc = str_split_1(id_nc,
                                                                                         ","), lon_cent = lon_cent, lat_cent = lat_cent) %>%
                dplyr::mutate(across(colnames(.), ~as.numeric(.x))) %>%
                left_join(data_var_ref, by = "id_nc",
                          relationship = "many-to-many") %>%
                dplyr::mutate(dist = abs(lat_cent -
                                           lat) + abs(lon_cent - lon)) %>%
                dplyr::select(id,
                              t, id_nc, dist, all_of(Predictor.name)) %>%
                distinct()

              if (run_mean_SDspace) {
                outMsd <- as.data.table(outM)
                rm(outM)
                if (all(c("SDspace", "mean") %in% pred.type)) {
                  outMsd <- outMsd[, c(list(id_nc = id_nc[which.min(dist)]),
                                       setNames(fmean(.SD, na.rm = TRUE),
                                                paste0(Predictor.name, ".mean")),
                                       setNames(fsd(.SD, na.rm = TRUE),
                                                paste0(Predictor.name, ".SDspace"))),
                                   by = list(id, t), .SDcols = Predictor.name]
                }
                else if ("mean" %in% pred.type) {
                  outMsd <- outMsd[, c(list(id_nc = id_nc[which.min(dist)]),
                                       setNames(fmean(.SD, na.rm = TRUE),
                                                paste0(Predictor.name, ".mean"))),
                                   by = list(id, t), .SDcols = Predictor.name]
                }
                else if ("SDspace" %in% pred.type) {
                  outMsd <- outMsd[, c(list(id_nc = id_nc[which.min(dist)]),
                                       setNames(fsd(.SD, na.rm = TRUE),
                                                paste0(Predictor.name, ".SDspace"))),
                                   by = .(id, t), .SDcols = Predictor.name]
                }
                outM <- outMsd %>% as_tibble()
                colnames(outM)[str_detect(colnames(outM),
                                          fixed(".mean")) | str_detect(colnames(outM),
                                                                       fixed(".SDspace"))] <- paste0(colnames(outM)[str_detect(colnames(outM),
                                                                                                                               fixed(".mean")) | str_detect(colnames(outM),
                                                                                                                                                            fixed(".SDspace"))], "_", pixel.radius,
                                                                                                     "p")
                outM <- outM %>% left_join(data_var_ref %>%
                                             dplyr::select(id_nc, t, all_of(Predictor.name)),
                                           by = c("t", "id_nc")) %>% dplyr::select(-id_nc)
              }
              outM <- as.data.table(outM)
              cols1 <- colnames(outM)[stringr::str_detect(colnames(outM),
                                                          fixed(".mean")) | stringr::str_detect(colnames(outM),
                                                                                                fixed(".SDspace"))]
              cols2 <- Predictor.name
              outfinal <- map(all_time.period, function(time.period) {
                if (time.period == 1) {
                  outf <- outM[which(outM$t %in% (max(all_time.period) -
                                                    time.period + 1):numtimes), -c("t")]
                  colnames(outf)[colnames(outf) %in%
                                   cols1] <- paste0(cols1, ".mean")
                  colnames(outf)[colnames(outf) %in%
                                   cols2] <- paste0(cols2, ".center")
                }
                else {
                  if (any(c("SDspace", "mean") %in%
                          pred.type) & (any(c("center", "SDtime") %in%
                                            pred.type) & pixel.radius == all_pixel.radius[1])) {
                    res1 <- outM[t %in% (max(all_time.period) -
                                           time.period + 1):numtimes, c(setNames(lapply(.SD,
                                                                                        mean, na.rm = TRUE), paste0(cols1,
                                                                                                                    ".mean")), setNames(lapply(.SD,
                                                                                                                                               sd, na.rm = TRUE), paste0(cols1,
                                                                                                                                                                         ".SDtime"))), by = id, .SDcols = cols1]
                    res2 <- outM[t %in% (max(all_time.period) -
                                           time.period + 1):numtimes, c(setNames(lapply(.SD,
                                                                                        mean, na.rm = TRUE), paste0(cols2,
                                                                                                                    ".center")), setNames(lapply(.SD,
                                                                                                                                                 sd, na.rm = TRUE), paste0(cols2,
                                                                                                                                                                           ".SDtime"))), by = id, .SDcols = cols2]
                    outf <- merge(res1, res2, by = "id")
                    rm(res1)
                    rm(res2)
                  }
                  else if (any(c("SDspace", "mean") %in%
                               pred.type)) {
                    outf <- outM[t %in% (max(all_time.period) -
                                           time.period + 1):numtimes, c(setNames(lapply(.SD,
                                                                                        mean, na.rm = TRUE), paste0(cols1,
                                                                                                                    ".mean")), setNames(lapply(.SD,
                                                                                                                                               sd, na.rm = TRUE), paste0(cols1,
                                                                                                                                                                         ".SDtime"))), by = id, .SDcols = cols1]
                  }
                  else if (any(c("center", "SDtime") %in%
                               pred.type)) {
                    outf <- outM[t %in% (max(all_time.period) -
                                           time.period + 1):numtimes, c(setNames(lapply(.SD,
                                                                                        mean, na.rm = TRUE), paste0(cols2,
                                                                                                                    ".center")), setNames(lapply(.SD,
                                                                                                                                                 sd, na.rm = TRUE), paste0(cols2,
                                                                                                                                                                           ".SDtime"))), by = id, .SDcols = cols2]
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
                  out_final <- out_final %>% left_join(outfinal[[i]] %>%
                                                         dplyr::select(id, all_of(colnames(outfinal[[i]])[!(colnames(outfinal[[i]]) %in%
                                                                                                           colnames(out_final))])),
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
    cat("\n\n")
  }
  invisible("Finished")
}


