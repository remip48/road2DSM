#' Calculate SDspace in addition to return center and SDtime over a prediction grid
#'
#' Can be used on files from extract_grid (prediction grids, that were possibly upscaled by upscale_grid or created by creategrid). Calculate SDspace only on columns that contain ".center", as returned by extract_nc. Return all columns from the initial grid.
#' return SDspace variables are called as paste0(centers, ".SDspace_", pixel.radius, "p") where center is the name of the variable center column (e.g. SST.center_1d)
#'
#' @param all_pixel.radius
#' @param n_cores
#' @param outfile
#' @param file_directory
#' @param writing_directory
#' @param variables variables to calculate SDspace. Let null if you want all variables with .center to be calculated.
#' @param crs_meter
#' @param SDspace_radius km or pixel, depending on the unit of values you put in all_pixel.radius (if SDspace is used)
#'
#' @return
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @export
#'
#' @examples
upscale_SDspace <- function (file_directory, writing_directory, all_pixel.radius, variables = NULL,
                             crs_meter = 3035, SDspace_radius = "pixel",
                             n_cores = NULL, outfile = "log.txt")
{

  list_files <- list.files(file_directory)

  list_files <- list_files[str_detect(list_files, fixed(".rds")) & str_detect(list_files, fixed("prediction_grid_"))]

  static_grid <- read_sf(paste0(file_directory, "/prediction_static_grid.shp"))

  write_sf(static_grid,
           paste0(writing_directory, "/prediction_static_grid.shp"))

  rm(static_grid)

  if (is.null(n_cores)) {
    n_cores <- detectCores() * 2/4
  }
  cat("Parallel processing with", n_cores, "cores will be used. Check that computed can handle it.\n")

  cl <- makeCluster(n_cores, outfile = outfile)
  registerDoParallel(cl)

  # c("purrr",
  #   "dplyr", "timeDate", "ncdf4", "stringr", "data.table",
  #   "collapse")

  try({
    all_run <- foreach(lf = list_files, .packages = c("stringr", "dplyr", "data.table", "collapse", "purrr"),
                       .noexport = ls()[!(ls() %in% c("file_directory",
                                                      "crs_meter", "SDspace_radius",
                                                      "writing_directory", "all_pixel.radius", "variables"))]) %dopar%
      {
        x <- str_split_1(lf, "_")
        file_to_save <- paste0(paste(x[-length(x)], collapse = "_"), "_SDspace_", last(x))

        if (!file.exists(file.path(writing_directory, file_to_save))) {

          data <- readRDS(file.path(file_directory, lf))

          data_var_ref <- data %>%
            dplyr::select(id, X, Y, lon_cent, lat_cent, all_of(colnames(.)[str_detect(colnames(.), "center") |
                                                                             (str_detect(colnames(.), "SDtime") & !(str_detect(colnames(.), "mean")
                                                                                                                    | str_detect(colnames(.), "SDspace")))])) %>%
            dplyr::mutate(id_nc = id)

          centers <- colnames(data_var_ref)[str_detect(colnames(data_var_ref), "center")]

          if (any(!is.null(variables))) {
            centers <- na.omit(do.call("c", lapply(centers, function(c) {
              ifelse(any(str_detect(c, fixed(paste0(variables, ".")))),
                     c,
                     NA)
            })))
          }

          listX <- unique(data_var_ref$X)
          listY <- unique(data_var_ref$Y)

          final <- lapply(all_pixel.radius, function(pixel.radius) {
            resY <- mean(sort(listX)[-1] - sort(listX)[-length(listX)],
                         na.rm = T) * (pixel.radius + 0.5)
            resX <- mean(sort(listY)[-1] - sort(listY)[-length(listY)],
                         na.rm = T) * (pixel.radius + 0.5)

            if (SDspace_radius == "pixel") {
              outM <- map_dfr(unique(data_var_ref$X),
                              function(l) {
                                # print(l)
                                data.var_refX <- data_var_ref[abs(l - data_var_ref$X) <= resX, ]
                                return(data_var_ref[data_var_ref$X == l, ] %>%
                                         group_by(id) %>%
                                         group_map(~{
                                           c(id = .y$id,
                                             id_nc = paste(data.var_refX$id_nc[abs(.x$Y - data.var_refX$Y) <= resY],
                                                           collapse = ","))
                                         }) %>%
                                         bind_rows())
                              })
            } else if (SDspace_radius == "km") {
              outM <- map_dfr(unique(data_var_ref$X),
                              function(l) {
                                # print(l)
                                data.var_refX <- data_var_ref[abs(l - data_var_ref$X) <= pixel.radius, ]
                                return(data_var_ref[data_var_ref$X == l, ] %>%
                                         group_by(id) %>%
                                         group_map(~{
                                           c(id = .y$id,
                                             id_nc = paste(data.var_refX$id_nc[sqrt((.x$Y - data.var_refX$Y)^2 +
                                                                                      (.x$X - data.var_refX$X)^2) <= pixel.radius],
                                                           collapse = ","))
                                         }) %>%
                                         bind_rows())
                              })
            } else {
              stop("SDspace_radius is different from 'km' and 'pixel': please choose one of them even if SDspace is not calculated.")
            }
             %>% group_by(id) %>%
              dplyr::reframe(id_nc = str_split_1(id_nc, ",")) %>%
              dplyr::mutate(across(colnames(.), ~ as.numeric(.x))) %>%
              left_join(data_var_ref %>%
                          dplyr::select(id_nc, all_of(centers)), by = "id_nc",
                        relationship = "many-to-many") %>%
              distinct()

            setDT(outM)

            outM <- outM[, lapply(.SD, sd, na.rm = TRUE),
                         by = id, .SDcols = centers]

            setnames(outM, old = centers, new = paste0(centers, ".SDspace_", pixel.radius, "p"))

            outM <- outM %>%
              as_tibble()

            # outM <- data_var_ref %>%
            #   left_join(outM, by = c("id")) %>%
            #   dplyr::select(-id_nc)

            return(outM)
          })

          out <- data# %>%
          # dplyr::select(-all_of(colnames(.)[str_detect(colnames(.), "center") |
          #                                     str_detect(colnames(.), "SDspace") |
          #                                     str_detect(colnames(.), "SDtime") |
          #                                     str_detect(colnames(.), "mean")]))

          for (i in 1:length(final)) {
            if (any(!is.null(final[[i]]))) {
              out <- out %>% left_join(final[[i]] %>%
                                         dplyr::select(id, all_of(colnames(final[[i]])[!(colnames(final[[i]]) %in%
                                                                                           colnames(out))])), by = "id")
            }
          }

          saveRDS(out,
                  file.path(writing_directory, file_to_save))

        }
        return(NULL)
      }
  })

  try(stopCluster(cl))
  gc()

  # return(NULL)
  invisible("Finished")
}

