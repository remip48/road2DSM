#' Calculate SDspace in addition to return center and SDtime over a prediction grid
#'
#' Can be used on files from extract_grid (prediction grids, that were possibly upscaled by upscale_grid or created by creategrid). Calculate SDspace only on columns that contain ".center", as returned by extract_nc. Return all columns from the initial grid.
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
upscale_SDspace <- function (file_directory, writing_directory, all_pixel.radius,
                            n_cores = NULL, outfile = "log.txt")
{

  list_files <- list.files(file_directory)

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
                       .noexport = ls()[!(ls() %in% c("file_directory", "writing_directory", "all_pixel.radius"))]) %dopar%
      {
        x <- str_split_1(lf, "_")
        file_to_save <- paste0(paste(x[-length(x)], collapse = "_"), "_SDspace_", last(x))

        if (!file.exists(file.path(writing_directory, file_to_save))) {

          data <- readRDS(file.path(file_directory, lf))

          data.var_ref <- data %>%
            dplyr::select(id, X, Y, lon_cent, lat_cent, all_of(colnames(.)[str_detect(colnames(.), "center") |
                                                                             (str_detect(colnames(.), "SDtime") & !(str_detect(colnames(.), "mean")
                                                                                                                    | str_detect(colnames(.), "SDspace")))])) %>%
            dplyr::mutate(id_nc = id)

          centers <- colnames(data.var_ref)[str_detect(colnames(data.var_ref), "center")]

          listX <- unique(data.var_ref$X)
          listY <- unique(data.var_ref$Y)

          final <- lapply(all_pixel.radius, function(pixel.radius) {
            resY <- mean(sort(listX)[-1] - sort(listX)[-length(listX)],
                         na.rm = T) * (pixel.radius + 0.5)
            resX <- mean(sort(listY)[-1] - sort(listY)[-length(listY)],
                         na.rm = T) * (pixel.radius + 0.5)

            outM <- map_dfr(unique(data.var_ref$X),
                            function(l) {
                              # print(l)
                              data.var_refX <- data.var_ref[abs(l - data.var_ref$X) <= resX, ]
                              return(data.var_ref[data.var_ref$X == l, ] %>%
                                       group_by(id) %>%
                                       group_map(~{
                                         c(id = .y$id,
                                           id_nc = paste(data.var_refX$id_nc[abs(.x$Y - data.var_refX$Y) <= resY],
                                                         collapse = ","))
                                       }) %>%
                                       bind_rows())
                            }) %>% group_by(id) %>%
              dplyr::reframe(id_nc = str_split_1(id_nc, ",")) %>%
              dplyr::mutate(across(colnames(.), ~as.numeric(.x))) %>%
              left_join(data.var_ref %>%
                          dplyr::select(id_nc, all_of(centers)), by = "id_nc",
                        relationship = "many-to-many")

            setDT(outM)

            outM <- outM[, lapply(.SD, sd, na.rm = TRUE),
                         by = id, .SDcols = centers]

            setnames(outM, old = centers, new = paste0(centers, ".SDspace_", pixel.radius, "p"))

            outM <- outM %>%
              as_tibble()

            # outM <- data.var_ref %>%
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

