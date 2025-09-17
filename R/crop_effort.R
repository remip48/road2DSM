#' Crop effort to the grid extracted from nc files
#'
#' The grid extracted from nc file is supposed to have been extracted with extract_nc
#'
#' @param effort
#' @param file_set_directory
#' @param keep
#'
#' @return
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @export
#'
#' @examples
crop_effort <- function(effort,
                        file_set_directory, # path to extracted nc files, so that .shp grids are loaded directly. Let NULL if you want to use only study_area
                        keep = "all_that_intersect", # or "only_fully_in"
                        study_area = NULL # additional sf file to filter
                        ) {

  if (!is.null(file_set_directory)) {
    cat("Load sf grids from", file_set_directory, " to filter effort with keep =", keep, ".\n")
    listf <- list.files(file_set_directory)

    grid_NS_BS <- map_dfr(listf[str_detect(listf, fixed(".shp")) & str_detect(listf, "file_set_")], function(f) {
      read_sf(file.path(file_set_directory, f)) %>%
        dplyr::mutate(file = f) %>%
        dplyr::select(file)
    }) %>%
      st_transform(crs = 3035) %>%
      group_by() %>%
      dplyr::summarise(do_union = F) %>%
      st_cast("MULTIPOLYGON")

    if (keep == "all_that_intersect") {
      labels <- effort %>%
        st_make_valid() %>%
        dplyr::filter(c(st_intersects(.,
                                      grid_NS_BS %>%
                                        st_make_valid(),
                                      sparse = F))) %>%
        st_drop_geometry() %>%
        pull(label)

      out <- effort %>%
        dplyr::filter(label %in% labels)

    } else if (keep == "only_fully_in") {
      labels <- effort %>%
        st_make_valid() %>%
        st_difference(.,
                      grid_NS_BS %>%
                        st_make_valid()) %>%
        st_drop_geometry() %>%
        pull(label)

      out <- effort %>%
        dplyr::filter(!(label %in% labels))
    }
  } else {
    cat("No sf grids from NC files were set in file_set_directory.\n")
  }

  if (all(!is.null(study_area))) {
    cat("Filter effort in study_area with keep =", keep, ".\n")

    study_area <- study_area  %>%
      st_transform(crs = 3035) %>%
      group_by() %>%
      dplyr::summarise(do_union = F) %>%
      st_cast("MULTIPOLYGON") %>%
      st_make_valid()

    if (keep == "all_that_intersect") {
      labels <- out %>%
        st_make_valid() %>%
        dplyr::filter(c(st_intersects(.,
                                      study_area,
                                      sparse = F))) %>%
        st_drop_geometry() %>%
        pull(label)

      out <- out %>%
        dplyr::filter(label %in% labels)

    } else if (keep == "only_fully_in") {
      labels <- out %>%
        st_make_valid() %>%
        st_difference(.,
                      study_area) %>%
        st_drop_geometry() %>%
        pull(label)

      out <- out %>%
        dplyr::filter(!(label %in% labels))
    }
  } else {
      cat("No additional study_area was provided to filter effort.\n")
  }

  cat ("Done!")

  return(out)
}
