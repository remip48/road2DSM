#' Title
#'
#' @param effort
#' @param file_set_directory
#' @param keep
#'
#' @return
#' @export
#'
#' @examples
crop_effort <- function(effort,
                        file_set_directory,
                        keep = "all_that_intersect" # or "only_fully_in"
                        ) {
  listf <- list.files(file_set_directory)

  grid_NS_BS <- map_dfr(listf[str_detect(listf, fixed(".shp")) & str_detect(listf, "file_set_")], function(f) {
    read_sf(f) %>%
      dplyr::mutate(file = f) %>%
      dplyr::select(f)
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
                      st_make_valid())
      st_drop_geometry() %>%
      pull(label)

      out <- effort %>%
        dplyr::filter(!(label %in% labels))
  }

  return(out)
}
