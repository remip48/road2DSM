#' Create sf grid based on effort bbox
#'
#' Use effort and addition input bbox to create a grid - that can be crop to the study area or the countries, as wished.
#' Effort must be an sf object, and grid will be created in the same crs than it. "resolution" must therefore also be in the same unit than this crs.
#'
#' @param effort
#' @param resolution
#' @param sea
#' @param country
#'
#' @return
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @export
#'
#' @examples
create_grid <- function(effort,
                         resolution,
                         bbox_grid = NULL,
                         sea = NULL,
                         country = NULL) {
  # pkg::fun(stars)
  # pkg::fun(sf)
  # pkg::fun(raster)

  bb <- st_bbox(effort)
  bb[1] <- bb[1] - resolution/2
  bb[2] <- bb[2] - resolution/2

  cat("If bbox_grid = NULL, grid will automatically be extended from st_bbox(effort) ± resolution/2.\n")

  if (all(!is.null(bbox_grid))) {
    if (length(bbox_grid) != 4) {
      cat("bbox_grid needs 4 elements: xmin, ymin, xmax, ymax. If you want to use the bbox_grid of the effort for one of them, just let NA in it respective value.\n")
    }

    for (i in 1:4) {
      if (!is.na(bbox_grid[i])) {
        bb[i] <- bbox_grid[i]
      }
    }
  }

  # bb
  # bb[1] <- floor(bb[1] / 1000) * 1000
  # bb[2] <- floor(bb[2] / 1000) * 1000
  bb[3] <- bb[1] + ceiling((bb[3] - bb[1]) / resolution) * resolution
  bb[4] <- bb[2] + ceiling((bb[4] - bb[2]) / resolution) * resolution
  # bb
  # (bb[3] - bb[1]) / resolution
  # (bb[4] - bb[2]) / resolution

  listX <- seq(bb[1] + resolution / 2, bb[3] - resolution / 2, resolution)
  listY <- seq(bb[2] + resolution / 2, bb[4] - resolution / 2, resolution)

  grid5km <- raster::raster(vals = 1,
                            xmn = bb[1], ymn = bb[2],
                            xmx = bb[3], ymx = bb[4],
                            resolution = c(resolution, resolution),
                            crs = st_crs(effort)$proj4string
  ) %>%
    st_as_stars() %>%
    st_as_sf() %>%
    st_transform(crs = st_crs(effort))

  grid5km <- grid5km %>%
    dplyr::mutate(X = st_coordinates(st_centroid(.))[,1],
                  Y = st_coordinates(st_centroid(.))[,2]) %>%
    group_by(X, Y) %>%
    dplyr::mutate(X = listX[which.min(abs(X - listX))],
                  Y = listY[which.min(abs(Y - listY))]) %>%
    ungroup()

  cc <- grid5km %>%
    sf::st_drop_geometry() %>%
    sf::st_as_sf(coords = c("X", "Y"), crs = sf::st_crs(grid)) %>%
    sf::st_transform(crs = 4326) %>%
    dplyr::mutate(lon_cent = st_coordinates(.)[,1],
                  lat_cent = st_coordinates(.)[,2])

  if (all(!is.null(sea))) {
    sea <- sea %>%
      st_transform(crs = st_crs(effort)) %>%
      group_by() %>%
      dplyr::summarise(do_union = F) %>%
      st_cast("MULTIPOLYGON") %>%
      # st_simplify() %>%
      st_make_valid()

    grid5km <- grid5km %>%
      # dplyr::filter(c(st_intersects(.,
      #                               sea,
      #                               sparse = F))) %>%
      st_intersection(.,
                      sea) %>%
      st_make_valid()

    cat("Coordinates written in the output files are coordinates of the initial centroid (before intersecting with sea).\n")
  }

  if (all(!is.null(country))) {
    country <- country %>%
      st_transform(crs = st_crs(effort)) %>%
      group_by() %>%
      dplyr::summarise(do_union = F) %>%
      st_cast("MULTIPOLYGON") %>%
      # st_simplify() %>%
      st_make_valid()

    grid5km <- grid5km %>%
      # dplyr::filter(!c(st_intersects(.,
      #                                country,
      #                               sparse = F))) %>%
      st_difference(.,
                    country) %>%
      st_make_valid()

    cat("Coordinates written in the output files are coordinates of the initial centroid (before intersecting with country).\n")
  }

  grid5km <- grid5km %>%
    dplyr::mutate(area_km2 = units::drop_units(st_area(.)) / 10^6)

  grid5km <- grid5km %>%
    dplyr::mutate(id = 1:n(),
                  lon_cent = cc$lon_cent,
                  lat_cent = cc$lat_cent
    ) %>%
    dplyr::select(id, X, Y, lon_cent, lat_cent) %>%
    st_cast()

  return(grid5km)
}
