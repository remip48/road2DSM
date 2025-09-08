#' Upscale input grid to higher resolution
#'
#' Can be used on the nc file grid created by extract_nc. Resolution must be in the same unit than the input grid, for which the crs will be used to create the new grid.
#'
#' @param grid
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
upscale_grid <- function(grid,
                         resolution) {
  # pkg::fun(stars)
  # pkg::fun(sf)
  # pkg::fun(raster)

  bb <- st_bbox(grid)

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
                            crs = st_crs(grid)$proj4string
  ) %>%
    st_as_stars() %>%
    st_as_sf() %>%
    st_transform(crs = st_crs(grid))

  allgrid <- grid %>%
    group_by() %>%
    dplyr::summarise(do_union = F) %>%
    st_cast("MULTIPOLYGON") %>%
    # st_simplify() %>%
    st_make_valid()

  grid5km <- grid5km %>%
    dplyr::filter(c(st_intersects(.,
                                  allgrid,
                                  sparse = F))) %>%
    # st_intersection(.,
    #                 allgrid) %>%
    st_make_valid()

  grid5km <- grid5km %>%
    dplyr::mutate(area_km2 = units::drop_units(st_area(.)) / 10^6)

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

  grid5km <- grid5km %>%
    dplyr::mutate(id = 1:n(),
                  lon_cent = cc$lon_cent,
                  lat_cent = cc$lat_cent
    ) %>%
    dplyr::select(id, X, Y, lon_cent, lat_cent) %>%
    st_cast()

  return(grid5km)
}
