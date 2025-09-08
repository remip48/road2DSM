#' create soap smooth object for DSM
#'
#' @param n
#' @param target_study_area
#' @param data
#' @param keep_point
#' @param n_try
#'
#' @return
#' @export
#'
#' @examples
create_soap <- function(n,
                        target_study_area, # sf polygon object
                        data, # sf object
                        keep_point = 1, ## simplification of the study area : more polygons and more points in each polygons are in the output,
                        ## higher will be the model uncertainties in some areas of low effort coverage. If effort is high
                        ## in every locations, should be okay
                        ### on my part I had to find a compromise between the buffer and the ratio of simplify to not cut a part of the study area
                        ## and to diminush at maximum the number of polygons & points per polygons
                        n_try = 25  # number of tries to optimize the location of knots, in order to get the lower number of knots possible inside
                        ) {

  cat("Function will place", n, "x", n, "knots on X,Y (ETRS3035) and will retain only those contained within the soap boundaries ",
      paste0("(maximum possible number of knots: ", n*n, ").\n"))
  # lapply(c("sf", "tidyverse", "rmapshaper",
  # ), # pelaDSM also ?
  # pkg::fun, character.only = TRUE
  # )
  # sjmisc

  # ###### define the parameter for the function of creating the soap smooth in etrs 3035
  # n = 20 ## number of knots on X and Y (so same distance between each knot)
  # n_try = 50 ## number of tries to optimize the location of knots, in order to get the higher/lower number of knots possible inside
  # keep_point = 1 ## simplification of the study area : more polygons and more points in each polygons are in the output,
  # ## higher will be the model uncertainties in some areas of low effort coverage. If effort is high
  # ## in every locations, should be okay
  # ### on my part I had to find a compromise between the buffer and the ratio of simplify to not cut a part of the study area
  # ## and to diminush at maximum the number of polygons & points per polygons
  # method = "min"

  ####
  target_study_area <- target_study_area %>%
    st_transform(crs = 3035) %>%
    st_simplify()

  data <- data %>%
    st_transform(crs = 3035) %>%
    st_centroid() %>%
    mutate(x = st_coordinates(.)[,1],
           y = st_coordinates(.)[,2])

  xmin  <- min(data$x)
  xmax  <- max(data$x)
  ymin  <- min(data$y)
  ymax  <- max(data$y)

  res <- ifelse((xmax - xmin) > (ymax - ymin), ymax - ymin, xmax - xmin) / (n-1) ## what is the distance between knots to use

  list_point <- as.data.frame(cbind(x = rep(xmin + (0:ceiling((xmax-xmin)/res)) * res, ceiling((ymax-ymin)/res) + 1),
                                    y = ymin + rep(0:ceiling((ymax-ymin)/res), each = ceiling((xmax-xmin)/res) + 1) * res)) ## initial locations of the knots

  adjust <- cbind(xtry = rep(0:(ceiling(sqrt(n_try))-1) * res / ceiling(sqrt(n_try)), ceiling(sqrt(n_try))),
                  ytry = rep(0:(ceiling(sqrt(n_try))-1), each = ceiling(sqrt(n_try))) * res / ceiling(sqrt(n_try))) ## possible locations for each knots

  which <- apply(adjust, 1, function(try) {
    point_in <- list_point + cbind(xtry = rep(try[1], nrow(list_point)),
                                   ytry = rep(try[2], nrow(list_point))) %>%
      as.data.frame()

    point_in <- point_in %>%
      st_as_sf(coords = c("x", "y"), crs = 3035) %>%
      st_intersection(target_study_area)

    return(nrow(point_in)
    )
  })

  # if (method == "max") {
  #   which <- which.max(which)
  # } else if (method == "min") {
    which <- which.min(which)
  # }

  knots <- list_point + cbind(x = rep(adjust[which, 1], nrow(list_point)),
                              y = rep(adjust[which, 2], nrow(list_point))) %>%
    as.data.frame() ## create results of the location of the knots

  zone_simp <- target_study_area %>%
    ms_simplify(keep = keep_point, keep_shapes=T) ## simplify sf object of the study area

  missing_points <- data %>%
    dplyr::filter(!c(st_intersects(.,
                                   zone_simp,
                                   sparse = FALSE)))

  if (nrow(missing_points) > 0) {
    # Distance from missing points to polygon
    max_dist <- max(st_distance(missing_points, zone_simp)) / 1000

    # Buffer polygon outward by max_dist
    zone_simp <- zone_simp %>%
      st_buffer(units::set_units(max_dist, "km"))
  }

  line <- as.data.frame(zone_simp %>%
                          st_cast("MULTILINESTRING") %>%
                          st_coordinates()) ## extract coordinates

  bnd <- lapply(unique(line$L1), function(l) {
    temp <- line %>%
      filter(L1 == l)

    out <- list(X = temp$X,
                Y = temp$Y)
    return(out)
  })

  knots <- knots %>%
    st_as_sf(coords = c("x", "y"), crs = 3035) %>%
    st_intersection(zone_simp %>%
                      st_buffer(units::set_units(-res / 10, m))) %>%
    st_coordinates(.) %>%
    as.data.frame()

  print(ggplot() +
          geom_sf(data = target_study_area, fill = "alice_blue", color = "midnightblue") +
          geom_sf(data = data, fill = "orange", color = "orange", size = .5) +
          geom_path(data = map_dfr(1:length(bnd), function(b) {
            return(data.frame(X = bnd[[b]]$X,
                              Y = bnd[[b]]$Y) %>%
                     dplyr::mutate(id = b))
          }), aes(x = X, y = Y, group = id), color = "black") +
          geom_point(data = knots, aes(x = X, y = Y), color = "red", size = 2) +
          labs(title = "Soap boundaries and knots, overlaying the study area and data:") +
          theme_bw() +
          theme(panel.background = element_rect(fill = "white")))

  return(list(knots = knots, bnd = bnd))

}
