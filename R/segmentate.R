#' Segmentate GPS data
#'
#' from line-transect surveys.
#'
#' @param effort
#' @param LEG
#' @param LENGTH
#' @param TOLERANCE
#' @param MIN_SEG_LENGTH
#' @param CRS
#' @param bad_conditions
#' @param conditions
#' @param MINROW
#' @param MAXDIST
#' @param COORD_COL
#' @param SEED
#' @param sf_linestring
#' @param plot
#'
#' @return
#' @export
#'
#' @examples
segmentate <- function(effort,
                       LEG="Leg",
                       LENGTH=10,
                       TOLERANCE=0,
                       MIN_SEG_LENGTH=2.5,
                       CRS=4326,
                       bad_conditions=c("ll", "lx", "xl", "xx", "pp", "lp", "pl", "xp", "px",""),
                       conditions="subj",
                       MINROW=1,
                       MAXDIST=2.5,
                       COORD_COL=c('lon','lat'),
                       SEED=1234,
                       sf_linestring=F,
                       plot = T) {

  # pkg::fun(patchwork)

  ##--------------------------------------------------------------------------------------------------------
  ## SCRIPT : transect segmentize
  ##
  ## Authors : Rémi Pigeault
  ## Last update : 2022-02-01
  ## R version 4.0.1 (2020-06-27) -- "See Things Now"
  ## Copyright (C) 2020 The R Foundation for Statistical Computing
  ## Platform: x86_64-w64-mingw32/x64 (64-bit)
  ##--------------------------------------------------------------------------------------------------------

  ## Notes :  TRANSECTCOL are columns used to define the effort portions to be segmented (named Leg)

  # EXTRA is an extra column for TRANSECTCOL to create Legs

  # spatial-lines are created with defined Legs, DATA rows must therefore be in the order collected.

  # Non-effort data during a Leg (e.g. when flying above islands) must be indicated in the TRANSECTCOL or EXTRA !!BEFORE!! the function
  # to not create Line on this non-effort, even if there are no points for this non-effort part. A column "part" can for example be created
  # to specify the different Legs part separated by a non-effort portion.
  # Can be found especially with "distance" column == 0 (new start of effort). See script "TransectSegmentize" for data standardization.

  # MAXDIST divide Legs into sub-Legs according to the distance between parts : if distance between 2 parts > MAXDIST, the two parts are
  # segmented separately.

  # write COORD_COL in the appropriate order : X, or longitude ; Y, or latitude

  # LENGTH, MIN_SEG_LENGTH & MAXDIST must be in kilometers

  # TOLERANCE can be used to add variability within segments length (TOL = .1 add for example 10% of variability)

  # SEED is useful for reproductive research : for a single seed, random length-different segments are always the same and lengths identical.

  # Check if any TRANSECTCOL = "" or NA in dataset firstly, as well as the summary of DISTANCECOL to detect preliminary errors (e.g.
  # unexpected high distance) before the function.

  # sf_linestring: if TRUE, data is treated as an sf object of type linestring, and COORD_COL are ignored, as well as CRS and MINROW. If FALSE, segmented as individual points of gps.

  # parallelize: to use several cores, via doParallel

  # core_to_use: either numeric values for the number of core to use, or by default NA (number of cores - 1). If an error is returned due to memory, decrease the number of core used

  # return the same dataframe than initially but with Leg, SubLeg (differentiating parts of a same Leg separated by a gap), SubLeg (Legs parts being separated by a gap > MAXDIST), Seg (segment of the SubLeg) and label (idientifying Leg, SubLeg and Seg in the same time)
  # if sf_linestring = T, return also distance (distance in km of the segment)

  .segmentate<-function(DATA,LEG="Leg",LENGTH=5,TOLERANCE=0,MIN_SEG_LENGTH=2.5,CRS=4326,bad_conditions=c("ll", "lx", "xl", "xx", "pp", "lp", "pl", "xp", "px",""),
                        conditions="subj",MINROW=1,MAXDIST=2.5,COORD_COL=c('lon','lat'),SEED=1234,sf_linestring=F){

    # load_packages <- function(p) {
    #   if (p %in% rownames(installed.packages())) {
    #     # pkg::fun(p, character.only = TRUE)
    #   } else {
    #     install.packages(p)
    #     # pkg::fun(p, character.only = TRUE)
    #   }
    # }
    # lapply(c("sf", "dplyr", "lwgeom", "purrr"), load_packages)
    # if (parallelize) {load_packages("doParallel")}

    cat("Data is ordered by Leg/time\n")

    ## Leg = initial Legs to be segmented
    ## SubLeg = section of Legs separated by bad conditions or distance = 0 in data, defined by "Subpart" variable. Allow to separate parts of Legs separated by a long distance
    ## SubLeg = Leg taken into account gap > MAXDIST, calculated from SubLeg. Each SubLeg is segmented individually.
    ## Seg : segments on each SubLeg.

    # for_SubLeg <- DATA[, TRANSECTCOL] ## create column for SubLeg
    # if (any(class(for_SubLeg) == "sf")) {
    #   for_SubLeg <- for_SubLeg %>%
    #     st_drop_geometry()
    # }
    # Legs <- apply(for_SubLeg, 1, function(t) {
    #   x <- paste(t, collapse = "_")
    #   return(x)
    # }) %>%
    #   unname()

    if ("Leg" %in% colnames(DATA) & LEG != "Leg") {
      DATA <- DATA %>%
        rename(Leg_init = Leg)

      DATA <- DATA %>%
        mutate(Leg = DATA %>%
                 pull(LEG))
    } else if (!("Leg" %in% colnames(DATA)) & LEG != "Leg") {
      DATA <- DATA %>%
        mutate(Leg = DATA %>%
                 pull(LEG))
    }

    DATA <- DATA %>%
      arrange(Leg)

    to_remove <- which(DATA %>%
                         pull(conditions) %in% bad_conditions)
    if (length(to_remove) > 0) {
      tot <- DATA[-to_remove, ] %>%
        pull(Leg) %>%
        unique()
    } else {
      tot <- DATA %>%
        pull(Leg) %>%
        unique()
    }

    options(dplyr.summarise.inform = F)

    if (!sf_linestring) { ## segment the dataframe by creating a sf linestring object and assigning points to closest segment

      if ("sf" %in% class(DATA)) {
        stop("Enter a non-sf object or select sf_linestring = TRUE")
      }

      cat("\nEstimating Leg parts (non-interrupted effort) \n")
      pb <- txtProgressBar(min = 0, max = length(tot), style = 3)

      DATA_legpart <- map(tot, function(l) { ## for optimisation, applied to each Leg
        # print(l)
        setTxtProgressBar(pb, match(l, tot))

        DATA_filt_leg <- DATA %>%
          filter(Leg == l)

        DATA_filt_leg <- DATA_filt_leg %>%
          dplyr::mutate(idd_data = 1:n())

        to_remove <- which(DATA_filt_leg %>%
                             pull(conditions) %in% bad_conditions)
        if (length(to_remove) > 0) {
          DATA_filt_leg <- DATA_filt_leg[-to_remove,]
        }

        DATA_filt_leg[, c("tlon", "tlat")] <- DATA_filt_leg[, COORD_COL]

        ## identify SubLegs
        ll <- DATA_filt_leg$idd_data[-1] - DATA_filt_leg$idd_data[-nrow(DATA_filt_leg)]
        ll <- unique(c(1,
                       which(ll > 1) + 1,
                       nrow(DATA_filt_leg) + 1))
        # ll <- unique(c(1,
        #                which(DATA_filt_leg[, conditions] %in% bad_conditions),
        #                nrow(DATA_filt_leg) + 1))

        if (length(ll) > 1) {
          DATA_filt_leg <- map_dfr(2:length(ll), function(i) {
            id <- ll[i-1]:(ll[i] - 1)
            out <- DATA_filt_leg[id,]
            out$LegPart <- rep(paste0("P", formatC(i - 1, width = 4, format = "d", flag = "0")), length(id))
            return(out)
          })
        } else if (length(l) == 1) {
          DATA_filt_leg$LegPart <- paste0("P", formatC(1, width = 4, format = "d", flag = "0"))
        }

        DATA_filt_leg <- map_dfr(unique(DATA_filt_leg$LegPart), function(tp) {

          temp <- DATA_filt_leg %>%
            filter(LegPart == tp) %>%
            st_as_sf(coords = COORD_COL, crs = CRS) %>%
            st_transform(crs = 3035)

          if (nrow(temp) > 1) {
            tempsf <- st_distance(temp[1:(nrow(temp) - 1),], temp[2:nrow(temp),], by_element = T) %>%
              units::drop_units()

            ll <- c(1,
                    which(tempsf > MAXDIST*1000) + 1,
                    nrow(temp) + 1)

            temp <- map_dfr(2:length(ll), function(i) {
              id <- ll[i-1]:(ll[i] - 1)
              out <- temp[id,]
              out$LegPart <- paste(out$LegPart, i-1, sep = "_")

              return(out)
            })
          } else {
            temp$LegPart[1] <- paste(temp$LegPart[1], 1, sep = "_")
          }

          return(temp)

        })

        return(DATA_filt_leg)
      })

      cat("\nSegmenting sub-legs: parts separated by distance >", MAXDIST, "km (MAXDIST)", "\n")
      DATA_sf_line <- map(DATA_legpart, function(DATA_filt_leg) {
        ## segmentize the Legs, distances must be provided in meters
        setTxtProgressBar(pb, match(DATA_filt_leg$Leg[1], tot))
        # print(unique(DATA_filt_leg$Leg))

        DATA_out <- .segmentize(sldf = DATA_filt_leg, seed = SEED, TOLERANCE = TOLERANCE,
                                step = LENGTH*1000, MINLENGTH = MIN_SEG_LENGTH*1000, MAX_DIST = MAXDIST*1000,
                                sf_linestring = sf_linestring, MINROW = MINROW)

        DATA_out$init_data <- DATA_filt_leg
        return(DATA_out)
      })

      cat("\nAssigning data points to the segments\n")
      DATA_seg <- map_dfr(DATA_sf_line, function(DATA_out) {

        DATA_filt_leg <- DATA_out$init_data
        DATA_line_seg <- DATA_out$dataframe
        missing_data <- DATA_out$missing

        l <- DATA_filt_leg$Leg[1]
        setTxtProgressBar(pb, match(l, tot))

        if (!is.null(DATA_line_seg) & nrow(DATA_line_seg) > 0 & length(DATA_line_seg) > 0) {

          temp_data <- DATA_filt_leg %>%
            filter(Leg == l & !(idd_data %in% missing_data$idd_data))

          min_seg <- temp_data %>%  ## assigning now data points to the closest Segment (if not further than MAXDIST)
            st_distance(DATA_line_seg) %>%
            units::drop_units() %>%
            as.matrix()

          seg_match <- apply(min_seg, 1, function(x) {
            out <- which.min(x)
            if (x[out]/1000 > MAXDIST) {out = integer(0)}
            return(ifelse(length(out) == 0, NA, out))
          })

          temp_data <- temp_data %>%
            st_transform(crs = CRS)

          missing_data <- missing_data %>%
            dplyr::mutate(Seg = NA, SubLeg = NA)

          out <- temp_data %>%
            dplyr::mutate(Seg = DATA_line_seg$Seg[seg_match],
                          SubLeg = DATA_line_seg$SubLeg[seg_match],
                          reason = "KEEP",
                          distance_to_closest_LegPart = 0) %>%
            st_drop_geometry() %>%
            rbind(missing_data)
        } else {
          out <- missing_data %>%
            dplyr::mutate(Seg = NA, SubLeg = NA)
        }

        out[, COORD_COL] = out[, c("tlon", "tlat")]

        out <- out %>%
          dplyr::select(!c("tlon","tlat","idd_data"))

        return(out)
      })

    } else if (sf_linestring) {

      if (!("sf" %in% class(DATA))) {
        stop("Enter a sf object or select sf_linestring = FALSE")
      }

      if (length(to_remove) > 0) {
        DATA <- DATA[-to_remove,]
      }

      cat("\nSegmenting\n")
      pb <- txtProgressBar(min = 0, max = length(tot), style = 3)

      DATA_seg <- map_dfr(tot, function(l) {
        # print(l)
        setTxtProgressBar(pb, match(l, tot))

        DATA_filt_leg <- DATA %>%
          filter(Leg == l)

        DATA_filt_leg <- DATA_filt_leg %>%
          dplyr::mutate(LegPart = paste0("P", formatC(1:n(), flag = "0", width = 4)))

        ## segmentize the Legs, distances must be provided in meters
        DATA_out <- .segmentize(sldf = DATA_filt_leg, seed = SEED, TOLERANCE = TOLERANCE,
                                step = LENGTH*1000, MINLENGTH = MIN_SEG_LENGTH*1000, MAX_DIST = MAXDIST*1000,
                                sf_linestring = sf_linestring)

        out <- DATA_out$dataframe
        missing_data <- DATA_out$missing %>%
          dplyr::mutate(Seg = NA, SubLeg = NA, distance = NA)

        if (!is.null(out)) {
          out <- out %>%
            dplyr::mutate(distance = units::drop_units(st_length(out) / 1000),
                          reason = "KEEP") %>%
            rbind(missing_data)
        } else {
          out <- missing_data %>%
            dplyr::mutate(distance = 0)
        }

        return(out)
      })
    }

    errors <- DATA_seg %>%
      filter(reason != "KEEP")

    DATA_seg <- DATA_seg %>%
      filter(reason == "KEEP") %>%
      dplyr::select(-reason)

    if (!sf_linestring) {
      DATA_seg <- DATA_seg %>%
        dplyr::select(-distance_to_closest_LegPart)
    }

    ## defining label of the Segment (Leg + SubLeg + Seg)
    DATA_seg <- DATA_seg %>%
      dplyr::mutate(label = paste(Leg, "_P", SubLeg, "_S",formatC(Seg, width = 4, format = "d", flag = "0"),sep=''))

    cat("\nAggregating data\n\n")

    logM <- c()
    for (reas in unique(errors$reason)) {
      temp <- errors %>%
        filter(reason == reas)

      warning(nrow(temp), " data discarded due to: ", reas, ".",
              ifelse(reas == "n_point LegPart < MINROW & dist_other_LegPart > MAXDIST",
                     paste("\nDistance to closest LegPart:",
                           paste(temp %>%
                                   arrange(distance_to_closest_LegPart) %>%
                                   pull(distance_to_closest_LegPart) %>%
                                   round(., 0) %>%
                                   unique(), collapse = ", "),
                           "(meters). Increase MAXDIST to integrate them."),
                     ""))

      logM <- c(logM, paste0(nrow(temp), " data discarded due to: ", reas, ".",
                             ifelse(reas == "n_point LegPart < MINROW & dist_other_LegPart > MAXDIST",
                                    paste("Distance to closest LegPart:",
                                          paste(temp %>%
                                                  arrange(distance_to_closest_LegPart) %>%
                                                  pull(distance_to_closest_LegPart) %>%
                                                  round(., 0) %>%
                                                  unique(), collapse = ", "),
                                          "(meters). Increase MAXDIST to integrate them."),
                                    "")))
    }

    if (sf_linestring) {
      out <- list(data_segmented = DATA_seg,
                  data_discarded = errors,
                  log_message = logM)
    } else {
      DATA_sf <- map_dfr(DATA_sf_line, function(d) {

        out <- d$dataframe

        if (!is.null(out) & nrow(out) > 0 & length(out) > 0) {
          out <- out %>%
            st_transform(crs = CRS) %>%
            dplyr::mutate(label = paste(Leg, "_P", SubLeg, "_S",formatC(Seg, width = 4, format = "d", flag = "0"),sep='')) %>%
            dplyr::select(Leg, LegPart, SubLeg, Seg, label)
        } else {
          out <- NULL
        }

        return(out)
      })

      out <- list(data_segmented = DATA_seg,
                  sf_line_segmented = DATA_sf,
                  data_discarded = errors,
                  log_message = logM)
    }

    return(out)
  }

  ## to segmentize an sf linestring object
  .segmentize <- function(sldf, # spatial line dataframe as sf object to segmentate : apply this to individual lines to be segmentate
                          step = 5e3,
                          MINLENGTH = 2.5e3,
                          seed = 1234,
                          TOLERANCE = 0,
                          MAX_DIST = 5e3,
                          # COORD_COL = c("lon","lat"),
                          # CRS = 4326,
                          MINROW = 1,
                          sf_linestring = sf_linestring
  ) {
    # for reproducibility
    set.seed(seed)
    # function to create relative step lengths of a given line with variability around the mean step length
    partition <- function(l, # length of the leg
                          step, # segment length
                          tolerance = 0, # variability in %, 0% by default
                          minlength
    ) {
      ### how many segments?
      n <- ifelse((l/step - floor(l / step)) < minlength / step,
                  floor(l / step),
                  ceiling(l / step))
      if(n %in% c(0:1)) {
        r <- c(0, 1)
      } else {
        ### store results
        r <- numeric(n)
        ### use logit scale
        logit <- function(p) { log(p / (1 - p)) }
        m <- logit(step / l)
        s <- (logit(min(0.99, (1 + tolerance) * step / l)) - logit(max(0.01, (1 - tolerance) * step / l))) / 4
        ### draws the steps
        r <- plogis(m + s * rnorm(n))
        ### randomize
        prop <- sample.int(n, size = 1)
        ### ensure that r sums to 1
        r[prop] <- 1 - sum(r[-prop])

        if (any(r <= minlength/l) & length(r) > 2) {
          min <- r[base::which.min(r)]
          r <- r[-base::which.min(r)]
          r[base::which.min(r)] <- r[base::which.min(r)] + min
        }
        if (any(r > (step + minlength)/l) & length(r) > 1) {
          rmax <- base::which.max(r)
          max <- r[rmax]
          if (rmax == 1) {
            r <- c(max / 2, max / 2, r[(base::which.max(r) + 1):length(r)])
          } else if (rmax == length(r)) {
            r <- c(r[1:(base::which.max(r) - 1)], max / 2, max / 2)
          } else {
            r <- c(r[1:(base::which.max(r) - 1)], max / 2, max / 2, r[(base::which.max(r) + 1):length(r)])
          }
        }
        r <- c(0, cumsum(r))
        r[length(r)] <- 1
      }
      ### wrap-up
      return(r)
    }

    ## transform in spatial object
    if (!sf_linestring) {

      sd <- sldf %>%
        dplyr::mutate(LegParti = LegPart)

      # tt <- ggplot() +
      #   geom_sf(data = out %>%
      #             st_transform(crs = 4326), aes(color = factor(Seg)), linewidth = 2) +
      #   scale_color_manual(values = rep(viridis::viridis(2), ceiling(nrow(out)))) +
      #   geom_point(data = sldf, aes(x = tlon, y = tlat), size = .1, color = "red")
      #
      # png("test.png")
      # print(tt)
      # dev.off()

      inf_minrow <- names(which(table(sd$LegPart) <= MINROW))
      if (length(inf_minrow) > 0 & length(unique(sd$LegPart)) == 1) {

        out_infminrow <- sd %>%
          st_drop_geometry() %>%
          dplyr::mutate(reason = "n_point Leg < MINROW",
                        distance_to_closest_LegPart = NA) %>%
          dplyr::select(-LegParti)

        sd <- sd %>%  ## correct the LegPart with not enough rows
          dplyr::select(-c(LegParti, tlon, tlat))

        sdf <- sd %>% ## create sf lines for each legs
          group_by(LegPart) %>%
          dplyr::summarise(Leg = unique(Leg),
                           do_union = FALSE) %>%
          st_cast("LINESTRING")

        sdf <- sdf[0,]

      } else {

        if (length(inf_minrow) > 0) { ## to assign to closest LegPart if number of points < MINROW and dist < MAXDIST
          sd <- .assign_LegPart_infminrow(sd, inf_minrow, MINROW, MAX_DIST)
        } else {
          sd$distance_to_closest_LegPart <- 0
        }

        out_infminrow <- sd %>%
          st_drop_geometry() %>%
          filter(LegPart == "NA_LegPart") %>%
          dplyr::mutate(LegPart = LegParti,
                        reason = ifelse(is.na(distance_to_closest_LegPart),
                                        "n_point LegPart < MINROW & next+previous LegPart have (n_point < MINROW & dist_other_LegPart > MAXDIST)",
                                        "n_point LegPart < MINROW & dist_other_LegPart > MAXDIST")) %>%
          dplyr::select(-LegParti)

        sd <- sd %>%  ## correct the LegPart with not enough rows
          filter(LegPart != "NA_LegPart") %>%
          dplyr::select(-c(LegParti, tlon, tlat))

        sdf <- sd %>% ## create sf lines for each legs
          group_by(LegPart) %>%
          dplyr::summarise(Leg = unique(Leg),
                           do_union = FALSE) %>%
          st_cast("LINESTRING")
      }

      ## filter LegPart with Length = 0 meters
      # length_sdf <- units::drop_units(st_length(sdf))
      #
      # sdf <- sdf %>%
      #   filter(length_sdf != 0)

    } else {
      sd <- sldf %>%
        st_transform(crs = 3035)

      ## filter LegPart with Length = 0 meters
      length_sdf <- units::drop_units(st_length(sd))

      out_infminrow <- sd %>%
        filter(length_sdf == 0) %>%
        dplyr::mutate(reason = "Length_LegPart = 0 meters")

      sd <- sd %>%
        filter(length_sdf != 0)

      sdf <- sd
    }

    it = 1 ## create SubLegs, separated by distance > MAXDIST and merging the LegParts closer
    sub_leg = rep(1, nrow(sdf)) ## carefull here, must be applied each individual leg independtly, Otheriwse, to be changed here
    if(nrow(sdf) > 1) {
      for (i in 2:length(unique(sdf$LegPart))) {
        if (units::drop_units(st_distance(sd[first(which(sd$LegPart == unique(sdf$LegPart)[i])), ],
                                          sd[last(which(sd$LegPart == unique(sdf$LegPart)[i-1])), ])) > MAX_DIST) {
          it = it+1}
        sub_leg[i] = it
      }
    }
    sdf <- sdf %>%
      dplyr::mutate(SubLeg = sub_leg)

    Legseg <- map_dfr(unique(sdf$SubLeg), function(s) { ## segmenting each SubLeg
      temp <- sdf %>%
        filter(SubLeg == s)

      ## create proportions of SubLeg where cut
      length_temp <- units::drop_units(st_length(temp))
      seg <- partition(l = sum(length_temp), step = step, tolerance = TOLERANCE, minlength = MINLENGTH)
      seg[length(seg)] <- 1

      if (nrow(temp) == 1) {
        out <- map_dfr(2:length(seg), function(x) { ## temp : line to partition, x : vector with the relative step lengths
          return(temp %>%
                   st_linesubstring(seg[x-1], seg[x]))
        }) %>%
          dplyr::mutate(Seg = 1:n())
      } else {
        prop <- c(0,cumsum(length_temp / sum(length_temp)))
        prop[length(prop)] <- 1
        out <- map_dfr(1:nrow(temp), function(i) { ## Segmenting for each line (LegPart) of the LegPart
          temp1 <- temp[i,]

          n_int <- seg[which(seg > prop[i] & seg < prop[i+1])]
          n_inf <- last(which(seg <= prop[i]))

          SubLeg_seg <- na.omit(c(prop[i], n_int, prop[i+1]) - prop[i]) / (prop[i+1] - prop[i])
          SubLeg_seg[length(SubLeg_seg)] <- 1

          seg_out <- map_dfr(2:length(SubLeg_seg), function(x) { ## temp : line to partition, x : vector with the relative step lengths
            return(temp1 %>%
                     st_linesubstring(SubLeg_seg[x-1], SubLeg_seg[x]))
          }) %>%
            dplyr::mutate(Seg = c(n_inf, which(seg > prop[i] & seg < prop[i+1])))

          return(seg_out)

        })
      }

      return(out)
    })

    return(list(dataframe = Legseg,
                missing = out_infminrow))
  }

  .assign_LegPart_infminrow <- function(sd, inf_minrow, MINROW, MAX_DIST) { ## to rename LegPart with n_point < MINROW and distance to one successive LegPart < MAXDIST

    sd <- sd %>%
      mutate(iLegPart = LegPart) %>%
      .sub_function_assign_LegPart(., inf_minrow, MINROW, MAX_DIST)

    inf_minrow <- names(which(table(sd$LegPart) <= MINROW))
    if (length(inf_minrow) > 0) { ## to replace be sure that replacement of LegPart is finished (if case consective LegPart have both been replaced)
      sd <- .sub_function_assign_LegPart(sd, inf_minrow, MINROW, MAX_DIST)
    }

    inf_minrow <- names(which(table(sd$LegPart[sd$LegPart != "NA_LegPart"]) <= MINROW))

    if (length(inf_minrow) > 0) {
      to_rename <- sd %>%
        filter(LegPart %in% inf_minrow & iLegPart == LegPart) %>%
        pull(LegPart)

      if (length(to_rename) > 0) {
        for (i in 1:(length(to_rename)/2)) {
          sd <- sd %>%
            dplyr::mutate(LegPart = ifelse(LegPart == to_rename[i*2], to_rename[i*2 - 1], LegPart))
        }
      }
    }

    sd <- sd %>%
      dplyr::select(-iLegPart)

    return(sd)
  }

  .sub_function_assign_LegPart <- function(sd, inf_minrow, MINROW, MAX_DIST) { ## subfunction for .assign_LegPart_infminrow

    list_LegPart <- sd %>%
      pull(LegPart) %>%
      unique()

    new_LegPart <- map_dfr(inf_minrow, function(lp) {
      if (lp == first(list_LegPart)) {
        prev_next_LegPart <- sd %>%
          filter(LegPart == list_LegPart[2])
        prev_next_LegPart <- prev_next_LegPart[1, ]
      } else if (lp == last(list_LegPart)) {
        prev_next_LegPart <- sd %>%
          filter(LegPart == list_LegPart[length(list_LegPart) - 1])
        prev_next_LegPart <- prev_next_LegPart[nrow(prev_next_LegPart), ]
      } else {
        prev_LegPart <- sd %>%
          filter(LegPart == list_LegPart[which(list_LegPart == lp) - 1])
        prev_LegPart <- prev_LegPart[nrow(prev_LegPart), ]

        next_LegPart <- sd %>%
          filter(LegPart == list_LegPart[which(list_LegPart == lp) + 1])
        next_LegPart <- next_LegPart[1, ]

        prev_next_LegPart <- rbind(prev_LegPart, next_LegPart)
      }

      dist <- st_distance(sd %>%
                            filter(LegPart == lp),
                          prev_next_LegPart) %>%
        as.data.frame() %>%
        units::drop_units()

      if (any(dist <= MAX_DIST)) {
        out <- data.frame(LegPart = prev_next_LegPart$LegPart[which.min(dist)],
                          distance_to_closest_LegPart = min(dist))
        return(out)
      } else {
        out <- data.frame(LegPart = "NA_LegPart",
                          distance_to_closest_LegPart = min(dist))
        return(out)
      }
    })

    sd <- sd %>%
      # filter(LegPart %in% inf_minrow) %>%
      dplyr::mutate(LegPart = ifelse(LegPart %in% inf_minrow, new_LegPart$LegPart[match(LegPart, inf_minrow)], LegPart),
                    distance_to_closest_LegPart = ifelse(LegPart == "NA_LegPart",
                                                         new_LegPart$distance_to_closest_LegPart[match(LegPart, inf_minrow)],
                                                         0))

    return(sd)
  }

  cat("Leg is created based on unique Transect (EFF_Track) & date\n\n")
  effort <- effort %>%
    dplyr::mutate(
      Leg = paste(EFF_Track, date, sep = "_")
    ) %>%
    arrange(Leg, datetime)

  effort_seg <- .segmentate(DATA = effort,
                            LEG=LEG,
                            LENGTH=LENGTH,
                            TOLERANCE=TOLERANCE,
                            MIN_SEG_LENGTH=MIN_SEG_LENGTH,
                            CRS=CRS,
                            bad_conditions=bad_conditions,
                            conditions=conditions,
                            MINROW=MINROW,
                            MAXDIST=MAXDIST,
                            COORD_COL=COORD_COL,
                            SEED=SEED,
                            sf_linestring=sf_linestring)

  effort_seg_ref <- effort_seg

  effort_sf_line <- effort_seg_ref$sf_line_segmented

  effort_seg <- effort_seg_ref$data_segmented

  # summary(effort_seg$distance)

  ## set to 0 the distances between GPS points that are too large
  effort_seg_corr <- effort_seg %>%
    mutate(effort_km2 = ifelse(distance > MAXDIST, 0, effort_km2),
           distance = ifelse(distance > MAXDIST, 0, distance))

  effort_seg_corr <- effort_seg_corr %>%
    arrange(Leg, datetime)

  ## group the segments by label (1 label = 1 single segment).
  effort <- effort_sf_line %>%
    group_by(label) %>%
    dplyr::summarise(#Leg = unique(Leg),
      # LegPart = paste(unique(LegPart), collapse = ", "), ## LegPart, SugLeg and Seg are columns used internally for segmenting.
      # SubLeg = unique(SubLeg),
      # Seg = unique(Seg),
      do_union = F) %>%
    st_cast("MULTILINESTRING")

  # effort <- effort %>%
  #   filter(units::drop_units(st_length(.)) > 0)

  coords_sf <- effort %>%
    st_transform(crs = 3035) %>%
    st_centroid() %>%
    st_transform(crs = 4326) %>%
    st_coordinates()

  effort <- effort %>%
    mutate(lon_cent = coords_sf[,1],
           lat_cent = coords_sf[,2],
           sf_length = units::drop_units(st_length(.)) / 1000)

  # summary(effort$sf_length)
  # hist(effort$sf_length)

  # any(is.na(effort_seg_corr$sub_sig[effort_seg_corr$species == "ppho"])) # if one sub_sig is NA for ppho, it must be corrected!

  ## Keep only sightings made in good and moderate conditions.
  all_bad_conditions <- do.call("c", lapply(bad_conditions, function(b) {
    strsplit(b, "")[[1]]
  }))

  effort_seg_corr <- effort_seg_corr %>%
    dplyr::mutate(calves = ifelse(!(sub_sig %in% all_bad_conditions),
                                  calves,
                                  0),
                  podsize = ifelse(!(sub_sig %in% all_bad_conditions),
                                   podsize,
                                   0),
                  SIG_number = ifelse(!(sub_sig %in% all_bad_conditions),
                                      SIG_number,
                                      0))

  info <- effort_seg_corr %>%
    group_by(label) %>%
    dplyr::summarise(distance = sum(distance, na.rm=T),
                     effort_km2 = sum(effort_km2, na.rm=T),
                     n_point = n(),
                     date = unique(date, na.rm=T),
                     year = unique(year(datetime), na.rm=T),
                     month = unique(month(datetime), na.rm=T),
                     day = unique(day(datetime), na.rm=T),
                     surveyID = paste(unique(surveyID, na.rm=T), collapse = "_"),
                     transectID = unique(EFF_Track, na.rm=T),
                     time_start = first(datetime),
                     time_end = last(datetime),
                     sea = paste(unique(sea, na.rm=T), collapse = ","),
                     com_EFF = paste(unique(com_EFF, na.rm=T), collapse = "_"),
                     project = paste(unique(project, na.rm=T), collapse = "_"),
                     source = "ITAW") %>%
    arrange(label)

  effort <- effort %>%
    left_join(info, by = "label")

  if (plot) {
    p1 <- ggplot() +
      scale_color_manual(values = rep(viridis::viridis(3), ceiling(length(effort_sf_line$label) / 3))) +
      # geom_point(data = effort_seg_ref$data_segmented, aes(x = lon, y = lat), color = "red", size = 2) +
      # geom_point(data = effort, aes(x = lon, y = lat), color = "purple") +
      geom_sf(data = effort_sf_line, aes(color = label), show.legend = F, linewidth = 1.5) +
      labs(title = paste0(LENGTH, " km-target segments:")) +
      theme_bw() +
      theme(panel.background = element_rect(fill = "aliceblue"))

    p2 <- ggplot() +
      geom_histogram(data = effort, aes(x = distance), fill = "midnightblue") +
      labs(title = "Histogram of length from segmented Legs:") +
      scale_x_continuous(name = "Length (km)")
      theme_bw()

    print((p1 * p2))
  }

  for (sp in unique(na.omit(effort_seg_corr$species))) {
    if (!(str_remove_all(sp, fixed(" ")) %in% "")) {
      obs <- effort_seg_corr %>%
        filter(species == sp)

      obs <- obs %>%
        group_by(label) %>%
        dplyr::summarise(!!paste0("n_", sp) := sum(SIG_number, na.rm = T),
                         !!paste0("n_group_", sp) := length(which(SIG_number > 0)),
                         !!paste0("n_calves_", sp) := sum(calves, na.rm = T))

      effort <- effort %>%
        left_join(obs, by = "label") %>%
        as.data.frame() %>%
        dplyr::mutate(!!paste0("n_", sp) := ifelse(is.na(get(paste0("n_", sp))), 0, get(paste0("n_", sp))),
                      !!paste0("n_group_", sp) := ifelse(is.na(get(paste0("n_group_", sp))), 0, get(paste0("n_group_", sp))),
                      !!paste0("n_calves_", sp) := ifelse(is.na(get(paste0("n_calves_", sp))), 0, get(paste0("n_calves_", sp))))
    }
  }

  effort <- effort %>%
    st_cast("MULTILINESTRING") %>%
    st_sf()

  return(effort)

}
