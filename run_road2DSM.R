library(road2DSM)
library(dplyr)

rm(list = ls())
gc()

setwd(WorkDir <- "C:/Documents/run_road2DSM")

### aggregation of the GPS effort data

### segmentate effort
segments <- segmentate(effort = data,
                       LENGTH=10, # length in kilometers of segments
                       TOLERANCE=.1, # tolerance in segment lengths
                       MIN_SEG_LENGTH=5, # minimum allowed segment length, if possible to achieve
                       CRS=4326,
                       bad_conditions=c("ll", "xx", "pp", "lx", "xl", "lp", "pl", "xp", "px",""),
                       conditions="subj", # column where conditions are found
                       MINROW=1,
                       MAXDIST=1, # if Legs are distant by more than 1 km, treat them separately as different transects
                       COORD_COL=c('lon','lat'),
                       SEED=1234,
                       sf_linestring = F)

### let's say you want to predict on every day from the first to the last annual surveyed date:
library(timeDate)
list_dates <- segments %>%
  arrange(date) %>%
  group_by(year) %>%
  dplyr::reframe(dates = as.character(timeSequence(from = min(date), to = max(date), by = "day"))) %>%
  pull(date)

if (!dir.exists(file.path(WorkDir, "NC_files_directory"))) {
  dir.create(file.path(WorkDir, "NC_files_directory"))
}
### extract variables from nc files: let's say you want to extract on a period 2016-2024 which is actually covered by 2 NC files. And that you have Physiographic
# and biological NC files set. On the Physiographic, you want to extract the center, SDtime and SDspace, while only the center and SDtime on the biological NC files.
# Let's say you also want to calculate EKE based on the eastern_current and western_current:
extract_nc(nc.path = file.path(WorkDir, "NC_files_directory"), # all your nc files must be in this directory
            list_variable = c("SST", "eastern_current", "western_current", "NPPV"),
            nc_files = data.frame(file = c("NC_file_PHY_period_1_2016-2020.nc",
                                           "NC_file_PHY_period_2_2021-2024.nc",
                                           "NC_file_BIO_period_1_2016-2020.nc",
                                           "NC_file_BIO_period_2_2021-2024.nc"),
                                  date_start = c("2016-01-01",
                                                 "2021-01-01",
                                                 "2016-01-01",
                                                 "2021-01-01"),
                                  date_end = c("2020-12-31",
                                               "2024-12-31",
                                               "2020-12-31",
                                               "2024-12-31"),
                                  type = c("center,SDtime,SDspace",
                                           "center,SDtime,SDspace",
                                           "center,SDtime",
                                           "center,SDtime"),
                                  expr = c("EKE = sqrt(eastern_current^2 + western_current^2)",
                                           "EKE = sqrt(eastern_current^2 + western_current^2)")),
            all_pixel.radius = c(1, 2), # to extract the SD space over 1 and 2 pixels
            all_time.period = c(1, 8, 15), # to extract the center, SDspace and SDtime over 1 day, 8 days and 15 days. SDtime_1day will be NA as it is over only 1 value, but will be filtered later.
            dates = list_dates,
            n_cores = 8,
            outfile = file.path(WorkDir, "log.txt"))

### let's create a prediction grid with the same extent than the nc files you have provided:
list_grids <- list.files(file.path(WorkDir, "NC_files_directory"))
list_grids <- list_grids[str_detect(list_grids, fixed(".shp")) & str_detect(list_grids, fixed("file_set_"))]

grid_NC <- map_dfr(list_grids, function(f) {
  read_sf(file.path(WorkDir, "NC_files_directory", f)) %>%
    dplyr::mutate(file = f) %>%
    dplyr::select(file)
}) %>%
  st_transform(crs = 3035) %>%
  group_by() %>%
  dplyr::summarise(do_union = F) %>%
  st_cast("MULTIPOLYGON")

grid5km <- upscale_grid(grid_NC,
                         resolution = 5000)

### Let's crop your effort to the area covered by your nc files:
segments <- crop_effort(segments,
                        file_set_directory = file.path(WorkDir, "NC_files_directory"),
                        keep = "all_that_intersect", # or "only_fully_in"
                        study_area = NULL # additional sf file to filter
)

### now you want to extract all variables on the prediction grids (except eastern and western current as now we want only the EKE) and store them in a new folder
if (!dir.exists(file.path(WorkDir, "prediction_grids"))) {
  dir.create(file.path(WorkDir, "prediction_grids"))
}

## load bathy and coastline
bathy <- terra::raster(bathy_raster_file)
coatline <- read_sf(coastline_sf_file)

extract_grid(grid = grid5km,
             variable = c("SST", "EKE", "NPPV"),
             file_set_directory = file.path(WorkDir, "NC_files_directory"),
             writing_directory = file.path(WorkDir, "prediction_grids"),
             rasters = list(bathy = bathy),
             distance_to = list(distance_to_coast = coatline),
             dates = list_dates,
             n_cores = 5,
             outfile = file.path(WorkDir, "log.txt")
)

########################
### OPTIONAL STEP!!! Let's say you are also interested in the SDspace but when considering a larger resolution than the nc files (let's say 2x2km): now you want the SDspace considering cells of 5 km to get more mesoscale patterns:
# this step is definitevly not necessary but interesting for exploration of data.
if (!dir.exists(file.path(WorkDir, "prediction_grids", "with_upscale_SDspace", "file_set_1"))) {
  dir.create(file.path(WorkDir, "prediction_grids", "with_upscale_SDspace", "file_set_1"))
}
upscale_SDspace(file_directory = file.path(WorkDir, "prediction_grids"),
               writing_directory = file.path(WorkDir, "prediction_grids", "with_upscale_SDspace", "file_set_1"), # we are reproducing here the folder structure created by extract_nc
               # this will be useful to extract covariates on effort after
               all_pixel.radius = c(1, 2),
               n_cores = 8)

# to reproduce exactly the folder structure created by extract_nc, we need to have the corresponding sf grid of file_set_1 in the parent directory:
write_sf(grid5km,
         file.path(WorkDir, "prediction_grids", "with_upscale_SDspace", "file_set_1.shp"))
########################

# If we used upscale_SDspace:
predgrid_directory <- file.path(WorkDir, "prediction_grids", "with_upscale_SDspace", "file_set_1")
file_set_directory <- file.path(WorkDir, "prediction_grids", "with_upscale_SDspace")
# Otherwise:
predgrid_directory <- file.path(WorkDir, "prediction_grids")
file_set_directory <- file.path(WorkDir, "NC_files_directory")

### Now we can extract all the variables on the effort:
# first we create a buffer around the segments
effort <- segments %>%
  st_transform(crs = 3035) %>%
  st_buffer(units::set_units(2.5, km))

### and now we extract:
effort_extracted <- extract_effort(effort = effort,
                                   variable = c("SST", "EKE", "NPPV"),
                                   rasters = list(bathy = bathy),
                                   distance_to = list(distance_to_coast = coatline),
                                   n_cores = 10,
                                   file_set_directory = file_set_directory)

### now that all the files are ready, we can fit models:
calibration_data <- effort_extracted %>%
  st_drop_geometry()

run_models <- run_all_DSM(segdata_obs = calibration_data,
                          response = "n_ppho",
                          predictors = c("SST.center_1d", "SST.SDspace.mean_1d", "SST.SDtime_8d",
                                         "EKE.center_1d", "EKE.SDspace_1d",
                                         "NPPV.center_1d", "NPPV.SDtime_8d"),
                          nb_max_pred = 3,
                          nb_min_pred = 5,
                          fit_all_once = T,
                          ncores = 8,
                          complexity = 10,
                          use_loo = T,
                          intermediate_model_save = file.path(WorkDir, "Model/models.rds"),
                          k = 5, # number of best models to return
                          offset_effort = "effort_km2",
                          load_saved_models = F
)

### you can now compare the models, run predictions and estimate abundances:
# In the first step
if (!dir.exists(file.path(WorkDir, "predictions"))) {
  dir.create(file.path(WorkDir, "predictions"))
}

model_comparison(run_models = run_models,
                 seg_data = calibration_data, # segments used for run_all_DSM. Should not be an sf object.
                 variable = c("SST.center_1d", "SST.SDspace.mean_1d", "SST.SDtime_8d",
                              "EKE.center_1d", "EKE.SDspace_1d",
                              "NPPV.center_1d", "NPPV.SDtime_8d"),
                 effort_column = "effort_km2",
                 version_preds = "2025-09-01", # version of predictions. You can put a new value if you want to rerun all the function.
                 output_file = paste0("2025-09-01", "_model_comparison"), # without extension, will be the name of html file
                 log1p_trans = NULL,
                 grid_folder = predgrid_directory, # prediction grids
                 prediction_folder = file.path(WorkDir, "predictions"), # folder to store the model predictions
                 correct_bias = T,
                 save_results_bias_corrected = F,
                 use_threshold = T,
                 quantile_mgcv_fixed = "mgcv", # "mgcv", "fixed" or quantile
                 threshold = 1.5, ## either a multiplication factor to use with max(mgcv_density) if quantile_mgcv_fixed = "mgcv",
                 # or a (fixed) maximum overall allowed density if quantile_mgcv_fixed = "fixed",
                 # or a quantile (0 to 1) to use on the density values if quantile_mgcv_fixed = "quantile"
                 breaks_plot = c(-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.25, 1.5, 1000),
                 labels_plot = c("0.00 - 0.10", "0.11 - 0.20", "0.21 - 0.30", "0.31 - 0.40", "0.41 - 0.50",
                                 "0.51 - 0.75", "0.76 - 1.00", "1.01 - 1.25", "1.26 - 1.5", "> 1.50"),
                 corr_groupsize = NULL, # correction factor for groupsize, if groupsize is used, to multiply the predicted densities.
                 response = "n_ppho", # should contain "group" if it is modelling the number of groups rather than of individuals. But please
                 # use n_SpeciesCode for the number of individuals, or n_group_SpeciesCode for the number of groups!
                 subspecies, # in case response contains several species. The groupsize estimate will then account for all the subspecies
                 filter_year_month_not_in = "0000-00", # year and month that should not be used for prediction
                 run_all = F,
                 study_area = NULL, # in case you want to predict only on a part of your prediction grid
                 block_file = NULL, # add here the sf file containing your sub-blocks for your area if you want to use groupsizes as response (containing the column Name for the sub-names), otherwise let NULL
                 data_file = NULL, # add here the intial GPS points to calculate group sizes if needed and print observations. Must contain the column AU for which groupsizes are averaged per AU, and the column Platform (aerial or ship).
                 sub_area_analysis_file = NULL,
                 n_cores = NULL, # path for sf file containing sub-areas if you want to investigate abundance per sub-area in addition to globally. Must contain the column Name for each sub-area.
                 outfile = "log.txt",
                 save_posterior_distribution = F)

### you can also run the gap analysis:
if (!dir.exists(file.path(WorkDir, "gap_analysis"))) {
  dir.create(file.path(WorkDir, "gap_analysis"))
}
gap_analysis(seg_data = calibration_data, # segments used for run_all_DSM. Should not be an sf object.
             variable = c("SST.center_1d", "SST.SDtime_8d",
                          "EKE.center_1d",
                          "NPPV.center_1d"), # let's say you want only with variables from the "best model"
             crs = 3035,
             grid_resolution = 5000,
             version_preds = "2025-09-01",
             output_file = paste0("2025-09-01", "_gap_analysis"), # without extension, will be the name of html file
             grid_folder = predgrid_directory,
             save_results_dsmextra = file.path(WorkDir, "gap_analysis"), # put NULL if you dont want to save results. Path where to save results.
             filter_year_month_not_in = "0000-00", # year and month that should not be used for gap_analysis
             n_cores = NULL,
             outfile = file.path(WorkDir, "gap_analysis", "log.txt"),
             run_all = F,
             study_area = NULL # in case you want to predict only on a part of your prediction grid
             )
