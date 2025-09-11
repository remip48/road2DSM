
# LAST UPDATE:

## **11 September 2025** <br>

**Please update if your installation of road2DSM is prior to this
revision.**

The full pipeline can be followed in the script **run_road2DSM.R** from
the main folder, where the package is explained step by step. Below is
the more detailed explanation of the functions.

### Important

The functions are adapted to the type of data used and created by the
R-package. If you want to use external grids or effort, please adapt the
format to the one created by the functions below.

### Functions

-   **segmentate**: use GPS effort (**PorpAbundance** format) to cut
    segments. **label** is the unique identifier of the segment in the
    output. **effort_km2** and **distance** are the effectively covered
    area and length of the segment based on the input data. For each
    unique code in species, 3 columns are created in the ouput:
    **n_code** (number of individuals based on SIG_number column),
    **n_group_code** (number of sightings) and **n_calves_code** (based
    on the calves column).

-   **extract_nc**: will extract automatically the center, SDspace,
    SDtime and mean of your variables in your nc files. You can define
    for each file if you want all of the 4 variable types or, for
    example, dont want to extract the mean and SDtime. You can extract
    variables over different periods and different pixel radius. Results
    are saved in a folder with the nc grid saved as sf.

-   **upscale_grid** / **create_grid**: upscale an input grid to the
    desired resolution, or create a new grid based on the extent you
    provide. The unique identifier for each cell is **id.** In
    **create_grid**, you can provide the country or sea (i.e. study
    area) shapefile to automatically intersect your grid with these
    ones.

-   **upscale_SDspace**: is used to calculate SDspace over larger
    spatial grid resolution than the nc file (i.e. study mesoscale
    events rather than small scale events). The input data are the grids
    created by **extract_grid** (see below) and is based on
    **upscale_grid** / **create_grid**.

-   **crop_effort**: will load all sf grids starting with **file_set\_**
    in your directory to combine them in a single sf POLYGON file. Your
    effort is then cropped to this joined POLYGON file.

-   **extract_effort**: will extract all covariates extracted from
    **extract_nc** in your effort (must be an sf POLYGON object). You
    can also extract the covariates from **upscale_SDspace** by
    reproducing the folder structure created by **extract_nc**: this
    means store the grids from **upscale_SDspace** in a folder
    **file_set_X**, where X is a random and unique number, and add in
    the parent directory the corresponding sf grid (that was used to run
    **upscale_SDspace**) under the name **file_set_X.shp**.

-   **extract_grid**: extract all covariates from **extract_nc** on the
    desired grid (e.g. created by **upscale_grid** / **create_grid**).
    Based on this output, **upscale_SDspace** can be run. Then,
    **extract_grid** can be used again to extract covariates from
    **upscale_SDspace** on the prediction grid. In this case, the folder
    structure of **extract_nc** must be reproduce to extract results of
    **upscale_SDspace** as explained in **extract_effort**.

-   **gap_analysis**: run **dsmextra** R-package to evaluate
    extrapolations, nearby data and most-influencal covariates on
    extrapolations in a all-in-one function. Will write a markdown with
    results.

-   **run_all_DSM** and **backward_selection**: model fitting. The first
    fit all combinations of possible covariates (flexible function), the
    second is running the backward selection (not running yet!).

-   **model_comparison**: compare the models fitted above with different
    statistics. Plot splines, predictions, and estimate abundances. You
    can define sub-areas to calculate abundances as well. Correct bias
    from predictions, and write final results. Will write a markdown
    with results.
