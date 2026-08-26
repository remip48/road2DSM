#' Title
#'
#' @param ynew
#' @param y
#'
#' @return
#' @export
#'
#' @examples
rescale2 <- function (ynew, y = NULL)
{
  # assert_that(is.numeric(ynew) && is.vector(ynew) && length(ynew) >=
  #               3)
  # assert_that((is.numeric(y) && is.vector(y) && length(y) >=
  #                3) || is.null(y))
  if (is.null(y)) {
    out <- (ynew - mean(ynew, na.rm = TRUE))/(sd(ynew, na.rm = TRUE))
  }
  else {
    out <- (ynew - mean(y, na.rm = TRUE))/(sd(y, na.rm = TRUE))
  }
  return(out)
}
