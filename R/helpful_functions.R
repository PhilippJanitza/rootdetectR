#' @title Calculate Standard Error
#' @description This function calculates the standard error for a numeric-vector x. NA values will be removed.
#' @param x a numeric vector
#' @return numeric value; standard error of vector x
#' @examples
#' x <- 1:20
#' se(x)
#'
#' @export
se <- function(x) sd(x, na.rm = T) / sqrt(length(na.omit(x)))


