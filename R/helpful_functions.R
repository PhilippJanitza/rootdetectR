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



#' @title Remove Outlier
#' @description This function removes outlier from a vector and returns either a vector without the outliers or replaces the outliers with NA.
#' Outliers are defined as 1.5 x IQR. Caution: NAs already present in the input vector will be removed.
#' @param x numeric vector
#' @param fill.na logical; If TRUE all outliers present in x will be replaced by NA. If False all outliers will be deleted.
#' @returns numeric vector; without outliers or outliers replaced by NA
#' @examples
#' # get some example vector
#' root_test_norm <- norm_10mm_standard(root_output)
#' test_vector <- root_test_norm[root_test_norm$Label == 'weitar1_1;28',]
#'
#' # just delete outliers
#' rm_outlier(test_vector, fill_na = F)
#'
#' # replace outliers with NA
#' rm_outlier(test_vector, fill_na = T)
#'
#' @export
rm_outlier <- function(x, fill_na = F) {

  low_border <- quantile(x, 0.25, na.rm = T) - (IQR(x, na.rm = T) * 1.5)
  high_border <- quantile(x, 0.75, na.rm = T) + (IQR(x, na.rm = T) * 1.5)

  if(fill_na){

    # outliers will be changed to NA
    x[which(x < low_border | x > high_border)] <- NA

  } else {
    # outliers will be just deleted
    x <- x[which(x > low_border & x < high_border)]

  }
  return(x)
}


