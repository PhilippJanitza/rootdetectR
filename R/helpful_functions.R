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
#' @return numeric vector; without outliers or outliers replaced by NA
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



#' @title Create subsets by Factor2
#' @description This function splits Normalized Rootdetection Output in subsets for each Factor2 control and treatment pair.
#' To be functional Factor2 must contain a control and at least two other conditions (treatments)
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @param control string; name of the Factor2 control condition
#' @return list of data.frames; each subset will be in an seperated data.frame stored in a list of data.frames
#' @example
#' # will come in next version with an data.frame containing multiple Factor2
#'
#' @export
subset_fac2 <- function(root_norm, control = '20'){

  root_norm$Factor2 <- as.factor(root_norm$Factor2)

  if(length(levels(root_norm$Factor2)) <= 2){

    warning('There are only 2 levels of Factor2 - no subsetting by Factor2 possible!')

    #else create subtables
  } else {

    # variable containing only Factor2 control
    root_norm_control <- root_norm[root_norm$Factor2 == control,]

    # create empty list to put subsets in
    dfl <- list()


    # loop over all levels in Factor2 that are not control
    for (i in levels(root_norm$Factor2)[levels(root_norm$Factor2) != control]) {

      # get subset for each i
      sub <- root_norm[root_norm$Factor2 == i,]
      #bind root_norm_control with sub_temp
      sub_control <- rbind(sub, root_norm_control)

      # paste in list and name
      name <- paste(unique(root_norm_control$Factor2), '_and_', i, sep = '')
      dfl[[name]] <- sub_control

    }
  }
  return(dfl)
}

