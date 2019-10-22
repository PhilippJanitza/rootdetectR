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
#' @examples
#' # create normalized data.set with multiple Factor2
#' root_norm_multfac2 <- norm_cust_standard(root_output_multfac2, label_standard = '20mm', standard_length_mm = '20')
#'
#' # create subsets for each control_treatment pair
#' subset_fac2(root_norm_multfac2, control = '20')
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



#' @title Create rootdetectR Compatible Matrix from TukeyHSD output
#' @description The function will create a rootdetectR compatible Matrix from TukeyHSD output. The returened matrix can be used to retrieve plotting letters or can be used directly inside the plotting functions (plot_rel, plot_abs).
#' @param tukeyHSD_output matrix or data.frame; Output from the TukeyHSD() function.
#' @return matrix; containing the p-values from TukeyHSD() in an appropriate format for other rootdetectR functions.
#' @examples
#' # conduct ANOVA and TukeyHSD
#' root_norm <- norm_10mm_standard(root_output)
#' aov_all_vs_all <- aov(LengthMM ~ Factor1 * Factor2, data = root_norm)
#' tuk <- TukeyHSD(aov_all_vs_all, ordered = FALSE)$`Factor1:Factor2`
#'
#' # use Tukey Output to retrieve matrix
#' tukey_to_matrix(tuk)
#'
#' @export
tukey_to_matrix <- function(tukeyHSD_output) {

  if(!is.data.frame(tukeyHSD_output)){

    tukeyHSD_output <- as.data.frame(tukeyHSD_output)

  }

  temp <- data.frame(name = rownames(tukeyHSD_output), p.val = tukeyHSD_output$`p adj`)

  temp_new <- tidyr::separate(temp, 'name', into = c('V1', 'V2'), sep = '-')
  labs <- sort(unique(c(temp_new$V1, temp_new$V2)))
  nr_labs <- length(labs)
  # create empty matrix
  mat <- matrix(NA, nrow = nr_labs, ncol = nr_labs)
  colnames(mat) <- labs
  rownames(mat) <- labs


  for (j in 1:(nr_labs - 1)) {
    for (k in (j + 1):nr_labs) {

      # get p-values and put them into the matrix
      idx <- which(paste(labs[j], '-', labs[k], sep = '') == temp$name)
      if (length(idx) == 0) {
        idx <- which(paste(labs[k], '-', labs[j], sep = '') == temp$name)
      }
      if (length(idx) != 0) {
        mat[j, k] <- temp[idx, 2]
      }
    }
  }
  return(mat)
}



#' @title Create Significance Letter encoding
#' @description The function takes a matrix as input (from rootdetectR ANOVA analysis or from tukey_to_matrix function) and returns a table with significance letters.
#' @param tukmatrix matrix; Output from ANOVA functions (onefacaof_fac1, onefacaov_fac2, twofacaov, interaction_twofacaov)
#' @return data.frame; containing Label and the corresponding significance letter
#' @examples
#' # conduct ANOVA and TukeyHSD
#' root_norm <- norm_10mm_standard(root_output)
#' aov_all_vs_all <- aov(LengthMM ~ Factor1 * Factor2, data = root_norm)
#' tuk <- TukeyHSD(aov_all_vs_all, ordered = FALSE)$`Factor1:Factor2`
#'
#' # use Tukey Output to retrieve matrix
#' tukey_to_matrix(tuk)
#'
#' # get significance_letters
#' get_sig_letters(tuk)
#'
#'
#' @export
get_sig_letters <- function(tukmatrix){

  # get Letters for twofacaov output
  mat_names <- character()
  mat_values <- numeric()
  # loop over matrix and get names + values
  for (j in 1:(length(row.names(tukmatrix)) - 1)) {
    for (k in (j + 1):length(colnames(tukmatrix))) {
      v <- tukmatrix[j, k]
      t <- paste(row.names(tukmatrix)[j],
                 colnames(tukmatrix)[k], sep = "-")
      mat_names <- c(mat_names, t)
      mat_values <- c(mat_values, v)
    }
  }

  # combine names + values
  names(mat_values) <- mat_names
  # get df with letters and replace : with label delim!!
  letters <-
    data.frame(multcompView::multcompLetters(mat_values)["Letters"])

  return(letters)

}
