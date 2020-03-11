#' @title Calculate Standard Error
#' @description This function calculates the standard error for a numeric-vector x. NA values will be removed.
#' @param x numeric vector
#' @return numeric value; standard error of vector x
#' @examples
#' x <- 1:20
#' se(x)
#' @export
se <- function(x) sd(x, na.rm = T) / sqrt(length(na.omit(x)))



#' @title Remove Outlier
#' @description This function removes outlier from a vector or replaced outliers by NA.
#' According to box plots the function defines outliers as 1.5 x IQR. Caution: NAs already present in the input vector will be removed first.
#' @param x numeric vector
#' @param fill.na logical; If TRUE all outliers present in x will be replaced by NA. If FALSE all outliers will be deleted.
#' @return numeric vector; without outliers or outliers replaced by NA
#' @examples
#' # get some example vector
#' root_norm <- norm_10mm_standard(root_output)
#' x <- root_norm[root_norm$Label == 'weitar1_1;28',]
#'
#' # delete outliers
#' rm_outlier(x, fill_na = F)
#'
#' # replace outliers with NA
#' rm_outlier(x, fill_na = T)
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



#' @title Create Subsets Of Normalized Rootdetection Standard Data Set By Factor2
#' @description This function creates subsets for each grouping variable 2 (Factor2) control and treatment pair of a normalized Rootdetection data set.
#' (At least two treatment conditions are required)
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param control string; name of the grouping variable 2 (Factor2) control condition
#' @return list of data.frames; each subset will be in an separated data.frame stored in a list
#' @examples
#' # create normalized data.set with multiple Factor2
#' root_norm_multfac2 <- norm_cust_standard(root_output_multfac2, label_standard = '20mm', standard_length_mm = '20')
#'
#' # create subsets for each control_treatment pair
#' subset_fac2(root_norm_multfac2, control = '20')
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



#' @title Create Matrix From TukeyHSD Output
#' @description The function will create a rootdetectR compatible Matrix from TukeyHSD output.
#' The returened matrix can be used to retrieve plotting letters with the get_sig_letters function from rootdetectR.
#' @param tukeyHSD_output matrix or data.frame; Output from the TukeyHSD() function.
#' @return matrix; containing p-values from TukeyHSD() in an appropriate format for other rootdetectR functions.
#' @examples
#' # conduct ANOVA and TukeyHSD
#' root_norm <- norm_10mm_standard(root_output)
#' aov_all_vs_all <- aov(LengthMM ~ Factor1 * Factor2, data = root_norm)
#' tuk <- TukeyHSD(aov_all_vs_all, ordered = FALSE)$`Factor1:Factor2`
#'
#' # use Tukey Output to retrieve matrix
#' tukey_to_matrix(tuk)
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



#' @title Create Significance Letter Encoding
#' @description The function takes a matrix as input (from rootdetectR ANOVA analysis or from tukey_to_matrix function) and returns a table with significance letters.
#' @param tukmatrix matrix; Output from ANOVA functions (onefacaof_fac1, onefacaov_fac2, twofacaov, interaction_twofacaov) or tukey_to_matrix function
#' @return data.frame; containing Label and the corresponding significance letter
#' @examples
#' # conduct ANOVA and TukeyHSD
#' root_norm <- norm_10mm_standard(root_output)
#' aov_all_vs_all <- aov(LengthMM ~ Factor1 * Factor2, data = root_norm)
#' tuk <- TukeyHSD(aov_all_vs_all, ordered = FALSE)$`Factor1:Factor2`
#'
#' # use Tukey Output to retrieve an appropriate matrix
#' tuk_mat <- tukey_to_matrix(tuk)
#'
#' # obtain significance_letters
#' get_sig_letters(tuk_mat)
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



#' @title Sort Grouping Variable Of Normalized Rootdetection Data Set
#' @description The function will sort the grouping variable (Label) according to Controls for plotting.
#' By defining control_fac1 and control_fac2 you can set the corresponding grouping variable (label) to the first place in order.
#' _Caution: The column that should be ordered must be organized as a factor._
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param label_delim character; defines how Factor1 and Factor2 are separated in Label
#' @param col_label string; name of the column carring the grouping variable (Label)
#' @param control_fac1 string; name of the control condition of grouping variable1 (Factor1)
#' @param control_fac2 string; name of the control condition of grouping variable2 (Factor2)
#' @return data.frame; with sorted col_label
#' @examples
#' # obtain normalized Rootdetection data.frame
#' root_norm <- norm_10mm_standard(root_output)
#'
#' # check the default sorting of the grouping variable (Labels)
#' # sorted alphabetically by default
#' levels(root_norm$Label)
#'
#' # change order of levels by setting yucOx and 28 as controls
#' sort_label(root_norm, control_fac1 = 'yucOx', control_fac2 = '28')
#' @export
sort_label <- function(root_norm, label_delim = ';', col_label = 'Label', control_fac1, control_fac2){


  colnames(root_norm)[colnames(root_norm) == col_label] <- 'Label'
  root_norm$Label <- as.factor(root_norm$Label)

  # get control positions
  con_fac1_position <- grep(paste(control_fac1, label_delim, sep = ''), levels(root_norm$Label))
  # get all other positions and put controls in first place
  sort_fac1 <- c(levels(root_norm$Label)[con_fac1_position], levels(root_norm$Label)[-con_fac1_position])

  fac1 <- unique(sub(paste(label_delim, ".*", sep = ''), '', sort_fac1))
  fac2 <- unique(sub(paste(".*", label_delim, sep = ''), '', sort_fac1))

  sorted_fac1_fac2 <- character(0)

  for (i in fac1) {

    temp <- sort_fac1[grep(paste(i, label_delim, sep = ''), sort_fac1)]
    con_fac2_position <- grep(control_fac2, unique(sub(paste(".*", label_delim, sep = ''), '', temp)))
    sort_fac2 <- c(temp[con_fac2_position], temp[-con_fac2_position])
    sorted_fac1_fac2 <- c(sorted_fac1_fac2, sort_fac2)
  }


  root_norm$Label <- factor(root_norm$Label,levels = sorted_fac1_fac2)

  root_norm$Label <- as.factor(root_norm$Label)
  colnames(root_norm)[colnames(root_norm) == 'Label'] <- col_label


  return(root_norm)


}



#' @title Detach All Elements
#' @description This function will detach all elements except th base R packages:
#' .GlobalEnv', 'package:stats', 'package:graphics', 'package:grDevices', 'package:utils',
#' 'package:datasets', 'package:methods', 'Autoloads', 'package:base', 'tools:rstudio'
#' _Caution: This function will also detach rootdetectR_
#' @param except vector; elements that should be excluded from detaching
#' @examples
#' # detach all elements in search()
#' detach_all()
#'
#' # detach all element except rootdetectR
#' detach_all(except = 'rootdetectR')
#'
#' @export
detach_all <- function(except){


  # create vector with base packages
  base_pkg <- c('.GlobalEnv', 'package:stats', 'package:graphics', 'package:grDevices', 'package:utils', 'package:datasets', 'package:methods', 'Autoloads', 'package:base', 'tools:rstudio')

  if(!missing(except)){
    base_pkg <- c(base_pkg, except)
  }

  # which are not in base packages
  att_pkg <- search()[which(!search() %in% base_pkg)]

  # remove all that are not in base list
  for(e in att_pkg){

    detach(name = e, character.only = T)

  }

}



#' @title Change Names of Elements In Specific Column
#' @description The function will rename elements in a specific column by vectors "old_name" and "new_name".
#' Simply give the function a vector with old names that should be renamed by a vector with new names.
#' _Be carefull that elements in new_name have the rigth position_ (e.g. the second element in "old_name" will be replaced by the second element in "new_name").
#' @param df data.frame; input data.frame
#' @param colname character; name of the column storing the elements to be renamed
#' @param old_name vector; elements that should be renamed in column specified by 'colname'
#' @param new_name vector; new names of the elements specified by "old_name"
#' @return data.frame; with renamed elements in "colname"
#' @examples
#' rename_element(root_output, colname = 'Label', oldname = c('Col_0;28', 'Col_0;20'), new_name = c('Col0_28', 'Col0_20'))
#' @export
rename_element <- function(df, colname, old_name, new_name){

  # check if colname is factor -> transform to character
  if(is.factor(df[,colname])){
    df[,colname] <- as.character(df[,colname])
  }

  # check if all old_names matches
  if(all(old_name %in% df[,colname])){
    # check if old_name and new_name have same length
    if(length(old_name) == length(new_name)){

      for(i in 1:length(old_name)){

        df[,colname][df[,colname] == old_name[i]] <- new_name[i]

      }

      df[,colname] <- as.factor(df[,colname])
      return(df)

    } else {
      stop('"old_name" must have the same length as "new_name"')
    }
  } else {
    stop(c('element/s ', paste(noquote(old_name[which(!old_name %in% df[,colname])]),collapse=','), ' do not exist in ', colname))
  }
}
