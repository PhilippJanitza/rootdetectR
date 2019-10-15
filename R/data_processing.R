#' @title Test Input data.frame for Rootdetection Output standard
#' @description The function tests wether a given object fulfils all criteria of a standard rootdetection output. These standards are:
#' - The object has to be a data.frame
#' - The columns Label, LengthMM and LengthPx needs to be present
#' - LengthPx needs to store numeric values
#' - A length standard needs to be present (default = 10mm) - can be changed by argument length_standard
#' - hyphens (-) in Labels are not allowed
#' @param root_output data.frame; output *.csv (list) from Rootdetection
#' @return logical; TRUE or FALSE
#' @examples
#' is_root_output(root_output)
#'
#' @export
is_root_output <- function(root_output, length_standard = '10mm') {

  # get name of the input variable
  object_name <- deparse(substitute(root_output))

  # check if input is a data.frame
  dataframe_bool <- is.data.frame(root_output)
  if (!dataframe_bool) {
    stop(paste(object_name, ' is not a data.frame', sep = ''), call. = FALSE)
  }

  # check if column Label is present
  col_Label <- 'Label' %in% colnames(root_output)
  if (!col_Label) {
    stop(paste(object_name, ' is missing the column Label', sep = ''), call. = FALSE)
  }
  # check if column LengthMM is present
  col_LMM <- 'LengthMM' %in% colnames(root_output)
  if (!col_LMM) {
    stop(paste(object_name, ' is missing the column LengthMM', sep = ''), call. = FALSE)
  }
  # check if column LengthPx is present
  col_LPx <- 'LengthPx' %in% colnames(root_output)
  if (!col_LPx) {
    stop(paste(object_name, ' is missing the column LengthPx', sep = ''), call. = FALSE)
  }

  # check if LengthPx is numeric
  pxnum <- is.numeric(root_output$LengthPx)
  if (!pxnum) {
    stop('The column LenghtPx needs to store numeric values!', call. = FALSE)
  }

  # check if length standard is present
  std <- length_standard %in% root_output$Label
  if (!std) {
    stop(paste('Length standard', length_standard, 'is missing! Is the dataset already normalized or do you use a customized lenght standard?'))
  }

  # check if Labels contain hyphens
  hyph <- any(grepl("-", root_output$Label))
  if (hyph == T) {
    stop("Labels containing '-' are not allowed!", call. = TRUE)
  }

  # check if all testings are passed
  if (all(c(dataframe_bool, col_Label, col_LMM, col_LPx, pxnum, std)) && !all(c(hyph))) {
    return(TRUE)
  } else {
    stop('The present input object does not fulfill the
             RootdetectionOutput standard.', call. = FALSE)
  }
}




#' @title Test Input data.frame for normalized Rootdetection standard
#' @description The function tests wether a given object fulfils all criteria of a normalized rootdetection dataset. These standards are:
#' - The object has to be a data.frame
#' - The columns Label, LengthMM and LengthPx, Factor1 and Factor2 needs to be present
#' - LengthMM needs to store numeric values
#' - No length standard must be present
#' - hyphens (-) in Labels are not allowed
#' @param root_output data.frame; output *.csv (list) from Rootdetection
#' @return logical; TRUE or FALSE
#' @examples
#' is_root_output(root_output)
#'
#' @export
is_root_norm <- function(root_norm) {

  # get name of the input variable
  object_name <- deparse(substitute(root_norm))

  # check if input is a data.frame
  dataframe_bool <- is.data.frame(root_norm)
  if (!dataframe_bool) {
    stop(paste(object_name, ' is not a data.frame', sep = ''), call. = FALSE)
  }

  # check if column Label is present
  col_Label <- 'Label' %in% colnames(root_norm)
  if (!col_Label) {
    stop(paste(object_name, ' is missing the column Label', sep = ''), call. = FALSE)
  }
  # check if column LengthMM is present
  col_LMM <- 'LengthMM' %in% colnames(root_norm)
  if (!col_LMM) {
    stop(paste(object_name, ' is missing the column LengthMM', sep = ''), call. = FALSE)
  }
  # check if column LengthPx is present
  col_LPx <- 'LengthPx' %in% colnames(root_norm)
  if (!col_LPx) {
    stop(paste(object_name, ' is missing the column LengthPx', sep = ''), call. = FALSE)
  }
  # check if column Factor1 is present
  col_LPx <- 'Factor1' %in% colnames(root_norm)
  if (!col_LPx) {
    stop(paste(object_name, ' is missing the column Factor1', sep = ''), call. = FALSE)
  }
  # check if column Factor2 is present
  col_LPx <- 'Factor2' %in% colnames(root_norm)
  if (!col_LPx) {
    stop(paste(object_name, ' is missing the column Factor2', sep = ''), call. = FALSE)
  }

  # check if LengthPx is numeric
  mmnum <- is.numeric(root_norm$LengthMM)
  if (!mmnum) {
    stop('The column LenghtPx needs to store numeric values!', call. = FALSE)
  }

  # check if length standard is present
  std <- '10mm' %in% root_norm$Label
  if (std) {
    stop(paste('Length standard 10mm is still present! Use norm_10mm_standard to normalize your dataset?'))
  }

  # check if Labels contain hyphens
  hyph <- any(grepl("-", root_norm$Label))
  if (hyph == T) {
    stop("Labels containing '-' are not allowed!", call. = TRUE)
  }

  # check if all testings are passed
  if (all(c(dataframe_bool, col_Label, col_LMM, col_LPx, mmnum)) && !all(c(hyph, std))) {
    return(TRUE)
  } else {
    stop('The present input object does not fulfill the
             RootdetectionOutput standard.', call. = FALSE)
  }
}



#' @title Calculate LengthMM from LengthPx for Rootdetection standard
#' @description This function takes 10mm standard values to calculate LengthMM from LengthPx
#' @param root_output data.frame; *.csv output from Rootdetection containing 10mm values
#' @param sort logical; if TRUE data.frame is sorted and Label splitted in Factor1 and Factor2
#' @param label_delim character; defines how Factor1 and Factor2 are seperated in Label
#' @return data.frame; containing normalized length values
#' @examples
#' norm_10mm_standard(root_output)
#' norm_10mm_standard(root_output, sort = FALSE)
#'
#' @export
norm_10mm_standard <- function(root_output, sort = TRUE, label_delim = ";") {

    # features --> change name of 10mm to something else choose length that
    # should be used as normalisation (other than 10mm)
    # calc 10mm
    standard_mm <- subset(root_output, Label == "10mm")
    standard_mm_mean <- mean(standard_mm$LengthPx)
    root_output$LengthMM <- root_output$LengthPx / (standard_mm_mean / 10)
    # delete 10mm
    root_output <- subset(root_output, Label != "10mm")
    root_output$Label <- droplevels(root_output$Label)

    # if sort is TRUE order the columns and devide labels according to
    if (sort == TRUE) {
        root_output <- tidyr::separate(data = root_output, col = Label,
                                into = c("Factor1", "Factor2"),
                                sep = label_delim, remove = F)
        root_output <- root_output[, c("Nr", "Filename", "RootNr", "Label",
                                       "Factor1", "Factor2", "LengthPx",
                                       "LengthMM")]

    }

    # return data.frame
    return(root_output)

}



#' @title Calculate LengthMM from LengthPx for Rootdetection standard
#' @description This function takes customized standard values to calculate LengthMM from LengthPx. The label_standard must match the label in the data.frame. In addition the mm value of the standard must be provided.
#' @param root_output data.frame; *.csv output from Rootdetection containing 10mm values
#' @param sort logical; if TRUE data.frame is sorted and Label splitted in Factor1 and Factor2
#' @param label_delim character; defines how Factor1 and Factor2 are seperated in Label
#' @param label_standard string; defines  how the standard is labled in the data.frame
#' @param standard_length_mm numeric; defines the length (in mm) of the standard
#' @return data.frame; containing normalized length values
#' @examples
#' # use data_set with standard 10mm label
#' # add label_standard and standard_length_mm
#'
#' norm_cust_standard(root_output, label_delim = ';', label_standard = '10mm', standard_length_mm = '10')
#'
#' @export
norm_cust_standard <- function(root_output, sort = TRUE, label_delim = ";", label_standard = '10mm', standard_length_mm = '10') {

  # features --> change name of 10mm to something else choose length that
  # should be used as normalisation (other than 10mm)
  # calc 10mm
  standard_mm <- subset(root_output, Label == label_standard)
  standard_mm_mean <- mean(standard_mm$LengthPx)
  root_output$LengthMM <- root_output$LengthPx / (standard_mm_mean / as.numeric(standard_length_mm))
  # delete 10mm
  root_output <- subset(root_output, Label != label_standard)
  root_output$Label <- droplevels(root_output$Label)

  # if sort is TRUE order the columns and devide labels according to
  if (sort == TRUE) {
    root_output <- tidyr::separate(data = root_output, col = Label,
                                   into = c("Factor1", "Factor2"),
                                   sep = label_delim, remove = F)
    root_output <- root_output[, c("Nr", "Filename", "RootNr", "Label",
                                   "Factor1", "Factor2", "LengthPx",
                                   "LengthMM")]

  }

  # return data.frame
  return(root_output)

}



#' @title Calculate Summary Statistics for Rootdetection standard
#' @description This function calculates some summary statistics for a LengthMM normalized Rootdetection output. Calculated values are: sample size (n), median, mean, standard deviation (sd), standard error (se)
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @return data.frame; containing summary statistics for each Factor1 - Factor2 combination
#' @examples
#' root_test_norm <- norm_10mm_standard(root_output)
#' summary_stat(root_test_norm)
#'
#' @export
summary_stat <- function(root_norm) {

    # wählen was man alles haben möchte??
    # check if data.frame is already normalized (containing no 10mm values)
    sum <- plyr::ddply(root_norm, plyr::.(Factor1, Factor2), plyr::summarize,
                 n = length(LengthMM), median = median(LengthMM),
                 mean = mean(LengthMM), sd = sd(LengthMM),
                 se = se(LengthMM))
    return(sum)

}



#' @title Normality Test for Rootdetection standard
#' @description The function performes a normality test for each Factor1 Factor 2 combination in Rootdetecion standard data.frame. Until now only Shapiro-Wilk test is implemented.
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @return data.frame; containing p-values for each Factor1 Facto2 combinations.
#' @examples
#' root_test_norm <- norm_10mm_standard(root_output)
#' normality_test(root_test_norm)
#'
#' @export
normality_test <- function(root_norm) {

  norm_table <- data.frame(Label = character(0), p.value = numeric(0))
  for (lev in 1:length(levels(root_norm$Label))) {

    # create subset contain only the data for one label
    hist_sub <- subset(root_norm, Label == levels(root_norm$Label)[lev])
    temp <-   data.frame(Label = levels(root_norm$Label)[lev], p.value = shapiro.test(hist_sub$LengthMM)$p.value)
    norm_table <- rbind(norm_table, temp)

  }
    return(norm_table)
}



#' @title Create Relative data from Rootdetection standard
#' @description The function takes the control values of Factor2 and calculates the relative data for the other Factor2
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @return data.frame; containing data for each Factor2 relative to Factor2 control
#' @examples
#' root_test_norm <- norm_10mm_standard(root_output)
#' rel_data(root_test_norm)
#'
#' @export
rel_data <- function(root_norm, control = "20") {
    # create subset containing only mock (control) data
    # (name of the control was assigned in the beginning)
    rel_table_mock <- subset(root_norm, Factor2 == control)
    # calc median for all levels of Factor 1 an save to LengthMM_median_control
    rel_table_mock_median <- plyr::ddply(rel_table_mock, plyr::.(Factor1),
                                  plyr::summarize, LengthMM_median_control =
                                  median(LengthMM))
    # merge tables
    rel_table_merge <- merge(root_norm, rel_table_mock_median, by = "Factor1")
    # new subset without the mock data
    rel_table <- subset(rel_table_merge, Factor2 != control)
    # calculate relative values for the new 'treatment-only' table
    rel_table$relative_value <-
      (100 * rel_table$LengthMM) / rel_table$LengthMM_median_control
    rel_root_table <- rel_table[, c("Label", "Factor1", "Factor2", "RootNr",
                                    "relative_value")]

    return(rel_root_table)
}


#' @title Remove Outlier from Rootdetection standard
#' @description The function removes or replaces outliers for every Factor1 Factor2 (Label) combination. Outliers are defined as 1.5 x IQR.
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @param fill_na logical; If TRUE all outliers present in x will be replaced by NA. If False all outliers will be deleted.
#' @return data.frame; containing data without outliers or outliers replaced by NA
#' @examples
#' root_test_norm <- norm_10mm_standard(root_output)
#'
#' # transform outliers to NA
#' rm_outlier_df(root_test_norm, fill_na = T)
#' # remove outliers
#' rm_outlier_df(root_test_norm, fill_na = F)
#'
#' @export
rm_outlier_df <- function(root_norm, fill_na = F) {
  # remove all NAs before from LengthMM
  root_norm <- root_norm[!is.na(root_norm$LengthMM),]

  # create empty data-frame
  final <- root_norm[0,]

  # loop over labels in data.frame
  for(name in unique(root_norm$Label)){
    #extract subsets and replace outliers with NA
    temp <- root_norm[root_norm$Label == name,]
    temp$LengthMM <- rm_outlier(temp$LengthMM, fill_na = T)

    final <- rbind(final, temp)
  }

  if(fill_na == F){
    final <- final[!is.na(final$LengthMM),]
  }
  return(final)
}
