#' @title Test Input data.frame for Rootdetection standard
#' @description The function tests wether a given data.frame has all important columns and the corresponding data types.
#' @param root_output data.frame; output *.csv from Rootdetection
#' @return logical; TRUE or FALSE
#' @examples
#' is.RootdetectionOutput(root_output)
#'
#' @export
is.root_detection_output <- function(root_output) {

    # sammelkommentar um weitere Abfragen zu sammeln sind 0 Werte vorhanden?

    d.f_bool <- is.data.frame(root_output)
    if (!d.f_bool) {
        stop("RootdetectionOutput is not a data.frame", call. = FALSE)
    }
    cols <- all(c("Label", "LengthMM", "LengthPx") %in% colnames(root_output))
    if (!cols) {
        stop("data.frame is missing one or more of the important columns
             Label, LenghtPx and LengthMM - check colnames!", call. = FALSE)
    }
    pxnum <- is.numeric(root_output$LengthPx)
    if (!pxnum) {
        stop("The column LenghtPx needs to store numeric values", call. = FALSE)
    }
    # print messages als warnings??
    std <- "10mm" %in% root_output$Label
    if (!std) {
        warning("10mm standard is missing")
        mmnum <- is.numeric(root_output$LengthMM)
        pxmm <- all(root_output$LengthMM == root_output$LengthPx)
        if (mmnum == T && pxmm == F) {
            warning("Are the LengthMM values already normalized
                    with a standard?")
        } else {
            stop("LengthMM seems to be not normalized with a standard!")
        }
    }

    hyph <- any(grepl("-", root_output$Label))
    if (hyph == T) {
        stop("Labels containing '-' are not allowed!", call. = TRUE)
    }

    labna <- any(is.na(root_output$Label))
    if (labna == T) {
        stop("Labels containing NA!", call. = TRUE)
    }

    pxna <- any(is.na(root_output$LengthPx))
    if (pxna == T) {
        stop("LengthPx containing NA!", call. = TRUE)
    }

    pxzero <- any(root_output$LengthPx == 0)
    if (pxzero == T) {
        warning("LengthPx containg one or more 0")
    }


    if (all(c(d.f_bool, cols, pxnum)) && !all(hyph, labna, pxna)) {
        return(TRUE)
    } else {
        stop("The present input object does not fulfill the
             RootdetectionOutput standard.", call. = FALSE)
    }
}



#' @title Calculate LengthMM from LengthPx for Rootdetection standard
#' @description This function takes 10mm standard values to calculate LengthMM from LengthPx
#' @param root_output data.frame; *.csv output from Rootdetection containing 10mm values
#' @param sort logical; if TRUE data.frame is sorted and Label splitted in Factor1 and Factor2
#' @param label_delim character; defin how Factor1 and Factor2 are seperated in Label
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
