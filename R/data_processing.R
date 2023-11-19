#' @title Test Input Data.frame For Rootdetection Output
#' @description The function tests whether a given object fulfills all criteria of a standard Rootdetection output. These standards are:
#'
#' - The object has to be a data.frame
#'
#' - The columns Label and LengthPx needs to be present
#'
#' - LengthPx needs to store numeric values
#'
#' - A length standard needs to be present (default = 10mm) - can be changed by argument length_standard
#'
#' - hyphens (-) in Labels are not allowed
#' @param root_output data.frame; output *.csv (list) from Rootdetection
#' @param length_standard string; set name of the standard
#' @return logical; TRUE or FALSE
#' @examples
#' is_root_output(root_output)
#' @export
is_root_output <- function(root_output, length_standard = "10mm") {
  # get name of the input variable
  object_name <- deparse(substitute(root_output))

  # check if input is a data.frame
  dataframe_bool <- is.data.frame(root_output)
  if (!dataframe_bool) {
    stop(paste(object_name, " is not a data.frame", sep = ""), call. = FALSE)
  }

  # check if column Label is present
  col_Label <- "Label" %in% colnames(root_output)
  if (!col_Label) {
    stop(paste(object_name, " is missing the column Label", sep = ""), call. = FALSE)
  }
  # check if column LengthPx is present
  col_LPx <- "LengthPx" %in% colnames(root_output)
  if (!col_LPx) {
    stop(paste(object_name, " is missing the column LengthPx", sep = ""), call. = FALSE)
  }

  # check if LengthPx is numeric
  pxnum <- is.numeric(root_output$LengthPx)
  if (!pxnum) {
    stop("The column LenghtPx needs to store numeric values!", call. = FALSE)
  }

  # check if length standard is present
  std <- length_standard %in% root_output$Label
  if (!std) {
    stop(paste("Length standard", length_standard, "is missing! Is the dataset already normalized or do you use a customized lenght standard?"))
  }

  # check if Labels contain hyphens
  hyph <- any(grepl("-", root_output$Label))
  if (hyph == T) {
    stop("Labels containing '-' are not allowed!", call. = TRUE)
  }

  # check if all testings are passed
  if (all(c(dataframe_bool, col_Label, col_LPx, pxnum, std)) && !all(c(hyph))) {
    return(TRUE)
  } else {
    stop("The present input object does not fulfill the
             RootdetectionOutput standard.", call. = FALSE)
  }
}



#' @title Test Input Data.frame For Normalized Rootdetection Data Set
#' @description The function tests whether a given object fulfills all criteria of a normalized Rootdetection data set. These standards are:
#'
#' - The object has to be a data.frame
#'
#' - The columns Label, LengthMM and LengthPx, Factor1 and Factor2 needs to be present
#'
#' - LengthMM needs to store numeric values
#'
#' - No length standard must be present
#'
#' - hyphens (-) in Labels are not allowed
#' @param root_norm data.frame; normalized Rootdetection data set
#' @return logical; TRUE or FALSE
#' @examples
#' # normalize root_output
#' root_norm <- norm_10mm_standard(root_output)
#'
#' # check data set for normalization
#' is_root_norm(root_norm)
#' @export
is_root_norm <- function(root_norm) {
  # get name of the input variable
  object_name <- deparse(substitute(root_norm))

  # check if input is a data.frame
  dataframe_bool <- is.data.frame(root_norm)
  if (!dataframe_bool) {
    stop(paste(object_name, " is not a data.frame", sep = ""), call. = FALSE)
  }

  # check if column Label is present
  col_Label <- "Label" %in% colnames(root_norm)
  if (!col_Label) {
    stop(paste(object_name, " is missing the column Label", sep = ""), call. = FALSE)
  }
  # check if column LengthMM is present
  col_LMM <- "LengthMM" %in% colnames(root_norm)
  if (!col_LMM) {
    stop(paste(object_name, " is missing the column LengthMM", sep = ""), call. = FALSE)
  }
  # check if column LengthPx is present
  col_LPx <- "LengthPx" %in% colnames(root_norm)
  if (!col_LPx) {
    stop(paste(object_name, " is missing the column LengthPx", sep = ""), call. = FALSE)
  }
  # check if column Factor1 is present
  col_LPx <- "Factor1" %in% colnames(root_norm)
  if (!col_LPx) {
    stop(paste(object_name, " is missing the column Factor1", sep = ""), call. = FALSE)
  }
  # check if column Factor2 is present
  col_LPx <- "Factor2" %in% colnames(root_norm)
  if (!col_LPx) {
    stop(paste(object_name, " is missing the column Factor2", sep = ""), call. = FALSE)
  }

  # check if LengthPx is numeric
  mmnum <- is.numeric(root_norm$LengthMM)
  if (!mmnum) {
    stop("The column LenghtPx needs to store numeric values!", call. = FALSE)
  }

  # check if length standard is present
  std <- "10mm" %in% root_norm$Label
  if (std) {
    stop(paste("Length standard 10mm is still present! Use norm_10mm_standard to normalize your dataset?"))
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
    stop("The present input object does not fulfill the
             RootdetectionOutput standard.", call. = FALSE)
  }
}



#' @title Normalize Rootdetection Output According To Standard
#' @description This function transforms Rootdetection output (*.csv as list) to a normalized Rootdetection data set.
#' Therefore the function calculates the length in mm (LengthMM) from length in pixel (LengthPx).
#' The standard measurements has to be labeled as '10mm' and has to store pixel values representing 10 mm.
#' @param root_output data.frame; *.csv output from Rootdetection containing standard values ('10mm')
#' @param split logical; if TRUE data.frame is sorted and Label is split in Factor1 and Factor2 (necessary for standard procedure)
#' @param label_delim character; defines how Factor1 and Factor2 are separated in Label
#' @return data.frame; Rootdetection data set,containing normalized length values
#' @examples
#' norm_10mm_standard(root_output)
#' norm_10mm_standard(root_output, split = FALSE)
#' @export
norm_10mm_standard <- function(root_output, split = TRUE, label_delim = ";") {
  # transform Label to factor
  root_output$Label <- as.factor(root_output$Label)

  # features --> change name of 10mm to something else choose length that
  # should be used as normalisation (other than 10mm)
  # calc 10mm
  standard_mm <- root_output[root_output$Label == "10mm", ]
  standard_mm_mean <- mean(standard_mm$LengthPx)
  root_output$LengthMM <- root_output$LengthPx / (standard_mm_mean / 10)
  # delete 10mm
  root_output <- root_output[root_output$Label != "10mm", ]
  root_output$Label <- droplevels(root_output$Label)

  # if split is TRUE order the columns and devide labels according to
  if (split == TRUE) {
    root_output <- tidyr::separate(
      data = root_output, col = "Label",
      into = c("Factor1", "Factor2"),
      sep = label_delim, remove = F
    )
    root_output <- root_output[, c(
      "Nr", "Filename", "RootNr", "Label",
      "Factor1", "Factor2", "LengthPx",
      "LengthMM"
    )]
  }

  # return data.frame
  return(root_output)
}



#' @title Normalize Rootdetection Output According To Customized Standard
#' @description This function transforms Rootdetection output (*.csv as list) to a normalized Rootdetection data set.
#' Compared to norm_10mm_standard function you can adjust several parameters defining the standard to obtain a normalized Rootdection data set.
#' Therefore label_standard must match the label of your length standard measurements in the data.frame. In addition the exact mm value of the length standard must be provided.
#' @param root_output data.frame; *.csv output from Rootdetection containing length standard values
#' @param split logical; if TRUE data.frame is sorted and label is split in Factor1 and Factor2 (necessary for standard procedure)
#' @param label_delim character; defines how Factor1 and Factor2 are separated in Label
#' @param col_label string; name of the column that should be used as grouping variable (Label)
#' @param col_value string; name of the column containing the measured values in pixel (depending variable) (LengthMM)
#' @param label_standard string; defines  how the length standard is labeled in col_label
#' @param standard_length_mm numeric; defines the length (in mm) of the standard
#' @return data.frame; Rootdetection data set, containing normalized length values
#' @examples
#' # use data_set with standard called '20mm' and 20 mm measured
#' norm_cust_standard(root_output_multfac2,
#'   label_delim = ";",
#'   col_label = "Label", col_value = "LengthPx",
#'   label_standard = "20mm", standard_length_mm = "20"
#' )
#' @export
norm_cust_standard <- function(root_output,
                               split = TRUE,
                               label_delim = ";",
                               col_label = "Label",
                               col_value = "LengthPx",
                               label_standard = "10mm",
                               standard_length_mm = "10") {
  colnames(root_output)[colnames(root_output) == col_label] <- "Label"
  colnames(root_output)[colnames(root_output) == col_value] <- "LengthPx"
  root_output$Label <- as.factor(root_output$Label)

  # calc 10mm
  standard_mm <- root_output[root_output$Label == label_standard, ]
  standard_mm_mean <- mean(standard_mm$LengthPx)
  root_output$LengthMM <- root_output$LengthPx / (standard_mm_mean / as.numeric(standard_length_mm))
  # delete 10mm
  root_output <- root_output[root_output$Label != label_standard, ]
  root_output$Label <- droplevels(root_output$Label)

  # if split is TRUE order the columns and devide labels according to
  if (split == TRUE) {
    root_output <- tidyr::separate(
      data = root_output, col = "Label",
      into = c("Factor1", "Factor2"),
      sep = label_delim, remove = F
    )
    root_output <- root_output[, c(
      "Nr", "Filename", "RootNr", "Label",
      "Factor1", "Factor2", "LengthPx",
      "LengthMM"
    )]
  }

  # return data.frame
  return(root_output)
}


#' @title Visualize and Inspect Rootdetection Data
#' @importFrom magrittr "%>%"
#' @description The function produces a summary of all lines and conditions of a normalized rootdetection dataset.
#' It visualizes a summary in a matrix styled plot and gives a description of data in the R console.
#' @param root_norm data.frame;  normalized Rootdetection data set
#' @param plot logical; if TRUE summary plot is produced
#' @param output logical; if TRUE summary is printed in R console
#' @return plot and/or console output
#' @examples
#' # normalize root_output
#' root_norm <- norm_10mm_standard(root_output)
#'
#' # data summary
#' inspect_root_norm(root_norm)
#' @export
inspect_root_norm <- function(root_norm, plot = TRUE, output = TRUE) {
  if (is_root_norm(root_norm) == F) {
    stop("Input does not fit the criteria - check your input with is_root_norm() function!")
  }

  rs_sum <- rootdetectR::summary_stat(root_norm)

  rs_sum$quality[rs_sum$n < 5] <- "bad"
  rs_sum$quality[rs_sum$n >= 5 & rs_sum$n < 10] <- "fair"
  rs_sum$quality[rs_sum$n >= 10] <- "good"

  if (output == T) {
    cat(paste("Factor1 levels present:", paste(as.character(unique(rs_sum$Factor1)), collapse = ", "), "\n"))
    cat(paste("Factor2 levels present:", paste(as.character(unique(rs_sum$Factor2)), collapse = ", "), "\n\n"))

    fair <- rs_sum[which(rs_sum$n < 10 & rs_sum$n >= 5), ]
    if (nrow(fair) != 0) {
      cat("Lines with n 5-10: ")
      cat(paste(as.character(fair$Factor1), ";", as.character(fair$Factor2), "\n", sep = ""))
    }

    bad <- rs_sum[which(rs_sum$n < 5 & rs_sum$n >= 1), ]
    if (nrow(bad) != 0) {
      cat("Lines with n 1-5: ")
      cat(paste(as.character(bad$Factor1), ";", as.character(bad$Factor2), "\n", sep = ""))
    }


    missing <- rs_sum[which(rs_sum$n == 0), ]
    if (nrow(missing) != 0) {
      cat("Missing Lines (n =0): ")
      cat(paste(as.character(missing$Factor1), ";", as.character(missing$Factor2), sep = ""))
    }
  }

  if (plot == T) {
    p <- ggplot2::ggplot(rs_sum, ggplot2::aes_string(y = "Factor1", x = "Factor2")) +
      ggplot2::geom_tile(ggplot2::aes_string(fill = "quality"), colour = "black") +
      ggplot2::geom_text(ggplot2::aes_string(label = "n")) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::scale_fill_manual(values = c("bad" = "red", "fair" = "yellow", "good" = "green"))

    return(p)
  }
}



#' @title Calculate Summary Statistics for Normalized Rootdetection Data Set
#' @description This function calculates sample size (n), median, mean, standard deviation (sd) and standard error (se) according to grouping variable/s
#' (Factor1 and Factor2) and a dependent variable (LengthMM).
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param col_grouping string or vector; name of the column/s that should be used as grouping variable/s (Label or c(Factor1, Factor2))
#' @param col_value string; name of the column containing values (dependent variable) (LengthMM)
#' @return data.frame; containing summary statistics for each Factor1 - Factor2 combination
#' @examples
#' # normalize Rootdetection Output
#' root_norm <- norm_10mm_standard(root_output)
#'
#' # calculate summary statistic of standard Rootdetection data set
#' summary_stat(root_norm)
#'
#' # calculate summary statistic grouped by Factor 1 only
#' summary_stat(root_norm, col_grouping = "Factor1", col_value = "LengthMM")
#' @export
summary_stat <- function(root_norm, col_grouping = c("Factor1", "Factor2"), col_value = "LengthMM") {
  if (!col_value %in% colnames(root_norm)) {
    stop(paste0("Column ", col_value, " not found in data.frame"))
  }

  colnames(root_norm)[colnames(root_norm) == col_value] <- "kleeche"
  kleeche <- NULL
  sum <- plyr::ddply(root_norm, col_grouping, plyr::summarize,
    n = length(kleeche), median = stats::median(kleeche),
    mean = mean(kleeche), sd = stats::sd(kleeche),
    se = se(kleeche)
  )
  return(sum)
}



#' @title Normality Test For Normalized Rootdetection Data Set
#' @description The function performs a Shapiro-Wilk test of normality of dependent variable (LengthMM) for each grouping variable (Label).
#' A data.frame containing p-values obtained from the Shapiro-Wilk test for each grouping variable (Label) is returned.
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param col_grouping string; name of the column that should be used as grouping variable (Label)
#' @param col_value string; name of the column containing values (dependent variable) (LengthMM)
#' @return data.frame; containing p-values for each grouping variable
#' @examples
#' # normalize Rootdetection Output
#' root_norm <- norm_10mm_standard(root_output)
#'
#' normality_test(root_norm)
#' @export
normality_test <- function(root_norm, col_grouping = "Label", col_value = "LengthMM") {
  # create empty dataframe
  norm_table <- data.frame(Label = character(0), p.value = numeric(0))

  # rename column in root_norm according to grouping
  colnames(root_norm)[colnames(root_norm) == col_grouping] <- "Label"
  colnames(root_norm)[colnames(root_norm) == col_value] <- "LengthMM"

  root_norm$Label <- as.factor(root_norm$Label)

  for (lev in 1:length(levels(root_norm$Label))) {
    # create subset contain only the data for one label
    hist_sub <- root_norm[root_norm$Label == levels(root_norm$Label)[lev], ]
    temp <- data.frame(Label = levels(root_norm$Label)[lev], p.value = stats::shapiro.test(hist_sub$LengthMM)$p.value)
    norm_table <- rbind(norm_table, temp)
  }

  colnames(norm_table)[1] <- col_grouping
  return(norm_table)
}



#' @title Calculate Relative Data Of Normalized Rootdetection Data Set
#' @description
#' The input must be a Normalized Rootdetection Data Set with columns Factor1, Factor2 and LengthMM present.
#' The function calculates values of grouping variable 2 (Factor2) treatments relative to the median of grouping variable 2 (Factor2) control.
#' Therefore for each grouping variable 1 (Factor1) the median of grouping variable 2 (Factor2) control is used to calculate relative values.
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param control string; name of the grouping variable 2 (Factor2) control condition
#' @return data.frame; containing data for each grouping variable 1 (Factor1) and grouping variable 2 (Factor2) treatments relative to grouping variable 2 (Factor2) control
#' @examples
#' # normalize Rootdetection Output
#' root_norm <- norm_10mm_standard(root_output)
#'
#' rel_data(root_norm)
#' @export
rel_data <- function(root_norm, control = "20") {
  # create subset containing only control data
  # (name of the control was assigned in the beginning)
  rel_table_mock <- root_norm[root_norm$Factor2 == control, ]
  # calc median for all levels of Factor 1 an save to LengthMM_median_control
  # rel_table_mock_median <- plyr::ddply(rel_table_mock, plyr::.(Factor1),
  #   plyr::summarize,
  #   LengthMM_median_control =
  #     stats::median(LengthMM)
  # )
  rel_table_mock_median <- median_rel(rel_table_mock)

  # merge tables
  rel_table_merge <- merge(root_norm, rel_table_mock_median, by = "Factor1")
  # new subset without the mock data
  rel_table <- rel_table_merge[rel_table_merge$Factor2 != control, ]
  # calculate relative values for the new 'treatment-only' table
  rel_table$relative_value <-
    (100 * rel_table$LengthMM) / rel_table$LengthMM_median_control
  rel_root_table <- rel_table[, c(
    "Label", "Factor1", "Factor2", "RootNr",
    "relative_value"
  )]

  return(rel_root_table)
}



#' @title Remove Outlier from Normalized Rootdetection Data Set
#' @description The function removes or replaces outliers for every grouping variable (Label) present in a Normalized Rootdetection Data Set.
#' According to outlier definition in box plots the function defines outliers as 1.5 x IQR.
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param col_grouping string; name of the column that should be used as grouping variable (Label)
#' @param col_value string; name of the column containing values (dependent variable) (LengthMM)
#' @param fill_na logical; If TRUE all outliers present in col_value will be replaced by NA. If FALSE all outliers will be removed.
#' @return data.frame; containing data without outliers or outliers replaced by NA
#' @examples
#' # normalize Rootdetection Output
#' root_norm <- norm_10mm_standard(root_output)
#'
#' # transform outliers to NA
#' rm_outlier_df(root_norm, fill_na = TRUE)
#' # remove outliers
#' rm_outlier_df(root_norm, fill_na = FALSE)
#' @export
rm_outlier_df <- function(root_norm, col_grouping = "Label", col_value = "LengthMM", fill_na = F) {
  colnames(root_norm)[colnames(root_norm) == col_grouping] <- "Label"
  colnames(root_norm)[colnames(root_norm) == col_value] <- "LengthMM"

  # remove all NAs before from LengthMM
  root_norm <- root_norm[!is.na(root_norm$LengthMM), ]

  # create empty data-frame
  final <- root_norm[0, ]

  # loop over labels in data.frame
  for (name in unique(root_norm$Label)) {
    # extract subsets and replace outliers with NA
    temp <- root_norm[root_norm$Label == name, ]
    temp$LengthMM <- rm_outlier(temp$LengthMM, fill_na = T)

    final <- rbind(final, temp)
  }

  if (fill_na == F) {
    final <- final[!is.na(final$LengthMM), ]
  }

  colnames(final)[colnames(final) == "Label"] <- col_grouping
  colnames(final)[colnames(final) == "LengthMM"] <- col_value

  return(final)
}
