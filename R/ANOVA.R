#' @title One-way ANOVA Over Grouping Variable 1
#' @description The function performs a one-way ANOVA and Tukey post-hoc test over grouping variable 1 (Factor1). The function iterates over a second grouping variable (Factor2) if given.
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param col_grouping1 string; name of the column that should be used as grouping variable 1 (Factor1)
#' @param col_grouping2 string; name of the column that should be used as grouping variable 2 (Factor2) to iterate over, can be set to NULL if not existing
#' @param col_value  string; name of the column containing values (dependent variable) (LengthMM)
#' @param summary_plots logical; If TRUE summary (plot(aov)) will be plotted
#' @param draw_out logical; If TRUE a matrix containing p-values is plotted in a pdf file
#' @param file_base string; file base name of the pdf output (is needed if draw_out = TRUE)
#' @param axis_label_size numeric; font size of axis labels in pdf file (if draw_out = TRUE)
#' @param p_value_size numeric; font size of p-values in pdf file (if draw_out = TRUE)
#' @return list; data.frames containing p-values of one-way ANOVA and Tukey post-hoc over grouping variable 1 (Factor1)
#' @examples
#' ### Usage Standard Rootdetection ###
#'
#' # get data.frame containg p-values from one-way ANOVA and Tukey post-hoc test over grouping variable 1
#' root_norm <- norm_10mm_standard(root_output)
#' onefacaov_fac1(root_norm, draw_out = F)
#'
#' # plot p-value matrix as pdf output
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' onefacaov_fac1(root_test_norm, draw_out = T, file_base = "1fac_ANOVA_factor1")
#' # function creates two pdf files: 1fac_ANOVA_factor1_20.pdf and 1fac_ANOVA_factor1_28.pdf
#'
#'
#' ### Usage for table containing only a single grouping variable ###
#'
#' # produce example data set containing only a single grouping variable
#' root_single_var <- root_test_norm[root_test_norm$Factor2 == "20", ]
#' root_single_var$Factor2 <- NULL
#' # rename some columns
#' colnames(root_single_var)[colnames(root_single_var) == "LengthMM"] <- "length"
#' colnames(root_single_var)[colnames(root_single_var) == "Factor1"] <- "lines"
#'
#' # use onefacaov_fac1 to conduct one-way ANOVA
#' onefacaov_fac1(root_single_var, col_grouping1 = "lines", col_grouping2 = NULL, col_value = "length")
#' @export
onefacaov_fac1 <- function(root_norm,
                           col_grouping1 = "Factor1", col_grouping2 = "Factor2",
                           col_value = "LengthMM",
                           summary_plots = F,
                           draw_out = F,
                           file_base = "1fac_ANOVA_factor1",
                           axis_label_size = 0.7,
                           p_value_size = 0.8) {
  # col_grouping2 = NULL --> no variable to loop over then:
  if (is.null(col_grouping2)) {
    # rename columns for grouping and dependend vars to match the rootdetection standard
    colnames(root_norm)[colnames(root_norm) == col_grouping1] <- "Factor1"
    colnames(root_norm)[colnames(root_norm) == col_value] <- "LengthMM"

    # prequisites
    fac1 <- levels(as.factor(root_norm$Factor1))
    nr_fac1 <- length(fac1)
    # matlist <- as.list(rep(NA, nr_fac2))

    # ANOVA
    lm_aov <- lm(LengthMM ~ Factor1, root_norm, na.action = na.omit)
    pval.aov <- anova(aov(lm_aov))[1, 5]

    # Tukey
    tuk <- as.data.frame(TukeyHSD(aov(lm_aov), ordered = FALSE)$Factor1)
    tuk <- cbind(rownames(tuk), tuk)
    colnames(tuk)[1] <- "comparison"

    # get tukey comparison and put them in matrix
    mat <- matrix(NA, nrow = nr_fac1, ncol = nr_fac1)
    colnames(mat) <- fac1
    rownames(mat) <- fac1
    rn <- rownames(tuk)
    for (j in 1:(nr_fac1 - 1)) {
      for (k in (j + 1):nr_fac1) {
        idx <- which(paste(fac1[k], "-", fac1[j], sep = "") == rn)
        if (length(idx) != 0) {
          mat[j, k] <- tuk[idx, 5]
        }
      }
    }

    if (draw_out) {
      pdf(file = paste(file_base, ".pdf", sep = ""))

      # Draw the Tuky-p-value-Matrix
      col <- matrix("black", nrow = nr_fac1, ncol = nr_fac1)
      col[lower.tri(col, diag = TRUE)] <- "white"
      image(t(mat[nr_fac1:1, ]),
        col = c("red", "white"),
        breaks = c(0, 0.05, 1), axes = FALSE
      )
      title(main = paste(fac2[tp], ", pval.aov = ",
        round(pval.aov, 5),
        sep = ""
      ))
      s <- seq(from = 0, to = 1, length = nr_fac1)
      axis(1, at = s, labels = fac1, las = 2, cex.axis = axis_label_size)
      axis(2, at = s, labels = rev(fac1), las = 2, cex.axis = axis_label_size)
      abline(
        h = c(s - (s[2] - s[1]) / 2, s[nr_fac1] + (s[2] - s[1]) / 2),
        v = c(s - (s[2] - s[1]) / 2, s[nr_fac1] + (s[2] - s[1]) / 2),
        col = "grey"
      )
      sg <- expand.grid(rev(s), s)
      text(sg[, 2], sg[, 1], format(as.vector(mat),
        digits = 2, nsmall = 3,
        scientific = T
      )[1:nr_fac1^2],
      col = as.vector(col), cex = p_value_size
      )

      if (summary_plots) {
        plot(aov(lm_aov))
      }

      dev.off()
    }


    return(mat)
  } else {
    # there is a second grouping var to loop over then (col_grouping2 != NULL)

    # rename columns for grouping and dependend vars to match the rootdetection standard
    colnames(root_norm)[colnames(root_norm) == col_grouping1] <- "Factor1"
    colnames(root_norm)[colnames(root_norm) == col_grouping2] <- "Factor2"
    colnames(root_norm)[colnames(root_norm) == col_value] <- "LengthMM"

    # prequisites
    fac2 <- levels(as.factor(root_norm$Factor2))
    fac1 <- levels(as.factor(root_norm$Factor1))
    nr_fac2 <- length(fac2)
    nr_fac1 <- length(fac1)
    matlist <- as.list(rep(NA, nr_fac2))

    for (tp in 1:nr_fac2) {
      # create subdataset with only one level of Factor2
      data_aov <-
        root_norm[which(root_norm$Factor2 == fac2[tp]), ]

      # ANOVA
      lm_aov <- lm(LengthMM ~ Factor1, data_aov, na.action = na.omit)
      pval.aov <- anova(aov(lm_aov))[1, 5]

      # Tukey
      tuk <- as.data.frame(TukeyHSD(aov(lm_aov), ordered = FALSE)$Factor1)
      tuk <- cbind(rownames(tuk), tuk)
      colnames(tuk)[1] <- "comparison"

      # get tukey comparison and put them in matrix
      mat <- matrix(NA, nrow = nr_fac1, ncol = nr_fac1)
      colnames(mat) <- fac1
      rownames(mat) <- fac1
      rn <- rownames(tuk)
      for (j in 1:(nr_fac1 - 1)) {
        for (k in (j + 1):nr_fac1) {
          idx <- which(paste(fac1[k], "-", fac1[j], sep = "") == rn)
          if (length(idx) != 0) {
            mat[j, k] <- tuk[idx, 5]
          }
        }
      }

      # save matrix in matlist
      matlist[[tp]] <- mat

      if (draw_out) {
        pdf(file = paste(file_base, "_", fac2[tp], ".pdf", sep = ""))

        # Draw the Tuky-p-value-Matrix
        col <- matrix("black", nrow = nr_fac1, ncol = nr_fac1)
        col[lower.tri(col, diag = TRUE)] <- "white"
        image(t(mat[nr_fac1:1, ]),
          col = c("red", "white"),
          breaks = c(0, 0.05, 1), axes = FALSE
        )
        title(main = paste(fac2[tp], ", pval.aov = ",
          round(pval.aov, 5),
          sep = ""
        ))
        s <- seq(from = 0, to = 1, length = nr_fac1)
        axis(1, at = s, labels = fac1, las = 2, cex.axis = axis_label_size)
        axis(2, at = s, labels = rev(fac1), las = 2, cex.axis = axis_label_size)
        abline(
          h = c(s - (s[2] - s[1]) / 2, s[nr_fac1] + (s[2] - s[1]) / 2),
          v = c(s - (s[2] - s[1]) / 2, s[nr_fac1] + (s[2] - s[1]) / 2),
          col = "grey"
        )
        sg <- expand.grid(rev(s), s)
        text(sg[, 2], sg[, 1], format(as.vector(mat),
          digits = 2, nsmall = 3,
          scientific = TRUE
        )[1:nr_fac1^2],
        col = as.vector(col), cex = p_value_size
        )

        if (summary_plots) {
          plot(aov(lm_aov))
        }

        dev.off()
      }
    }

    names(matlist) <- unique(root_norm$Factor2)
    return(matlist)
  }
}



#' @title One-way ANOVA Over Grouping Variable 2 Treatments
#' @description The function performs a one-way ANOVA and Tukey post-hoc test over grouping variable 2 (Factor2).
#' The function takes each control + treatment combination of grouping variable 2 and perfoms a one-way ANOVA iterating over grouping Variable 1 (Factor1).
#' The p-values are adjusted using Benjamini-Hochberg correction.
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param col_grouping1 string; name of the column that should be used as grouping variable 1 (Factor1)
#' @param col_grouping2 string; name of the column that should be used as grouping variable 2 (Factor2)
#' @param col_value  string; name of the column containing values (dependent variable) (LengthMM)
#' @param control string; name of the grouping variable 2 (Factor2) control condition
#' @param draw_out logical; If TRUE a matrix containing p-values is plotted in a pdf file
#' @param file_base string; file base name of the pdf output (is needed if draw_out = TRUE)
#' @param axis_label_size numeric; font size of axis labels in pdf file (if draw_out = TRUE)
#' @param p_value_size numeric; font size of p-values in pdf file (if draw_out = TRUE)
#' @return list; data.frames containing p-values for one-way ANOVA and Tukey post-hoc over grouping variable 2 (Factor2)
#' @examples
#' # get data.frame containg p-values for one-way ANOVA and Tukey post-hoc over grouping variable 2 (Factor2)
#' root_norm <- norm_10mm_standard(root_output)
#' onefacaov_fac2(root_norm, control = "20", draw_out = F)
#'
#' # plot p-value matrix as pdf output
#' onefacaov_fac2(root_test_norm, control = "20", draw_out = T, file_base = "1fac_ANOVA_factor2")
#' # function creates a pdf file 1fac_ANOVA_factor2_28.pdf
#' @export
onefacaov_fac2 <- function(root_norm,
                           col_grouping1 = "Factor1", col_grouping2 = "Factor2",
                           col_value = "LengthMM",
                           control = "20",
                           draw_out = F,
                           file_base = "1fac_ANOVA_factor2",
                           axis_label_size = 0.7,
                           p_value_size = 0.8) {
  colnames(root_norm)[colnames(root_norm) == col_grouping1] <- "Factor1"
  colnames(root_norm)[colnames(root_norm) == col_grouping2] <- "Factor2"
  colnames(root_norm)[colnames(root_norm) == col_value] <- "LengthMM"

  fac2 <- levels(as.factor(root_norm$Factor2))
  fac1 <- levels(as.factor(root_norm$Factor1))
  nr_fac1 <- length(fac1)
  con_position <- which(unique(root_norm$Factor2) == control)
  treat_position <- which(unique(root_norm$Factor2) != control)

  # create empty list to put in mats for every pair
  matl <- list()
  # loop over treat positions
  for (tp in treat_position) {
    mat <- matrix(NA, nrow = 1, ncol = nr_fac1)
    rownames(mat) <- c("pval")
    colnames(mat) <- fac1

    # loop over factor1
    for (i in 1:nr_fac1) {
      # for each factor1: get control and one single treatment
      d <- root_norm[which(root_norm$Factor1 == fac1[i] &
        (root_norm$Factor2 ==
          fac2[con_position] |
          root_norm$Factor2 ==
            fac2[tp])), ]
      # ANOVA
      lm_all <- lm(LengthMM ~ Factor2, d, na.action = na.omit)
      aov_all <- aov(lm_all)
      mat[1, i] <- anova(aov_all)[1, 5]
    }

    # adjust pval
    mat[1, ] <- multtest::mt.rawp2adjp(mat[1, ], proc = "BH")$adjp[, 2]

    if (draw_out) {
      # Visualize the p-values
      col <- matrix("black", nrow = 1, ncol = nr_fac1)
      col[lower.tri(col, diag = TRUE)] <- "white"
      # filename
      pdf(file = paste(file_base, "_", fac2[tp], ".pdf", sep = ""), width = 6, height = 2.5)
      image(t(mat),
        col = c("red", "white"), breaks = c(0, 0.05, 1),
        axes = FALSE
      )
      title(main = paste(fac2[con_position], "vs.", fac2[tp]))
      s <- seq(from = 0, to = 1, length = nr_fac1)
      axis(1, at = s, labels = fac1, las = 2, cex.axis = axis_label_size)
      axis(2, at = s[1], labels = "pval", las = 2, cex.axis = axis_label_size)
      abline(v = c(s - (s[2] - s[1]) / 2, s[nr_fac1] +
        (s[2] - s[1]) / 2), col = "grey")
      text(s, s[1], format(as.vector(mat),
        digits = 2, nsmall = 3,
        scientific = TRUE
      )[1:nr_fac1^2],
      cex = p_value_size
      )
      dev.off()
    }

    # add mat to list with name
    name <- fac2[tp]
    matl[[name]] <- mat
  }
  return(matl)
}



#' @title Two-Way ANOVA Over Grouping Variable 1 And 2
#' @description The function performs a two way ANOVA and Tukey post-hoc test for normalized Rootdetection standard over grouping variable 1 and 2 (Factor1 and Factor2).
#' Caution: In the actual version it is necessary to provide a Label column (grouping variable 1 and 2 separated by label delimiter (label_delim)) in addition to grouping variable 1 and 2 necessary.
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param col_grouping1 string; column name of the first grouping variable (Factor1)
#' @param col_grouping2 string; column name of the second grouping variable (Factor2)
#' @param col_value string; name of the column containing values (dependent variable) (LengthMM)
#' @param col_label string; column name of label (combining grouping variable 1 and 2 separated by delimiter)
#' @param label_delim character; defines how Factor1 and Factor2 are separated in Label
#' @param summary_plots logical; If TRUE summary plots were printed (plot(aov))
#' @param draw_out logical; If TRUE Matrix containing p-values is plotted in pdf file
#' @param file string; file name of the pdf output (is needed if draw_out = T)
#' @param axis_label_size numeric; font size of axis labels
#' @param p_value_size numeric; font size of the p-values printed in pdf file
#' @return list; data.frames containing p-values for one-way ANOVA and Tukey post-hoc over grouping variable 1 and 2 (Factor1 and Factor2)
#' @examples
#' # get data.frame containing p-values for two-way ANOVA and Tukey post-hoc over grouping variable 1 and 2 (Factor1 and Factor2)
#'
#' root_norm <- norm_10mm_standard(root_output)
#' twofacaov(root_norm, label_delim = ";", draw_out = F)
#'
#' # get data.frame and plot as pdf output
#'
#' root_norm <- norm_10mm_standard(root_output)
#' twofacaov(root_norm, label_delim = ";", draw_out = T, file = "2fac_ANOVA_all_vs_all.pdf")
#' #' # function creates a pdf file 2fac_ANOVA_all_vs_all.pdf
#' @export
twofacaov <- function(root_norm,
                      col_grouping1 = "Factor1", col_grouping2 = "Factor2",
                      col_value = "LengthMM",
                      col_label = "Label",
                      label_delim = ";",
                      summary_plots = F,
                      draw_out = F,
                      file = "2fac_ANOVA_all_vs_all.pdf",
                      axis_label_size = 0.7,
                      p_value_size = 0.8) {
  colnames(root_norm)[colnames(root_norm) == col_grouping1] <- "Factor1"
  colnames(root_norm)[colnames(root_norm) == col_grouping2] <- "Factor2"
  colnames(root_norm)[colnames(root_norm) == col_value] <- "LengthMM"
  colnames(root_norm)[colnames(root_norm) == col_label] <- "Label"

  # model that compares all vs all
  aov_all_vs_all <- aov(LengthMM ~ Factor1 * Factor2, data = root_norm)

  # tukey post hoc test
  tuk <- as.data.frame(TukeyHSD(aov_all_vs_all,
    ordered = FALSE
  )$`Factor1:Factor2`)
  tuk <- cbind(rownames(tuk), tuk)
  # which and how many labels are there
  labs <- levels(root_norm$Label)
  nr_labs <- length(labs)
  # create empty matrix
  mat <- matrix(NA, nrow = nr_labs, ncol = nr_labs)
  colnames(mat) <- labs
  rownames(mat) <- labs
  rn <- rownames(tuk)
  for (j in 1:(nr_labs - 1)) {
    for (k in (j + 1):nr_labs) {
      # get p-values and put them into the matrix
      idx <- which(paste(gsub(label_delim, ":", labs[k]), "-",
        gsub(label_delim, ":", labs[j]),
        sep = ""
      ) ==
        rn)
      if (length(idx) == 0) {
        idx <- which(paste(gsub(label_delim, ":", labs[j]), "-",
          gsub(label_delim, ":", labs[k]),
          sep = ""
        )
        == rn)
      }
      #
      if (length(idx) != 0) {
        mat[j, k] <- tuk[idx, 5]
      }
    }
  }
  if (draw_out) {
    # Draw the Tuky-p-value-Matrix
    col <- matrix("black", nrow = nr_labs, ncol = nr_labs)
    col[lower.tri(col, diag = TRUE)] <- "white"
    pdf(file = file)
    image(t(mat[nr_labs:1, ]),
      col = c("red", "white"),
      breaks = c(0, 0.05, 1), axes = FALSE
    )
    title(main = "ANOVA all vs all")
    s <- seq(from = 0, to = 1, length = nr_labs)
    axis(1, at = s, labels = labs, las = 2, cex.axis = axis_label_size)
    axis(2, at = s, labels = rev(labs), las = 2, cex.axis = axis_label_size)
    abline(
      h = c(s - (s[2] - s[1]) / 2, s[nr_labs] + (s[2] - s[1]) / 2), v =
        c(s - (s[2] - s[1]) / 2, s[nr_labs] + (s[2] - s[1]) / 2),
      col = "grey"
    )
    sg <- expand.grid(rev(s), s)
    text(sg[, 2], sg[, 1], format(as.vector(mat),
      digits = 2, nsmall = 3,
      scientific = TRUE
    )[1:nr_labs^2],
    col = as.vector(col), cex = p_value_size
    )

    # summary plots
    if (summary_plots) {
      plot(aov_all_vs_all)
    }

    dev.off()
  }
  return(mat)
}



#' @title Pairwise Two-Way ANOVA Of Interaction (Treatment) Effects
#' @description The function performs a pairwise two-way ANOVA of interaction effects of Factor2 control to Factor2 treatment for every Factor1.
#' The p-values are adjusted using Benjamini-Hochberg correction procedure.
#' #' Caution: In the actual version it is necessary to provide a Label column (grouping variable 1 and 2 separated by label delimiter (label_delim)) in addition to grouping variable 1 and 2 necessary.
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param col_grouping1 string; column name of the first grouping variable (Factor1)
#' @param col_grouping2 string; column name of the second grouping variable (Factor2)
#' @param col_value string; name of the column containing values (dependent variable) (LengthMM)
#' @param col_label string; column name of label (combining grouping variable 1 and 2 separated by delimiter)
#' @param label_delim character; defines how Factor1 and Factor2 are separated in Label
#' @param control string; name of the grouping variable 2 (Factor2) control condition
#' @param label_delim character; defines how Factor1 and Factor2 are separated in Label
#' @param draw_out logical; If TRUE Matrix containing p-values is plotted in pdf file
#' @param file_base string; file name of the pdf output (is needed if draw_out = T)
#' @param axis_label_size numeric; font size of axis labels
#' @param p_value_size numeric; font size of the p-values printed in pdf file
#' @return list; matrices containg p-values for pairwise two way ANOVA for each Factor1 per Factor2 control treatment effect
#' @examples
#' # get data.frame containing p-values for pairwise two-way ANOVA of interaction effects
#'
#' root_norm <- norm_10mm_standard(root_output)
#' interaction_twofacaov(root_norm, control = "20", label_delim = ";", draw_out = F)
#'
#' # get data.frame and plot as pdf output
#'
#' root_norm <- norm_10mm_standard(root_output)
#' interaction_twofacaov(root_norm, control = "20", label_delim = ";", draw_out = T, file_base = "2fac_ANOVA_BH_corrected")
#' # function creates a pdf file 2fac_ANOVA_BH_corrected_28.pdf
#' @export
interaction_twofacaov <- function(root_norm,
                                  col_grouping1 = "Factor1", col_grouping2 = "Factor2",
                                  col_value = "LengthMM",
                                  control = "20",
                                  label_delim = ";",
                                  draw_out = F,
                                  file_base = "2fac_ANOVA_BH_corrected",
                                  axis_label_size = 0.7,
                                  p_value_size = 0.8) {
  colnames(root_norm)[colnames(root_norm) == col_grouping1] <- "Factor1"
  colnames(root_norm)[colnames(root_norm) == col_grouping2] <- "Factor2"
  colnames(root_norm)[colnames(root_norm) == col_value] <- "LengthMM"


  root_norm$Factor1 <- as.factor(root_norm$Factor1)
  root_norm$Factor2 <- as.factor(root_norm$Factor2)

  # test if values with zero are present
  if (nrow(root_norm[root_norm$LengthMM == 0, ]) >= 1) {
    stop("There are measurements with LengthMM = 0 present")
  }

  # get all levels of factor2
  fac2 <- levels(as.factor(root_norm$Factor2))
  # get level position of control and treatment of Factor2
  con_position <- which(unique(root_norm$Factor2) == control)
  treat_position <- which(unique(root_norm$Factor2) != control)
  # get all levels of Factor1
  fac1 <- levels(root_norm$Factor1)
  # hoe many levels are there in Factor1
  nr_fac1 <- length(fac1)


  matl <- list()

  for (tp in treat_position) {
    # create empty matrix with NA / as factor1 names to rows and colums
    mat <- matrix(NA, nrow = nr_fac1, ncol = nr_fac1)
    rownames(mat) <- fac1
    colnames(mat) <- fac1
    # loop over all pairs of Factor1 i = 1 - length(Factor1-1)

    for (i in 1:(nr_fac1 - 1)) {
      for (j in (i + 1):nr_fac1) {
        # create data.frame with 2 ecotypes | 1. control 2. from loops
        d <- root_norm[which((root_norm$Factor1 ==
          fac1[i] | root_norm$Factor1
        == fac1[j]) &
          (root_norm$Factor2
          == fac2[con_position] |
            root_norm$Factor2 ==
              fac2[tp])), ]
        # definiere linear model for ANOVA
        lm_all <- lm(log2(LengthMM) ~ Factor1 * Factor2,
          d,
          na.action = na.omit
        )
        # get matrix with all the p-values!!
        mat[i, j] <- anova(aov(lm_all))[3, 5]
      }
    }

    # adjust p-values
    if (nr_fac1 > 2) {
      pvals <- mat[upper.tri(mat)]
      mat_adj <- matrix(NA, nrow = nr_fac1, ncol = nr_fac1)
      rownames(mat_adj) <- fac1
      colnames(mat_adj) <- fac1
      adjp <- multtest::mt.rawp2adjp(pvals, proc = "BH")
      idx <- order(adjp$index)
      mat_adj[upper.tri(mat_adj)] <- adjp$adjp[idx, "BH"]
      mat <- mat_adj
    }

    if (draw_out) {
      # Visualisiation of P-Value-Matrix font color
      col <- matrix("black", nrow = nr_fac1, ncol = nr_fac1)
      col[lower.tri(col, diag = TRUE)] <- "white"
      # filename for pdf
      pdf(file <- paste(file_base, "_", control, "_vs_", fac2[tp],
        ".pdf",
        sep = ""
      ))
      # draw matrix draw white and red cells
      image(t(mat[nr_fac1:1, ]),
        col = c("red", "white"),
        breaks = c(0, 0.05, 1), axes = FALSE
      )
      title(main = paste("treatment effect ", control, " vs ", fac2[tp]))
      # axis label
      s <- seq(from = 0, to = 1, length = nr_fac1)
      axis(1, at = s, labels = fac1, las = 1, cex.axis = axis_label_size)
      axis(2, at = s, labels = rev(fac1), las = 1, cex.axis = axis_label_size)
      # grid lines
      abline(
        h = c(s - (s[2] - s[1]) / 2, s[nr_fac1] + (s[2] - s[1]) / 2),
        v = c(s - (s[2] - s[1]) / 2, s[nr_fac1] + (s[2] - s[1]) / 2),
        col = "grey"
      )
      # print p-Werte into it
      sg <- expand.grid(rev(s), s)
      text(sg[, 2], sg[, 1], cex = p_value_size, format(as.vector(mat),
        digits = 1, nsmall = 3,
        scientific = TRUE
      )[1:nr_fac1^2], col = as.vector(col))
      dev.off()
    }

    # Generate dataframes with p-values for each condition (for each loop)
    # and put it in list matl
    name <- fac2[tp]
    matl[[name]] <- mat
  }
  return(matl)
}
