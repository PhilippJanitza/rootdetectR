#' @title 1 Factorial ANOVA Over Factor1
#' @description The function performes a one way ANOVA over Factor1. The function iterates over Factor 2.
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @param draw_out logical; If TRUE Matrix containg p-values is plotted in pdf file
#' @param file_base string; file base name of the pdf output (is needed if draw_out = T)
#' @return list; data.frames containg p-values for one way ANOVA over Factor 1
#' @examples
#' # get data.frame containg p-values for one way ANOVA over Factor 1
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' onefacaov_fac1(root_test_norm, draw_out = F)
#'
#' # get data.frame and plot as pdf output
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' onefacaov_fac1(root_test_norm, draw_out = T, file_base = '1fac_ANOVA_factor1')
#' # create two pdf files: 1fac_ANOVA_factor1_20.pdf and 1fac_ANOVA_factor1_28.pdf
#' @export
onefacaov_fac1 <- function(root_norm, draw_out = T,
                           file_base = "1fac_ANOVA_factor1") {

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
                #
                idx <- which(paste(fac1[k], "-", fac1[j], sep = "") == rn)
                #
                if (length(idx) != 0) {
                  mat[j, k] <- tuk[idx, 5]
                }
            }
        }
        matlist[[tp]] <- mat

        if (draw_out) {
            pdf(file = paste(file_base, "_", fac2[tp], ".pdf", sep = ""))
            # Draw the Tuky-p-value-Matrix
            col <- matrix("black", nrow = nr_fac1, ncol = nr_fac1)
            col[lower.tri(col, diag = TRUE)] <- "white"
            image(t(mat[nr_fac1:1, ]), col = c("red", "white"),
                  breaks = c(0, 0.05, 1), axes = FALSE)
            title(main = paste(fac2[tp], ", pval.aov = ",
                               round(pval.aov, 5),
                               sep = ""))
            s <- seq(from = 0, to = 1, length = nr_fac1)
            axis(1, at = s, labels = fac1, las = 2)
            axis(2, at = s, labels = rev(fac1), las = 2)
            abline(h = c(s - (s[2] - s[1]) / 2, s[nr_fac1] + (s[2] - s[1]) / 2),
                   v = c(s - (s[2] - s[1]) / 2, s[nr_fac1] + (s[2] - s[1]) / 2),
                   col = "grey")
            sg <- expand.grid(rev(s), s)
            text(sg[, 2], sg[, 1], format(c(as.vector(mat), 0.123456789),
                                          digits = 2, nsmall = 3,
                                          scientiffic = FALSE)[1:nr_fac1 ^ 2],
                                          col = as.vector(col), cex = 0.8)
            dev.off()
        }
    }

    names(matlist) <- unique(root_norm$Factor2)
    return(matlist)
}



#' @title 1 Factorial ANOVA Over Factor2
#' @description The function performes a one way ANOVA over Factor2. The function takes each control and treatment combination of Factor2 and perfomes a one way ANOVA.
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @param control string; name of the Factor2 control condition
#' @param draw_out logical; If TRUE Matrix containg p-values is plotted in pdf file
#' @param file_base string; file base name of the pdf output (is needed if draw_out = T)
#' @return list; data.frames containg p-values for one way ANOVA over Factor2
#' @examples
#' # get data.frame containg p-values for one way ANOVA over Factor 2
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' onefacaov_fac2(root_test_norm, control = '20', draw_out = F)
#'
#' # get data.frame and plot as pdf output
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' onefacaov_fac2(root_test_norm, control = '20', draw_out = T, file_base = '1fac_ANOVA_factor2')
#' # create a pdf file 1fac_ANOVA_factor2_28.pdf
#' @export
# 1 factorial ANOVA over factor2
onefacaov_fac2 <- function(root_norm, control = '20',  draw_out = T,
                           file_base = "1fac_ANOVA_factor2") {

    # only last mat is returned --> list of dfs!!!

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

            # add mat to list with name
            name <- fac2[tp]
            matl[[name]] <- mat

            if (draw_out) {
                # Visualize the p-values
                col <- matrix("black", nrow = 1, ncol = nr_fac1)
                col[lower.tri(col, diag = TRUE)] <- "white"
                # filename
                pdf(file = paste(file_base, "_", fac2[tp], ".pdf", sep = ""),
                    width = 6, height = 2.5)
                image(t(mat), col = c("red", "white"), breaks = c(0, 0.05, 1),
                      axes = FALSE)
                title(main = fac2[tp])
                s <- seq(from = 0, to = 1, length = nr_fac1)
                axis(1, at = s, labels = fac1, las = 2)
                axis(2, at = s[1], labels = "pval", las = 2)
                abline(v = c(s - (s[2] - s[1]) / 2, s[nr_fac1] +
                               (s[2] - s[1]) / 2), col = "grey")
                text(s, s[1], round(as.vector(mat), 3), cex = 0.8)
                dev.off()
            }
        }
    }
    return(matl)
}



#' @title Two Way ANOVA
#' @description The function performes a two way ANOVA for Rootdetection standard
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @param label_delim character; defin how Factor1 and Factor2 are seperated in Label
#' @param draw_out logical; If TRUE Matrix containg p-values is plotted in pdf file
#' @param file string; file name of the pdf output (is needed if draw_out = T)
#' @return list; data.frames containg p-values for one way ANOVA over Factor2
#' @examples
#' # get data.frame containg p-values for one way ANOVA over Factor 2
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' twofacaov(root_test_norm, label_delim = ';', draw_out = F)
#'
#' # get data.frame and plot as pdf output
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' twofacaov(root_test_norm, label_delim = ';', draw_out = T, file = '2fac_ANOVA_all_vs_all.pdf')
#' @export
twofacaov <- function(root_norm, label_delim = ";", draw_out = T,
                      file = "2fac_ANOVA_all_vs_all.pdf") {

    # model that compares all vs all
    aov_all_vs_all <- aov(LengthMM ~ Factor1 * Factor2, data = root_norm)
    # tukey post hoc test
    tuk <- as.data.frame(TukeyHSD(aov_all_vs_all,
                                  ordered = FALSE)$`Factor1:Factor2`)
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
                               gsub(label_delim, ":", labs[j]), sep = "") ==
                               rn)
            if (length(idx) == 0) {
                idx <- which(paste(gsub(label_delim, ":", labs[j]), "-",
                                   gsub(label_delim, ":", labs[k]), sep = "")
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
        image(t(mat[nr_labs:1, ]), col = c("red", "white"),
              breaks = c(0, 0.05, 1), axes = FALSE)
        title(main = "ANOVA all vs all")
        s <- seq(from = 0, to = 1, length = nr_labs)
        axis(1, at = s, labels = labs, las = 2, cex.axis = 0.6)
        axis(2, at = s, labels = rev(labs), las = 2, cex.axis = 0.6)
        abline(h = c(s - (s[2] - s[1]) / 2, s[nr_labs] + (s[2] - s[1]) / 2), v =
                 c(s - (s[2] - s[1]) / 2, s[nr_labs] + (s[2] - s[1]) / 2),
                 col = "grey")
        sg <- expand.grid(rev(s), s)
        text(sg[, 2], sg[, 1], format(c(as.vector(mat), 0.123456789),
                                      digits = 2, nsmall = 3,
                                      scientiffic = FALSE)[1:nr_labs ^ 2],
                                      col = as.vector(col), cex = 0.8)
        dev.off()
    }
    return(mat)
}


#' @title Pairwise Two Way ANOVA of Treatment Effects
#' @description The function performes a two way ANOVA. It compares the treatment effect of Factor2 control to Factor2 treatment for every Factor1 (pairwise). P-values are corrected using Benjamini-Hochberg procedure.
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @param control string; name of the Factor2 control condition
#' @param label_delim character; defin how Factor1 and Factor2 are seperated in Label
#' @param draw_out logical; If TRUE Matrix containg p-values is plotted in pdf file
#' @param file_base string; file name of the pdf output (is needed if draw_out = T)
#' @return list; data.frames containg p-values for pairwise two way ANOVA for each Factor1 per Factor2 control treatment effect
#' @examples
#' # get data.frame containg p-values for pairwise two way ANOVA
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' pairwise_2facaov(root_test_norm, control = '20', label_delim = ';', draw_out = F)
#'
#' # get data.frame and plot as pdf output
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' twofacaov(root_test_norm, control = '20', label_delim = ';', draw_out = T, file = '2fac_ANOVA_BH_corrected')
#' @export
pairwise_2facaov <- function(root_norm, control = "20", label_delim = ";",
                             draw_out = T, file_base =
                             "2fac_ANOVA_BH_corrected") {

    root_norm$Factor1 <- as.factor(root_norm$Factor1)
    root_norm$Factor2 <- as.factor(root_norm$Factor2)


    # get all levels of factor2
    fac2 <- levels(as.factor(root_norm$Factor2))
    # get level position of control and treatment of Factor2
    con_position <- which(unique(root_norm$Factor2) == control)
    treat_position <- which(unique(root_norm$Factor2) != control)
    # get all levels of Factor1
    fac1 <- levels(root_norm$Factor1)
    # hoe many levels are there in Factor1
    nr_fac1 <- length(fac1)


    dfl <- list()

    for (tp in treat_position) {
        # create empty matrix with NA / as factor1 names to rows and colums
        mat <- matrix(NA, nrow = nr_fac1, ncol = nr_fac1)
        rownames(mat) <- fac1
        colnames(mat) <- fac1
        # loop over all pairs of Factor1 i = 1 - length(Factor1-1)

        for (i in 1:(nr_fac1 - 1)) {
            for (j in (i + 1):nr_fac1) {
                # create data.frame with 2 ecotypes | 1. control 2. from loops
                d <- root_norm[which( (root_norm$Factor1 ==
                                           fac1[i] | root_norm$Factor1
                                           == fac1[j]) &
                                           (root_norm$Factor2
                                           == fac2[con_position] |
                                           root_norm$Factor2 ==
                                           fac2[tp])), ]
                # definiere linear model for ANOVA
                lm_all <- lm(log2(LengthMM) ~ Factor1 * Factor2,
                             d, na.action = na.omit)
                # get matrix with all the p-values!!
                mat[i, j] <- anova(aov(lm_all))[3, 5]
            }
        }

        # adjust p-values
        pvals <- mat[upper.tri(mat)]
        mat_adj <- matrix(NA, nrow = nr_fac1, ncol = nr_fac1)
        rownames(mat_adj) <- fac1
        colnames(mat_adj) <- fac1
        adjp <- multtest::mt.rawp2adjp(pvals, proc = "BH")
        idx <- order(adjp$index)
        mat_adj[upper.tri(mat_adj)] <- adjp$adjp[idx, "BH"]
        mat <- mat_adj
        if (draw_out) {
            # Visualisiation of P-Value-Matrix font color
            col <- matrix("black", nrow = nr_fac1, ncol = nr_fac1)
            col[lower.tri(col, diag = TRUE)] <- "white"
            # filename for pdf
            pdf(file <- paste(file_base, '_', control, "_vs_", fac2[tp],
                              ".pdf", sep = ""))
            # draw matrix draw white and red cells
            image(t(mat[nr_fac1:1, ]), col = c("red", "white"),
                  breaks = c(0, 0.05, 1), axes = FALSE)
            title(main = paste("treatment effect ", control, " vs ", fac2[tp]))
            # axis label
            s <- seq(from = 0, to = 1, length = nr_fac1)
            axis(1, at = s, labels = fac1, las = 1, cex.axis = 0.8)
            axis(2, at = s, labels = rev(fac1), las = 1, cex.axis = 0.8)
            # grid lines
            abline(h = c(s - (s[2] - s[1]) / 2, s[nr_fac1] + (s[2] - s[1]) / 2),
                   v = c(s - (s[2] - s[1]) / 2, s[nr_fac1] + (s[2] - s[1]) / 2),
                   col = "grey")
            # print p-Werte into it
            sg <- expand.grid(rev(s), s)
            text(sg[, 2], sg[, 1], cex = 0.8, format(c(as.vector(mat),
                 0.123456789), digits = 1, nsmall = 3,
                 scientiffic = FALSE)[1:nr_fac1 ^ 2], col = as.vector(col))
            dev.off()
        }

        # Generate dataframes with p-values for each condition (for each loop)
        # and put it in list dfl

        # get names from row and col
        rowcol <- expand.grid(rownames(mat), colnames(mat))
        # take only the upper.tri
        labs <- rowcol[as.vector(upper.tri(mat)), ]
        # add condition to labels
        labs[, 1] <- paste(labs[, 1], unique(root_norm$Factor2)[tp],
                           sep = label_delim)
        labs[, 2] <- paste(labs[, 2], unique(root_norm$Factor2)[tp],
                           sep = label_delim)
        # bind with values
        df <- cbind(labs, mat[upper.tri(mat)])
        # order data frame
        df <- df[, c(2, 1, 3)]
        # combine cells to new cell called name
        df <- tidyr::unite(df, names, c(Var2, Var1), sep = "-", remove = TRUE)
        colnames(df)[2] <- "p.value"

        name <- fac2[tp]
        dfl[[name]] <- df
    }
    return(dfl)
}
