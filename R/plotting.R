#' @title Plotting histograms and performing normailty tests
#' @description The function performes a shapiro-wilk normality test for each Factor1 Factor 2 combination in Rootdetecion standard data.frame and plots a histogram for each containing the p-value from shaprio-wilk test and colour coding - green = normaly distributed; red = not normaly distributed.
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @param draw_out logical; TRUE = prints plot in pdf file, False = no pdf output
#' @param file string; filename of the pdf output
#' @return list; containg the p-values of shapiro-wilk test for each element
#' @examples
#' #Only produce list of plots:
#'
#' root_test_norm <- norm_10mm_standard(root_test)
#' plot_hist(root_test_norm, draw_out = F)
#' #To view the first plot in R:
#' plot_hist[[1]]
#'
#' #Produce plots and print as pdf
#'
#' root_test_norm <- norm_10mm_standard(root_test)
#' plot_hist(root_test_norm, draw_out = T, file = 'test_plot_histograms.pdf')
#' #To view the first plot in R:
#' plot_hist[[1]]
#'
#' @export
plot_hist <- function(root_norm, draw_out = T,
                      file = "data_distribution.pdf") {
    # other tests than shapiro-wilk
    plot_list <- list()

    # make sure Factor1 and Factor 2 are factors
    root_norm$Factor1 <- as.factor(root_norm$Factor1)
    root_norm$Factor2 <- as.factor(root_norm$Factor2)

    for (lev in 1:length(levels(root_norm$Label))) {

        # create subset contain only the data for one label
        hist_sub <- subset(root_norm, Label == levels(root_norm$Label)[lev])
        hist_sub$Factor1 <- droplevels(hist_sub$Factor1)
        # shapiro-wilk test (> 0.05 significant normally distributed =
        # green else (< 0.05) = failure --> red)
        shap <- shapiro.test(hist_sub$LengthMM)
        if (shap$p.value >= 0.05) {

            brx <- pretty(range(hist_sub$LengthMM),
                          n = nclass.Sturges(hist_sub$LengthMM), min.n = 1)

            p <- ggplot2::ggplot(hist_sub, ggplot2::aes_string(x =
                          hist_sub$LengthMM)) +
                 ggplot2::geom_histogram(color = "black", fill = "green",
                          breaks = brx) +
                 ggplot2::theme_classic() +
                 ggplot2::scale_x_continuous("LengthMM") +
                 ggplot2::annotate("text", x = Inf, y = Inf,
                          label = round(shap$p.value, digits = 4), hjust = 1,
                          vjust = 1, size = 3) +
                 ggplot2::ggtitle(paste(unique(hist_sub$Factor1),
                         unique(hist_sub$Factor2, sep = " ")))

        } else {

            brx <- pretty(range(hist_sub$LengthMM),
                          n = nclass.Sturges(hist_sub$LengthMM), min.n = 1)

            p <- ggplot2::ggplot(hist_sub, ggplot2::aes_string(x =
                          hist_sub$LengthMM)) +
                 ggplot2::geom_histogram(color = "black", fill = "red", breaks = brx) +
                 ggplot2::theme_classic() +
                 ggplot2::scale_x_continuous("LengthMM") +
                 ggplot2::annotate("text", x = Inf, y = Inf,
                          label = round(shap$p.value, digits = 4), hjust = 1,
                          vjust = 1, size = 3) +
                 ggplot2::ggtitle(paste(unique(hist_sub$Factor1),
                         unique(hist_sub$Factor2, sep = " ")))
        }
        plot_list[[lev]] <- p
    }

    if (draw_out == T) {
        pdf(file)
        gridExtra::grid.arrange(grobs = plot_list, ncol = 3, nrow = 4)
        dev.off()
    }
    return(plot_list)
}



#' @title Plotting absolute data of Rootdetection standard
#' @description Absolute data are plotted as boxplot or box-jitter-plot combination. Signifcances can be illustrated in the plot.
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @param plot_significance logical; if TRUE significances will be drawn as letters
#' @param twofacaov_output data.frame; Output of twofacaov (needed if plot_significance = T)
#' @param label_delim character; defin how Factor1 and Factor2 are seperated in Label
#' @param type string; "box" = will produce Boxplot, 'jitter' = will produce combination of box and jitter plot
#' @param plot_colours vector; provide colours for the boxplot - colour vector must have the same length than Factor2
#' @param letter_height numeric; defines the position of the significance letters
#' @return plot; box or box-jitter-plot of the absolute data
#' @examples
#' # Plot without significance
#'
#' root_test_norm <- norm_10mm_standard(root_test)
#' # boxplot
#' plot_abs(root_test_norm, plot_significance = F, type = 'box')
#' # jitterplot
#' plot_abs(root_test_norm, plot_significance = F, type = 'jitter')
#'
#' # Plot with siginficance letters
#' root_test_norm <- norm_10mm_standard(root_test)
#' root_stat <- twofacaov(root_test_norm, draw_out = F)
#' # boxplot with statistics
#' plot_abs(root_test_norm, plot_significance = T, twofacaov_output = root_stat, type = 'box')
#'
#' # provide own colours (ggplot2-colours or hex-code) - must have the same length than level(root_test_norm$Factor2)
#' root_test_norm <- norm_10mm_standard(root_test)
#' plot_abs(root_test_norm, plot_significance = F, type = 'jitter', plot_colours = c('dodgerblue3', 'firebrick2'))
#'
#' @export
plot_abs <- function(root_norm, plot_significance = F,
                     twofacov_output, label_delim = ";", type = "jitter",
                     plot_colours, letter_height = 2) {

    # kann hart gekürzt werden indem stat_n einfach angefügt wird ohne den
    # ganzen ggplot2 Befehl neu auszuführen! Style und so!
    # if length(plot_colours) != length(levels(root_norm$Factor2)) --> error!
    # braucht man wirklich label_delim als input???? ist ja schon mal getrennt wurden

    if (plot_significance) {

        # get Letters for twofacaov output
        mat_names <- character()
        mat_values <- numeric()
        # loop over matrix and get names + values
        for (j in 1:(length(row.names(twofacov_output)) - 1)) {
            for (k in (j + 1):length(colnames(twofacov_output))) {
                v <- twofacov_output[j, k]
                t <- paste(row.names(twofacov_output)[j],
                           colnames(twofacov_output)[k], sep = "-")
                mat_names <- c(mat_names, t)
                mat_values <- c(mat_values, v)
            }
        }

        # combine names + values
        names(mat_values) <- mat_names
        # get df with letters and replace : with label delim!!
        letters <-
          data.frame(multcompView::multcompLetters(mat_values)["Letters"])
        letters$Label <- gsub(":", label_delim, rownames(letters))


        # letter 1/5 of highest value
        y_coord <- plyr::ddply(root_norm, plyr::.(Label), plyr::summarize,
                         y = fivenum(LengthMM)[5])
        y_coord$y <- y_coord$y + letter_height
        plot_letters <- merge(letters, y_coord, by = "Label")

        if (type == "jitter") {


            # Input von oben : Farben, Größen
            abs_plot <- ggplot2::ggplot() +
              ggplot2::geom_boxplot(data = root_norm,
                                    ggplot2::aes(x = root_norm$Label,
                                    y = root_norm$LengthMM),
                                    outlier.shape = NA) +
              ggplot2::geom_jitter(data = root_norm,
                            ggplot2::aes(x = root_norm$Label,
                            y = root_norm$LengthMM,
                            colour = as.factor(root_norm$Factor2)),
                            position = ggplot2::position_jitter(0.2,
                            seed = 1)) +
              ggplot2::geom_text(data = plot_letters,
                            ggplot2::aes(x = plot_letters$Label,
                            y = plot_letters$y,
                            label = plot_letters$Letters)) +
              ggplot2::labs(colour = "condition") +
              ggplot2::theme_classic() +
              ggplot2::scale_y_continuous(name = "hypocotyl length [mm]") +
              ggplot2::scale_x_discrete(name = "") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1, vjust = 1), legend.title =
                ggplot2::element_text(size = 16), legend.text =
                ggplot2::element_text(size = 12)) +
              EnvStats::stat_n_text(data = root_norm,
                ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM),
                angle = 90, size = 3)

        } else if (type == "box") {

            abs_plot <- ggplot2::ggplot() +
              ggplot2::geom_boxplot(data = root_norm,
                ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM,
                fill = as.factor(root_norm$Factor2))) +
              ggplot2::geom_text(data = plot_letters,
                ggplot2::aes(x = plot_letters$Label, y = plot_letters$y,
                label = plot_letters$Letters)) +
              ggplot2::labs(fill = "condition") +
              ggplot2::theme_classic() +
              ggplot2::scale_y_continuous(name = "hypocotyl length [mm]") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1, vjust = 1)) +
              EnvStats::stat_n_text(data = root_norm,
                ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM),
                angle = 90, size = 3)

        } else {
            stop("For type only \"box\" or \"jitter\" are yet implemented")
        }

    } else {

        if (type == "jitter") {

            # Input von oben : Farben, Größen
            abs_plot <- ggplot2::ggplot() +
              ggplot2::geom_boxplot(data = root_norm,
                ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM),
                outlier.shape = NA) +
              ggplot2::geom_jitter(data = root_norm,
                ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM,
                colour = as.factor(root_norm$Factor2)), position =
                ggplot2::position_jitter(0.2, seed = 1)) +
              ggplot2::labs(colour = "condition") +
              ggplot2::theme_classic() +
              ggplot2::scale_y_continuous(name = "hypocotyl length [mm]") +
              ggplot2::scale_x_discrete(name = "") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1, vjust = 1), legend.title =
                ggplot2::element_text(size = 16), legend.text =
                ggplot2::element_text(size = 12)) +
              EnvStats::stat_n_text(data = root_norm,
                ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM),
                angle = 90, size = 3)

        } else if (type == "box") {

            abs_plot <- ggplot2::ggplot() +
              ggplot2::geom_boxplot(data = root_norm,
                ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM,
                fill = as.factor(root_norm$Factor2))) +
              ggplot2::labs(fill = "condition") +
              ggplot2::theme_classic() +
              ggplot2::scale_y_continuous(name = "hypocotyl length [mm]") +
              ggplot2::scale_x_discrete(name = "") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1, vjust = 1)) +
              EnvStats::stat_n_text(data = root_norm,
                 ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM),
                 angle = 90, size = 3)

        } else {
            stop("For type only \"box\" or \"jitter\" are yet implemented")
        }

    }

    if (!missing(plot_colours)) {

        abs_plot <- abs_plot +
          ggplot2::scale_fill_manual(values = plot_colours) +
          ggplot2::scale_colour_manual(values = plot_colours)
    }

    return(abs_plot)
}



#' @title Plotting relative data of Rootdetection standard
#' @description Relative data are plotted as boxplot or box-jitter-plot combination. If Significances should be illustrated multiple plots are generated. Each Factor2 control Factor2 treatment combination will produce a plot.
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @param plot_significance logical; if TRUE significances will be drawn as letters
#' @param pairwise_2facaov_output data.frame; Output of pairwise_2facaov (needed if plot_significance = T)
#' @param control string; name of the Factor2 control condition
#' @param type string; "box" = will produce Boxplot, 'jitter' = will produce combination of box and jitter plot
#' @param letter_height numeric; defines the position of the significance letters
#' @return plot or list of plots; box or box-jitter-plot of relative data, if significance = T and multiple Factor2 treatments a list of boxplots will be generated
#' @examples
#' # Plot without significance
#'
#' root_test_norm <- norm_10mm_standard(root_test)
#' # boxplot
#' plot_rel(root_test_norm, plot_significance = F, control = '20', type = 'box')
#' # jitterplot
#' plot_rel(root_test_norm, plot_significance = F, control = '20', type = 'jitter')
#'
#' # Plot with siginficance letters
#'
#' root_test_norm <- norm_10mm_standard(root_test)
#' root_stat <- pairwise_2facaov(root_test_norm, draw_out = F, label_delim = ';')
#' # boxplot with statistics
#' plot_rel(root_test_norm, plot_significance = T, pairwise_2facaov_output = root_stat, control = '20', type = 'box')
#'
#' @export
plot_rel <- function(root_norm, plot_significance = T,  pairwise_2facaov_output,
                     control = "20", type = "jitter", letter_height = 10) {

    rel_root_norm <- rootdetectR::rel_data(root_norm, control = control)

    if (plot_significance == F) {

        if (type == "jitter") {

            relative_plot <- ggplot2::ggplot() +
              ggplot2::geom_boxplot(data = rel_root_norm,
                ggplot2::aes(x = rel_root_norm$Label, y =
                rel_root_norm$relative_value), outlier.shape = NA) +
              ggplot2::geom_jitter(data = rel_root_norm,
                ggplot2::aes(x = rel_root_norm$Label, y =
                rel_root_norm$relative_value, colour =
                as.factor(rel_root_norm$Label)), position =
                ggplot2::position_jitter(0.2, seed = 1)) +
              ggplot2::labs(colour = "ecotypes") +
              ggplot2::theme_classic() +
              ggplot2::scale_y_continuous(name = "% growth") +
              ggplot2::scale_x_discrete(name = "") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1, vjust = 1), legend.title =
                ggplot2::element_text(size = 16), legend.text =
                ggplot2::element_text(size = 12))

        } else if (type == "box") {

            relative_plot <- ggplot2::ggplot() +
              ggplot2::geom_boxplot(data = rel_root_norm, ggplot2::aes(x =
                rel_root_norm$Label, y = rel_root_norm$relative_value,
                fill = as.factor(rel_root_norm$Label))) +
              ggplot2::labs(fill = "ecotypes") +
              ggplot2::theme_classic() +
              ggplot2::scale_y_continuous(name = "% growth") +
              ggplot2::scale_x_discrete(name = "") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1, vjust = 1), legend.title =
                ggplot2::element_text(size = 16), legend.text =
                ggplot2::element_text(size = 12))

        } else {
            stop("For type only box and jitter are allowed")
        }


    } else if (plot_significance == T) {

        if (type == "jitter") {
            # loop over length elements of 2facaov output
            relative_plot <- list()
            for (i in 1:length(pairwise_2facaov_output)) {

                # get labels for every table
                labs <-
                  data.frame(multcompView::multcompLetters(
                  structure(pairwise_2facaov_output[[i]][, 2],
                  names =
                 as.character(pairwise_2facaov_output[[i]][, 1])))["Letters"])
                labs$Label <- rownames(labs)
                y_coord <- plyr::ddply(rel_root_norm, plyr::.(Label),
                             plyr::summarize, y = fivenum(relative_value)[5])
                y_coord$y <- y_coord$y + letter_height
                plot_letter <- merge(labs, y_coord, by = "Label")


                rel_sub <- subset(rel_root_norm,
                                  Factor2 == names(pairwise_2facaov_output)[i])
                rel_sub$Factor2 <- as.factor(rel_sub$Factor2)
                rel_sub$Factor2 <- droplevels(rel_sub$Factor2)

                # use annotate instead of geom_text() to use within a loop
                relative_plot_temp <- ggplot2::ggplot() +
                  ggplot2::geom_boxplot(data = rel_sub,
                    ggplot2::aes_string(x = rel_sub$Label, y =
                    rel_sub$relative_value), outlier.shape = NA) +
                  ggplot2::geom_jitter(data = rel_sub,
                    ggplot2::aes_string(x = rel_sub$Label, y =
                    rel_sub$relative_value, colour = as.factor(rel_sub$Label)),
                    position = ggplot2::position_jitter(0.2, seed = 1)) +
                  ggplot2::labs(colour = "ecotypes") +
                  ggplot2::theme_classic() +
                  ggplot2::scale_y_continuous(name = "% growth") +
                  ggplot2::scale_x_discrete(name = "") +
                  ggplot2::theme(axis.text.x = ggplot2::element_text(angle =
                    45, hjust = 1, vjust = 1), legend.title =
                    ggplot2::element_text(size = 16), legend.text =
                    ggplot2::element_text(size = 12)) +
                  ggplot2::annotate("text", x = plot_letter$Label,
                    y = plot_letter$y, label = plot_letter$Letters)

                relative_plot[[i]] <- relative_plot_temp
            }

        } else if (type == "box") {

            # loop over length elements of output
            relative_plot <- list()
            for (i in 1:length(pairwise_2facaov_output)) {

                # get labels for every table
                labs <- data.frame(multcompView::multcompLetters(
                  structure(pairwise_2facaov_output[[i]][, 2],
                  names =
                  as.character(pairwise_2facaov_output[[i]][, 1])))["Letters"])
                 labs$Label <- rownames(labs)
                y_coord <- plyr::ddply(rel_root_norm, plyr::.(Label),
                  plyr::summarize, y = fivenum(relative_value)[5] + 10)
                plot_letter <- merge(labs, y_coord, by = "Label")


                rel_sub <- subset(rel_root_norm, Factor2 ==
                                    names(pairwise_2facaov_output)[i])
                rel_sub$Factor2 <- as.factor(rel_sub$Factor2)
                rel_sub$Factor2 <- droplevels(rel_sub$Factor2)

                # use annotate instead of geom_text() to use within a loop
                relative_plot_temp <- ggplot2::ggplot() +
                  ggplot2::geom_boxplot(data = rel_sub, ggplot2::aes_string(x =
                    rel_sub$Label, y = rel_sub$relative_value, fill =
                    rel_sub$Label)) +
                  ggplot2::labs(fill = "ecotypes") +
                  ggplot2::theme_classic() +
                  ggplot2::scale_y_continuous(name = "% growth") +
                  ggplot2::scale_x_discrete(name = "") +
                  ggplot2::theme(axis.text.x = ggplot2::element_text(angle =
                    45, hjust = 1, vjust = 1), legend.title =
                    ggplot2::element_text(size = 16), legend.text =
                    ggplot2::element_text(size = 12)) +
                  ggplot2::annotate("text", x = plot_letter$Label,
                          y = plot_letter$y, label = plot_letter$Letters)

                relative_plot[[i]] <- relative_plot_temp
            }
        } else {
            stop("For type only box and jitter are allowed")
        }

    } else {
        stop("plot_significance can only be TRUE or FALSE")
    }

    return(relative_plot)
}
