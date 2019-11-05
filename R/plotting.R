#' @title Plotting Histograms
#' @description The function produces histogram plots for each grouping varibale 1 and 2 combinations (Factor1 and Factor2).
#' In addition a Shapiro-Wilk test of normaility is performed and p-value is plotted in the subtititle of the plots.
#' Also there is colour coding - green = normaly distributed according to Shapiro-Wilk test; red = not normaly distributed according to Shapiro-Wilk test
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param draw_out logical; TRUE = prints plot in pdf file, FALSE = no pdf output
#' @param file string; filename of the pdf output
#' @return list; containg the p-values of shapiro-wilk test for each element
#' @examples
#' # produce list of plots:
#'
#' root_norm <- norm_10mm_standard(root_output)
#' plot_hist(root_norm, draw_out = F)
#' # view the first plot in R:
#' plot_hist[[1]]
#'
#' # produce list of plots and print as pdf
#' root_norm <- norm_10mm_standard(root_output)
#' plot_hist(root_norm, draw_out = T, file = 'test_plot_histograms.pdf')
#' # view the first plot in R:
#' plot_hist[[1]]
#' @export
plot_hist <- function(root_norm, draw_out = F,
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

        #calculate Shapiro-Wilk
        shap <- shapiro.test(hist_sub$LengthMM)

        #calculate breaks
        brx <- pretty(range(hist_sub$LengthMM),
                      n = nclass.Sturges(hist_sub$LengthMM), min.n = 1)

        p <- ggplot2::ggplot(hist_sub, ggplot2::aes_string(x =
                                                               hist_sub$LengthMM)) +
            {if(shap$p.value >= 0.05)ggplot2::geom_histogram(color = "black", fill = "green",
                                                             breaks = brx)} +
            {if(shap$p.value < 0.05)ggplot2::geom_histogram(color = "black", fill = "red",
                                                            breaks = brx)} +
            ggplot2::theme_classic() +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
                           plot.subtitle =  ggplot2::element_text(hjust = 0.5, size = 8)) +
            ggplot2::scale_x_continuous("LengthMM") +
            ggplot2::ggtitle(label = paste(unique(hist_sub$Factor1),
                                           unique(hist_sub$Factor2, sep = " ")),
                             subtitle =  paste('p.val =', round(shap$p.value, digits = 4)))

        plot_list[[lev]] <- p
    }

    if (draw_out == T) {
        pdf(file)

        fp <- length(plot_list) %/% 12
        x <- 1
        if(fp > 0){
            for(i in 1:fp){


                gridExtra::grid.arrange(grobs = plot_list[x:(x + 11)], ncol = 3, nrow = 4)
                x <- x + 12

            }
        }

        lp <- length(plot_list) %% 12
        if(lp >= 1){

            gridExtra::grid.arrange(grobs = plot_list[x:(x + (lp - 1))], ncol = 3, nrow = 4)

        }


        dev.off()
    }
    return(plot_list)
}



#' @title Plotting Absolute Data Of Normalized Rootdetection Data Set
#' @description Absolute data are plotted as box plot, box and jitter plot combination or violin plot.
#' Signifcances can be illustrated by letter encoding. There are a lot of possibilities to adjust the visualization.
#' @param root_norm data.frame; normalized Rootdetection data set
#' @param label_delim character; define how Factor1 and Factor2 are separated in Label
#' @param type string; "box" = will produce box plot, 'jitter' = will produce combination of box and jitter plot,  'violin' = will produce violin plot
#' @param size_jitter_dot number; defines the size of the dots in jitter plots
#' @param plot_n logical; if TRUE sample size (n) will be plotted
#' @param plot_colours vector; provide colours for the plot - color vector must have the same length than grouping variable 2 (Factor2)
#' @param plot_title string; sets a plot title
#' @param size_plot_title numeric; defines size of the plot title
#' @param y_label string; axes description of y-axes
#' @param x_label string; axes description of x-axes
#' @param legend_label string; title of legend
#' @param size_legend_title numeric; defines size of legend title
#' @param size_legend_text numeric; defines size of legend text
#' @param width_lines numeric; defines the line width of boxes or violins
#' @param width_axis numeric; defines the line width of the axis
#' @param size_x_axes numeric; defines size of x-axes labels
#' @param size_y_axes numeric; defines size of y-axes labels
#' @param size_n = numeric; defines size of n plotted (if plot_n = TRUE)
#' @param plot_significance logical; if TRUE significance will be drawn as letters
#' @param height_letter numeric; defines the position of the significance letters
#' @param size_letter numeric; defines size of the significance letters
#' @param angle_letter numeric, defines angle of significance letters
#' @return plot; box, box jitter or violin plot of the absolute data
#' @examples
#' # plot without significance letters
#'
#' root_norm <- norm_10mm_standard(root_output)
#' # box plot
#' plot_abs(root_norm, type = 'box')
#' # box and jitter plot
#' plot_abs(root_norm, type = 'jitter')
#' # violin plot
#' plot_abs(root_norm, type = 'violin')
#' # provide own colours (standard r colours or hex-code) - must have the same length than grouping variable 2 level(root_test_norm$Factor2)
#' plot_abs(root_norm, plot_colours = c('dodgerblue3', 'firebrick2'))
#'
#' # plot with significance letters
#'
#' root_norm <- norm_10mm_standard(root_output)
#' plot_abs(root_norm, plot_significance = T)
#'
#' # all customizable plotting parameters have a default value
#'  plot_abs_short(root_test_norm,
#'                 label_delim = ';',
#'                 type = 'jitter', size_jitter_dot = 2,
#'                 plot_n = T,
#'                 plot_colours = c('cornflowerblue', 'coral2'),
#'                 plot_title = 'My Plot', size_plot_title = 14,
#'                 y_label = 'length mm', x_label = '',
#'                 legend_label = 'condition', size_legend_title = 12, size_legend_text = 10,
#'                 width_lines = 0.5, width_axis = 0.5,
#'                 size_x_axes = 9, size_y_axes = 9,
#'                 size_n = 3,
#'                 plot_significance = T, height_letter = 5, size_letter = 5, angle_letter = 0)
#' @export
plot_abs <- function(root_norm,
                     label_delim = ";",
                     type = "jitter",
                     size_jitter_dot = 2,
                     plot_n = T,
                     plot_colours,
                     plot_title,
                     size_plot_title = 14,
                     y_label = 'hypocotyl length [mm]',
                     x_label = '',
                     legend_label = 'condition',
                     size_legend_title = 12,
                     size_legend_text = 10,
                     width_lines = 0.5,
                     width_axis = 0.5,
                     size_x_axes = 9,
                     size_y_axes = 9,
                     size_n = 3,
                     plot_significance = F,
                     height_letter = 2,
                     size_letter = 5,
                     angle_letter = 0) {

    ############## first chunk creating letters from statistics

    # create table with letters only if plot_significance = T
    if (plot_significance) {

        # conduct ANOVA
        twofacov_output <- rootdetectR::twofacaov(root_norm, label_delim = label_delim)

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
        y_coord$y <- y_coord$y + height_letter
        plot_letters <- merge(letters, y_coord, by = "Label")
    }


    abs_plot <- ggplot2::ggplot() +
        {if(type == 'jitter')ggplot2::geom_boxplot(data = root_norm,
                                                   ggplot2::aes(x = root_norm$Label,
                                                                y = root_norm$LengthMM),
                                                   outlier.shape = NA, lwd = width_lines)} +
        {if(type == 'jitter')ggplot2::geom_jitter(data = root_norm,
                                                  ggplot2::aes(x = root_norm$Label,
                                                               y = root_norm$LengthMM,
                                                               colour = as.factor(root_norm$Factor2)),
                                                  position = ggplot2::position_jitter(0.2,
                                                                                      seed = 1), size = size_jitter_dot)} +
        {if(type == 'jitter')ggplot2::labs(colour = legend_label)} +
        {if(type == 'box')ggplot2::geom_boxplot(data = root_norm,
                                                ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM,
                                                             fill = as.factor(root_norm$Factor2)), lwd = width_lines)} +
        {if(type == 'box')ggplot2::labs(fill = legend_label)} +
        {if(type == 'violin')ggplot2::geom_violin(data = root_norm,
                                                   ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM,
                                                                fill = as.factor(root_norm$Factor2)), lwd = width_lines)} +
        {if(type == 'violin')ggplot2::labs(fill = legend_label)} +
        {if(plot_significance)ggplot2::geom_text(data = plot_letters,
                                                 ggplot2::aes(x = plot_letters$Label,
                                                              y = plot_letters$y,
                                                              label = plot_letters$Letters,
                                                              angle = angle_letter),
                                                 size = size_letter)} +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(name = y_label) +
        ggplot2::scale_x_discrete(name = x_label) +
        {if(!missing(plot_title))ggplot2::ggtitle(label = plot_title)} +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                           hjust = 1, vjust = 1,
                                                           colour = 'black',
                                                           size = size_x_axes),
                       axis.text.y = ggplot2::element_text(colour = 'black',
                                                           size= size_y_axes),
                       axis.line = ggplot2::element_line(size = width_axis),
                       legend.title = ggplot2::element_text(size = size_legend_title),
                       legend.text = ggplot2::element_text(size = size_legend_text),
                       plot.title = ggplot2::element_text(hjust = 0.5, size = size_plot_title)) +
        {if(plot_n)EnvStats::stat_n_text(data = root_norm,
                                         ggplot2::aes(x = root_norm$Label, y = root_norm$LengthMM),
                                         angle = 90, size = size_n)} +
        {if(!missing(plot_colours))ggplot2::scale_fill_manual(values = plot_colours)} +
        {if(!missing(plot_colours))ggplot2::scale_colour_manual(values = plot_colours)}

    return(abs_plot)
}













#' @title Plotting Relative Data Of Normalized Rootdetection Data Set
#' @description Relative data is plotted as boxplot or box-jitter-plot combination. If Significances should be illustrated multiple plots are generated. Each Factor2 control Factor2 treatment combination will produce a plot.
#' @param root_norm data.frame; LengthMM normalized output from Rootdetection containing NO 10mm values
#' @param label_delim character; define how Factor1 and Factor2 are seperated in Label
#' @param control string; name of the Factor2 control condition
#' @param type string; "box" = will produce Boxplot, 'jitter' = will produce combination of box and jitter plot, 'violin' = will produce violin plot
#' @param size_jitter_dot number; defines the size of the dots in jitter plots
#' @param plot_colours vector; provide colours for the boxplot - colour vector must have the same length than Factor1
#' @param plot_title string; defines the plot title
#' @param size_plot_title numeric; defines size of the plot title
#' @param y_label string; axes description of y-axes
#' @param x_label string; axes description of x-axes
#' @param legend_label string; title of legend
#' @param size_legend_title numeric; defines size of legend title
#' @param size_legend_text numeric; defines size of legend text
#' @param width_lines numeric; defines the line width of boxes or violins
#' @param width_axis numeric; defines the line width of the axis
#' @param size_x_axes numeric; defines size of x-axes labels
#' @param size_y_axes numeric; defines size of y-axes labels
#' @param plot_significance logical; if TRUE significances will be drawn as letters
#' @param height_letter numeric; defines the position of the significance letters
#' @param size_letter numeric; defines size of the significance letters
#' @param angle_letter numeric, defines angle of significance letters
#' @return plot or list of plots; box or box-jitter-plot of relative data, if significance = T and multiple Factor2 treatments a list of boxplots will be generated
#' @examples
#' # Plot without significance
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' # boxplot
#' plot_rel(root_test_norm, plot_significance = F, control = '20', type = 'box')
#' # jitterplot
#' plot_rel(root_test_norm, plot_significance = F, control = '20', type = 'jitter')
#'
#' # Plot with siginficance letters
#'
#' root_test_norm <- norm_10mm_standard(root_output)
#' root_stat <- interaction_twofacaov(root_test_norm, draw_out = F, label_delim = ';')
#' # boxplot with statistics
#' plot_rel(root_test_norm, plot_significance = T, interaction_twofacaov_output = root_stat, control = '20', type = 'box')
#'
#' # all customizable plotting parameters have a default value
#' plot_rel(root_test_norm,
#'          label_delim = ';',
#'          control = 20,
#'          type = "jitter",
#'          plot_colours = c('blue', 'red', 'orange', 'green'),
#'          y_label = '% growth',
#'          x_label = '',
#'          legend_label = 'Label',
#'          size_legend_title = 12,
#'          size_legend_text = 10,
#'          size_x_axes = 9,
#'          size_y_axes = 9,
#'          plot_significance = T,
#'          interaction_twofacaov_output = root_stat,
#'          height_letter = 10,
#'          size_letter = 5,
#'          angle_letter = 0)
#'
#' @export
plot_rel <- function(root_norm,
                     label_delim = ';',
                     control = 20,
                     type = "jitter",
                     size_jitter_dot = 2,
                     plot_colours,
                     plot_title,
                     size_plot_title = 14,
                     y_label = 'hypocotyl growth [%]',
                     x_label = '',
                     legend_label = 'Label',
                     size_legend_title = 12,
                     size_legend_text = 10,
                     width_lines = 0.5,
                     width_axis = 0.5,
                     size_x_axes = 9,
                     size_y_axes = 9,
                     plot_significance = F,
                     height_letter = 10,
                     size_letter = 5,
                     angle_letter = 0) {

    rel_root_norm <- rootdetectR::rel_data(root_norm, control = control)

    if (plot_significance == T) {

        # conduct ANOVA
        interaction_twofacaov_output <- rootdetectR::interaction_twofacaov(root_norm,
                                                                           control = control,
                                                                           label_delim = label_delim)

        relative_plot <- list()

        for (i in 1:length(interaction_twofacaov_output)) {

            # create rel subsets
            rel_sub <- subset(rel_root_norm,
                              Factor2 == names(interaction_twofacaov_output)[i])
            rel_sub$Factor2 <- as.factor(rel_sub$Factor2)
            rel_sub$Factor2 <- droplevels(rel_sub$Factor2)



            # create subset of ANOVA matrices
            sub_aov <- interaction_twofacaov_output[[i]]

            # get Letters from ANOVA subset
            mat_names <- character()
            mat_values <- numeric()
            # loop over matrix and get names + values
            for (j in 1:(length(row.names(sub_aov)) - 1)) {
                for (k in (j + 1):length(colnames(sub_aov))) {
                    v <- sub_aov[j, k]
                    t <- paste(row.names(sub_aov)[j],
                               colnames(sub_aov)[k], sep = "-")
                    mat_names <- c(mat_names, t)
                    mat_values <- c(mat_values, v)
                }
            }

            # combine names + values
            names(mat_values) <- mat_names
            # get df with letters and replace : with label delim!!
            letters <-
                data.frame(multcompView::multcompLetters(mat_values)["Letters"])
            letters$Label <- paste(rownames(letters), names(interaction_twofacaov_output)[i], sep = label_delim)


            # letter 1/5 of highest value
            y_coord <- plyr::ddply(rel_sub, plyr::.(Label), plyr::summarize,
                                   y = fivenum(relative_value)[5])
            y_coord$y <- y_coord$y + height_letter
            plot_letters <- merge(letters, y_coord, by = "Label")




            # use annotate instead of geom_text() to use within a loop
            relative_plot_temp <- ggplot2::ggplot() +
                {if(type == 'jitter')ggplot2::geom_boxplot(data = rel_sub,
                                                           ggplot2::aes_string(x = rel_sub$Label, y =
                                                                                   rel_sub$relative_value), outlier.shape = NA,
                                                           lwd = width_lines)} +
                {if(type == 'jitter')ggplot2::geom_jitter(data = rel_sub,
                                                          ggplot2::aes_string(x = rel_sub$Label, y =
                                                                                  rel_sub$relative_value, colour = as.factor(rel_sub$Label)),
                                                          position = ggplot2::position_jitter(0.2, seed = 1), size = size_jitter_dot)} +
                {if(type == 'jitter')ggplot2::labs(colour = legend_label)} +
                {if(type == 'box')ggplot2::geom_boxplot(data = rel_sub, ggplot2::aes_string(x =
                                                                                                rel_sub$Label, y = rel_sub$relative_value, fill =
                                                                                                rel_sub$Label),
                                                        lwd = width_lines)} +
                {if(type == 'box')ggplot2::labs(fill = legend_label)} +
                {if(type == 'violin')ggplot2::geom_violin(data = rel_sub,
                                                           ggplot2::aes_string(x = rel_sub$Label, y = rel_sub$relative_value,
                                                                               fill = as.factor(rel_sub$Label)), lwd = width_lines)} +
                {if(type == 'violin')ggplot2::labs(fill = legend_label)} +
                ggplot2::theme_classic() +
                ggplot2::scale_y_continuous(name = y_label) +
                ggplot2::scale_x_discrete(name = x_label) +
                {if(!missing(plot_title))ggplot2::ggtitle(label = plot_title)} +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle =
                                                                       45, hjust = 1, vjust = 1, colour = 'black', size = size_x_axes),
                               axis.text.y = ggplot2::element_text(colour = 'black', size = size_y_axes),
                               axis.line = ggplot2::element_line(size = width_axis),
                               legend.title = ggplot2::element_text(size = size_legend_title),
                               legend.text = ggplot2::element_text(size = size_legend_text),
                               plot.title = ggplot2::element_text(hjust = 0.5, size = size_plot_title)) +
                ggplot2::annotate("text", x = plot_letters$Label,
                                  y = plot_letters$y, label = plot_letters$Letters, size = size_letter,
                                  angle = angle_letter) +
                {if(!missing(plot_colours))ggplot2::scale_fill_manual(values = plot_colours)} +
                {if(!missing(plot_colours))ggplot2::scale_colour_manual(values = plot_colours)}

            relative_plot[[i]] <- relative_plot_temp
        }
    } else if(plot_significance == F) {

        relative_plot <- ggplot2::ggplot() +
            {if(type == 'jitter')ggplot2::geom_boxplot(data = rel_root_norm,
                                                       ggplot2::aes(x = rel_root_norm$Label, y =
                                                                        rel_root_norm$relative_value), outlier.shape = NA,
                                                       lwd = width_lines)} +
            {if(type == 'jitter')ggplot2::geom_jitter(data = rel_root_norm,
                                                      ggplot2::aes(x = rel_root_norm$Label, y =
                                                                       rel_root_norm$relative_value, colour =
                                                                       as.factor(rel_root_norm$Label)), position =
                                                          ggplot2::position_jitter(0.2, seed = 1), size = size_jitter_dot)} +
            {if(type == 'jitter')ggplot2::labs(colour = legend_label)} +
            {if(type == 'box')ggplot2::geom_boxplot(data = rel_root_norm, ggplot2::aes(x =
                                                                                           rel_root_norm$Label, y = rel_root_norm$relative_value,
                                                                                       fill = as.factor(rel_root_norm$Label)),
                                                    lwd = width_lines)} +
            {if(type == 'box')ggplot2::labs(fill = legend_label)} +
            {if(type == 'violin')ggplot2::geom_violin(data = rel_root_norm,
                                                       ggplot2::aes(x = rel_root_norm$Label, y = rel_root_norm$relative_value,
                                                                    fill = as.factor(rel_root_norm$Label)), lwd = width_lines)} +
            {if(type == 'violin')ggplot2::labs(fill = legend_label)} +
            ggplot2::theme_classic() +
            ggplot2::scale_y_continuous(name = y_label) +
            ggplot2::scale_x_discrete(name = x_label) +
            {if(!missing(plot_title))ggplot2::ggtitle(label = plot_title)} +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                               hjust = 1, vjust = 1, colour = 'black', size = size_x_axes),
                           axis.text.y = ggplot2::element_text(colour = 'black', size = size_y_axes),
                           axis.line = ggplot2::element_line(size = width_axis),
                           legend.title = ggplot2::element_text(size = size_legend_title),
                           legend.text =  ggplot2::element_text(size = size_legend_text),
                           plot.title = ggplot2::element_text(hjust = 0.5, size = size_plot_title)) +
            {if(!missing(plot_colours))ggplot2::scale_fill_manual(values = plot_colours)} +
            {if(!missing(plot_colours))ggplot2::scale_colour_manual(values = plot_colours)}


    }
    return(relative_plot)
}
