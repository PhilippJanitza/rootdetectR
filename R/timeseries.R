#' @title Calculate growth rates of rootdetection time series data
#' @description
#' @param root_timeseries data.frame; output *.csv (list) from Rootdetection
#' @param remove_na logcial; TRUE or FALSE
#' @return
#' @examples
#' is_root_output(root_output)
#' @export
growth_rates <- function(root_timeseries, remove_na = F) {

  data_rate <- root_timeseries
  data_rate$id <- paste0(data_rate$RootNr, data_rate$Factor1, data_rate$Factor2)
  data_rate$Rate <- NA
  for(root in unique(data_rate$id)){

    # zwischenschritt df um sicherzustellen dass die Reihenfolge stimmt!
    df_t <- data_rate[data_rate$id == root,]
    df_t <- df_t[order(df_t$Timepoint),]

    data_rate[data_rate$id == root,]$Rate <- c(NA, diff(df_t$LengthMM))

  }

  if(remove_na) {
    data_rate <- data_rate[!is.na(data_rate$Rate),]
  }

  data_rate$id <- NULL
  return(data_rate)

}


#' @title
#' @description
#' @param root_timeseries data.frame; output *.csv (list) from Rootdetection
#' @param remove_na logcial; TRUE or FALSE
#' @return
#' @examples
#' is_root_output(root_output)
#' @export
plot_timeseries <- function(root_timeseries,
                            stats = "mean",
                            smooth = F,
                            smooth_method = "gam",
                            smooth_step = F,
                            disp = "sd",
                            error = "ribbon",
                            points = F,
                            individual = F,
                            x_breaks = NULL,
                            y_label = "growth rate",
                            x_label = "timepoint",
                            plot_title,
                            size_plot_title = 14,
                            plot_colours,
                            size_x_axes = 9,
                            size_y_axes = 9,
                            size_axis_title = 14,
                            width_axis = 0.5,
                            legend_title,
                            size_legend_title = 12,
                            size_legend_text = 10
) {

  # stats = mean or median
  if(!stats %in% c("mean", "median")) {
    warning(paste0("Statistical measure ", stats, " not found --> using mean instead"))
    stats <- "mean"
  }

  # disp = sd or se
  if(!disp %in% c("sd", "se")) {
    warning(paste0("Statistical dispersion ", disp, " not found --> using standard deviation instead"))
    disp <- "sd"
  }

  if(!smooth_method %in% c("gam", "lm", "glm", "gam", "loess")) {
    warning(paste0("Method ", smooth_method, " not found --> using gam instead"))
    smooth_method <- "gam"
  }



  sum_df <- summary_stat(root_timeseries, col_grouping = c("Factor1", "Factor2", "Timepoint"), col_value = "LengthMM")

  time_plot <- ggplot2::ggplot() +
    {
      if(error == "ribbon") {
        ggplot2::geom_ribbon(data = sum_df, ggplot2::aes(x = Timepoint, ymin = get(stats) - get(disp), ymax = get(stats) + get(disp), group = interaction(Factor1, Factor2), fill = Factor2), alpha = .2)
      } else if(error == "bar") {
        ggplot2::geom_linerange(data = sum_df, ggplot2::aes(x = Timepoint, ymin = get(stats) - get(disp), ymax = get(stats) + get(disp), group = interaction(Factor1, Factor2), color = Factor2), alpha = .6, na.rm = T)
      }
    } +
    {
      if(smooth) {
        ggplot2::geom_smooth(data = sum_df, ggplot2::aes(x = Timepoint, y = get(stats), group = interaction(Factor1, Factor2), color = Factor2, fill = Factor2), method = smooth_method, span = .3, alpha = .2)
      } else {
        ggplot2::geom_line(data = sum_df, ggplot2::aes(x = Timepoint, y = get(stats), group = interaction(Factor1, Factor2), color = Factor2))
      }
    } +
    {
      if(smooth_step){
        ggplot2::geom_step(data = sum_df, ggplot2::aes(x = Timepoint, y = get(stats), group = interaction(Factor1, Factor2), color = Factor2))
      }
    } +
    {
      if(points){
        ggplot2::geom_point(data = root_timeseries, ggplot2::aes(x = Timepoint, y = LengthMM, group = interaction(RootNr, Factor1, Factor2), color = Factor2), shape = 20, alpha = .2)
      }
    } +
    {
      if(individual) {
        ggplot2::geom_line(data = root_timeseries, ggplot2::aes(x = Timepoint, y = LengthMM, group = interaction(RootNr, Factor1, Factor2), color = Factor2), alpha = .2)
      }
    } +
    {
      if (!missing(plot_title)) ggplot2::ggtitle(label = plot_title)
    } +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        colour = "black",
        size = size_x_axes,
        angle = 45, hjust = 1, vjust = 1 # hier wieder raus und variablen dafür rein!!!
      ),
      axis.text.y = ggplot2::element_text(
        colour = "black",
        size = size_y_axes
      ),
      axis.title = ggplot2::element_text(size = size_axis_title),
      axis.line = ggplot2::element_line(linewidth = width_axis),
      legend.title = ggplot2::element_text(size = size_legend_title),
      legend.text = ggplot2::element_text(size = size_legend_text),
      plot.title = ggplot2::element_text(hjust = 0.5, size = size_plot_title)
    ) +
    ggplot2::scale_y_continuous(name = y_label) +
    ggplot2::scale_x_continuous(name = x_label, n.breaks = x_breaks) +
    {
      if (!missing(plot_colours)) ggplot2::scale_colour_manual(values = plot_colours)
    } +
    {
      if (!missing(plot_colours)) ggplot2::scale_fill_manual(values = plot_colours)
    } +
    {
      if (!missing(legend_title)) ggplot2::labs(color = legend_title, fill = legend_title)
    } +
    ggplot2::expand_limits(y = 0)

  return(time_plot)
}


#' @title
#' @description
#' @param root_timeseries data.frame; output *.csv (list) from Rootdetection
#' @param remove_na logcial; TRUE or FALSE
#' @return
#' @examples
#' is_root_output(root_output)
#' @export
plot_growthrates <- function(root_timeseries,
                             stats = "mean",
                             smooth = F,
                             smooth_method = "gam",
                             smooth_step = F,
                             disp = "sd",
                             error = "ribbon",
                             points = F,
                             individual = F,
                             x_breaks = NULL,
                             y_label = "growth rate",
                             x_label = "timepoint",
                             plot_title,
                             size_plot_title = 14,
                             plot_colours,
                             size_x_axes = 9,
                             size_y_axes = 9,
                             size_axis_title = 14,
                             width_axis = 0.5,
                             legend_title,
                             size_legend_title = 12,
                             size_legend_text = 10) {

  if(!stats %in% c("mean", "median")) {
    warning(paste0("Statistical measure ", stats, " not found --> using mean instead"))
    stats <- "mean"
  }

  if(!disp %in% c("sd", "se")) {
    warning(paste0("Statistical dispersion ", disp, " not found --> using standard deviation instead"))
    disp <- "sd"
  }

  rates <- growth_rates(root_timeseries, remove_na = T)
  rates_sum <- summary_stat(rates, col_grouping = c("Factor1", "Factor2", "Timepoint"), col_value = "Rate")


  ggplot2::ggplot() +
    {
      if(error == "ribbon") {
        ggplot2::geom_ribbon(data = rates_sum, ggplot2::aes(x = Timepoint, ymin = get(stats) - get(disp), ymax = get(stats) + get(disp), group = interaction(Factor1, Factor2), fill = Factor2), alpha = .2)
      } else if(error == "bar") {
        ggplot2::geom_linerange(data = rates_sum, ggplot2::aes(x = Timepoint, ymin = get(stats) - get(disp), ymax = get(stats) + get(disp), group = interaction(Factor1, Factor2), color = Factor2), alpha = .6, na.rm = T)
      }
    } +
    {
      if(smooth) {
        #geom_smooth (alpha=0.3, size=0, span=0.5) +
        #stat_smooth (geom="line", alpha=0.3, size=3, span=0.5)

        ggplot2::geom_smooth(data = rates_sum, ggplot2::aes(x = Timepoint, y = get(stats), group = interaction(Factor1, Factor2), color = Factor2, fill = Factor2), na.rm = T, method = smooth_method, span = 0.3, alpha = .2)
      } else {
        ggplot2::geom_line(data = rates_sum, ggplot2::aes(x = Timepoint, y = get(stats), group = interaction(Factor1, Factor2), color = Factor2), na.rm = T)
      }
    } +
    {
      if(smooth_step){
        ggplot2::geom_step(data = rates_sum, ggplot2::aes(x = Timepoint, y = get(stats), group = interaction(Factor1, Factor2), color = Factor2))
      }
    } +
    {
      if(points){
        geom_point(data = rates, aes(x = Timepoint, y = Rate, group = interaction(RootNr, Factor1, Factor2), color = Factor2), shape = 20, alpha = .2)
      }
    } +
    {
      if(individual) {
        geom_line(data = rates, aes(x = Timepoint, y = Rate, group = interaction(RootNr, Factor1, Factor2), color = Factor2), alpha = .2)
      }
    } +
    {
      if (!missing(plot_title)) ggplot2::ggtitle(label = plot_title)
    } +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        colour = "black",
        size = size_x_axes,
        angle = 45, hjust = 1, vjust = 1 # hier wieder raus und variablen dafür rein!!!
      ),
      axis.text.y = ggplot2::element_text(
        colour = "black",
        size = size_y_axes
      ),
      axis.title = ggplot2::element_text(size = size_axis_title),
      axis.line = ggplot2::element_line(linewidth = width_axis),
      legend.title = ggplot2::element_text(size = size_legend_title),
      legend.text = ggplot2::element_text(size = size_legend_text),
      plot.title = ggplot2::element_text(hjust = 0.5, size = size_plot_title)
    ) +
    ggplot2::scale_y_continuous(name = y_label) +
    ggplot2::scale_x_continuous(name = x_label, n.breaks = x_breaks) +
    {
      if (!missing(plot_colours)) ggplot2::scale_colour_manual(values = plot_colours)
    } +
    {
      if (!missing(plot_colours)) ggplot2::scale_fill_manual(values = plot_colours)
    } +
    {
      if (!missing(legend_title)) labs(color = legend_title, fill = legend_title)
    } +
    ggplot2::expand_limits(y = 0)


}


#' @title
#' @description
#' @param root_timeseries data.frame; output *.csv (list) from Rootdetection
#' @param remove_na logcial; TRUE or FALSE
#' @return
#' @examples
#' is_root_output(root_output)
#' @export
extract_timepoint <- function(root_timeseries, tp) {

  if(!tp %in% unique(root_timeseries$Timepoint)){
    stop(paste0("timepoint ", tp, " not in timeseries data"))
  }


  df <- root_timeseries[root_timeseries$Timepoint == tp,]
  df$Timepoint <- NULL

  df$Label <- gsub(paste0(";",tp),"", as.character(df$Label))

  return(df)

}

