% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ANOVA.R
\name{onefacaov_fac1}
\alias{onefacaov_fac1}
\title{One-way ANOVA Over Grouping Variable 1}
\usage{
onefacaov_fac1(
  root_norm,
  col_grouping1 = "Factor1",
  col_grouping2 = "Factor2",
  col_value = "LengthMM",
  summary_plots = F,
  draw_out = F,
  file_base = "1fac_ANOVA_factor1",
  axis_label_size = 0.7,
  p_value_size = 0.8
)
}
\arguments{
\item{root_norm}{data.frame; normalized Rootdetection data set}

\item{col_grouping1}{string; name of the column that should be used as grouping variable 1 (Factor1)}

\item{col_grouping2}{string; name of the column that should be used as grouping variable 2 (Factor2)
to iterate over, can be set to NULL if not existing}

\item{col_value}{string; name of the column containing values (dependent variable) (LengthMM)}

\item{summary_plots}{logical; If TRUE summary (plot(aov)) will be plotted}

\item{draw_out}{logical; If TRUE a matrix containing p-values is plotted in a pdf file}

\item{file_base}{string; file base name of the pdf output (is needed if draw_out = TRUE)}

\item{axis_label_size}{numeric; font size of axis labels in pdf file (if draw_out = TRUE)}

\item{p_value_size}{numeric; font size of p-values in pdf file (if draw_out = TRUE)}
}
\value{
list; data.frames containing p-values of one-way ANOVA and Tukey post-hoc
}
\description{
The function performs a one-way ANOVA and Tukey post-hoc test over grouping variable 1 (Factor1).
The function iterates over a second grouping variable (Factor2) if given.
}
\examples{
### Usage Standard Rootdetection ###

# get data.frame containg p-values
root_norm <- norm_10mm_standard(root_output)
onefacaov_fac1(root_norm, draw_out = FALSE)

# plot p-value matrix as pdf output

root_test_norm <- norm_10mm_standard(root_output)
onefacaov_fac1(root_test_norm, draw_out = TRUE, file_base = "1fac_ANOVA_factor1")
# function creates two pdf files: 1fac_ANOVA_factor1_20.pdf and 1fac_ANOVA_factor1_28.pdf


### Usage for table containing only a single grouping variable ###

# produce example data set containing only a single grouping variable
root_single_var <- root_test_norm[root_test_norm$Factor2 == "20", ]
root_single_var$Factor2 <- NULL
# rename some columns
colnames(root_single_var)[colnames(root_single_var) == "LengthMM"] <- "length"
colnames(root_single_var)[colnames(root_single_var) == "Factor1"] <- "lines"

# use onefacaov_fac1 to conduct one-way ANOVA
onefacaov_fac1(root_single_var, col_grouping1 = "lines", col_grouping2 = NULL, col_value = "length")
}
