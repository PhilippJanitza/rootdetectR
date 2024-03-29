% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ANOVA.R
\name{twofacaov}
\alias{twofacaov}
\title{Two-Way ANOVA Over Grouping Variable 1 And 2}
\usage{
twofacaov(
  root_norm,
  col_grouping1 = "Factor1",
  col_grouping2 = "Factor2",
  col_value = "LengthMM",
  col_label = "Label",
  label_delim = ";",
  summary_plots = F,
  draw_out = F,
  file = "2fac_ANOVA_all_vs_all.pdf",
  axis_label_size = 0.7,
  p_value_size = 0.8
)
}
\arguments{
\item{root_norm}{data.frame; normalized Rootdetection data set}

\item{col_grouping1}{string; column name of the first grouping variable (Factor1)}

\item{col_grouping2}{string; column name of the second grouping variable (Factor2)}

\item{col_value}{string; name of the column containing values (dependent variable) (LengthMM)}

\item{col_label}{string; column name of label (combining grouping variable 1 and 2 separated by delimiter)}

\item{label_delim}{character; defines how Factor1 and Factor2 are separated in Label}

\item{summary_plots}{logical; If TRUE summary plots were printed (plot(aov))}

\item{draw_out}{logical; If TRUE Matrix containing p-values is plotted in pdf file}

\item{file}{string; file name of the pdf output (is needed if draw_out = T)}

\item{axis_label_size}{numeric; font size of axis labels}

\item{p_value_size}{numeric; font size of the p-values printed in pdf file}
}
\value{
list; data.frames containing p-values for one-way ANOVA and Tukey post-hoc over grouping
variable 1 and 2 (Factor1 and Factor2)
}
\description{
The function performs a two way ANOVA and Tukey post-hoc test for normalized Rootdetection
standard over grouping variable 1 and 2 (Factor1 and Factor2).
Caution: In the actual version it is necessary to provide a Label column
(grouping variable 1 and 2 separated by label delimiter (label_delim)) in addition to grouping
variable 1 and 2 necessary.
}
\examples{
# get data.frame containing p-values for two-way ANOVA and Tukey post-hoc

root_norm <- norm_10mm_standard(root_output)
twofacaov(root_norm, label_delim = ";", draw_out = FALSE)

# get data.frame and plot as pdf output

root_norm <- norm_10mm_standard(root_output)
twofacaov(root_norm, label_delim = ";", draw_out = TRUE, file = "2fac_ANOVA_all_vs_all.pdf")
#' # function creates a pdf file 2fac_ANOVA_all_vs_all.pdf
}
