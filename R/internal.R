#' calculate fivenum for y coordinae of plotting significant letters
#' @param root_norm normalized root dataset
#' @export
#' @noRd
fivenum_root <- function(root_norm, values = "LengthMM") {
  sum_df <- data.frame(
    Label = character(0),
    y = numeric(0)
  )

  for (e in unique(root_norm$Label)) {
    t_vals <- root_norm[root_norm$Label == e, ][, values]
    t_df <- data.frame(
      Label = e,
      y = stats::fivenum(t_vals)[5]
    )
    sum_df <- rbind(sum_df, t_df)
  }

  return(sum_df)
}

#' calculate fivenum for y coordinae of plotting significant letters
#' @param root_norm normalized root dataset
#' @export
#' @noRd
median_rel <- function(rel_table_mock) {
  sum_df <- data.frame(
    Factor1 = character(0),
    LengthMM_median_control = numeric(0)
  )

  for (e in unique(rel_table_mock$Factor1)) {
    t_vals <- rel_table_mock[rel_table_mock$Factor1 == e, ]$LengthMM
    t_df <- data.frame(
      Factor1 = e,
      LengthMM_median_control = stats::median(t_vals)
    )
    sum_df <- rbind(sum_df, t_df)
  }

  return(sum_df)
}
