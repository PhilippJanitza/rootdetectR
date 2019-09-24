#' @title Extract Color From Picture
#' @description The function extractes the n most prominent Colors of a given Picture (jpeg or png format). The picture must be loaded with readPNG() or readJPEG() from the packages 'png' or 'jpeg'.
#' To extract the most prominent colors kmeans clustering algorithm is used. The n-value must be provided by the user (standard is 6).
#' @param image matrix; must contain colour information for each pixel of the picture. (use readPNG or readJPEG functions to obtain the right data structure)
#' @param n numeric; Number of colors to be extracted from the picture
#' @param plot_out logical; If TRUE a plot of the picture containing only the extracted colors will be generated
#' @return data.frame; containing the extracted colors (r,g and b value and the corresponding hex code)
#' @examples
#' # extract the 4 most prominent colors from the example picture primroses
#'
#' extract_colors(primroses, n = 4)
#'
#' # show the picture with the reduced number of colors
#'
#'extract_colors(primroses, n = 4, plot_out = T)
#'
#' @export
extract_color <- function(image, n = 6, plot_out = F){

  # get dataframe from image
  df <- data.frame(red = matrix(image[,,1], ncol=1),
                   green = matrix(image[,,2], ncol=1),
                   blue = matrix(image[,,3], ncol=1))

  # compute the k-means clustering
  set.seed(42)
  K <- kmeans(df,n)
  df$label <- K$cluster

  # get the "center" coloring of K-means clustering
  colors <- data.frame(label = 1:nrow(K$centers),
                       R = K$centers[,"red"],
                       G = K$centers[,"green"],
                       B = K$centers[,"blue"])
  colors$rgb <- rgb(colors$R, colors$G, colors$B)

  if(plot_out == T){

    # merge color codes on to df
    # IMPORTANT: we must maintain the original order of the df after the merge!
    df$order <- 1:nrow(df)
    df <- merge(df, colors)
    df <- df[order(df$order),]
    df$order <- NULL

    # get mean color channel values for each row of the df.
    R_print <- matrix(df$R, nrow=dim(image)[1])
    G_print <- matrix(df$G, nrow=dim(image)[1])
    B_print <- matrix(df$B, nrow=dim(image)[1])

    # reconstitute the segmented image in the same shape as the input image
    image.segmented = array(dim=dim(image))
    image.segmented[,,1] = R_print
    image.segmented[,,2] = G_print
    image.segmented[,,3] = B_print

    # View the result
    grid::grid.raster(image.segmented)

  }

  return(colors)

}



#' @title Show Color Palette
#' @description Plotting a vector containing hex encoded colors. The function will return a plot showing all the colors that are stored in a vector.
#' @param x vector; must contain hex encoded colors
#' @return plot; showing the colors in an example plot
#' @examples
#' # extract the 4 most prominent colors from the example picture primroses
#'
#' primroses_colors <- extract_colors(primroses, n = 4)
#'
#' # show the color palette
#'
#' show_color(primroses_colors$rgb)
#'
#' # show the color palette for colorblind individuals
#' show_color(colors_wong)
#'
#' @export
show_color <- function(x){

  on.exit(par(mar = c(0.5, 0.5, 0.5, 0.5)))

  n <- length(x)
  image(1:n, 1, as.matrix(1:n), col = x,
        ylab = "", xlab = '', xaxt = "n", yaxt = "n", bty = "n")

}
