#' Rootdetection Output of an example study
#'
#' A dataset containing measured _Arabidopsis thaliana_ wild-type (Col-0) and mutant (tir1_1afb2_3, weitar1_1, yucOx) seedling roots measured at 20 °C and 28°C. Plants were grown on ATS plates for 4 days at 20°C and then stayed at 20°C or were shifted to 28°C. The columns labeled with '10mm' contains the length standard in colum LengthPx.
#'
#' @format A data frame with 125 rows of  6 variables:
#' \describe{
#'   \item{Nr}{Number given by Rootdetection}
#'   \item{Filename}{Filename of the photo}
#'   \item{RootNr}{Number given by Rootdetection}
#'   \item{Label}{Label describing the line}
#'   \item{LengthPx}{Length of measured root or Hypocotyl in Pixel}
#'   \item{LengthMM}{Length of measured root or Hypocotyl in mm}
#' }
"root_output"



#' Rootdetection Output of an example study with multiple Factor2
#'
#' A dataset containing measured _Arabidopsis thaliana_ wild-type (Col-0 and Ws2) and mutant (D31 and D64) seedling roots measured at 20 °C and 28°C. One part of the plants was grown at 20°C or 28°C for 8 days. Another batch was transferred after 4 days from 20°C to 28°C for another 4 days.
#' The columns labeled with '20mm' contains the length standard in colum LengthPx. In this dataset a different length standard is used (2cm are measured and are labled as 20mm).
#'
#' @format A data frame with 213 rows of  6 variables:
#' \describe{
#'   \item{Nr}{Number given by Rootdetection}
#'   \item{Filename}{Filename of the photo}
#'   \item{RootNr}{Number given by Rootdetection}
#'   \item{Label}{Label describing the line}
#'   \item{LengthPx}{Length of measured root or Hypocotyl in Pixel}
#'   \item{LengthMM}{Length of measured root or Hypocotyl in mm}
#' }
"root_output_multfac2"



#' Picture Of Colorful Primroses
#'
#' Picture of primroses loaded with 'jpeg' package. Picture made by 3268zauber~commonswiki and obtained from Wikimedia Commons.
#'
#' @format A matrix containing color information for each pixel.
#' @examples
#' # show picture in R using grid package
#'
#' grid.raster(primroses)
"primroses"


#' Color Palette For Colorblind Individuals.
#'
#' A String storing a color palette with 7 colors optimized for colorblind individuals in hex-format. The colors were proposed by Wong in following puclication:
#' Wong, B. (2011). Points of view: Color blindness. Nature Methods, 8(6), 441–441. doi:10.1038/nmeth.1618.
#'
#' @format A vector containg hex-codes for 7 colors
"colors_wong"


## IDEAS


# normality test erweitern --> über die Factoren??
# andere tests rein


# is.root... erweitern und verbessern
# eine Funktion mit 10mm und eine ohne??
# check dann in jede Funktion einbauen welche dieses als Input verwendet!


# ggplot2 producing figure captions??? mal versuchen!


# Label delim für twofacaov benötigt??
