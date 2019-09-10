#' Rootdetection Output of an example study
#'
#' A dataset containing measured _Arabidopsis thaliana_ wild-type (Col-0) and mutant (tir1_1afb2_3, weitar1_1, yucOx) seedling hypocotyls measured at 20 °C and 28°C. Plants were grown on ATS plates for 4 days at 20°C and then stayed at 20°C or were shifted to 28°C. The columns labeled with '10mm' contains the length standard in colum LengthPx.
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

## IDEAS

# normality test erweitern --> über die Factoren??
# andere tests rein

# ausreiser rauswerfen funktion

# is.root... erweitern und verbessern
# eine Funktion mit 10mm und eine ohne??
# check dann in jede Funktion einbauen welche dieses als Input verwendet!

# Plotting Parameter verfügbar machen
# Text an labs / title
# Schriftgrößen
# Letter height!
# ggplot2 producing figure captions??? mal versuchen!

# Label delim für twofacaov benötigt??
