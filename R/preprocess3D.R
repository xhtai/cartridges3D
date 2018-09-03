#' Perform all pre-processing steps
#'
#' Images must be x3p images containing 3D topographies. This function performs
#' all the pre-processing steps, consisting of \itemize{ \item{Automatically
#' selecting the breechface marks} \item{Leveling the image} \item{Removing
#' circular symmetry} \item{Filtering.} }
#'
#' @param fileName location of image file.
#'
#' @return A matrix of processed depth values. Non-breechface pixels are set to
#'   0.
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export

allPreprocess3D <- function(fileName) {
    readIn <- x3ptools::read_x3p(fileName)

    if (readIn$header.info$incrementY < 10^(-5)) { # ISU data is already in microns
        surfaceMat <- readIn$surface.matrix
        lateralResMicrons <- readIn$header.info$incrementY*10^6
    } else {
        surfaceMat <- readIn$surface.matrix/10^6
        lateralResMicrons <- readIn$header.info$incrementY
    }

    set.seed(0)
    testImage <- findPlaneRansac(surfaceMat)
    testImage <- levelBF3D(testImage)

    # now standardize
    testImage <- standardize3D(testImage, lateralResMicrons)

    testImage <- removeCircular3D(testImage)

    nonBF <- is.na(testImage)
    testImage[is.na(testImage)] <- 0
    testImage <- gaussianFilter3D(testImage, nonBF)

    return(testImage)
}
