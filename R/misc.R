#' Standardize 3D images
#'
#' Standardize topographies so that they have the same lateral resolution (6.25
#' microns), and center and crop images (NOTE: first find plane and then
#' standardize before further pre-processing) (TODO center and crop)
#'
#' @param surfaceMat matrix of depth values in microns
#' @param lateralResMicrons lateral resolution of image
#' @return Standardized image for further pre-processing.
#' @examples
#' \dontrun{
#'     testImage <- standardize3D(testImage, lateralResMicrons)
#' }
#' @export
#'
standardize3D <- function(surfaceMat, lateralResMicrons){
    # standardize the resolution to 6.25 which is the lowest resolution available
    # ret <- imager::imresize(imager::as.cimg(ret), scale = .5, interpolation = 5)[, , , ]
    # ret <- EBImage::resize(ret, w = floor(dim(ret)[1]/2), h = floor(dim(ret)[2]/2))
    scale <- lateralResMicrons/6.25 # make everything 6.25
	if (scale != 1) ret <- imager::resize(imager::as.cimg(surfaceMat), -scale * 100, -scale * 100, -100, interpolation_type = 5)[, , , ]
    return(ret)
}




#' Plot an image
#'
#' @param image matrix of pixel values to be plotted
#' @param type either \code{"original"} for images on the original scale (pixel
#'   values 0-255), \code{"residuals"} for residual pixel values (-255 to 255),
#'   or \code{"any"} for plotting using the range available (i.e. the minimum
#'   value is plotted as black and the maximum is white, with a linear scale for
#'   intermediate values). \code{"original"} is the default.
#' @param grayscale logical value indicating whether or not the grayscale is to
#'   be plotted. \code{FALSE} is the default.
#' @param main title for plots with grayscale. The default is an empty string.
#'
#' @details This function plots an image in the same orientation as the matrix
#' of pixel values, i.e. the top-left pixel of the image is the first entry of
#' the matrix, and the bottom-right pixel is the last entry.
#'
#' @examples
#' \dontrun{
#' plotImage(processedExample, type = "any")
#' }
#' @export


plotImage <- function(image, type = "original", grayscale = FALSE, main = "") {

    numRow <- nrow(image)
    numCol <- ncol(image)
    if (grayscale == FALSE) {
        par(mar = c(0, 0, 0, 0))
        plot(c(1, numCol), c(1, numRow), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", yaxs = "i", xaxs = "i")
        if (type == "original") {
            rasterImage(image/255, 1, 1, numCol, numRow, interpolate = FALSE)
        } else if (type == "residuals") {
            rasterImage((image + 255)/510, 1, 1, numCol, numRow, interpolate = FALSE)
        } else if (type == "any") {
            rasterImage((image - min(image, na.rm = TRUE))/(max(image, na.rm = TRUE) - min(image, na.rm = TRUE)), 1, 1, numCol, numRow, interpolate = FALSE)
        }
    } else if (grayscale == TRUE) {
        if (!requireNamespace("fields", quietly = TRUE)) {
            stop("Package fields required for grayscale. Please install it or plot with grayscale = FALSE", call. = FALSE)
        }
        par(mar = c( 2.1, 2.1, 6.1, 6.1))
        plot(c(1, numCol), c(1, numRow), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", yaxs = "i", xaxs = "i", main = main)
        if (type == "original") {
            rasterImage(as.raster(image/255), 1, 1, numCol, numRow, interpolate = FALSE)
            fields::image.plot(legend.only = TRUE, zlim = c(0, 255), col = gray(0:255/255))
        } else if (type == "residuals") {
            rasterImage(as.raster((image + 255)/510), 1, 1, numCol, numRow, interpolate = FALSE)
            fields::image.plot(legend.only = TRUE, zlim = c(-255, 255), col = gray(0:255/255), axis.args = list(at = seq(-250, 250, 50), labels = seq(-250, 250, 50), cex.axis = 0.8))
        } else if (type == "any") {
            rasterImage((image - min(image, na.rm = TRUE))/(max(image, na.rm = TRUE) - min(image, na.rm = TRUE)), 1, 1, numCol, numRow, interpolate = FALSE)
            fields::image.plot(legend.only = TRUE, zlim = c(min(image, na.rm = TRUE), max(image, na.rm = TRUE)), col = gray(0:255/255), cex.axis = 0.8)
        }
    }
}
