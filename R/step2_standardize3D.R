#' Standardize 3D images
#'
#' Standardize topographies so that they have the same lateral resolution (6.25
#' microns), and crop/pad images such that they are 701x701 pixel images (NOTE:
#' first find plane and level, and then only standardize before further
#' pre-processing)
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
	if (scale != 1) {
	    ret <- imager::resize(imager::as.cimg(surfaceMat), -scale * 100, -scale * 100, -100, interpolation_type = 5)[, , , ]
	} else {
	    ret <- imager::as.cimg(surfaceMat)
	}

	# pad with NAs
    dim1 <- dim(ret)
    padded <- matrix(NA, nrow = 701, ncol = 701)
    if (dim1[1] <= 701) {
        row_start <- floor((701 - dim1[1]) / 2) + 1
        row_end <- row_start + dim1[1] - 1
    } else {
        midpoint <- ceiling(dim1[1]/2)
        row_start <- midpoint - 350
        row_end <- midpoint + 350
        ret <- ret[row_start:row_end, ]
        row_start <- 1
        row_end <- 701
    }
    if (dim1[2] <= 701) {
        col_start <- floor((701 - dim1[2]) / 2) + 1
        col_end <- col_start + dim1[2] - 1
    } else {
        midpoint <- ceiling(dim1[2]/2)
        col_start <- midpoint - 350
        col_end <- midpoint + 350
        ret <- ret[, col_start:col_end]
        col_start <- 1
        col_end <- 701
    }
    padded[row_start:row_end, col_start:col_end] <- ret

    return(padded)
}
