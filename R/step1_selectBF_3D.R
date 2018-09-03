#' Finds plane of breechface marks using RANSAC
#'
#' Given input depths (in microns), find best-fitting plane using RANSAC. This
#' should be the plane that the breechface marks are on.
#'
#' @param inputDepths matrix of input depths in microns.
#' @return List object containing fitted plane and selected breechface marks.
#' @examples
#' \dontrun{
#'     set.seed(0)
#'     testImage <- findPlaneRansac(surfaceMat)
#' }
#' @export
#'
findPlaneRansac <- function(inputDepths) {
    maxInliers <- 0

    # sample from this
    locs <- data.frame(which(!is.na(inputDepths), arr.ind = TRUE))
    locs$depth <- inputDepths[!is.na(inputDepths)]

    for (iter in 1:75) { ########### vary this
        samp <- sample(nrow(locs), 3)
        fit <- lm(depth ~ row + col, data = locs[samp, ])
        preds <- predict(fit, locs)

        dists <- abs(preds - locs$depth)
        tmpInliers <- dists < .00001 # 10 microns ########### vary this -- based on range? IQR? sd? check Glock DW621US RP1

        if (sum(tmpInliers) > maxInliers) {
           maxInliers <- sum(tmpInliers)
           inliers <- tmpInliers
        }
    }

    # final coefs only computed using inliers
    fit <- lm(depth ~ row + col, data = locs[inliers, ])
    outImage <- matrix(NA, nrow = nrow(inputDepths), ncol = ncol(inputDepths))
    outImage[cbind(locs$row[inliers], locs$col[inliers])] <- locs$depth[inliers]

    ret <- list(fit = fit, outImage = outImage)
    return(ret)
}

#' Levels image
#'
#' Levels image using plane fitted in \code{findPlaneRansac}.
#'
#' @param ransacOutput list output from \code{findPlaneRansac}
#' @return Leveled image
#' @examples
#' \dontrun{
#'     testImage <- levelBF3D(testImage)
#' }
#' @export

levelBF3D <- function(ransacOutput) {
    fittedPlane <- ransacOutput$outImage
    fittedPlane[!is.na(fittedPlane)] <- ransacOutput$fit$fitted.values
    # then take residuals
    leveled <- ransacOutput$outImage - fittedPlane

    return(leveled)
}
