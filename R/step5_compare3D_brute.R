# ACKNOWLEDGEMENTS: the first two functions were adapted from MATLAB code
# provided by Joseph Roth of Michigan State University. The MATLAB code was
# worked on by Joseph Roth and Zach Richardson, as part of their Breech Face
# Ballistics work.

#' Compute correlation for two images
#'
#' Given two processed (and possibly rotated) images, compute the
#' cross-correlation function and return the largest correlation value and its
#' associated translation parameters. (This assumes images have mean zero)
#'
#' @param image1 first image matrix
#' @param image2 second image matrix. Note that the two matrices have to be the
#'   same size.
#'
#' @return A list with three items: 1) maximum correlation taking into account
#'   all possible translations, 2) dx, the horizontal translation corresponding
#'   to the maximum correlation, and 3) dy, the vertical translation
#'   corresponding to the maximum correlation.
#' @examples
#' \dontrun{
#' comparison(processedExample, processedExample2)
#' }
#'
#' @export

comparison <- function(image1, image2) {
    resp <- filterViaFFT(image1, image2) / (sqrt(sum(image1^2)) * sqrt(sum(image2^2)))
    corr <- max(resp)
    tmp <- which(resp == corr, arr.ind = TRUE)[1, ]
    d_offset <- floor(dim(image1)/2)
    dx <- tmp[["col"]] - d_offset[2] - 1
    dy <- -(tmp[["row"]] - d_offset[1] - 1)

    ret <- list(corr = corr, dx = dx, dy = dy)
    return(ret)
}

#' Compute maximum correlation for two images, taking into account a sequence of
#' rotation angles and all integer-valued translations
#'
#' Given two processed images, the images are first resized to 1/4 of their
#' original dimension. The smaller image is zero-padded to make the two images
#' the same size. Then the second image is rotated through all rotation angles
#' 2.5 degrees apart, and for each rotation we compute the maximum correlation
#' between the two images, considering all integer-valued translations in the x-
#' and y-dimensions. We then choose the rotation angle that results in the
#' largest correlation. We then rotate 0.5 degrees around this best rotation
#' angle, and again compute the maximum correlation. This serves as a similarity
#' metric for the two input images.
#'
#' @param image1 first image matrix
#' @param image2 second image matrix. The two images should be the same size.
#'   (update 9/28: input images are all 701 x 701)
#'
#' @return A list with four items: 1) maximum correlation taking into account
#'   many possible rotations and translations of the second image, 2) the
#'   corresponding dx, the horizontal translation, 3) dy, the vertical
#'   translation, and 4) theta, the corresponding rotation angle.
#' @examples
#' \dontrun{
#' calculateCCFmax(processedExample, processedExample2)
#' }
#'
#' @export

calculateCCFmaxNew <- function(image1, image2) {
    image1_small <- EBImage::resize(image1, w = floor(dim(image1)[1]/4), h = floor(dim(image1)[2]/4))
    image2_small <- EBImage::resize(image2, w = floor(dim(image2)[1]/4), h = floor(dim(image2)[2]/4))

    thetas <- seq(from = -177.5, to = 180, by = 2.5)
    allResults <- data.frame(thetas = thetas, dx = NA, dy = NA, corr = NA)

    for (i in 1:length(thetas)) {
        #if (i %% 10 == 0) cat(i, ", ")
        rotated <- cartridges3D:::bilinearInterpolation(image2_small, thetas[i])
        ## remove mean from image1_small and rotated
        mean1 <- mean(image1_small[image1_small != 0])
        image1_small[image1_small != 0] <- image1_small[image1_small != 0] - mean1
        mean2 <- mean(rotated[rotated != 0])
        rotated[rotated != 0] <- rotated[rotated != 0] - mean2
        ####### left with rounding error
        out <- cartridges3D:::comparison(image1_small, rotated)
        allResults[i, ] <- c(thetas[i], out$dx, out$dy, out$corr)
    }

    fineThetas <- seq(from = thetas[which.max(allResults$corr)] - 2, to = thetas[which.max(allResults$corr)] + 2, by = .5)
    fineResults <- data.frame(fineThetas = fineThetas, dx = NA, dy = NA, corr = NA)

    for (i in 1:length(fineThetas)) {
        #if (i %% 10 == 0) cat(i, ", ")
        rotated <- cartridges3D:::bilinearInterpolation(image2_small, fineThetas[i])
        ## remove mean from image1_small and rotated
        mean1 <- mean(image1_small[image1_small != 0])
        image1_small[image1_small != 0] <- image1_small[image1_small != 0] - mean1
        mean2 <- mean(rotated[rotated != 0])
        rotated[rotated != 0] <- rotated[rotated != 0] - mean2
        ####### left with rounding error
        out <- cartridges3D:::comparison(image1_small, rotated)
        fineResults[i, ] <- c(fineThetas[i], out$dx, out$dy, out$corr)
    }
    index <- which.max(fineResults$corr)
    ret <- list(corr = fineResults$corr[index], dx = fineResults$dx[index], dy = fineResults$dy[index], theta = fineResults$fineThetas[index])

    return(ret)
}


bilinearInterpolation <- function(inputMatrix, thetaDegs) {
    dimension <- nrow(inputMatrix) # assume square
    center <- floor(dimension/2) # floor in case dim is odd
    theta <- thetaDegs/180*pi
    cosT <- cos(theta)
    sinT <- sin(theta)
    minX <- 1 - center
    minY <- center - nrow(inputMatrix)
    maxX <- ncol(inputMatrix) - center
    maxY <- center - 1

    dataframe <- expand.grid(1:dimension, 1:dimension, KEEP.OUT.ATTRS = FALSE)
    names(dataframe) <- c("row", "col")
    dataframe$pixelValue <- inputMatrix[cbind(dataframe$row, dataframe$col)]

    dataframe$x <- dataframe$col - center
    dataframe$y <- center - dataframe$row
    dataframe$newX <- cosT*dataframe$x - sinT*dataframe$y
    dataframe$newY <- sinT*dataframe$x + cosT*dataframe$y
    dataframe$x1 <- floor(dataframe$newX)
    dataframe$x2 <- ceiling(dataframe$newX)
    dataframe$y1 <- floor(dataframe$newY)
    dataframe$y2 <- ceiling(dataframe$newY)
    dataframe$dx <- dataframe$newX - dataframe$x1
    dataframe$dy <- dataframe$newY - dataframe$y1

    indicesToFill <- which(!(dataframe$x1 < minX | dataframe$x2 >= maxX | dataframe$y1 < minY | dataframe$y2 >= maxY))
    dataframe$newGrid <- 0
    dataframe$newGrid[indicesToFill] <- (1 - dataframe$dy[indicesToFill])*((1 - dataframe$dx[indicesToFill]) * inputMatrix[cbind(center - dataframe$y1[indicesToFill], dataframe$x1[indicesToFill] + center)] + (dataframe$dx[indicesToFill]) * inputMatrix[cbind(center - dataframe$y1[indicesToFill], dataframe$x2[indicesToFill] + center)]) + dataframe$dy[indicesToFill] * ((1 - dataframe$dx[indicesToFill]) * inputMatrix[cbind(center - dataframe$y2[indicesToFill], dataframe$x1[indicesToFill] + center)] + (dataframe$dx[indicesToFill]) * inputMatrix[cbind(center - dataframe$y2[indicesToFill], dataframe$x2[indicesToFill] + center)])

    outputMatrix <- matrix(dataframe$newGrid, nrow = dimension)
    return(outputMatrix)
}

