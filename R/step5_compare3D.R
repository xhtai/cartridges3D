#' Align two images using the Lucas-Kanade algorithm
#'
#' Finds translation and rotation parameters to approximately minimize the
#' mean-squared error between two images.
#'
#' @param It image at time t, or image that is transformed
#' @param It1 image at time t+1, or image to align to
#'
#' @return homography transformation matrix \code{M}, transformation vector
#'   \code{p} (theta, dx, dy), number of iterations before convergence
#'   \code{iter}
#'
#' @examples
#' \dontrun{
#' # first run this:
#' }
#'
#' @export

LucasKanadeEuclidean <- function(It, It1) {
    # inverse compositional
	eps <- .0001
	p <- rep(0, 3)

	maxY <- dim(It1)[1]
	maxX <- dim(It1)[2]

	# get (x, y) coordinates going down column
	# X <- rep(1:maxX, each = maxY)
	# Y <- rep(1:maxY, maxX)
	# [X, Y] = meshgrid(1:maxX, 1:maxY);
	mgrid <- pracma::meshgrid(1:maxX, 1:maxY)

	# [Gx, Gy] = imgradientxy(It);
    tmpImager <- imager::as.cimg(It)
    Gy <- imager::get_gradient(tmpImager, axes = "x", scheme = 2)[[1]][, , , ] # scheme 2 for Sobel
    Gx <- imager::get_gradient(tmpImager, axes = "y", scheme = 2)[[1]][, , , ]

	# % multiply gradient with Jacobian
	A <- cbind(-c(mgrid$Y)*c(Gx) + c(mgrid$X)*c(Gy), c(Gx), c(Gy))

	H <- solve(t(A) %*% A)
	HAp <- H %*% t(A)

    for (iter in 1:10000) {
		tmpM <- matrix(c(cos(p[1]), -sin(p[1]), p[2], sin(p[1]), cos(p[1]), p[3], 0, 0, 1), byrow = TRUE, nrow = 3)
		tmp <- tmpM %*% rbind(c(mgrid$X), c(mgrid$Y), rep(1, maxX*maxY))

		Xq <- matrix(tmp[1, ], nrow = nrow(mgrid$X))
		Yq <- matrix(tmp[2, ], nrow = nrow(mgrid$X))
		warpedT1 <- pracma::interp2(1:maxX, 1:maxY, It1, Xq, Yq)
		warpedT1 <- matrix(warpedT1, nrow = nrow(mgrid$X))

		b <- warpedT1 - It

		# % only rows that are not NA
		b <- c(b)
		tmpHAp <- HAp[ , !is.na(b)] # % 3xD
		b <- b[!is.na(b)]

		deltaP <- tmpHAp %*% b # % 3x1
		# % update warp
		Winv <- solve(matrix(c(cos(deltaP[1]), -sin(deltaP[1]), deltaP[2], sin(deltaP[1]), cos(deltaP[1]), deltaP[3], 0, 0, 1), byrow = TRUE, nrow = 3))
		W <- tmpM %*% Winv

		theta <- acos(W[1, 1])
		if (W[2, 1] < 0) {
            theta <- 2*pi - theta
		}

		p <- c(theta, W[1, 3], W[2, 3])

		if ( sum(deltaP^2) < eps ) break
    }

	M <- matrix(c(cos(p[1]), -sin(p[1]), p[2], sin(p[1]), cos(p[1]), p[3], 0, 0, 1), nrow = 3, byrow = TRUE)

    ret <- list(M = M, p = p, iter = iter)
    return(ret)
}

#'Compare two images and return similarities
#'
#'Runs the Lucas-Kanade algorithm to find best alignment parameters and returns
#'correlation and mean-squared error
#'
#'@param I1 image that is transformed
#'@param I2 image to align to
#'
#'@return transformation vector \code{p} (theta, dx, dy), number of iterations
#'  before convergence, final correlation (only breechface area), MSE (only
#'  breechface area), overlapping number of pixels in breechface area
#' @examples
#' \dontrun{
#' # first run this:
#' }
#'
#'@export

compareTwo3D <- function(I1, I2) {

	tmpI1 <- imager::resize(imager::as.cimg(I1), -.1 * 100, -.1 * 100, -100, interpolation_type = 5)[, , , ]
	tmpI2 <- imager::resize(imager::as.cimg(I2), -.1 * 100, -.1 * 100, -100, interpolation_type = 5)[, , , ]

    # % pad if different sizes
	out <- padImages(tmpI1, tmpI2)
    I1 <- out$padded1
	I2 <- out$padded2

	# tStart = tic;
	out <- LucasKanadeEuclidean(I1, I2)
	# timing = toc(tStart);
    M <- out$M
    p <- out$p
    iter <- out$iter

	# % for plot
	maxY <- dim(I1)[1]
	maxX <- dim(I1)[2]
	# [X, Y] = meshgrid(1:maxX, 1:maxY);
	mgrid <- pracma::meshgrid(1:maxX, 1:maxY)

	tmp <- solve(M) %*% rbind(c(mgrid$X), c(mgrid$Y), rep(1, maxX*maxY)) #% M inverse instead of M to do the actual warp
	Xq <- matrix(tmp[1, ], nrow = nrow(mgrid$X))
	Yq <- matrix(tmp[2, ], nrow = nrow(mgrid$X))
	warpedT <- pracma::interp2(1:maxX, 1:maxY, I1, Xq, Yq)
	warpedT <- matrix(warpedT, nrow = nrow(mgrid$X))

	# % for cor, MSE and overlap
	tmp <- c(warpedT) - c(I2)
	# cor <- ccf(warpedT[!is.na(tmp)], I2[!is.na(tmp)], 0, 'correlation', plot = FALSE)
	# overlap <- sum(!is.na(tmp))
	cor <- ccf(warpedT[warpedT != 0 & I2 != 0], I2[warpedT != 0 & I2 != 0], 0, 'correlation', plot = FALSE, na.action = na.pass)
	overlap <- sum(warpedT != 0 & I2 != 0, na.rm = TRUE)
	# MSE <- sum((tmp)^2, na.rm = TRUE) / overlap
	MSE <- sum((tmp[warpedT != 0 & I2 != 0])^2, na.rm = TRUE) / overlap
	# note that 0's are not considered here but were considered in minimization in LucasKanadeEuclidean (i.e. LucasKanadeEuclidean will try to align nonBF areas)

	ret <- list(p = p, iter = iter, cor = cor, MSE = MSE, overlap = overlap)
	return(ret)
}

#'Zero pad pairs of images
#'
#' Given a pair of images, zero-pad smaller one so that they are the same size
#'
#'@param I1 first image
#'@param I2 second image
#'
#'@return list with two images
#' @examples
#' \dontrun{
#' # first run this:
#' }
#'
#'@export

padImages <- function(I1, I2) {

    dim1 <- dim(I1)
    dim2 <- dim(I2)
    maxRow <- max(dim1[1], dim2[1])
	maxCol <- max(dim1[2], dim2[2])
    padded1 <- matrix(0, nrow = maxRow, ncol = maxCol)
    padded2 <- matrix(0, nrow = maxRow, ncol = maxCol)

    row_start <- floor((maxRow - dim1[1]) / 2)
    col_start <- floor((maxCol - dim1[2]) / 2)
    padded1[(row_start + 1):(row_start + dim1[1]), (col_start + 1):(col_start + dim1[2])] <- I1

    row_start <- floor((maxRow - dim2[1]) / 2)
    col_start <- floor((maxCol - dim2[2]) / 2)
    padded2[(row_start + 1):(row_start + dim2[1]), (col_start + 1):(col_start + dim2[2])] <- I2

    ret <- list(padded1 = padded1, padded2 = padded2)
    return(ret)
}

