# ACKNOWLEDGEMENTS: all the code in this file was adapted from MATLAB code
# provided by Joseph Roth of Michigan State University. The MATLAB code was
# worked on by Joseph Roth and Zach Richardson, as part of their Breech Face
# Ballistics work.

#' Filter image.
#'
#' This function applies a low and high-pass Gaussian filter. The low-pass
#' filter uses an 11 x 11 matrix with \code{sigma = 6}, and the high-pass filter
#' uses a 149 x 149 matrix with \code{sigma = 44}. This is translated from code
#' by Joseph Roth and Zach Richardson.
#'
#' @param inputImage image matrix
#' @param nonBF binary matrix of the same dimensions as \code{inputImage}, where
#'   1s indicate the non-breechface region.
#'
#' @return The image matrix after filtering.
#'
#' @examples
#' \dontrun{
#' # first run this:
#' croppedExample <- cropBorders(removedExample, centeredExample$centeredPrimer)
#' # then
#' nonBF <- is.na(croppedExample)
#' processedExample <- gaussianFilter(inpaintedExample, nonBF)
#' }
#'
#' @export

gaussianFilter3D <- function(inputImage, nonBF) {
    short_F <- EBImage::makeBrush(size = 3*6 - 1, shape = "Gaussian", sigma = 3)
    long_F <- -EBImage::makeBrush(size = 24*6 - 1, shape = "Gaussian", sigma = 24)
    # short_F <- EBImage::makeBrush(size = 11, shape = "Gaussian", sigma = 6)
    # long_F <- -EBImage::makeBrush(size = 149, shape = "Gaussian", sigma = 44)
    long_F[ceiling(nrow(long_F)/2), ceiling(ncol(long_F)/2)] <- long_F[ceiling(nrow(long_F)/2), ceiling(ncol(long_F)/2)] + 1

    # % Amplifies information at a certain wavelength
    outputImage <- cartridges3D:::filterViaFFT(inputImage, long_F)
    # % handle boundary of circle case
    tmp_ <- cartridges3D:::filterViaFFT(nonBF, long_F)
    outputImage <- outputImage + tmp_ * inputImage
    outputImage[nonBF] <- 0

    # % Smooths image to remove noise
    # % TODO: Can try smoothing before amplifying.
    outputImage <- cartridges3D:::filterViaFFT(outputImage, short_F)
    outputImage[nonBF] <- 0

    return(outputImage)
}


filterViaFFT <- function(A, B) {
    # size of full filter
    m <- dim(A)
    n <- dim(B)
    x <- m + n - 1

    # pad images with 0 so that we do not have circular issues with FFT
    padA <- matrix(0, nrow = x[1], ncol = x[2])
    padB <- matrix(0, nrow = x[1], ncol = x[2])
    padA[1:m[1], 1:m[2]] <- A
    padB[1:n[1], 1:n[2]] <- B

    # Filter in frequency domain
    C <- circshift(fftshift( fft(fft(padA)*Conj(fft(padB)), inverse = TRUE)/(prod(x)) ), round2((n - m)/2, 0))

    half_n <- round2(n/2, 0)
    C <- C[half_n[1]:(half_n[1] + m[1] - 1), half_n[2]:(half_n[2] + m[2] - 1)]
    if (all.equal(c(Im(C)), rep(0, prod(dim(C)))) == FALSE) {
        stop("Non-zero imaginary part")
    }
    return(Re(C))
}

# to copy behavior of round() in matlab
# http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
round2 = function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
}

# http://stackoverflow.com/questions/38230794/how-to-write-fftshift-and-ifftshift-in-r
fftshift <- function(input_matrix) {

    rows <- dim(input_matrix)[1]
    cols <- dim(input_matrix)[2]

    swap_up_down <- function(input_matrix) {
        rows_half <- ceiling(rows/2)
        return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
    }

    swap_left_right <- function(input_matrix) {
        cols_half <- ceiling(cols/2)
        return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
    }

    input_matrix <- swap_up_down(input_matrix)
    return(swap_left_right(input_matrix))

}

# http://stackoverflow.com/questions/18791212/equivalent-to-numpy-roll-in-r
circshift <- function(x, vec) {
    dimx <- dim(x)
    # row first
    if (vec[1] != 0) {
        #out <- rbind(x[(dimx[1] - vec[1] + 1):dimx[1], ], x[1:(dimx[1] - vec[1]), ])
        tmp <- c(t(x))
        x <- matrix(c( tail(tmp, vec[1]*dimx[2]) , head(tmp, -vec[1]*dimx[2]) ), byrow = TRUE, nrow = dimx[1])
    }
    # col
    if (vec[2] != 0) {
        tmp <- c(x)
        x <- matrix(c( tail(tmp, vec[2]*dimx[1]) , head(tmp, -vec[2]*dimx[1]) ), nrow = dimx[1])
    }
    return(x)
}
