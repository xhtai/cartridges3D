% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step3_removeCircular3D.R
\name{statisticsByBasisFunction}
\alias{statisticsByBasisFunction}
\title{Get statistics for each basis function}
\usage{
statisticsByBasisFunction(ij, image, requestedStats)
}
\arguments{
\item{ij}{ij coordinates of the basis function. Input is typically one item
of the list of basis functions}

\item{image}{image that statistics are to be calculated for}

\item{requestedStats}{vector with each element being a requested statistic,
e.g. \code{c("coef","RSS")}. Possible options are \code{"numPixels"} for
the number of pixels in the input basis function, \code{"coef"} for the
coefficient as described above, \code{"RSS"} for the residual sum of
squares, \code{"sum"} for the sum of pixel values, and \code{"max"} and
\code{"min"} for the maximum and minimum pixel values.}
}
\value{
A list with the requested statistics in the order provided.
}
\description{
Get statistics on pixel values in each basis function: mean of the pixel
values (basis function coefficient), residual sum of squares using fitted
coefficient, sum of pixel values.
}
\examples{
basis3 <- getBasisFunctions(3)
stats <- statisticsByBasisFunction(basis3[[2]],
         matrix(1:9, nrow = 3), c("numPixels", "coef"))
}
