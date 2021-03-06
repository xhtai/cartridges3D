% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{plotImage}
\alias{plotImage}
\title{Plot an image}
\usage{
plotImage(image, type = "original", grayscale = FALSE, main = "")
}
\arguments{
\item{image}{matrix of pixel values to be plotted}

\item{type}{either \code{"original"} for images on the original scale (pixel
values 0-255), \code{"residuals"} for residual pixel values (-255 to 255),
or \code{"any"} for plotting using the range available (i.e. the minimum
value is plotted as black and the maximum is white, with a linear scale for
intermediate values). \code{"original"} is the default.}

\item{grayscale}{logical value indicating whether or not the grayscale is to
be plotted. \code{FALSE} is the default.}

\item{main}{title for plots with grayscale. The default is an empty string.}
}
\description{
Plot an image
}
\details{
This function plots an image in the same orientation as the matrix
of pixel values, i.e. the top-left pixel of the image is the first entry of
the matrix, and the bottom-right pixel is the last entry.
}
\examples{
\dontrun{
plotImage(processedExample, type = "any")
}
}
