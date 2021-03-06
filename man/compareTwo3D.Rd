% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step5_compare3D.R
\name{compareTwo3D}
\alias{compareTwo3D}
\title{Compare two images and return similarities}
\usage{
compareTwo3D(I1, I2)
}
\arguments{
\item{I1}{image that is transformed}

\item{I2}{image to align to}
}
\value{
transformation vector \code{p} (theta, dx, dy), number of iterations
 before convergence, final correlation (only breechface area), MSE (only
 breechface area), overlapping number of pixels in breechface area
}
\description{
Runs the Lucas-Kanade algorithm to find best alignment parameters and returns
correlation and mean-squared error
}
\examples{
\dontrun{
# first run this:
}

}
