% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step5_compare3D.R
\name{LucasKanadeEuclidean}
\alias{LucasKanadeEuclidean}
\title{Align two images using the Lucas-Kanade algorithm}
\usage{
LucasKanadeEuclidean(It, It1)
}
\arguments{
\item{It}{image at time t, or image that is transformed}

\item{It1}{image at time t+1, or image to align to}
}
\value{
homography transformation matrix \code{M}, transformation vector
  \code{p} (theta, dx, dy), number of iterations before convergence
  \code{iter}
}
\description{
Finds translation and rotation parameters to approximately minimize the
mean-squared error between two images.
}
\examples{
\dontrun{
# first run this:
}

}
