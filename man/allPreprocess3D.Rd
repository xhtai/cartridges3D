% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess3D.R
\name{allPreprocess3D}
\alias{allPreprocess3D}
\title{Perform all pre-processing steps}
\usage{
allPreprocess3D(fileName)
}
\arguments{
\item{fileName}{location of image file.}
}
\value{
A matrix of processed depth values. Non-breechface pixels are set to
  0.
}
\description{
Images must be x3p images containing 3D topographies. This function performs
all the pre-processing steps, consisting of \itemize{ \item{Automatically
selecting the breechface marks} \item{Leveling the image} \item{Removing
circular symmetry} \item{Filtering.} }
}
\examples{
\dontrun{
}

}
