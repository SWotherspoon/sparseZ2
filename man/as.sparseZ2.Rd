% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparseZ2.R
\name{as.sparseZ2}
\alias{as.sparseZ2}
\title{Convert Matrix to Sparse Format.}
\usage{
as.sparseZ2(mat)
}
\arguments{
\item{mat}{A matrix of 0s and 1s}
}
\value{
A sparse matrix representation
}
\description{
Convert a matrix to a \code{sparseZ2} matrix.
}
\details{
Converts a matrix to a \code{sparseZ2} format, with any non-zero values
treated as 1.
}
