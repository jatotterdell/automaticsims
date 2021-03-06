% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{prob_rank}
\alias{prob_rank}
\title{Rank rows of a matrix.}
\usage{
prob_rank(mat)
}
\arguments{
\item{mat}{A matrix of samples}
}
\value{
A matrix of rank probabilities. Each row is a column in mat, and each column reflects a rank.
}
\description{
Rank the rows of a matrix by decreasing value. So, rank 1 is maximum, rank K is minimum.
}
\examples{
prob_rank(mvnfast::rmvn(1e4, 1:3, diag(1, 3)))
}
