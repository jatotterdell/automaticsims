% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{prob_h_best}
\alias{prob_h_best}
\title{Probability best h according to rank r}
\usage{
prob_h_best(mat, r = order(expected_rank(mat)), delta = 0)
}
\arguments{
\item{mat}{A matrix of samples}

\item{r}{A ranking vector giving the order from smallest to largest}
}
\value{
A vector giving the probability the h ranked means are superior
to the other k-h ranked means
}
\description{
Probability that the h best populations ranked according to r
are the best h.
}
\examples{
# Expect high probability that two best identified
# Others should be 50/50
D <- mvnfast::rmvn(4e4, c(1,1,2,2), diag(0.1,4))
prob_h_best(D)
}
