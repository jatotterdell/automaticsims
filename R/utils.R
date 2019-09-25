#' Probability that each column in mat is the maximum of all columns
#' 
#' @param mat A matrix of samples
#' @return A vector of probabilities that each dimension is maximum
#' @examples 
#' prob_max(matrix(rnorm(1e4, 1:10), 1000, 10, byrow = T))
prob_max <- function(mat) {
  as.numeric(prop.table(table(factor(max.col(mat), levels = 1:ncol(mat)))))
}

#' Rank rows of a matrix.
#' 
#' Rank the rows of a matrix by decreasing value. So, rank 1 is maximum, rank K is minimum.
#' 
#' @param mat A matrix of samples
#' @return A matrix of rank probabilities. Each row is a column in mat, and each column reflects a rank.
#' @examples 
#' prob_rank(mvnfast::rmvn(1e4, 1:3, diag(1, 3)))
prob_rank <- function(mat) {
  k <- ncol(mat)
  apply(matrixStats::rowRanks(-mat), 2, matrixStats::binCounts, bx = 1:(k+1), right = F) / nrow(mat)
}

#' Probability that items in a are inferior to b
#' 
#' @param a A matrix or vector of items to compare to b
#' @param b A vector or matrix to which a is compared
#' @param delta The minimal effect bound
#' @return A vector of probabilities
#' @examples 
#' prob_inferior(mvnfast::rmvn(1e4, 1:3, diag(1, 3)), rnorm(1e4), 1)
#' prob_inferior(rnorm(1e4), mvnfast::rmvn(1e4, 1:3, diag(1, 3)), 1)
prob_inferior <- function(a, b, delta = 0) {
  return(colMeans(a < b + delta))
}

#' Probability that items in a are superior to b
#' 
#' @param a A matrix or vector of items to compare to b
#' @param b A vector or matrix to which a is compared
#' @param delta The minimal effect bound
#' @return A vector of probabilities
#' @examples 
#' prob_superior(mvnfast::rmvn(1e4, 1:3, diag(1, 3)), rnorm(1e4), 1)
#' prob_superior(rnorm(1e4), mvnfast::rmvn(1e4, 1:3, diag(1, 3)), 1)
prob_superior <- function(a, b, delta = 0) {
  return(colMeans(a > b + delta))
}

#' Probability superior
#' 
#' Probability that a is superior to all items in b
#' 
#' @param a A vector
#' @param b A matrix
#' @param delta The minimal effect bound
#' @return A numeric probability
prob_superior_all <- function (a, b, delta = 0) {
  return(mean(a > matrixStats::rowMaxs(b) + delta))
}


#' Probability superior
#' 
#' Probability that all in a are superior b
#' 
#' @param a A matrix
#' @param b A vector
#' @param delta The minimal effect bound
#' @return A numeric probability
#' @examples 
#' # All superior
#' prob_all_superior(matrix(rnorm(1e4, 1:10), 1000, 10, byrow = T), rnorm(1e4, 0), 0.5)
#' # All non-inferior
#' prob_all_superior(matrix(rnorm(1e4, 1:10), 1000, 10, byrow = T), rnorm(1e4, 0), -0.5)
prob_all_superior <- function (a, b, delta = 0) {
  return(mean(matrixStats::rowMins(a) > b + delta))
}


#' Probability each superior to all
#' 
#' Probability that each column in mat is superior to all other columns.
#' 
#' @param mat A matrix of samples
#' @param delta The minimal effect bound
#' @return A vector of probabilities that each dimension is maximum
#' @examples 
#' # Superiority
#' prob_each_superior_all(matrix(rnorm(1e4, 1:10), 1000, 10, byrow = T), 0.5)
#' # Non-inferiority
#' prob_each_superior_all(matrix(rnorm(1e4, 1:10), 1000, 10, byrow = T), -0.5)
prob_each_superior_all <- function(mat, delta = 0) {
  return(
    sapply(seq_len(ncol(mat)), 
           function(x) 
             prob_superior_all(mat[, x], mat[, -x], delta))
  ) 
}

#' Pairwise difference matrix
#' 
#' A function to generate pairwise difference matrix for unique combinations.
#' The function only generates the lower triangular differences, that is K(K-1)/2
#' where K is the number of columns.
#' 
#' @param mat A matrix of samples
#' @return A matrix of K(K-1)/2 pairwise differences
pairwise_diff <- function (mat) {
  pair_comp <- arrangements::combinations(ncol(mat), 2)
  pair_mat <- apply(pair_comp, 1, function(x) mat[, x[1]] - mat[, x[2]])
  colnames(pair_mat) <- apply(pair_comp - 1, 1, paste, collapse = "-")
  return(pair_mat)
}

#' Pairwise difference matrix
#' 
#' A function to generate pairwise difference matrix for all combinations.
#' The function generates all K(K-1) differences (excludes the diagonal)
#' where K is the number of columns.
#' 
#' @param mat A matrix of samples
#' @return A matrix of K(K-1) pairwise differences in both directions.
pairwise_diff_all <- function(mat, ...) {
  pair_comp <- arrangements::permutations(ncol(mat), 2, ...)
  pair_mat <- apply(pair_comp, 1, function(x) mat[, x[1]] - mat[, x[2]])
  colnames(pair_mat) <- apply(pair_comp - 1, 1, paste, collapse = "-")
  return(pair_mat)
}


#' Pairwise superiority
#' 
#' @param mat A matrix of samples
#' @return A vector of K(K-1)/2 pairwise probabiliy of superiority
pairwise_superiority <- function(mat) {
  pmat <- pairwise_diff(mat)
  apply(pmat, 2, function(x) mean(x > 0))
}


#' Pairwise superiority for both directions
#' 
#' @param mat A matrix of samples
#' @param delta The minimal effect of interest
#' @return A vector of K(K-1)/2 pairwise probabiliy of superiority
pairwise_superiority_all <- function(mat, delta = 0, ...) {
  pmat <- pairwise_diff_all(mat, ...)
  apply(pmat, 2, function(x) mean(x > delta))
}


#' Threshold sequence
#' 
#' @param a0 The starting threshold
#' @param a1 The ending threshold
#' @param r The change transformation
#' @param t The sequence length from start to end
#' @return A vector of length t + 1 giving the threshold sequence
thres_seq <- function(a0, a1, r, t) {
  a0 + (a1 - a0) * (0:t / t)^r
}


#' Fit VB model
#' 
#' @export
vb_mod <- function(y, n, ...) {
  mod <- varapproxr::vb_logistic_n(
    X_con, y, n, 
    mu0 = rep(0, ncol(X_con)), Sigma0 = diag(10, ncol(X_con)),
    mu_init = rep(0, ncol(X_con)), Sigma_init = diag(1, ncol(X_con)),
    alg = "sj", maxiter_jj = 100)
  return(list(
    mu = drop(X_con %*% mod$mu),
    Sigma = X_con %*% mod$Sigma %*% t(X_con),
    beta_mu = drop(Q %*% mod$mu),
    beta_Sigma  = Q %*% mod$Sigma %*% t(Q)
  ))
}
