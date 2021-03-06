#' Mass-weighted urn randomisation
#'
#' @param target_alloc The target allocation ratios
#' @param sample_size The number of allocations to generate
#' @param alpha Parameter to control imbalance between arms
#' @return A list detailing the mass-weighted-urn process.
#' @export
mass_weighted_urn_design <- function(
  target_alloc,
  sample_size,
  alpha = 4
) {
  arms <- length(target_alloc)
  prob_alloc <- target_alloc / sum(target_alloc)
  # Masses
  x <- matrix(0, sample_size + 1, arms)
  x[1, ] <- alpha * prob_alloc
  # Sample size
  n <- matrix(0, sample_size + 1, arms)
  # Random number
  y <- runif(sample_size)
  # Conditional selection probability
  p <- matrix(0, sample_size + 1, arms)
  # Imbalance
  d <- rep(0, sample_size)
  # Allocation Predictability
  g <- rep(0, sample_size + 1)
  # Treatment assignment
  trt <- rep(0, sample_size)
  
  imbalance_cap <- sqrt(sum(((alpha - 1)*(1 - prob_alloc) + (arms - 1))^2))
  
  for(i in 2:(sample_size + 1)) {
    # Update allocation probabilities
    p[i - 1, ] <- pmax(alpha * prob_alloc - n[i - 1, ] + (i - 1)*prob_alloc, 0)
    p[i - 1, ] <- p[i - 1, ] / sum(p[i - 1, ])
    trt[i-1] <- findInterval(y[i - 1], c(0, cumsum(p[i - 1, ])))
    # Update sample sizes
    n[i, ] <- n[i - 1, ]
    n[i, trt[i-1]] <- n[i, trt[i-1]] + 1
    # Update urn masses
    x[i, trt[i-1]] <- x[i - 1, trt[i-1]] - 1 + prob_alloc[trt[i-1]]
    x[i, -trt[i-1]] <- x[i - 1, -trt[i-1]] + prob_alloc[-trt[i-1]]
    # Calculate imbalance
    d[i - 1] <- sqrt(sum((n[i, ] - (i - 1)*prob_alloc)^2))
    # Calculate allocation predictability
    g[i] <- d[i - 1] / alpha
  }
  return(list(
    max_imbalance_bound = imbalance_cap,
    imbalance = d,
    alloc_predict = g,
    rand_num = y,
    trt = trt,
    mass = x,
    sample_size = n,
    selection_prob = p))
}


#' Find first element of vector satisfying condition
#' 
#' @param x Logical vector
#' @param v Value of NA
#' @return First element satisfying condition
findfirst <- function(x, v = NA) {
  j <- which(x)
  if(length(j)) min(j) else v
}

#' Calculate the differential entropy of multivariate normal
#' 
#' @param S A square matrix
#' @return The log-determinant of S
mvn_entropy <- function(S) {
  res <- 0.5*(ncol(S)*(1 + log(2*pi)) + determinant(S)[[1]])
  attributes(res) <- NULL
  return(res)
}

#' Probability that each column in mat is the maximum of all columns
#' 
#' @param mat A matrix of samples
#' @return A vector of probabilities that each dimension is maximum
#' @examples 
#' prob_max(matrix(rnorm(1e4, 1:10), 1000, 10, byrow = T))
#' @export
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
#' @export
prob_rank <- function(mat) {
  k <- ncol(mat)
  out <- apply(matrixStats::rowRanks(-mat), 2, matrixStats::binCounts, bx = 1:(k+1), right = F) / nrow(mat)
  dimnames(out) <- list("rank" = 1:k, "arm" = 0:(k - 1))
  return(out)
}


collapse_prob_rank <- function(mat) {
  pair_rank <- arrangements::permutations(ncol(mat), 2, replace = T)
  pair_rank[, 1] <- pair_rank[, 1] - 1
  rank_vec <- c(mat)
  names(rank_vec) <- apply(pair_rank, 1, paste, collapse = "-")
  return(rank_vec)
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
#' @export
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
#' @export
prob_superior <- function(a, b, delta = 0) {
  return(matrixStats::colMeans2(a > b + delta))
}

#' Probability superior to all
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


#' Probability all superior
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
#' @export
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
#' @param delta Reference value of interest (negative for noninferiority)
#' @return A vector of K(K-1)/2 pairwise probabiliy of superiority
#' @export
pairwise_superiority <- function(mat, delta = 0) {
  pmat <- pairwise_diff(mat)
  apply(pmat, 2, function(x) mean(x > delta))
}


#' Pairwise superiority for both directions
#' 
#' Pairwise superiority for all K*(K-1) comparisons rather than
#' the reduced K(K-1)/2 one-way comparisons.
#' 
#' @param mat A matrix of samples
#' @param delta Reference value of interest (negative for noninferiority)
#' @return A vector of K(K-1) pairwise probabiliy of superiority
#' @export
pairwise_superiority_all <- function(mat, delta = 0, ...) {
  pmat <- pairwise_diff_all(mat, ...)
  apply(pmat, 2, function(x) mean(x > delta))
}


#' Probability equivalent
#' 
#' Probability all arms in mat are equivalent within delta
#' @param mat A matrix of samples
#' @param delta Equivalence bound
#' @export
prob_all_equivalent <- function(mat, delta) {
  a <- rep(1/ncol(mat), ncol(mat))
  mean(apply(sweep(mat, 1, drop(mat %*% a), `-`), 1,
             function(x) all(abs(x) <= delta)))
}


#' Expected rank
#' 
#' Calculate expected ranks according to posterior draws.
#' 
#' @param mat A matrix of samples
#' @return A vector of size ncol(mat) giving expected rank for each column
#' @export
expected_rank <- function(mat) {
  matrixStats::colMeans2(matrixStats::rowRanks(mat))
}


#' Probability best h according to rank r
#' 
#' Probability that the h best populations ranked according to r
#' are the best h.
#' 
#' @param mat A matrix of samples
#' @param r A ranking vector giving the order from smallest to largest
#' @return A vector giving the probability the h ranked means are superior
#' to the other k-h ranked means
#' @examples
#' # Expect high probability that two best identified
#' # Others should be 50/50
#' D <- mvnfast::rmvn(4e4, c(1,1,2,2), diag(0.1,4))
#' prob_h_best(D)
#' @export
prob_h_best <- function(mat, r = order(expected_rank(mat)), delta = 0) {
  k <- ncol(mat)
  r_mat <- mat[, r]
  c_r_max <- matrixStats::rowCummaxs(r_mat)
  c_r_min <- matrixStats::rowCummins(r_mat[, k:1])[, k:1]
  setNames(c(1, matrixStats::colMeans2(c_r_min[, -1] - c_r_max[, -k] > delta)), k:1)
}

#' Probability that the h best populations ranked according to ord
#' are the best h.
#' 
#' @param mat A matrix of samples
#' @param r A ranking vector giving the order from smallest to largest
#' @param delta An indifference zone value
#' @return A vector giving the probability of indifference (equivalence) 
#' between the h ranked means
#' @examples
#' # Expect high probability that two best identified
#' # Others should be 50/50
#' D <- mvnfast::rmvn(4e4, c(1,1,2,2), diag(0.1,4))
#' rbind(prob_h_best(D), prob_h_indiff(D, delta = 1))
#' @export
prob_h_indiff <- function(mat, r = order(expected_rank(mat)), delta) {
  k <- ncol(mat)
  r_mat <- mat[, r][, k:1]
  setNames(rev(matrixStats::colMeans2(
    matrixStats::rowCummaxs(r_mat) - matrixStats::rowCummins(r_mat) < delta)), k:1)
}


#' Threshold sequence
#' 
#' @export
#' @param a0 The starting threshold
#' @param a1 The ending threshold
#' @param r The change transformation
#' @param t The sequence length from start to end
#' @return A vector of length t + 1 giving the threshold sequence
thres_seq <- function(a0, a1, r, t) {
  if(t == 0) return(a1)
  return(a0 + (a1 - a0) * (0:t / t)^r)
}


#' Fit VB model
#' 
#' @param y Respones
#' @param n Sample size
#' 
#' @export
vb_mod <- function(y, n, ...) {
  mod <- varapproxr::vb_logistic_n(
    X_con, y, n, 
    mu0 = rep(0, ncol(X_con)),
    mu_init = rep(0, ncol(X_con)), Sigma_init = diag(1, ncol(X_con)),
    alg = "sj", maxiter_jj = 100, ...)
  return(list(
    mod_mu = mod$mu,
    mod_Sigma = mod$Sigma,
    mu = drop(X_con %*% mod$mu),
    Sigma = X_con %*% mod$Sigma %*% t(X_con),
    beta_mu = drop(Q %*% mod$mu),
    beta_Sigma  = Q %*% mod$Sigma %*% t(Q)
  ))
}


#' Fit VB model independent
#' 
#' @param y Respones
#' @param n Sample size
#' 
#' @export
vb_mod_ind <- function(y, n, ...) {
  X <- cbind(1, rbind(0, diag(1, 12)))
  mod <- varapproxr::vb_logistic_n(
    X, y, n, 
    mu0 = rep(0, ncol(X)),
    mu_init = rep(0, ncol(X)), Sigma_init = diag(1, ncol(X)),
    alg = "sj", maxiter_jj = 100, ...)
  return(list(
    beta_mu = drop(mod$mu),
    beta_Sigma = mod$Sigma,
    mu = drop(X %*% mod$mu),
    Sigma = X %*% mod$Sigma %*% t(X)
  ))
}


#' Fit VB model with no control
#' 
#' @export
vb_mod_trt <- function(y, n, ...) {
  X_red <- X_con[-1, ][, -1]
  Q_red <- Q[-1, ][, -1]
  mod <- varapproxr::vb_logistic_n(
    X_red, y, n, 
    mu0 = rep(0, ncol(X_red)), Sigma0 = diag(10, ncol(X_red)),
    mu_init = rep(0, ncol(X_red)), Sigma_init = diag(1, ncol(X_red)),
    alg = "sj", maxiter_jj = 100)
  return(list(
    mod_mu = mod$mu,
    mod_Sigma = mod$Sigma,
    mu = drop(X_red %*% mod$mu),
    Sigma = X_red %*% mod$Sigma %*% t(X_red),
    beta_mu = drop(Q %*% mod$mu),
    beta_Sigma  = Q %*% mod$Sigma %*% t(Q)
  ))
}

#' Named list
#' 
#' @export
nlist <- function (...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) 
    FALSE
  else nzchar(names(out))
  if (all(has_name)) 
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  }
  else {
    names(out)[!has_name] <- nms[!has_name]
  }
  return(out)
}
