# library(varapproxr)
# library(automaticsims)
# library(matrixStats)
# 
# Sigma0 <- diag(1.5, ncol(X_con))
# mu <- c(rep(1, 7), rep(1.5, 6))
# p <- plogis(mu) / sum(plogis(mu))
# n <- rmultinom(1, 10000, p)[, 1]
# y <- rbinom(13, n, plogis(mu))
# mod <- automaticsims::vb_mod(y, n, Sigma0 = Sigma0)
# m <- mod$mu
# v <- diag(mod$Sigma)
# draws <- mvnfast::rmvn(1e5, m, sigma = mod$Sigma)
# P_draws <- plogis(draws)
# P_maxes <- rowMaxs(P_draws)
# w <- automaticsims:::prob_max(P_draws[, -1])
# 
# # Non-inferiority to the maximum
# active <- which(w > 0.01) + 1
# automaticsims:::prob_all_superior(P_draws[, active], P_maxes, -0.025)
# automaticsims:::prob_all_superior(P_draws[, active], P_draws[, 1], 0)
# colMeans2(sweep(P_draws[, active], 1, rowMaxs(P_draws) - 0.025, ">="))
# colMeans2(sweep(P_draws[, active], 1, P_draws[, 1], ">"))
# 
# P_maxs_a <- rowMaxs(P_draws[, active])
# P_mins_a <- rowMins(P_draws[, active])
# P_maxs_i <- rowMaxs(P_draws[, !active][, -1])
# P_mins_i <- rowMins(P_draws[, !active][, -1])
# mean(P_mins_a - P_maxs_a > -.025)
# mean(P_mins_a - P_draws[, 1] > 0)
# mean(P_mins_a - P_maxs_i > 0)
# automaticsims:::prob_all_equivalent(P_draws[, active], 0.025)
# 
# plot(density(P_maxs_a), xlim = c(0.5, 1))
# lines(density(P_mins_a))
# 
# # Density of the max(mu_1,...,mu_13)
# x <- seq(-5, 5, length.out = 10000)
# r <- sapply(1:13, function(b) sapply(x, function(a) dnorm(a, m[b], sqrt(v[b]))*prod(pnorm(a, m[-b], sqrt(v[-b])))))
# hist(rowMaxs(qlogis(P_draws)), freq = F, breaks = 50, xlim = c(0, 4), col = "grey50", border = "grey50")
# lines(x, rowSums(r), col = "blue", lwd = 2)
# for(a in active) {
#   curve(dnorm(x, m[a], sqrt(v[a])), col = "red", add = T, n = 1000)
# }
# 
# 
# P_noninf <- sweep(P_draws[, active], 1, P_maxes, "-")
# gather(as.data.frame(P_noninf)) %>% ggplot(., aes(value)) + facet_wrap( ~ key) + stat_ecdf() + geom_vline(xintercept = -0.025)
# 
# automaticsims:::prob_rank(P_draws)
# 
# plot(density(P_max))
# lines(density(P_draws[, 13]), col = 'red')
# 
# P_max <- sapply(1:nrow(P_draws), function(i) P_draws[i, sample.int(13, 1, replace = T, w)])
# colMeans2(sweep(P_draws, 1, P_max - 0.025, ">="))


#' @export
gen_potential_outcomes <- function(n, p) {
  k <- length(p)
  matrix(rbinom(n*k, 1, p), n, k, byrow = T,
         dimnames = list('id' = 1:n, 'o' = paste0('y_', 1:k)))
}

#' @export
gen_allocations <- function(n, a) {
  k <- length(a)
  sample.int(k, n, prob = a, replace = T)
}

#' @export
run_trial <- function(
  tdm,
  n_seq,
  use_cum_p_max = TRUE,
  brar_h = function(i) 0
) {

  n_cohort <- diff(c(0, n_seq))
  n_interim <- length(n_seq)
  n_arms <- ncol(tdm)
  n_max <- nrow(tdm)
  
  dnames1 <- list("interim" = 0:n_interim, "arm" = seq_len(n_arms))
  dnames2 <- list("interim" = 1:n_interim, "arm" = seq_len(n_arms))
  dnames3 <- list("interim" = 1:n_interim, "arm" = seq_len(n_arms)[-1])
  
  # Storage
  alloc_prob        <- matrix(rep(1 / n_arms, n_arms), n_interim + 1, n_arms, dimnames = dnames1)
  p_max_all         <- matrix(0, n_interim, n_arms, dimnames = dnames2)
  p_max             <- matrix(0, n_interim, n_arms - 1, dimnames = dnames3)
  p_beat_ctrl       <- matrix(0, n_interim, n_arms - 1, dimnames = dnames3)
  p_worst_beat_ctrl <- matrix(0, n_interim, 1, dimnames = list("interim" = 1:n_interim, "p_worst_beat_ctrl"))
  p_worst_beat_drop <- matrix(0, n_interim, 1, dimnames = list("interim" = 1:n_interim, "p_worst_beat_drop"))
  p_worst_ninf_best <- matrix(0, n_interim, 1, dimnames = list("interim" = 1:n_interim, "p_worst_ninf_best"))
  
  active             <- matrix(1, n_interim, n_arms - 1, dimnames = dnames3)
  arm_above_cumthres <- matrix(1, n_interim, n_arms - 1, dimnames = dnames3)
  rnk                <- matrix(0, n_interim, n_arms - 1, dimnames = dnames3)
  best               <- matrix(0, n_interim, 1, dimnames = list("interim" = 1:n_interim, "best"))
  
  mu_mean <- matrix(0, n_interim, n_arms, dimnames = dnames2)
  mu_var <- matrix(0, n_interim, n_arms, dimnames = dnames2)
  mu_hdi <- array(0, c(n_interim, n_arms, 2), dimnames = c(dnames2, list("bound" = c("lo", "hi"))))
  
  y_cum <- matrix(0, n_interim, n_arms, dimnames = dnames2)
  n_cum <- matrix(0, n_interim, n_arms, dimnames = dnames2)
  alloc <- factor(vector("integer", n_max), levels = seq_len(n_arms))
  dat   <- tibble::tibble(id = seq_len(n_max), a = alloc, y = NA)
  
  for(i in seq_len(n_interim)) {
    alloc_id <- (c(0, n_seq)[i]+ 1):n_seq[i] # who to allocate
    anlys_id <- seq_len(n_seq[i])
    alloc[alloc_id] <- gen_allocations(n_cohort[i], alloc_prob[i, ]) # what do they get
    dat[alloc_id, -1] <- cbind(alloc[alloc_id], tdm[cbind(seq_along(alloc[alloc_id]), alloc[alloc_id])]) # fill in the data
    dat_agg <- dat[alloc_id, ] %>% group_by(a, .drop=F) %>% summarise(y = sum(y), n = n())
    if(i == 1) {
      n_cum[i, ] <- dat_agg[["n"]]
      y_cum[i, ] <- dat_agg[["y"]]
    } else {
      n_cum[i, ] <- n_cum[i - 1, ] + dat_agg[["n"]]
      y_cum[i, ] <- y_cum[i - 1, ] + dat_agg[["y"]]
    }
    
    # Model and posterior summaries
    mod <- vb_mod(y_cum[i, ], n_cum[i, ])
    mu_mean[i, ] <- mod$mu
    mu_var[i, ] <- diag(mod$Sigma)
    mu_hdi[i, , ] <- mu_mean[i, ] + cbind(-sqrt(mu_var[i, ]), sqrt(mu_var[i, ]))
    
    # Posterior samples
    draws <- mvnfast::rmvn(1e4, mod$mu, sigma = mod$Sigma)
    pi_draws <- plogis(draws)
    beta_draws <- draws %*% X_con_inv_t_Q_t
    max_draws <- matrixStats::rowMaxs(draws)
    
    p_max_all[i, ] <- prob_max(draws)
    p_max[i, ] <- prob_max(draws[, -1])
    p_beat_ctrl[i, ] <- prob_superior(draws[, -1], draws[, 1], 0)
    rnk[i, ] <- rank(1 - p_max[i, ], ties = "random")
    best[i, ] <- which(rnk[i, ] == 1) + 1
    arm_above_cumthres[i, as.integer(names(which(cumsum(p_max[i, ][order(rnk[i, ])]) > 0.95)[-1])) - 1] <- 0
    
    # Identify active arms and update allocation probabilities according to rules
    if(use_cum_p_max) {
      active[i, ] <- arm_above_cumthres[i, ]
    } else {
      active[i, ] <- as.numeric(p_max[i, ] > kappa_act[i] & p_beat_ctrl[i, ] > 1 - kappa_ctr[i])   
    }
    w <- (p_max[i, ] / (n_cum[i, -1] + 1))^brar_h(i)
    w[!active[i, ]] <- 0
    alloc_prob[i + 1, -1] <- (1 - alloc_prob[i, 1]) * w / sum(w)   
    
    p_worst_beat_ctrl[i, ] <- mean(matrixStats::rowMins(draws[, which(arm_above_cumthres[i, ] == 1) + 1, drop = F]) - draws[, 1] > 0)
    p_worst_beat_drop[i, ] <- mean(matrixStats::rowMins(draws[, which(arm_above_cumthres[i, ] == 1) + 1, drop = F]) - 
                                     matrixStats::rowMaxs(draws[, which(arm_above_cumthres[i, ] == 0) + 1, drop = F]) > 0)
    p_worst_ninf_best[i, ] <- prob_all_superior(
      pi_draws[, setdiff(which(active[i, ] == 1) + 1, best[i, ]), drop = F],
      pi_draws[, best[i, ]], 
      -0.025)
    
  }
  
  alloc_prob <- alloc_prob[-1, ]
  
  return(
    nlist(
      n_cum, y_cum, mu_mean, mu_var, mu_hdi,
      p_max, p_max_all, p_beat_ctrl, p_worst_beat_ctrl, p_worst_beat_drop, p_worst_ninf_best,
      alloc_prob, active, rnk, best
    )
  )
}
