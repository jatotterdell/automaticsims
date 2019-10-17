#' Run an indifference trial with interims
#'
#' An indifference stopping...
#'
#' @export
run_indifference_trial <- function(
  id,
  mu,
  delta_indiff = 0.25,
  delta_ctrl = delta_indiff,
  kappa_sup = 0.95,
  kappa_ind = 0.75,
  kappa_ctr = 0.95,
  ctrl_alloc = 1/13
) {
  
  # Setup
  mod_draws <- 1e4
  N <- 10000
  P <- 13
  K <- 20
  M <- N / K
  arm_labs <- sprintf("%02d", 0:12)
  mes_labs <- paste0("m", 1:4)
  tim_labs <- paste0("t", 1:3)
  
  p <- setNames(rep(1/P, P), arm_labs)
  p[1] <- ctrl_alloc
  p[-1] <- (1 - p[1])/(P - 1)
  
  n <- setNames(rep(0, P), arm_labs)
  y <- setNames(rep(0, P), arm_labs)
  
  for(k in 1:K) {
    # Generate data
    n_new <- c(rmultinom(1, M, p))
    y_new <- rbinom(P, n_new, plogis(mu)) 
    n <- n + n_new
    y <- y + y_new
    
    # Fit model
    mod <- vb_mod(y, n)
    m <- setNames(mod$mu, arm_labs)
    v <- setNames(diag(mod$Sigma), arm_labs)
    
    # Compute posterior quantities
    draws <- mvnfast::rmvn(mod_draws, m, sigma = mod$Sigma)
    beta_draws <- draws %*% X_con_inv_t_Q_t
    
    p_max <- setNames(prob_max(draws[, -1]), arm_labs[-1])
    e_rank <- setNames(expected_rank(draws[, -1]), arm_labs[-1])
    o_rank <- order(e_rank)
    best <- which.max(e_rank)
    p_h_best <- prob_h_best(draws[, -1], o_rank)
    p_h_idif <- prob_h_indiff(draws[, -1], o_rank, delta = delta_indiff)
    p_beat_ctrl <- setNames(prob_superior(draws[, -1], draws[, 1], delta_ctrl), arm_labs[-1])
    
    w <- sqrt(p_max * v[-1] / n[-1])
    p[-1] <- (1 - p[-1])* w / sum(w)
    
    found_best <- any(p_h_best > kappa_sup & p_h_idif > kappa_ind)
    all_indiff <- any(p_h_idif[1] > kappa_ind)
    if(found_best) {
      which_best <- o_rank[(P-1):max(which(p_h_best > kappa_sup & p_h_idif > kappa_ind))]
      in_best <- rep(FALSE, P - 1)
      in_best[which_best] <- TRUE
    } else {
      which_best <- NA
      in_best <- rep(0, P)
    }
    beat_ctrl <- p_beat_ctrl > kappa_ctr
    lose_ctrl <- p_beat_ctrl < 1 - kappa_ctr
    lose <- all(lose_ctrl)
    best_beat_ctrl <- all(beat_ctrl[best])
    all_best_beat_ctrl <- ifelse(found_best, all(beat_ctrl[in_best]), 0)
    
    if((found_best & best_beat_ctrl) | lose | all_indiff) break
  }
  
  return(list(
    id = id,
    mu = mu,
    delta_indiff = delta_indiff,
    delta_ctrl = delta_ctrl,
    kappa_sup = kappa_sup,
    kappa_ind = kappa_ind,
    kappa_ctr = kappa_ctr,
    ctrl_alloc = ctrl_alloc,
    inteirm = k,
    stopped = as.numeric(k < K),
    found_best = found_best,
    best_beat_ctrl = best_beat_ctrl,
    all_best_beat_ctrl = all_best_beat_ctrl,
    all_indiff = all_indiff,
    lose = lose,
    p = p,
    n = n,
    y = y,
    m = m,
    v = v,
    p_max = p_max,
    e_rank = e_rank,
    best = best,
    in_best = in_best,
    p_h_best = p_h_best,
    p_h_idif = p_h_idif,
    p_beat_ctrl = p_beat_ctrl,
    beat_ctrl = beat_ctrl
  ))
}