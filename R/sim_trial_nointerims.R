#' Run a non-inferiority trial
#' 
#' Run a single trial with no interims.
#' 
#' @param id The trial ID
#' @param mu The true mean, must be length 13
#' @param delta The superiority margin
#' @param kappa_lo The starting threshold to deactivate poor arms
#' @param kappa_hi The starting threshold for superiority
#' @param kappa_no The starting threshold for non-inferiority
#' @param kappa_ctrl The starting threshold for non-inferiority
#' @param ind_comp_ctrl Should deactivation be based on P(max trt arms) or P(max all arms)
#'
#' @return A list of trial quantities
#'
#' @export
run_a_noninf_trial_nointerims <- function(
  id,
  mu,
  delta,
  kappa_lo = 0.01,
  kappa_hi = 0.75,
  kappa_no = 0.5,
  kappa_ctrl = 0.95,
  ctrl_alloc = 1/13,
  ind_comp_ctrl = FALSE
) {
  
  # Setup
  N <- 10000
  P <- 13
  arm_labs <- sprintf("%02d", 0:12)
  mes_labs <- paste0("m", 1:4)
  tim_labs <- paste0("t", 1:3)

  p <- setNames(rep(1/P, P), arm_labs)
  p[1] <- ctrl_alloc
  p[-1] <- (1 - p[1])/(P - 1)

  p_max_mes <- setNames(rep(0, P), 1:4)
  p_max_tim <- setNames(rep(0, P), 1:3)
  p_noninf <- 0
  p_best_beat_inactive <- 0
  p_best_beat_inactive_ctrl <- 0
  best <- 0
  is_sup <- setNames(rep(0, P), arm_labs)
  active <- setNames(rep(0, P - 1), arm_labs[-1])
  
  # Generate data
  n <- setNames(table(factor(sample.int(P, N, replace = TRUE, prob = p), levels = 1:P)), arm_labs)
  y <- setNames(rbinom(P, n, plogis(mu)), arm_labs)

  # Fit model
  mod <- vb_mod(y, n)
  m <- setNames(mod$mu, arm_labs)
  v <- setNames(diag(mod$Sigma), arm_labs)
    
  # Compute posterior quantities
  draws <- mvnfast::rmvn(1e4, m, sigma = mod$Sigma)
  beta_draws <- draws %*% X_con_inv_t_Q_t
  p_sup <- setNames(prob_each_superior_all(draws, delta), arm_labs)
  p_sup_trt <- setNames(prob_each_superior_all(draws[, -1], delta), arm_labs[-1])
  p_max_all <- setNames(prob_max(draws), arm_labs)
  p_max <- setNames(prob_max(draws[, -1]), arm_labs[-1])
  
  p_sup_mes <- prob_each_superior_all(beta_draws[, 3:6], delta)
  p_sup_tim <- prob_each_superior_all(beta_draws[, 7:9], delta)
  p_max_mes <- prob_max(beta_draws[, 3:6])
  p_max_tim <- prob_max(beta_draws[, 7:9])
  p_beat_ctrl <- setNames(prob_superior(draws[, -1], draws[, 1], delta), arm_labs[-1])
  p_rank <- prob_rank(draws)
  e_rank <- expect_rank(p_rank)
  p_rank_cum <- matrixStats::colCumsums(p_rank)
  dimnames(p_rank_cum) <- list("rank" = 1:P, "arm" = 0:(P - 1))
  best <- unname(which.max(p_max))
    
  # Deactivation rule
  if(!ind_comp_ctrl) {
    active <- setNames(as.numeric(p_sup[-1] > kappa_lo), arm_labs[-1])
    is_sup <- setNames(p_sup > kappa_hi, arm_labs)
    superior <- any(is_sup)
  } else {
    active <- setNames(as.numeric(p_sup_trt > kappa_lo & p_beat_ctrl > 1 - kappa_ctrl), arm_labs[-1])
    is_sup <- setNames(c(FALSE, p_sup_trt > kappa_hi & p_beat_ctrl > kappa_ctrl), arm_labs)
    superior <- any(is_sup)
  }
  if(sum(active) > 1) {
    # Probability all active noninferior to superior by delta
    p_noninf <- prob_all_superior(
      draws[, -1][, active & !(1:(P-1) == best), drop = F],
      draws[, -1][, best], 
      -delta)
    # Probability best active superior to all active and control
    p_best_beat_inactive_ctrl <- prob_superior_all(
      draws[, -1][, best], 
      draws[, c(1, which(active == 0) + 1), drop = F],
      delta)
    p_best_beat_inactive <- prob_superior_all(
      draws[, -1][, best], 
      draws[, which(active == 0) + 1, drop = F],
      delta)
  }
  noninferior <- any(p_noninf > kappa_no & p_best_beat_inactive > kappa_hi)
  nonsuperior <- all(!active)
  stopped <- superior | noninferior | nonsuperior
  
  return(
    list(
      id = id,
      mu = mu,
      delta = delta,
      kappa_lo = kappa_lo,
      kappa_hi = kappa_hi,
      kappa_no = kappa_no,
      stopped = stopped,
      superior = superior,
      noninferior = noninferior,
      nonsuperior = nonsuperior,
      p = p,
      n = n,
      y = y,
      m = m,
      v = v,
      p_sup = p_sup,
      p_sup_trt = p_sup_trt,
      p_max_all = p_max_all,
      p_max = p_max,
      p_max_mes = p_max_mes,
      p_max_tim = p_max_tim,
      p_beat_ctrl = p_beat_ctrl,
      p_rank = collapse_prob_rank(p_rank),
      e_rank = e_rank,
      p_rank_cum = collapse_prob_rank(p_rank_cum),
      p_noninf = p_noninf,
      p_best_beat_inactive_ctrl = p_best_beat_inactive_ctrl,
      p_best_beat_inactive = p_best_beat_inactive,
      p_sup_mes = p_sup_mes,
      p_sup_tim = p_sup_tim,
      p_sup_pairwise = pairwise_superiority_all(draws, delta, replace = TRUE),
      p_sup_mes_pairwise = pairwise_superiority_all(beta_draws[, 3:6], delta, replace = TRUE),
      p_sup_tim_pairwise = pairwise_superiority_all(beta_draws[, 7:9], delta, replace = TRUE),
      active = active,
      best = best,
      sup = is_sup,
      beat_ctrl = p_beat_ctrl > kappa_hi
    )
  )
}
