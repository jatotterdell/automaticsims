#' Run a non-inferiority trial
#' 
#' Run a single trial using non-inferiority stopping.
#' Stopping for noninferiority requires that:
#' 
#' @param id The trial ID
#' @param mu The true mean, must be length 13
#' @param delta The superiority margin
#' @param kappa_lo_0 The starting threshold to deactivate poor arms
#' @param kappa_lo_1 The final threshold to deactivate poor arms
#' @param kappa_hi_0 The starting threshold for superiority
#' @param kappa_hi_1 The final threshold for superiority
#' @param kappa_no_0 The starting threshold for non-inferiority
#' @param kappa_no_1 The final threshold for non-inferiority
#' @param brar Use Response Adaptive Randomisation?
#' @param allocate_inactive Continue to allocate subjects to arms meeting the inactive threshold?
#' @param return_all Return all value, or only those from final analysis?
#' @param ind_comp_ctrl Should deactivation be based on P(max trt arms) or P(max all arms)
#'
#' @return A list of trial quantities
#'
#' @export
run_a_noninf_trial <- function(
  id,
  mu,
  delta_sup = 0.1,
  delta_ctr = delta_sup,
  delta_noninf = delta_sup,
  delta_pair = delta_sup,
  kappa_act_0 = 0.01,
  kappa_act_1 = 0.05,
  kappa_sup_0 = 0.95,
  kappa_sup_1 = 0.80,
  kappa_ctr_0 = 0.95,
  kappa_ctr_1 = 0.95,
  kappa_noninf_0 = 0.60,
  kappa_noninf_1 = 0.60,
  kappa_nonsup_0 = 0.05,
  kappa_nonsup_1 = 0.05,
  brar = FALSE,
  allocate_inactive = FALSE,
  return_all = FALSE,
  ind_comp_ctrl = FALSE,
  ctrl_alloc = 1/13
) {
  
  # Setup
  N <- 10000
  K <- 20
  M <- 500
  P <- 13
  arm_labs <- sprintf("%02d", 0:12)
  mes_labs <- paste0("m", 1:4)
  tim_labs <- paste0("t", 1:3)
  
  kappa_act <- thres_seq(kappa_act_0, kappa_act_1, 1/2, K - 1)
  kappa_sup <- thres_seq(kappa_sup_0, kappa_sup_1, 1/2, K - 1)
  kappa_ctr <- thres_seq(kappa_ctr_0, kappa_ctr_1, 1/2, K - 1)
  kappa_noninf <- thres_seq(kappa_noninf_0, kappa_noninf_1, 1/2, K - 1)
  kappa_nonsup <- thres_seq(kappa_nonsup_0, kappa_nonsup_1, 1/2, K - 1)
  
  p <- matrix(1/P, K + 1, P, dimnames = list("interim" = 0:K, "arm" = arm_labs))
  p[, 1] <- ctrl_alloc
  p[, -1] <- (1 - p[, 1])/12
  n <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  y <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  m <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  v <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  
  p_sup <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  p_sup_trt <- matrix(0, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  p_max_all <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  p_max <- matrix(0, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  p_max_mes <- matrix(0, K, 4, dimnames = list("interim" = 1:K, "mes" = 1:4))
  p_max_tim <- matrix(0, K, 3, dimnames = list("interim" = 1:K, "tim" = 1:3))
  p_beat_ctrl <- matrix(0, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  p_noninf <- rep(0, K)
  p_best_beat_inactive <- rep(0, K)
  
  best <- rep(0, K)
  is_sup <- setNames(rep(0, P), arm_labs)
  active <- matrix(1, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  
  for(i in 1:20) {
    # Generate Data
    n_new <- table(factor(sample.int(P, M, replace = TRUE, prob = p[i, ]), levels = 1:P))
    if(i == 1) {
      n[i, ] <- n_new
      y[i, ] <- rbinom(P, n_new, plogis(mu))
    } else {
      n[i, ] <- n[i - 1, ] + n_new
      y[i, ] <- y[i - 1, ] + rbinom(P, n_new, plogis(mu))
    }
    
    # Fit model
    mod <- vb_mod(y[i, ], n[i, ])
    m[i, ] <- mod$mu
    v[i, ] <- diag(mod$Sigma)
    
    # Compute posterior quantities
    draws <- mvnfast::rmvn(1e4, m[i, ], sigma = mod$Sigma)
    beta_draws <- draws %*% X_con_inv_t_Q_t
    p_sup[i, ] <- prob_each_superior_all(draws, delta_sup)
    p_sup_trt[i, ] <- prob_each_superior_all(draws[, -1], delta_sup)
    p_max_all[i, ] <- prob_max(draws)
    p_max[i, ] <- prob_max(draws[, -1])
    p_max_mes[i, ] <- prob_max(beta_draws[, 3:6])
    p_max_tim[i, ] <- prob_max(beta_draws[, 7:9])
    p_beat_ctrl[i, ] <- prob_superior(draws[, -1], draws[, 1], delta_ctr)
    best[i] <- unname(which.max(p_max[i, ]))
    
    # Deactivation rule
    if(!ind_comp_ctrl) {
      # active[i, ] <- as.numeric(p_sup[i, -1] > kappa_lo[i])
      active[i, ] <- as.numeric(p_max[i, ] > kappa_act[i] & p_beat_ctrl[i, ] > 1 - kappa_ctr[i])
      is_sup <- p_sup[i, ] > kappa_sup[i]
      superior <- any(is_sup)
    } else {
      active[i, ] <- as.numeric(p_max[i, ] > kappa_act[i] & p_beat_ctrl[i, ] > 1 - kappa_ctr[i])  
      is_sup <- c("00" = FALSE, p_max[i, ] > kappa_sup[i] & p_beat_ctrl[i, ] > kappa_ctr[i])
      superior <- any(is_sup)
    }
    if(sum(active[i, ]) > 1) {
      # Probability all active noninferior to superior by delta
      p_noninf[i] <- prob_all_superior(
        draws[, -1][, active[i, ] & !(1:(P-1) == best[i]), drop = F],
        draws[, -1][, best[i]], 
        -delta_noninf)
      # Probability best active superior to all active and control
      p_best_beat_inactive[i] <- prob_superior_all(
        draws[, -1][, best[i]], 
        draws[, c(1, which(active[i, ] == 0) + 1), drop = F],
        delta_sup)
    }
    noninferior <- any(p_noninf[i] > kappa_noninf[i] & p_best_beat_inactive[i] > kappa_sup[i])
    # nonsuperior <- all(!active[i, ])
    nonsuperior <- max(p_sup_trt[i, ]) < kappa_nonsup[i]
    lose <- all(!active[i, ])
    stopped <- superior | noninferior | nonsuperior | lose
    if(brar) {
      if(!allocate_inactive) {
        w <- sqrt(p_max[i, ] * v[i, -1] / (n[i, -1] + 1))
        w[!active[i, ]] <- 0
        p[i + 1, -1] <- (1 - p[i, 1]) * w / sum(w)   
      } else {
        w <- sqrt(p_max[i, ] * v[i, -1] / (n[i, -1] + 1))
        p[i + 1, -1] <- (1 - p[i, 1]) * w / sum(w)
      }
    } else{
      if(!allocate_inactive) {
        p[i + 1, -1] <- (1 - p[i, 1]) * active[i, ] / sum(active[i, ])    
      } else {
        p[i + 1, ] <- p[i, ]
      }
    }
    if(stopped) break
  }
  
  if(return_all) {ret_seq <- seq_len(i)} else {ret_seq <- i}  
  
  return(
    list(
      id = id,
      mu = mu,
      delta_sup = delta_sup,
      delta_ctr = delta_ctr,
      delta_noninf = delta_noninf,
      delta_pair = delta_pair,
      kappa_act_0 = kappa_act_0,
      kappa_act_1 = kappa_act_1,
      kappa_sup_0 = kappa_sup_0,
      kappa_sup_1 = kappa_sup_1,
      kappa_ctr_0 = kappa_ctr_0,
      kappa_ctr_1 = kappa_ctr_1,
      kappa_noninf_0 = kappa_noninf_0,
      kappa_noninf_1 = kappa_noninf_1,
      kappa_nonsup_0 = kappa_nonsup_0,
      kappa_nonsup_1 = kappa_nonsup_1,
      brar = brar,
      allocate_inactive = allocate_inactive,
      interim = i,
      stopped = stopped,
      superior = superior,
      noninferior = noninferior,
      nonsuperior = nonsuperior,
      lose = lose,
      p = p[ret_seq, ],
      n = n[ret_seq, ],
      y = y[ret_seq, ],
      m = m[ret_seq, ],
      v = v[ret_seq, ],
      p_sup = p_sup[ret_seq, ],
      p_sup_trt = p_sup_trt[ret_seq, ],
      p_max_all = p_max_all[ret_seq, ],
      p_max = p_max[ret_seq, ],
      p_max_mes = p_max_mes[ret_seq, ],
      p_max_tim = p_max_tim[ret_seq, ],
      p_beat_ctrl = p_beat_ctrl[ret_seq, ],
      p_noninf = p_noninf[ret_seq],
      p_best_beat_inactive = p_best_beat_inactive[ret_seq],
      p_sup_pairwise = pairwise_superiority_all(draws, delta_pair, replace = TRUE),
      p_sup_mes_pairwise = pairwise_superiority_all(beta_draws[, 3:6], delta_pair, replace = TRUE),
      p_sup_tim_pairwise = pairwise_superiority_all(beta_draws[, 7:9], delta_pair, replace = TRUE),
      active = active[ret_seq, ],
      best = best[ret_seq],
      sup = is_sup,
      beat_ctrl = p_beat_ctrl[i, ] > kappa_ctr[i]
    )
  )
}
