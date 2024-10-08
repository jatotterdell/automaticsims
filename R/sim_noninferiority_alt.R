#' Run a non-inferiority trial
#' 
#' As sim_noninferiority.R, except
#' rather than stopping for superiority by
#' delta, always base superiority on prob(max).
#' 
#' @param id The trial ID
#' @param mu The true mean, must be length 13
#' @param delta_sup The superiority margin
#' @param delta_ctr The superiority margin for control comparisons
#' @param delta_noninf The noninferiority margin
#' @param delta_pair The margin for pairwise comparisons
#' @param kappa_act_0 The starting threshold to deactivate poor arms
#' @param kappa_act_1 The final threshold to deactivate poor arms
#' @param kappa_sup_0 The starting threshold for superiority
#' @param kappa_sup_1 The final threshold for superiority
#' @param kappa_ctr_0 The starting threshold for beating control
#' @param kappa_ctr_1 THe final threshold for beating control
#' @param kappa_noninf_0 The starting threshold for non-inferiority
#' @param kappa_noninf_1 The final threshold for non-inferiority
#' @param kappa_nonsup_0 The starting threshold for non-superiority
#' @param kappa_nonsup_1 The final threshold for non-superiority
#' @param brar Use Response Adaptive Randomisation?
#' @param allocate_inactive Continue to allocate subjects to arms meeting the inactive threshold?
#' @param return_all Return all value, or only those from final analysis?
#' @param ind_comp_ctrl Should deactivation be based on independent 
#' P(max trt arms) and P(beat ctrl) or just global P(max all arms)
#'
#' @return A list of trial quantities
#'
#' @export
run_a_noninf_trial_alt <- function(
  id,
  mu,
  delta_sup = 0.1,
  delta_ctr = 0,
  delta_noninf = 0.1,
  delta_pair = 0,
  kappa_act_0 = 0.01,
  kappa_act_1 = 0.05,
  kappa_act_r = 0.5,
  kappa_sup_0 = 0.95,
  kappa_sup_1 = 0.80,
  kappa_sup_r = 0.5,
  kappa_ctr_0 = 0.95,
  kappa_ctr_1 = 0.95,
  kappa_ctr_r = 0.5,
  kappa_noninf_0 = 0.60,
  kappa_noninf_1 = 0.60,
  kappa_noninf_r = 0.5,
  kappa_nonsup_0 = 0.1,
  kappa_nonsup_1 = 0.1,
  kappa_nonsup_r = 0.5,
  brar = TRUE,
  active_pmax = TRUE,
  allocate_inactive = FALSE,
  return_all = FALSE,
  ind_comp_ctrl = FALSE,
  ctrl_alloc = 1/13,
  Nseq = seq(500, 10000, length.out = 20)
) {
  
  # Setup
  which_best   <- which(mu == max(mu)) - 1
  which_indiff <- which(mu >= max(mu) - delta_sup) - 1
  N <- max(Nseq)
  K <- length(Nseq)
  M <- diff(c(0, Nseq))
  P <- 13
  arm_labs <- sprintf("%02d", 0:12)
  mes_labs <- paste0("m", 1:4)
  tim_labs <- paste0("t", 1:3)
  
  kappa_act <- thres_seq(kappa_act_0, kappa_act_1, kappa_act_r, K - 1)
  kappa_sup <- thres_seq(kappa_sup_0, kappa_sup_1, kappa_sup_r, K - 1)
  kappa_ctr <- thres_seq(kappa_ctr_0, kappa_ctr_1, kappa_ctr_r, K - 1)
  kappa_noninf <- thres_seq(kappa_noninf_0, kappa_noninf_1, kappa_noninf_r, K - 1)
  kappa_nonsup <- thres_seq(kappa_nonsup_0, kappa_nonsup_1, kappa_nonsup_r, K - 1)
  
  p <- matrix(1/P, K + 1, P, dimnames = list("interim" = 0:K, "arm" = arm_labs))
  p[, 1] <- ctrl_alloc
  p[, -1] <- (1 - p[, 1])/12
  n <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  y <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  m <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  v <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  mu_maxes <- matrix(0, K, 5, dimnames = list("interim" = 1:K, "val" = c("mu", "sig", "lo", "hi", "wi")))
  
  p_sup <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  p_sup_trt <- matrix(0, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  p_max_all <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  p_max <- matrix(0, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  p_max_mes <- matrix(0, K, 4, dimnames = list("interim" = 1:K, "mes" = 1:4))
  p_max_tim <- matrix(0, K, 3, dimnames = list("interim" = 1:K, "tim" = 1:3))
  p_beat_ctrl <- matrix(0, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  p_noninf <- rep(0, K)
  p_best_beat_inactive <- rep(0, K)
  p_ctrl_beat_all_trt <- rep(0, K)
  
  best_trt <- rep(0, K)
  best <- rep(0, K)
  is_sup <- setNames(rep(0, P), arm_labs)
  active <- matrix(1, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  
  for(i in 1:K) {
    # Generate Data
    n_new <- table(factor(sample.int(P, M[i], replace = TRUE, prob = p[i, ]), levels = 1:P))
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
    max_draws <- matrixStats::rowMaxs(draws)
    
    mu_maxes[i, 1:4] <- c(mean(max_draws), var(max_draws), setNames(quantile(max_draws, prob = c(0.025, 0.975)), c("lo", "hi")))
    mu_maxes[i, 5] <- mu_maxes[i, 4] - mu_maxes[i, 3]
    
    p_sup[i, ] <- prob_each_superior_all(draws, delta_sup)
    p_sup_trt[i, ] <- prob_each_superior_all(draws[, -1], delta_sup)
    p_max_all[i, ] <- prob_max(draws)
    p_max[i, ] <- prob_max(draws[, -1])
    p_max_mes[i, ] <- prob_max(beta_draws[, 3:6])
    p_max_tim[i, ] <- prob_max(beta_draws[, 7:9])
    p_beat_ctrl[i, ] <- prob_superior(draws[, -1], draws[, 1], delta_ctr)
    p_ctrl_beat_all_trt[i] <- prob_superior_all(draws[, 1], draws[, -1], delta_ctr)
    best_trt[i] <- unname(which.max(p_max[i, ]))
    best[i] <- unname(which.max(p_max_all[i, ])) - 1
    
    # Superiority rule
    if(!ind_comp_ctrl) {
      is_sup <- p_max_all[i, ] > kappa_sup[i]
      superior <- any(is_sup)
      
    } else {
      is_sup <- p_max[i, ] > kappa_sup[i] & p_beat_ctrl[i, ] > kappa_ctr[i]
      superior <- any(is_sup)
    }
    
    # Activation rule
    if(active_pmax) {
      if(!ind_comp_ctrl) {
        active[i, ] <- as.numeric(p_max_all[i, -1] > kappa_act[i])
      } else {
        active[i, ] <- as.numeric(p_max[i, ] > kappa_act[i] & p_beat_ctrl[i, ] > 1 - kappa_ctr[i]) 
      }
      
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
    } else {
      if(brar) {
        if(!allocate_inactive) {
          w <- sqrt(p_max[i, ] * v[i, -1] / (n[i, -1] + 1))
          w[!active[i, ]] <- 0
          p[i + 1, -1] <- (1 - p[i, 1]) * w / sum(w)
          active[i, ] <- p[i + 1, -1] > kappa_act[i]
          p[i + 1, -1][!active[i, ]] <- 0
          p[i + 1, -1] <- (1 - p[i, 1]) * p[i + 1, -1]/ sum(p[i + 1, -1])
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
    }


    
    if(sum(active[i, ]) > 1) {
      # Probability all active noninferior to superior by delta
      p_noninf[i] <- prob_all_superior(
        draws[, -1][, active[i, ] & !(1:(P-1) == best_trt[i]), drop = F],
        draws[, -1][, best_trt[i]], 
        -delta_noninf)
      # Probability best active superior to all active and control
      if(!ind_comp_ctrl) {
        p_best_beat_inactive[i] <- prob_superior_all(
          draws[, -1][, best_trt[i]], 
          draws[, c(1, which(active[i, ] == 0) + 1), drop = F], 0) 
      } else {
        p_best_beat_inactive[i] <- prob_superior_all(
          draws[, -1][, best_trt[i]], 
          draws[, which(active[i, ] == 0) + 1, drop = F], 0) 
      }

    }
    
    # Stopping flags
    if(!ind_comp_ctrl) {
      noninferior <- any(p_noninf[i] > kappa_noninf[i] & p_best_beat_inactive[i] > kappa_sup[i])
      nonsuperior <- max(p_sup[i, ]) < kappa_nonsup[i]
    } else {
      noninferior <- any(p_noninf[i] > kappa_noninf[i] & p_best_beat_inactive[i] > kappa_sup[i] & p_beat_ctrl[best_trt[i], ] > kappa_ctr[i])
      nonsuperior <- max(p_sup_trt[i, ]) < kappa_nonsup[i]
    }
    nonsuperior <- max(p_sup_trt[i, ]) < kappa_nonsup[i]
    lose <- all(p_beat_ctrl[i, ] < 1 - kappa_ctr[i]) # Everything worse than control so may as well stop
    stopped <- (superior | noninferior | nonsuperior | lose) & (i < K)
    
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
      n_first = Nseq[1],
      brar = brar,
      allocate_inactive = allocate_inactive,
      ctrl_alloc = ctrl_alloc,
      interim = i,
      N = Nseq[i],
      stopped = stopped,
      superior = superior,
      noninferior = noninferior,
      nonsuperior = nonsuperior,
      lose = lose,
      p = p[ret_seq + 1, ],
      n = n[ret_seq, ],
      y = y[ret_seq, ],
      m = m[ret_seq, ],
      v = v[ret_seq, ],
      mu_maxes = mu_maxes[ret_seq, ],
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
      p_sup_trt_ctr = mean(beta_draws[, 2] > 0),
      active = active[ret_seq, ],
      best_trt = best_trt[ret_seq],
      best = best[ret_seq],
      best_in_best = best[i] %in% which_best,
      best_in_indiff = best[i] %in% which_indiff,
      sup = is_sup,
      beat_ctrl = p_beat_ctrl[i, ] > kappa_ctr[i]
    )
  )
}
