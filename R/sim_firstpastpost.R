#' Run a first past the post trial
#' 
#' Run a single trial using where stopping occurs as soon as one
#' arm is superior to control.
#' 
#' @param id The trial ID
#' @param mu The true mean, must be length 13
#' @param delta The superiority margin
#' @param kappa_lo_0 The starting threshold to deactivate poor arms
#' @param kappa_lo_1 The final threshold to deactivate poor arms
#' @param kappa_hi_0 The starting threshold for superiority
#' @param kappa_hi_1 The final threshold for superiority
#' @param brar Use Response Adaptive Randomisation?
#' @param allocate_inactive Continue to allocate subjects to arms meeting the inactive threshold?
#' @param return_all Return all value, or only those from final analysis?
#'
#' @return A list of trial quantities
#'
#' @export
run_a_firstpastpost_trial <- function(
  id,
  mu,
  delta,
  kappa_lo_0 = 0.10,
  kappa_lo_1 = 0.10,
  kappa_hi_0 = 0.95,
  kappa_hi_1 = 0.95,
  brar = FALSE,
  allocate_inactive = FALSE,
  return_all = FALSE
) {
  
  # Setup
  N <- 10000
  K <- 20
  M <- 500
  P <- 13
  arm_labs <- sprintf("%02d", 0:12)
  mes_labs <- paste0("m", 1:4)
  tim_labs <- paste0("t", 1:3)
  
  kappa_lo <- thres_seq(kappa_lo_0, kappa_lo_1, 1/2, K - 1)
  kappa_hi <- thres_seq(kappa_hi_0, kappa_hi_1, 1/2, K - 1)
  
  p <- matrix(1/P, K + 1, P, dimnames = list("interim" = 0:K, "arm" = arm_labs))
  n <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  y <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  m <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  v <- matrix(0, K, P, dimnames = list("interim" = 1:K, "arm" = arm_labs))
  
  p_max <- matrix(0, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  p_max_mes <- matrix(0, K, 4, dimnames = list("interim" = 1:K, "mes" = 1:4))
  p_max_tim <- matrix(0, K, 3, dimnames = list("interim" = 1:K, "tim" = 1:3))
  p_beat_ctrl <- matrix(0, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  
  best <- rep(0, K)
  active <- matrix(1, K, P - 1, dimnames = list("interim" = 1:K, "arm" = arm_labs[-1]))
  
  for(i in 1:20) {
    n_new <- table(factor(sample.int(P, M, replace = TRUE, prob = p[i, ]), levels = 1:P))
    if(i == 1) {
      n[i, ] <- n_new
      y[i, ] <- rbinom(P, n_new, plogis(mu))
    } else {
      n[i, ] <- n[i - 1, ] + n_new
      y[i, ] <- y[i - 1, ] + rbinom(P, n_new, plogis(mu))
    }
    mod <- vb_mod(y[i, ], n[i, ])
    m[i, ] <- mod$mu
    v[i, ] <- diag(mod$Sigma)
    draws <- mvnfast::rmvn(1e4, m[i, ], sigma = mod$Sigma)
    beta_draws <- draws %*% X_con_inv_t_Q_t
    p_max[i, ] <- prob_max(draws[, -1])
    p_max_mes[i, ] <- prob_max(beta_draws[, 3:6])
    p_max_tim[i, ] <- prob_max(beta_draws[, 7:9])
    p_beat_ctrl[i, ] <- prob_superior(draws[, -1], draws[, 1], delta)
    best[i] <- unname(which.max(p_max[i, ]))
    active[i, ] <- as.numeric(p_beat_ctrl[i, ] > kappa_lo[i])
    superior <- any(p_beat_ctrl[i, ] > kappa_hi[i])
    futile <- all(!active[i, ])
    
    stopped <- superior | futile
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
      mu = mu,
      delta = delta,
      kappa_lo_0 = kappa_lo_0,
      kappa_lo_1 = kappa_lo_1,
      kappa_hi_0 = kappa_hi_0,
      kappa_hi_1 = kappa_hi_1,
      brar = brar,
      allocate_inactive = allocate_inactive,
      interim = i,
      stopped = stopped,
      superior = superior,
      futile = futile,
      p = p[ret_seq, ],
      n = n[ret_seq, ],
      y = y[ret_seq, ],
      m = m[ret_seq, ],
      v = v[ret_seq, ],
      p_max = p_max[ret_seq, ],
      p_max_mes = p_max_mes[ret_seq, ],
      p_max_tim = p_max_tim[ret_seq, ],
      p_beat_ctrl = p_beat_ctrl[ret_seq, ],
      p_sup_pairwise = pairwise_superiority_all(draws, 0, replace = TRUE),
      p_sup_mes_pairwise = pairwise_superiority_all(beta_draws[, 3:6], 0, replace = TRUE),
      p_sup_tim_pairwise = pairwise_superiority_all(beta_draws[, 7:9], 0, replace = TRUE),
      active = active[ret_seq, ],
      best = best[ret_seq],
      beat_ctrl = p_beat_ctrl[i, ] > kappa_hi[i]
    )
  )
}