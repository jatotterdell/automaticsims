#' Run automatic trial simulation
#' 
#' @param mu The true response rate in each arm
#' @param n_seq The interim sequence
#' @param delta_sup The reference value for meaningful difference
#' @param thres_sup Threshold for superiority
#' @param thres_inf Threshold for inferiority
#' @param thres_fut Threshold for futility
#' @param thres_nin Threshold for non-inferiority
#' @param thres_eff Threshold for effectiveness
#' @param thres_equ Threshold for equality
#' @param thres_ina Threshold for inadequacy
#' @param use_mwu Use mass-weighted-urn randomisation
#' @param mc_draws Number of Monte Carlo draws for posterior quantity calculation
#' @param rar_scale Scaling factor for RAR
#' @param rar_rule Which RAR rule to use
#' @param drop_rule Which rule to use for dropping control
#' @param ct_alloc What fixed allocation to control
#' @param drop_perm Permanently drop interventions?
#' @param enable_stopping Apply stopping rule or let trial run its course to max SS?
#' @export
run_automatic_trial <- function(
  mu,
  Xdes = 1,
  mu0 = c(qlogis(0.8), rep(0, 12)),
  Sigma0 = diag(c(2^2, rep(1, 12))),
  n_seq = seq(500, 10000, 500),
  delta_sup = log(1.1),
  thres_sup = function(t) 0.80,
  thres_inf = function(t) 1 - thres_sup(t),
  thres_fut = 0.95,
  thres_nin = 0.90,
  thres_eff = function(t) 0.99,
  thres_hrm = function(t) 0.01,
  thres_equ = 0.90,
  thres_ina = 0.99,
  thres_ctr = 0.99,
  use_mwu = TRUE,
  mc_draws = 1e4,
  rar_scale = 0.5,
  rar_min  = 0.05/12,
  rar_rule = 1,
  drop_rule = 1,
  stop_rule = 1,
  ct_alloc = sqrt(12) / (12 + sqrt(12)),
  drop_perm = TRUE,
  enable_stopping = FALSE,
  return_all = TRUE,
  gen_prior = FALSE,
  use_info_n = TRUE
) {
  
  if(gen_prior) {
    mu <- drop(mvnfast::rmvn(1, mu0, Sigma0))
  }
  
  # Truth
  which_eff   <- which(mu[-1] > mu[1])
  which_ineff <- setdiff(1:12, which_eff)
  which_sup   <- which(sapply(1:(length(mu)-1), function(a) all(mu[-1][a] >= mu[-1][-a])))
  which_inf   <- setdiff(1:12, which_sup)
  which_fut   <- which(sapply(1:(length(mu)-1), function(a) all(mu[-1][a] < mu[-1][-a] + delta_sup)))
  which_ina   <- which(mu[-1] < mu[-1] + delta_sup)
  which_ada   <- setdiff(1:12, which_ina)
  num_eff <- length(which_eff)
  num_sup <- length(which_sup)
  num_fut <- length(which_fut)
  num_ina <- length(which_ina)
  
  # Design
  if(Xdes == 1) {
    X <- X_con
  } else if(Xdes == 2) {
    X <- X_con_alt
  } else if(Xdes == 3) {
    X <- X_trt
  } else if(Xdes == 4) {
    X <- X_pol
  } else {
    X <- diag(1, 13)
  }
  unconT <- Q %*% pinv_X_con %*% X
  ctrT   <- cbind(-1, diag(1, 12))
  bstT   <- cbind(-1, diag(1, 11))
  In     <- diag(1, 12)
  
  
  n_max <- max(n_seq)
  n_int <- length(n_seq)
  n_new <- diff(c(0, n_seq))
  n_arms <- nrow(Xmap)
  arm_labs <- rownames(Xmap)
  mes_labs <- paste0("m", 1:4)
  tim_labs <- paste0("t", 1:3)
  par_labs  <- paste0("b", 0:12)
  beta_labs <- c("b_ctr", "b_trt", paste0("b_m", 1:4), paste0("b_t", 1:3), 
                 paste0(rep(paste0("b_m", 1:4), each = 3), rep(paste0("t", 1:3), times = 4)))
  dn_all <- list("interim" = 1:n_int, "arm" = arm_labs)
  dn_trt <- list("interim" = 1:n_int, "arm" = arm_labs[-1])
  dn_par  <- list("interim" = 1:n_int, "parameter" = par_labs)
  
  # Data
  n <- matrix(0, n_int, n_arms, dimnames = dn_all)
  y <- matrix(0, n_int, n_arms, dimnames = dn_all)
  
  # Parameters
  arm_mean  <- matrix(0, n_int, n_arms, dimnames = dn_all)
  arm_var   <- matrix(0, n_int, n_arms, dimnames = dn_all)
  
  par_mean  <- matrix(0, n_int, n_arms, dimnames = dn_par)
  par_var   <- matrix(0, n_int, n_arms, dimnames = dn_par)
  
  eff_mean  <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  eff_var   <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  
  nin_mean  <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  nin_var   <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  
  beta_mean <- matrix(0, n_int, 21, dimnames = list("interim" = 1:n_int, "par" = beta_labs))
  beta_var  <- matrix(0, n_int, 21, dimnames = list("interim" = 1:n_int, "par" = beta_labs))
  
  ctr_mean   <- matrix(0, n_int, nrow(Xmain), dimnames = list("interim" = 1:n_int, "contrast" = rownames(Xmain)))
  ctr_var    <- matrix(0, n_int, nrow(Xmain), dimnames = list("interim" = 1:n_int, "contrast" = rownames(Xmain)))
  
  p_ctr_gt0  <- matrix(0, n_int, nrow(Xmain), dimnames = list("interim" = 1:n_int, "contrast" = rownames(Xmain)))
  p_ctr_ltd  <- matrix(0, n_int, nrow(Xmain), dimnames = list("interim" = 1:n_int, "contrast" = rownames(Xmain)))
  p_ctr_equ  <- matrix(0, n_int, nrow(Xmain), dimnames = list("interim" = 1:n_int, "contrast" = rownames(Xmain)))
  is_ctr_gt0 <- matrix(FALSE, n_int, nrow(Xmain), dimnames = list("interim" = 1:n_int, "contrast" = rownames(Xmain)))
  is_ctr_lt0 <- matrix(FALSE, n_int, nrow(Xmain), dimnames = list("interim" = 1:n_int, "contrast" = rownames(Xmain)))
  is_ctr_ltd <- matrix(FALSE, n_int, nrow(Xmain), dimnames = list("interim" = 1:n_int, "contrast" = rownames(Xmain)))
  
  # Store probability quantities
  p_rand    <- matrix(1/n_arms, n_int + 1, n_arms, dimnames = list("interim" = 0:n_int, "arm" = arm_labs))
  p_rand[, 1] <- ct_alloc
  p_rand[, -1] <- (1 - ct_alloc) / (n_arms - 1)
  
  p_sup     <- matrix(0, n_int, n_arms, dimnames = dn_all)
  p_inf     <- matrix(0, n_int, n_arms, dimnames = dn_all)
  p_sup_trt <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  p_fut_trt <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  p_sup_act <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  p_eff     <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  p_equ     <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  p_ina     <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  p_nin     <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  # p_best_mes <- matrix(0, n_int, 4, dimnames = list("interim" = 1:n_int, "message" = mes_labs))
  # p_best_tim <- matrix(0, n_int, 3, dimnames = list("interim" = 1:n_int, "timing" = tim_labs))
  # p_in_best  <- matrix(0, n_int, ncol(Xmap), dimnames = list("interim" = 1:n_int, "treatment" = colnames(Xmap)))
  
  # Store flags
  is_act <- matrix(TRUE, n_int, n_arms, dimnames = dn_all)      # Active
  is_sup <- matrix(FALSE, n_int, n_arms, dimnames = dn_all) 
  is_inf <- matrix(FALSE, n_int, n_arms, dimnames = dn_all) 
  is_sup_trt <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Superior (best)
  is_sup_act <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Superior (best)
  is_inf_trt <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Superior (best)
  is_fut <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Futile (not best by enough)
  is_nin <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Non-inferior (close enough to best)
  is_ina <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Inadequate (not effective enough)
  is_eff <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Effective (better than control)
  is_hrm <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Harmful (worse than control)
  is_equ <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Equivalent (same as control)
  # is_crb <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Current best (currently best)
  
  stopped <- FALSE
  
  for(i in 1:n_int) {
    
    # Generate allocations
    if(use_mwu) {
      x <- factor(mass_weighted_urn_design(p_rand[i, ], n_new[i], alpha = 5)$trt, levels = 1:n_arms)
      x_agg <- table(x)
    } else {
      x <- factor(sample.int(n_arms, n_new[i], replace = TRUE, prob = p_rand[i, ]), levels = 1:n_arms)
      x_agg <- table(x)
    }
    
    # Generate data
    if(i == 1) {
      n[i, ] <- x_agg
      y[i, ] <- rbinom(n_arms, x_agg, plogis(mu))
    } else {
      n[i, ] <- n[i - 1, ] + x_agg
      y[i, ] <- y[i - 1, ] + rbinom(n_arms, x_agg, plogis(mu))
    }
    
    # Current information
    if(use_info_n) {
      info <- sum(n[i, ]) / n_max  
    } else {
      info <- i / n_int
    }
    
    
    # Fit model
    mod <- varapproxr::vb_logistic_n(
      X, y[i, ], n[i, ], mu0, Sigma0, 
      rep(0, 13), diag(1, 13), alg = "sj", maxiter_jj = 100)
    
    # Summarise parameters
    eta_mu      <- X %*% mod$mu
    eta_sigma   <- X %*% mod$Sigma %*% t(X)
    beta_mu     <- unconT %*% mod$mu
    beta_sigma  <- unconT %*% mod$Sigma %*% t(unconT)
    theta_mu    <- ctrT %*% eta_mu
    theta_sigma <- ctrT %*% eta_sigma %*% t(ctrT)
    
    par_mean[i, ]  <- drop(mod$mu)
    par_var[i, ]   <- diag(mod$Sigma)
    arm_mean[i, ]  <- drop(eta_mu)
    arm_var[i, ]   <- diag(eta_sigma)
    beta_mean[i, ] <- drop(beta_mu)
    beta_var[i, ]  <- diag(beta_sigma)
    eff_mean[i, ]  <- drop(theta_mu)
    eff_var[i, ]   <- diag(theta_sigma)
    ctr_mean[i, ]  <- drop(Xmain %*% beta_mu)
    ctr_var[i, ]   <- diag(Xmain %*% beta_sigma %*% t(Xmain))
    
    # Probabilities and decisions
    p_eff[i, ]  <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
    p_ina[i, ]  <- pnorm(delta_sup, eff_mean[i, ], sqrt(eff_var[i, ]))
    is_eff[i, ] <- p_eff[i, ] > thres_eff(info)
    is_hrm[i, ] <- p_eff[i, ] < thres_hrm(info)
    is_ina[i, ] <- p_ina[i, ] > thres_ina
    
    p_ctr_gt0[i, ]  <- 1 - pnorm(0, ctr_mean[i, ], sqrt(ctr_var[i, ]))
    p_ctr_ltd[i, ]  <- pnorm(delta_sup, ctr_mean[i, ], sqrt(ctr_var[i, ]))
    p_ctr_equ[i, ]  <- pnorm(delta_sup, ctr_mean[i, ], sqrt(ctr_var[i, ])) - pnorm(-delta_sup, ctr_mean[i, ], sqrt(ctr_var[i, ]))
    is_ctr_gt0[i, ] <- p_ctr_gt0[i, ] > thres_ctr
    is_ctr_lt0[i, ] <- p_ctr_gt0[i, ] < 1 - thres_ctr
    is_ctr_ltd[i, ] <- p_ctr_ltd[i, ] > thres_ina
    
    # Compute quantities requiring MC
    draws           <- mvnfast::rmvn(mc_draws, eta_mu, sigma = eta_sigma)
    p_sup[i, ]      <- prob_max(draws)
    p_sup_trt[i, ]  <- prob_max(draws[, -1])
    is_sup[i, ]     <- p_sup[i, ] > thres_sup(info)
    is_sup_trt[i, ] <- p_sup_trt[i, ] > thres_sup(info)
    
    
    if(i == 1) {
      p_sup_act[i, ]  <- p_sup_trt[i, ]
      is_inf[i, ]     <- p_sup[i, ] < thres_inf(info) / 13
      is_inf_trt[i, ] <- p_sup_act[i, ] < thres_inf(info) / 12   
    } else {
      p_sup_act[i, is_act[i - 1, -1]] <- prob_max(draws[, -1][, is_act[i - 1, -1], drop = F])
      is_inf[i, ]     <- p_sup[i, ] < thres_inf(info) / (13 - !is_act[i-1, 1])
      is_inf_trt[i, ] <- p_sup_act[i, ] < thres_inf(info) / (12 - drop_perm * sum(!is_act[i-1, -1])) # If permanent dropping, increase threshold
    }
    is_sup_act[i, ] <- p_sup_act[i, ] > thres_sup(info)

    # Relative to current best
    best               <- which.max(p_sup_trt[i, ])
    prmT               <- bstT %*% In[c(best, setdiff(1:12, best)), ]
    zeta_mu            <- prmT %*% eta_mu[-1]
    zeta_var           <- prmT %*% eta_sigma[-1, -1] %*% t(prmT)
    nin_mean[i, best]  <- 0
    nin_var[i, best]   <- 0
    nin_mean[i, -best] <- drop(zeta_mu)
    nin_var[i, -best]  <- diag(zeta_var)
    p_nin[i, ]         <- 1 - pnorm(-delta_sup, nin_mean[i, ], sqrt(nin_var[i, ]))
    is_nin[i, ]        <- p_nin[i, ] > thres_nin
    
    # Update active arms and allocations
    
    if(drop_rule == 1) { 
      # Only drop control if any effective or effective on average, drop intervention only if harmful
      is_act[i, 1]  <- !any(is_eff[i, ]) & !is_ctr_gt0[i, "b_trt"]
      is_act[i, -1] <- !is_hrm[i, ]
    } else if(drop_rule == 2) {
      # Only drop control if effective on average, drop intervention only if harmful
      is_act[i, 1]  <- !is_ctr_gt0[i, "b_trt"]
      is_act[i, -1] <- !is_hrm[i, ]
    } else if(drop_rule == 3) {
      # Never drop control, drop active intervention if inferior or harmful (original automatic dropping rule) 
      is_act[i, -1] <- !is_inf_trt[i, ] & !is_hrm[i, ] 
    } else if(drop_rule == 4) {
      # Never drop control, drop interventions if harmful 
      is_act[i, -1] <- !is_hrm[i, ]
    } else if(drop_rule == 5) { 
      # Drop any that are inferior including control
      is_act[i, ] <- !is_inf[i, ]
    } else if(drop_rule == 6) {
      # Only drop control if any effective or effective on average, drop intervention only if inadequate
      is_act[i, 1]  <- !any(is_eff[i, ]) & !is_ctr_gt0[i, "b_trt"]
      is_act[i, -1] <- !is_ina[i, ]
    } else if(drop_rule == 7) {
      # Only drop control if any effective or effective on average, drop intervention if harmful or inefrior
      is_act[i, 1]  <- !any(is_eff[i, ]) & !is_ctr_gt0[i, "b_trt"]
      is_act[i, -1] <- !is_inf_trt[i, ] & !is_hrm[i, ]    
    }
    
    if(i > 1) {
      # Always drop control permanently
      is_act[i, 1] <- as.logical(is_act[i, 1] * is_act[i - 1, 1])
      # Drop interventions permanently?
      if(drop_perm) {
        is_act[i, ] <- as.logical(is_act[i, ] * is_act[i - 1, ])  
      }
    }
    if(!is_act[i, 1]) p_rand[i + 1, 1] <- 0

    # Calculate active dependent quantities
    
    # Update allocations
    if(rar_rule == 1) {
      w <- (p_sup_trt[i, ])^rar_scale  
    } else if (rar_rule == 2) {
      w <- (p_sup_trt[i, ] * arm_var[i, -1] / (n[i, -1] + 1))^rar_scale  
    } else if (rar_rule == 3) {
      w <- (p_eff[i, ])^rar_scale
    } else if (rar_rule == 4) {
      w <- (p_eff[i, ] * eff_var[i, ] / (n[i, -1] + 1))^rar_scale
    } else if (rar_rule == 5) { # Set a minimum P(best) amount to be included in RAR, otherwise zero but not permanently dropped.
      w <- (p_sup_trt[i, ] * (p_sup_trt[i, ] > rar_min))^rar_scale  
    } else if (rar_rule == 6) {
      w <- (p_sup_trt[i, ] * (p_sup_trt[i, ] > rar_min) * arm_var[i, -1] / (n[i, -1] + 1))^rar_scale  
    }
    
    # Fix the allocation to zero if inactive
    w[!is_act[i, -1]] <- 0
    # If nothing is active re-activate control otherwise apply RAR
    if(all(!is_act[i, -1])) {
      p_rand[i + 1, 1] <- 1
      p_rand[i + 1, -1] <- 0
    } else {
      p_rand[i + 1, -1] <- (1 - p_rand[i + 1, 1]) * w / sum(w)  
    }
    
    # Assess stopping criteria
    # If 
    #  - any is best and control inactive
    #  - all are inactive but control
    #  - average effect no better than control
    any_best     <- any(is_sup_act[i, ] & !is_act[i, 1])
    all_inactive <- !any(is_act[i, -1])
    avg_ina      <- is_ctr_ltd[i, "b_trt"]
    avg_hrm      <- is_ctr_lt0[i, "b_trt"]
    if(enable_stopping) {
      if(stop_rule == 1) {
        # Stop if any single best, or all inactive, or on average inadequate
        stopped <- any_best | all_inactive | avg_ina    
      } else if(stop_rule == 2) {
        # Stop if any single best or all inactive
        stopped <- any_best | all_inactive   
      } else if(stop_rule == 3) {
        # Original stopping rule, any single best, or all non-superior, or all inactive
        p_fut_trt[i, ]  <- 1 - prob_each_superior_all(draws[, -1], delta_sup)
        all_fut         <- min(p_fut_trt[i, ]) > thres_fut
        stopped         <- any_best | all_fut | all_inactive
      } else if(stop_rule == 4) {
        # Stop if any single best, or all inactive, or on average harmful
        stopped <- any_best | all_inactive | avg_hrm     
      }
    }
    
    if(stopped) break
  }
  p_rand <- p_rand[-1, , drop = F]
  
  dec_sup_trt_at <- apply(is_sup_trt, 2, findfirst)
  dec_inf_at <- apply(is_inf, 2, findfirst)
  dec_eff_at <- apply(is_eff, 2, findfirst)
  dec_hrm_at <- apply(is_hrm, 2, findfirst)
  dec_ina_at <- apply(is_ina, 2, findfirst)
  dec_drp_at <- apply(!is_act, 2, findfirst)
  
  dec_ctr_gt0_at <- apply(is_ctr_gt0, 2, findfirst)
  dec_ctr_ltd_at <- apply(is_ctr_ltd, 2, findfirst)
  
  if(return_all) {
    idx <- 1:i  
  } else {
    idx <- i
  }
  
  # Arm specific quantities
  arm_quantities <- list(
    n          = n[idx, , drop = F],
    y          = y[idx, , drop = F],
    p_rand     = p_rand[idx, , drop = F],
    arm_mean   = arm_mean[idx, , drop = F], 
    arm_var    = arm_var[idx, , drop = F],
    eff_mean   = eff_mean[idx, , drop = F],
    eff_var    = eff_var[idx, , drop = F],
    nin_mean   = nin_mean[idx, , drop = F],
    nin_var    = nin_var[idx, , drop = F],
    p_sup      = p_sup[idx, , drop = F],
    p_sup_trt  = p_sup_trt[idx, , drop = F],
    p_sup_act  = p_sup_act[idx, , drop = F],
    p_eff      = p_eff[idx, , drop = F],
    p_ina      = p_ina[idx, , drop = F],
    p_nin      = p_nin[idx, , drop = F],
    is_sup     = is_sup[idx, , drop = F],
    is_sup_trt = is_sup_trt[idx, , drop = F],
    is_sup_act = is_sup_act[idx, , drop = F],
    is_inf_trt = is_inf_trt[idx, , drop = F],
    is_inf     = is_inf[idx, , drop = F],
    is_eff     = is_eff[idx, , drop = F],
    is_hrm     = is_hrm[idx, , drop = F],
    is_nin     = is_nin[idx, , drop = F],
    is_ina     = is_ina[idx, , drop = F],
    is_act     = is_act[idx, , drop = F]
  )
  if(stop_rule == 3) arm_quantities <- c(arm_quantities, list(p_fut_trt = p_fut_trt[idx, , drop = F]))
  
  # Contrast specific quantities
  ctr_quantities <- list(
    ctr_mean = ctr_mean[idx, , drop = F],
    ctr_var  = ctr_var[idx, , drop = F],
    p_ctr_gt0 = p_ctr_gt0[idx, , drop = F],
    p_ctr_ltd = p_ctr_ltd[idx, , drop = F],
    p_ctr_equ = p_ctr_equ[idx, , drop = F],
    is_ctr_gt0 = is_ctr_gt0[idx, , drop = F],
    is_ctr_ltd = is_ctr_ltd[idx, , drop = F]
  )
  
  final <- list(
    
    trial_quantitites = list(
      true_par  = setNames(drop(mu), colnames(par_mean)),
      true_eta  = setNames(drop(X_con %*% mu), colnames(arm_mean)),
      true_beta = setNames(drop(unconT %*% mu), colnames(beta_mean))
    ),
    
    result_quantities = list(
      stop_early   = stopped,
      stop_at      = ifelse(stopped, i, NA),
      any_best     = any_best,
      all_inactive = all_inactive,
      avg_ina      = avg_ina,
      all_fut      = all(is_fut[i, ]),
      true_sup     = any(is_sup[i, which_sup]),
      false_sup    = any(is_sup[i, which_inf]),
      false_inf    = any(is_inf[i, which_sup]),
      any_eff      = any(is_eff[i, ]),
      any_false_eff= any(is_eff[i, which_ineff]),
      any_true_eff = any(is_eff[i, which_eff]),
      all_true_eff = ifelse(num_eff > 0, all(is_eff[i, which_eff]), FALSE)
    ),
    
    arm_quantities = arm_quantities,
    ctr_quantities = ctr_quantities,
  
    par1_quantities = list(
      par_mean = par_mean[idx, , drop = F],
      par_var  = par_var[idx, , drop = F]
    ),
    
    par2_quantities = list(
      beta_mean = beta_mean[idx, , drop = F],
      beta_var  = beta_var[idx, , drop = F]
    ),
    
    dec1_quantities = loo::nlist(
      dec_sup_trt_at, 
      dec_inf_at, 
      dec_eff_at, 
      dec_hrm_at,
      dec_ina_at, 
      dec_drp_at
    ),
    
    dec2_quantities = loo::nlist(
      dec_ctr_gt0_at,
      dec_ctr_ltd_at
    )
    
  ) 
  return(final)
}


#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `run_automatic_trial` as a list
#' @param final Return data from final analysis or interims
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom purrr reduce
tibble_arm_quantities <- function(dat, final = TRUE, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["arm_quantities"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "arm", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "arm")) %>%
      dplyr::mutate(arm = forcats::fct_inorder(arm)) %>%
      dplyr::arrange(interim, arm)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial), interim = as.numeric(interim)) %>%
    dplyr::arrange(trial, interim)
}

#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `run_automatic_trial` as a list
#' @param final Return data from final analysis or interims
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom purrr reduce
tibble_treat_quantities <- function(dat, final = TRUE, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["treat_quantities"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "treatment", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "treatment")) %>%
      dplyr::mutate(treatment = forcats::fct_inorder(treatment)) %>%
      dplyr::arrange(interim, treatment)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial), interim = as.numeric(interim)) %>%
    dplyr::arrange(trial, interim)
}

#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `run_automatic_trial` as a list
#' @param final Return data from final analysis or interims
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom purrr reduce
tibble_ctr_quantities <- function(dat, final = TRUE, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["ctr_quantities"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "contrast", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "contrast")) %>%
      dplyr::mutate(contrast = forcats::fct_inorder(contrast)) %>%
      dplyr::arrange(interim, contrast)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial), interim = as.numeric(interim)) %>%
    dplyr::arrange(trial, interim)
}


#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `run_automatic_trial` as a list
#' @param final Return data from final analysis or interims
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom purrr reduce
tibble_par1_quantities <- function(dat, final = TRUE, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["par1_quantities"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "parameter", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "parameter")) %>%
      dplyr::mutate(parameter = forcats::fct_inorder(parameter)) %>%
      dplyr::arrange(interim, parameter)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial), interim = as.numeric(interim)) %>%
    dplyr::arrange(trial, interim)
}


#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `run_automatic_trial` as a list
#' @param final Return data from final analysis or interims
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom purrr reduce
tibble_par2_quantities <- function(dat, final = TRUE, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["par2_quantities"]]
    lapply(1:length(tq),
           function(x) {
             tidyr::gather(tidyr::as_tibble(tq[[x]], rownames = "interim"), "parameter", !!names(tq)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "parameter")) %>%
      dplyr::mutate(parameter = forcats::fct_inorder(parameter)) %>%
      dplyr::arrange(interim, parameter)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial), interim = as.numeric(interim)) %>%
    dplyr::arrange(trial, interim)
}


#' Group list of decision outcomes into a tibble
#'
#' @param dat The results of `ascot_trial2` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_dec1_quantities <- function(dat, ...) {
  dplyr::bind_rows(lapply(dat, function(i) {
    tq <- i[["dec1_quantities"]]
    tibble::enframe(tq) %>%
      tidyr::unnest_wider(value) %>%
      tidyr::gather(arm, value, -name) %>%
      tidyr::spread(name, value)
  }), .id = "trial")
}


#' Group list of decision outcomes into a tibble
#'
#' @param dat The results of `ascot_trial2` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_dec2_quantities <- function(dat, ...) {
  dplyr::bind_rows(lapply(dat, function(i) {
    tq <- i[["dec2_quantities"]]
    tibble::enframe(tq) %>%
      tidyr::unnest_wider(value) %>%
      tidyr::gather(contrast, value, -name) %>%
      tidyr::spread(name, value)
  }), .id = "trial")
}


#' Group list of trial outcomes into a tibble
#'
#' @param dat The results of `ascot_trial2` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_result_quantities <- function(dat, ...) {
  tibble::enframe(lapply(dat, function(x) x[["result_quantities"]]), name = "trial") %>%
    tidyr::unnest_wider(value)
}
