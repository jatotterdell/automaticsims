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
#' @export
run_automatic_trial <- function(
  mu,
  Sigma0 = diag(c(2^2, rep(1, 12))),
  n_seq = seq(500, 10000, 500),
  delta_sup = log(1.1),
  thres_sup = 0.85,
  thres_inf = 1 - thres_sup,
  thres_fut = 0.95,
  thres_nin = 0.90,
  thres_eff = 0.99,
  thres_equ = 0.90,
  thres_ina = 0.99,
  rar_scale = 0.5,
  use_mwu = TRUE,
  mc_draws = 1e4,
  ct_alloc = 0.25
) {
  
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
  
  n_max <- max(n_seq)
  n_int <- length(n_seq)
  n_new <- diff(c(0, n_seq))
  n_arms <- nrow(Xmap)
  arm_labs <- rownames(Xmap)
  mes_labs <- paste0("m", 1:4)
  tim_labs <- paste0("t", 1:3)
  beta_labs <- c("b_ctr", "b_trt", paste0("b_m", 1:4), paste0("b_t", 1:3), 
                 paste0(rep(paste0("b_m", 1:4), each = 3), rep(paste0("t", 1:3), times = 4)))
  dn_all <- list("interim" = 1:n_int, "arm" = arm_labs)
  dn_trt <- list("interim" = 1:n_int, "arm" = arm_labs[-1])
  
  # Data
  n <- matrix(0, n_int, n_arms, dimnames = dn_all)
  y <- matrix(0, n_int, n_arms, dimnames = dn_all)
  
  # Parameters
  arm_mean  <- matrix(0, n_int, n_arms, dimnames = dn_all)
  arm_var   <- matrix(0, n_int, n_arms, dimnames = dn_all)
  beta_mean <- matrix(0, n_int, 21, dimnames = list("interim" = 1:n_int, "par" = beta_labs))
  beta_var  <- matrix(0, n_int, 21, dimnames = list("interim" = 1:n_int, "par" = beta_labs))
  me_mean   <- matrix(0, n_int, 7, dimnames = list("interim" = 1:n_int, "treatment" = colnames(Xmap)[-1]))
  me_var    <- matrix(0, n_int, 7, dimnames = list("interim" = 1:n_int, "treatment" = colnames(Xmap)[-1]))
  
  # Store probability quantities
  p_rand    <- matrix(1/n_arms, n_int + 1, n_arms, dimnames = list("interim" = 0:n_int, "arm" = arm_labs))
  p_rand[, 1] <- ct_alloc
  p_rand[, -1] <- (1 - ct_alloc) / (n_arms - 1)
  
  p_sup     <- matrix(0, n_int, n_arms, dimnames = dn_all)
  p_inf     <- matrix(0, n_int, n_arms, dimnames = dn_all)
  p_sup_trt <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  p_fut     <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  p_eff     <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  p_equ     <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  p_ina     <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  # p_nin   <- matrix(0, n_int, n_arms - 1, dimnames = dn_trt)
  # p_best_mes <- matrix(0, n_int, 4, dimnames = list("interim" = 1:n_int, "message" = mes_labs))
  # p_best_tim <- matrix(0, n_int, 3, dimnames = list("interim" = 1:n_int, "timing" = tim_labs))
  # p_in_best  <- matrix(0, n_int, ncol(Xmap), dimnames = list("interim" = 1:n_int, "treatment" = colnames(Xmap)))
  
  # Store flags
  is_act <- matrix(TRUE, n_int, n_arms, dimnames = dn_all)      # Active
  is_sup <- matrix(TRUE, n_int, n_arms, dimnames = dn_all) 
  is_inf <- matrix(TRUE, n_int, n_arms, dimnames = dn_all) 
  is_sup_trt <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Superior (best)
  # is_inf <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Inferior (not best)
  # is_fut <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Futile (not best by enough)
  # is_nin <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Non-inferior (close enough to best)
  # is_ina <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Inadequate (not effective enough)
  is_eff <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Effective (better than control)
  # is_equ <- matrix(FALSE, n_int, n_arms - 1, dimnames = dn_trt) # Equivalent (same as control)
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
    
    # Fit model
    # mod <- vb_mod(y[i, ], n[i, ], Sigma0 = Sigma0)
    mod <- varapproxr::vb_logistic_n(
      X_con, y[i, ], n[i, ], rep(0, 13), Sigma0, 
      rep(0, 13), diag(1, 13), alg = "sj", maxiter_jj = 100)
    
    eta_mu     <- X_con %*% mod$mu
    eta_sigma  <- X_con %*% mod$Sigma %*% t(X_con)
    beta_mu    <- Q %*% mod$mu
    beta_sigma <- Q %*% mod$Sigma %*% t(Q)
    
    arm_mean[i, ]  <- drop(eta_mu)
    arm_var[i, ]   <- diag(eta_sigma)
    beta_mean[i, ] <- drop(beta_mu)
    beta_var[i, ]  <- diag(beta_sigma)
    me_mean[i, ]   <- drop(Xmain %*% beta_mu)
    me_var[i, ]    <- diag(Xmain %*% beta_sigma %*% t(Xmain))
    
    # Compute posterior draws and quantities
    draws      <- mvnfast::rmvn(mc_draws, eta_mu, sigma = eta_sigma)
    # beta_draws <- draws %*% X_con_inv_t_Q_t
    
    p_sup[i, ]     <- prob_max(draws)
    p_sup_trt[i, ] <- prob_max(draws[, -1])
    p_fut[i, ]     <- 1 - prob_each_superior_all(draws[, -1], delta_sup)
    best <- which.max(p_sup[i, ])

    # Relative to control
    diff_ctrl_draws <- sweep(draws[, -1], 1, draws[, 1])
    p_eff[i, ] <- colMeans(apply(diff_ctrl_draws, 2, function(x) x > 0))
    p_ina[i, ] <- 1 - colMeans(apply(diff_ctrl_draws, 2, function(x) x > delta_sup))
    p_equ[i, ] <- colMeans(apply(diff_ctrl_draws, 2, function(x) abs(x) < delta_sup))
    
    # Relative to current best
    # diff_best_draws <- sweep(draws[, -1], 1, draws[, best + 1])
    # p_in_best[i, ]  <- matrixStats::colMeans2((matrixStats::rowRanks(draws) == 13) %*% Xmap)
    # p_nin[i, ]      <- colMeans(apply(diff_best_draws, 2, function(x) x >= -delta_sup))
    # p_best_mes[i, ] <- prob_max(beta_draws[, 3:6])
    # p_best_tim[i, ] <- prob_max(beta_draws[, 7:9])
    
    # Check decisions
    is_sup[i, ]     <- p_sup[i, ] > thres_sup
    is_sup_trt[i, ] <- p_sup_trt[i, ] > thres_sup
    is_inf[i, ]     <- p_sup[i, ] < thres_inf / (13 - !is_act[i, 1])
    is_eff[i, ]     <- p_eff[i, ] > thres_eff
    is_act[i, ]     <- !is_inf[i, ]
    
    if(i > 1) {
      # Drop control permanently
      is_act[i, 1] <- is_act[i, 1] * is_act[i - 1, 1]
    }
    if(!is_act[i, 1]) p_rand[i + 1, 1] <- 0

    # Calculate active dependent quantities
    
    # Update allocations
    w <- (p_sup[i, -1] * arm_var[i, -1] / (n[i, -1] + 1))^rar_scale
    w[!is_act[i, -1]] <- 0
    if(all(!is_act[i, ])) {
      p_rand[i + 1, 1] <- 1
      p_rand[i + 1, -1] <- 0
    } else {
      p_rand[i + 1, -1] <- (1 - p_rand[i + 1, 1]) * w / sum(w)  
    }
    
    # Assess stopping criteria
    # If 
    #  - any is best and effective
    #  - all are inactive or all are futile
    #  - all active are non-inferior and best is effective
    any_best <- any(is_sup[i, ])
    # any_effective_best <- any(is_sup[i, ] & is_eff[i, is_crb[i, ]])
    all_inactive       <- !any(is_act[i, ])
    # all_futile         <- all(is_fut[i, ])
    # all_active_noninf  <- all(is_nin[i, is_act[i, ]] & is_eff[i, is_crb[i, ]])
    stopped <- any_best | all_inactive
    
    if(stopped) break
  }
  p_rand <- p_rand[-1, ]
  
  idx <- 1:i
  
  final <- list(
    
    result_quantities = list(
      stop_early = stopped,
      stop_at    = ifelse(stopped, i, NA),
      any_best = any_best,
      all_inactive = all_inactive,
      any_sup    = any(is_sup[i, ]),
      all_fut    = all(is_fut[i, ]),
      none_act   = !any(is_act[i, ]),
      true_sup   = any(is_sup[i, which_sup]),
      false_sup  = any(is_sup[i, which_inf]),
      false_inf  = any(is_inf[i, which_sup]),
      any_eff    = any(is_eff[i, ]),
      any_false_eff  = any(is_eff[i, which_ineff]),
      any_true_eff   = any(is_eff[i, which_eff]),
      all_true_eff   = ifelse(num_eff > 0, all(is_eff[i, which_eff]), FALSE)
      # any_ina = any(is_ina[i, ]),
      # any_false_ina = any(is_ina[i, which_ada])
    ),
    
    arm_quantities = list(
      n        = n[idx, , drop = F],
      y        = y[idx, , drop = F],
      p_rand   = p_rand[idx, , drop = F],
      arm_mean = arm_mean[idx, , drop = F], 
      arm_var  = arm_var[idx, , drop = F],
      p_sup  = p_sup[idx, , drop = F],
      p_sup_trt = p_sup_trt[idx, , drop = F],
      p_eff = p_eff[idx, , drop = F],
      is_sup = is_sup[idx, , drop = F],
      is_inf = is_inf[idx, , drop = F],
      is_eff = is_eff[idx, , drop = F],
      is_act = is_act[idx, , drop = F]
    ),
    
    treat_quantities = list(
      me_mean = me_mean[idx, , drop = F],
      me_var = me_var[idx, , drop = F]
    ),
    
    par_quantities = list(
      beta_mean = beta_mean[idx, , drop = F],
      beta_var  = beta_var[idx, , drop = F]
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
tibble_par_quantities <- function(dat, final = TRUE, ...) {
  dplyr::bind_rows(parallel::mclapply(dat, function(i) {
    tq <- i[["par_quantities"]]
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
#' @param dat The results of `ascot_trial2` as a list
#' @param ... Other arguments to `mclapply`
#' @export
#'
#' @importFrom dplyr %>%
tibble_result_quantities <- function(dat, ...) {
  tibble::enframe(lapply(dat, function(x) x[["result_quantities"]]), name = "trial") %>%
    tidyr::unnest_wider(value)
}
