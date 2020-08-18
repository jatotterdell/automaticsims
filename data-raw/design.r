mes <- 4
tim <- 3
arm <- mes*tim

# Design Matrix
Xmes <- kronecker(diag(1, mes), rep(1, tim))
colnames(Xmes) <- paste0("m", 1:4)
Xtim <- kronecker(rep(1, mes), diag(1, tim))
colnames(Xtim) <- paste0("t", 1:3)
Xarm <- diag(1, arm)
Xdes <- cbind(1, c(0, rep(1, arm)), rbind(0, Xmes), rbind(0, Xtim), rbind(0, Xarm))
colnames(Xdes) <- c("ctr", "trt", 
                 paste0("mes", 1:mes), paste0("tim", 1:tim), 
                 paste0(rep(paste0("mes", 1:mes), each = tim), 
                        rep(paste0("tim", 1:tim), times = mes)))

# Indicator for which combination involves which message/timing
Xmap <- cbind(0, rbind(0, cbind(Xmes, Xtim)))
colnames(Xmap)[1] <- "m0t0"
Xmap[1, 1] <- 1
rownames(Xmap) <- sapply(1:nrow(Xmap), function(a) paste(colnames(Xmap)[Xmap[a, ] == 1], collapse = ''))

# Main effect contrasts in terms of Xdes
Xmain <- rbind(
  c(0, 1, rep(0, 19)),
  cbind(0, 0, diag(1 + 1/3, 4) - matrix(1 / 3, 4, 4), matrix(0, 4, ncol(Xdes) - 6)),
  cbind(0, 0, 0, 0, 0, 0, diag(1 + 1/2, 3) - matrix(1 / 2, 3, 3), matrix(0, 3, ncol(Xdes) - 9)),
  c(0, 1, 1, rep(0, 18)),
  c(0, 1, 0, 1, rep(0, 17)),
  c(0, 1, 0, 0, 1, rep(0, 16)),
  c(0, 1, 0, 0, 0, 1, rep(0, 15)),
  c(0, 1, 0, 0, 0, 0, 1, rep(0, 14)),
  c(0, 1, 0, 0, 0, 0, 0, 1, rep(0, 13)),
  c(0, 1, 0, 0, 0, 0, 0, 0, 1, rep(0, 12))
)
colnames(Xmain) <- colnames(Xdes)
rownames(Xmain) <- c("b_trt", colnames(Xdes)[3:9], 
                     "m1_v_ctr", "m2_v_ctr", "m3_v_ctr", "m4_v_ctr", 
                     "t1_v_ctr", "t2_v_ctr", "t3_v_ctr")


# Constrained Design
Smes <- diag(1, mes) - 1/mes
Stim <- diag(1, tim) - 1/tim
Qmes <- eigen(Smes)$vector[, -mes]
Qtim <- eigen(Stim)$vector[, -tim]
Qarm <- kronecker(Qmes, Qtim)
Q <- as.matrix(Matrix::bdiag(1, 1, Qmes, Qtim, Qarm))
Xmes_con <- Xmes %*% Qmes
Xtim_con <- Xtim %*% Qtim
Xarm_con <- Xarm %*% Qarm

X_con <- Xdes %*% Q
X_con_inv <- solve(X_con)
X_con_inv_t <- t(X_con_inv)
Q_t <- t(Q)
X_con_inv_t_Q_t <- X_con_inv_t %*% Q_t
pinv_X_con <- MASS::ginv(X_con)

X_con_alt <- X_con
X_con_alt[, 1] <- c(1, rep(0, 12))
X_con_alt_inv <- solve(X_con_alt)
X_con_alt_inv_t <- t(X_con_alt_inv)

# Treatment coding
Q_trt <- as.matrix(Matrix::bdiag(1, 1, 
                                 contr.treatment(4), contr.treatment(3), 
                                 kronecker(contr.treatment(4), contr.treatment(3))))
X_trt <- Xdes %*% Q_trt

# Polynomial coding
Q_pol <- as.matrix(Matrix::bdiag(1, 1, 
                                 contr.poly(4), contr.poly(3), 
                                 kronecker(contr.poly(4), contr.poly(3))))
X_pol <- Xdes %*% Q_pol

# Indicator design
X_ind <- diag(1, 13)
C_ind <- matrix(
  c(
    c(1, rep(0, 12)),
    c(0, rep(1/12,12)),
    c(0, rep(1/3, 3), rep(0, 9)),
    c(rep(0, 4), rep(1/3, 3), rep(0, 6)),
    c(rep(0, 7), rep(1/3, 3), rep(0, 3)),
    c(rep(0, 10), rep(1/3, 3)),
    c(0, 1/4, 0, 0, 1/4, 0, 0, 1/4, 0, 0, 1/4, 0, 0),
    c(0, 0, 1/4, 0, 0, 1/4, 0, 0, 1/4, 0, 0, 1/4, 0),
    c(0, 0, 0, 1/4, 0, 0, 1/4, 0, 0, 1/4, 0, 0, 1/4)
  ),
  9, 13, byrow = T,
  dimnames = list("average" = c("b_ctr", "b_trt", "b_m1", "b_m2", "b_m3", "b_m4", "b_t1", "b_t2", "b_t3"))
)

# Prior
Sigma0 <- diag(c(10, 10, 10, 10, 10, 10, 10, rep(1, 6)))
Sigma0 <- diag(c(rep(1.5, 6), rep(1, 7)))

usethis::use_data(
  Xdes, Xmap, Xmain, 
  X_con, X_con_inv, X_con_inv_t, 
  X_con_alt, X_con_alt_inv, X_con_alt_inv_t,
  pinv_X_con,
  X_trt, Q_trt,
  X_pol, Q_pol,
  Q, Q_t, X_con_inv_t_Q_t, 
  
  Sigma0, 
  overwrite = TRUE, internal = TRUE)
