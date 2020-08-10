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
  cbind(0, 0, diag(1 + 1/3, 4) - matrix(1 / 3, 4, 4), matrix(0, 4, ncol(Xdes) - 6)),
  cbind(0, 0, 0, 0, 0, 0, diag(1 + 1/2, 3) - matrix(1 / 2, 3, 3), matrix(0, 3, ncol(Xdes) - 9))
)
colnames(Xmain) <- colnames(Xdes)
rownames(Xmain) <- colnames(Xdes)[3:9]


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
X_con <- as.matrix(cbind(1, rbind(0, cbind(1, Xmes_con, Xtim_con, Xarm_con))))
X_con_inv <- solve(X_con)
X_con_inv_t <- t(X_con_inv)
Q_t <- t(Q)
X_con_inv_t_Q_t <- X_con_inv_t %*% Q_t

# Prior
Sigma0 <- diag(c(10, 10, 10, 10, 10, 10, 10, rep(1, 6)))
Sigma0 <- diag(c(rep(1.5, 6), rep(1, 7)))

usethis::use_data(Xdes, Xmap, Xmain, X_con, X_con_inv, X_con_inv_t, Q, Q_t, X_con_inv_t_Q_t, Sigma0, overwrite = TRUE, internal = TRUE)
