mes <- 4
tim <- 3
arm <- mes*tim

# Design Matrix
Xmes <- kronecker(diag(1, mes), rep(1, tim))
Xtim <- kronecker(rep(1, mes), diag(1, tim))
Xarm <- diag(1, arm)
X <- cbind(1, c(0, rep(1, arm)), rbind(0, Xmes), rbind(0, Xtim), rbind(0, Xarm))
colnames(X) <- c("ctr", "trt", 
                 paste0("mes", 1:mes), paste0("tim", 1:tim), 
                 paste0(rep(paste0("mes", 1:mes), each = tim), 
                        rep(paste0("tim", 1:tim), times = mes)))

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

usethis::use_data(X_con, X_con_inv, X_con_inv_t, Q, Q_t, X_con_inv_t_Q_t, overwrite = TRUE, internal = TRUE)
