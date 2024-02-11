

create_sbw_constraints <- function(X, target, tol) {

  # balance constraint
  A1 <- t(X)
  l1 <- nrow(X) * (target - tol) 
  u1 <- nrow(X) * (target + tol) 

  # sum to n constraint
  A2 <- rep(1, nrow(X))
  l2 <- nrow(X)
  u2 <- nrow(X)

  # positivity constraint
  A3 <- Matrix::Diagonal(nrow(X))
  l3 <- numeric(nrow(X))
  u3 <- rep(Inf, nrow(X))

  # return(list(A = A1,
  #             l = l1,
  #             u = u1))

  return(list(A = rbind(A1, A2, A3),
              l = c(l1, l2, l3),
              u = c(u1, u2, u3)))
}

#' SBW implmenetation with different solver
#' @param X Matrix of covariates (standardized!)
#' @param  Z Vector of treatment assignment
#' @param tol Tolerance for imbalance
#' @return List with weights for all units (treated weights are zero) and imbalance vector
sbw_osqp <- function(X, Z, tol, ...) {

  # compute target
  target <- colMeans(X[Z == 1, ])

  # P matrix
  n0 <- sum(Z == 0)
  P <- Matrix::sparseMatrix(1:n0, 1:n0, x = rep(1, n0))

  consts <- create_sbw_constraints(X[Z == 0,, drop = F], target, tol)

  pars <- do.call(osqp::osqpSettings, list(...))

  sol <- osqp::solve_osqp(P = P, A = consts$A, l = consts$l, u = consts$u, pars = pars)
  wts <- numeric(nrow(X))
  wts[Z == 0] <- sol$x

  imbal <- c(t(sol$x) %*% X[Z == 0, ] / sum(Z == 0)) - target

  return(list(wts = wts, imbalance = imbal))

}
