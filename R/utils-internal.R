# Internal helpers — not exported

#' Build a quadrature grid over `[-6, 6]^D`
#'
#' @param D         Integer number of dimensions.
#' @param n_pts     Integer number of points per dimension.
#' @return List with components `grid` (matrix, nrow = `n_pts^D`, ncol = D) and
#'   `dx` (scalar step size).
#' @keywords internal
make_grid <- function(D, n_pts) {
  pts  <- seq(-6, 6, length.out = n_pts)
  grid <- as.matrix(expand.grid(replicate(D, pts, simplify = FALSE)))
  list(grid = grid, dx = 12 / (n_pts - 1L))
}

#' Evaluate the multivariate normal density on a grid
#'
#' @param grid       Numeric matrix (N x D).
#' @param theta_mean Numeric vector of length D.
#' @param theta_cov  Numeric D x D positive-definite matrix.
#' @return Numeric vector of length N.
#' @keywords internal
grid_density <- function(grid, theta_mean, theta_cov) {
  D <- ncol(grid)
  if (D == 1L) {
    stats::dnorm(grid[, 1L],
                 mean = theta_mean[1L],
                 sd   = sqrt(theta_cov[1L, 1L]))
  } else {
    mvtnorm::dmvnorm(grid, mean = theta_mean, sigma = theta_cov)
  }
}

#' Compute auxiliary integrals for M2PL Delta-method SEs (D-dimensional)
#'
#' Evaluates the four vector/scalar integrals (gamma, mu_vec, lambda_vec,
#' kappa) needed for the gradients of B_j, A_j (vector), and A_j (scalar)
#' w.r.t. item parameters `(a, d)` under a D-variate normal population.
#'
#' For D = 1 a 100-point 1-D quadrature is used; for D >= 2 a product grid
#' of 40 points per dimension is used.
#'
#' @param a          Numeric vector of length D: discrimination parameters.
#' @param d          Scalar intercept parameter.
#' @param theta_mean Latent mean vector (length D, default zeros).
#' @param theta_cov  Latent covariance matrix (D x D, default identity).
#' @param n_pts      Quadrature points per dimension (auto-selected if `NULL`).
#' @return Named list:
#'   * `gamma`      — scalar  integral of P * Q * f
#'   * `mu_vec`     — D-vector of integrals of P * Q * theta_k * f
#'   * `lambda_vec` — D-vector of integrals of P * Q * (1 - 2P) * theta_k * f
#'   * `kappa`      — scalar  integral of P * Q * (1 - 2P) * f
#' @keywords internal
aux_integrals <- function(a, d,
                          theta_mean = NULL, theta_cov = NULL,
                          n_pts = NULL) {
  D <- length(a)
  if (is.null(theta_mean)) theta_mean <- rep(0, D)
  if (is.null(theta_cov))  theta_cov  <- diag(D)
  if (is.null(n_pts))      n_pts      <- if (D == 1L) 100L else 40L

  g    <- make_grid(D, n_pts)
  dens <- grid_density(g$grid, theta_mean, theta_cov)

  eta <- as.numeric(g$grid %*% a) + d
  pp  <- stats::plogis(eta)
  pq  <- pp * (1 - pp)

  w     <- dens * (g$dx^D)
  pq_w  <- pq * w
  pql_w <- pq * (1 - 2 * pp) * w

  list(
    gamma      = sum(pq_w),
    mu_vec     = colSums(pq_w  * g$grid),
    lambda_vec = colSums(pql_w * g$grid),
    kappa      = sum(pql_w)
  )
}
