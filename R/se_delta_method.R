#' Delta-method SE for population difficulty (M2PL, any D)
#'
#' Computes the standard error of the population difficulty estimator \eqn{B_j}
#' via the multivariate Delta method. The gradient with respect to the item
#' parameter vector \eqn{\boldsymbol{\xi}_j = (a_1, \ldots, a_D, d)} is
#' \deqn{\nabla_{\boldsymbol{\xi}_j} B_j = (-\boldsymbol{\mu}, -\gamma)}
#' where \eqn{\gamma = \int P Q f\, d\boldsymbol{\theta}} and
#' \eqn{\boldsymbol{\mu}} is the D-vector
#' \eqn{\int P Q\, \boldsymbol{\theta}\, f\, d\boldsymbol{\theta}}.
#'
#' @param a          Numeric discrimination vector of length D.
#' @param d          Scalar intercept parameter.
#' @param vcov_j     \eqn{(D+1) \times (D+1)} covariance matrix for
#'   \eqn{(a_1, \ldots, a_D, d)}.
#' @param theta_mean Latent mean vector (length D). Defaults to zeros.
#' @param theta_cov  Latent covariance matrix (D x D). Defaults to identity.
#'
#' @return A non-negative scalar standard error.
#'
#' @seealso [se_pop_discrimination_scalar()], [se_pop_discrimination_vector()],
#'   [item_metrics()]
#' @export
se_pop_difficulty <- function(a, d, vcov_j,
                              theta_mean = NULL, theta_cov = NULL) {
  if (!is.numeric(a) || anyNA(a) || length(a) < 1L) {
    stop("`a` must be a non-empty numeric vector of discrimination parameters with no NA values.")
  }
  D <- length(a)
  if (!is.numeric(d) || length(d) != 1L || is.na(d)) {
    stop("`d` must be a single numeric intercept value.")
  }
  if (!is.numeric(vcov_j) || !is.matrix(vcov_j) || anyNA(vcov_j) ||
      nrow(vcov_j) != D + 1L || ncol(vcov_j) != D + 1L) {
    stop(sprintf("`vcov_j` must be a numeric %d x %d covariance matrix for (a_1, ..., a_D, d).",
                 D + 1L, D + 1L))
  }
  if (!is.null(theta_mean)) {
    if (!is.numeric(theta_mean) || anyNA(theta_mean) || length(theta_mean) != D) {
      stop(sprintf("`theta_mean` must be a numeric vector of length %d.", D))
    }
  }
  if (!is.null(theta_cov)) {
    if (!is.numeric(theta_cov) || !is.matrix(theta_cov) || anyNA(theta_cov) ||
        nrow(theta_cov) != D || ncol(theta_cov) != D) {
      stop(sprintf("`theta_cov` must be a numeric %d x %d covariance matrix.", D, D))
    }
  }
  ints <- aux_integrals(a, d, theta_mean, theta_cov)
  grad <- c(-ints$mu_vec, -ints$gamma)
  var_B <- as.numeric(crossprod(grad, vcov_j %*% grad))
  sqrt(max(var_B, 0))
}


#' Delta-method SEs for vector population discrimination (M2PL, any D)
#'
#' Computes component-wise standard errors for the vector population
#' discrimination \eqn{\mathbf{A}_j}. For M2PL,
#' \eqn{\mathbf{A}_j = \gamma_j \mathbf{a}_j} and the \eqn{D \times (D+1)}
#' Jacobian is
#' \deqn{J = \bigl[\mathbf{a} \otimes \boldsymbol{\lambda}^\top + \gamma\, I_D \;\big|\; \kappa\, \mathbf{a}\bigr].}
#'
#' @inheritParams se_pop_difficulty
#'
#' @return Numeric vector of length D: one SE per component of \eqn{\mathbf{A}_j}.
#'
#' @seealso [se_pop_difficulty()], [se_pop_discrimination_scalar()],
#'   [item_metrics()]
#' @export
se_pop_discrimination_vector <- function(a, d, vcov_j,
                                         theta_mean = NULL, theta_cov = NULL) {
  if (!is.numeric(a) || anyNA(a) || length(a) < 1L) {
    stop("`a` must be a non-empty numeric vector of discrimination parameters with no NA values.")
  }
  D <- length(a)
  if (!is.numeric(d) || length(d) != 1L || is.na(d)) {
    stop("`d` must be a single numeric intercept value.")
  }
  if (!is.numeric(vcov_j) || !is.matrix(vcov_j) || anyNA(vcov_j) ||
      nrow(vcov_j) != D + 1L || ncol(vcov_j) != D + 1L) {
    stop(sprintf("`vcov_j` must be a numeric %d x %d covariance matrix for (a_1, ..., a_D, d).",
                 D + 1L, D + 1L))
  }
  if (!is.null(theta_mean)) {
    if (!is.numeric(theta_mean) || anyNA(theta_mean) || length(theta_mean) != D) {
      stop(sprintf("`theta_mean` must be a numeric vector of length %d.", D))
    }
  }
  if (!is.null(theta_cov)) {
    if (!is.numeric(theta_cov) || !is.matrix(theta_cov) || anyNA(theta_cov) ||
        nrow(theta_cov) != D || ncol(theta_cov) != D) {
      stop(sprintf("`theta_cov` must be a numeric %d x %d covariance matrix.", D, D))
    }
  }
  ints <- aux_integrals(a, d, theta_mean, theta_cov)

  J      <- cbind(outer(a, ints$lambda_vec) + ints$gamma * diag(D),
                  ints$kappa * a)
  Sigma_A <- J %*% vcov_j %*% t(J)
  sqrt(pmax(diag(Sigma_A), 0))
}


#' Delta-method SE for scalar population discrimination (M2PL, any D)
#'
#' Computes the standard error of the scalar population discrimination
#' \eqn{A_j = \gamma_j \|\mathbf{a}_j\|_{\Sigma_\theta}} via the Delta method.
#' The gradient with respect to \eqn{\boldsymbol{\xi}_j = (\mathbf{a}, d)} is
#' \deqn{\nabla_a A_j = \frac{\|\mathbf{a}\|_\Sigma^2\,\boldsymbol{\lambda} + \gamma\,\boldsymbol{\Sigma}_\theta\,\mathbf{a}}{\|\mathbf{a}\|_\Sigma}, \qquad
#'       \nabla_d A_j = \|\mathbf{a}\|_\Sigma\,\kappa.}
#'
#' @inheritParams se_pop_difficulty
#'
#' @return A non-negative scalar SE, or `NA` if \eqn{\|\mathbf{a}\|_\Sigma < 10^{-10}}.
#'
#' @seealso [se_pop_difficulty()], [se_pop_discrimination_vector()],
#'   [item_metrics()]
#' @export
se_pop_discrimination_scalar <- function(a, d, vcov_j,
                                          theta_mean = NULL, theta_cov = NULL) {
  if (!is.numeric(a) || anyNA(a) || length(a) < 1L) {
    stop("`a` must be a non-empty numeric vector of discrimination parameters with no NA values.")
  }
  D <- length(a)
  if (!is.numeric(d) || length(d) != 1L || is.na(d)) {
    stop("`d` must be a single numeric intercept value.")
  }
  if (!is.numeric(vcov_j) || !is.matrix(vcov_j) || anyNA(vcov_j) ||
      nrow(vcov_j) != D + 1L || ncol(vcov_j) != D + 1L) {
    stop(sprintf("`vcov_j` must be a numeric %d x %d covariance matrix for (a_1, ..., a_D, d).",
                 D + 1L, D + 1L))
  }
  if (!is.null(theta_mean)) {
    if (!is.numeric(theta_mean) || anyNA(theta_mean) || length(theta_mean) != D) {
      stop(sprintf("`theta_mean` must be a numeric vector of length %d.", D))
    }
  }
  if (!is.null(theta_cov)) {
    if (!is.numeric(theta_cov) || !is.matrix(theta_cov) || anyNA(theta_cov) ||
        nrow(theta_cov) != D || ncol(theta_cov) != D) {
      stop(sprintf("`theta_cov` must be a numeric %d x %d covariance matrix.", D, D))
    }
  }
  if (is.null(theta_cov)) theta_cov <- diag(D)

  norm_a_sq <- as.numeric(crossprod(a, theta_cov %*% a))
  norm_a    <- sqrt(max(norm_a_sq, 0))
  if (norm_a < 1e-10) return(NA_real_)

  ints <- aux_integrals(a, d, theta_mean, theta_cov)

  grad_a <- (norm_a_sq * ints$lambda_vec + ints$gamma * (theta_cov %*% a)) /
              norm_a
  grad_d <- norm_a * ints$kappa
  grad   <- c(grad_a, grad_d)

  var_A <- as.numeric(crossprod(grad, vcov_j %*% grad))
  sqrt(max(var_A, 0))
}
