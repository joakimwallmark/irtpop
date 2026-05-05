#' Population difficulty of an IRT item
#'
#' Computes the population difficulty index
#' \deqn{B_j = \int \left[1 - \frac{E(X_j \mid \boldsymbol{\theta})}{\mathrm{max}_j}\right] f(\boldsymbol{\theta})\, d\boldsymbol{\theta}}
#' by numerical quadrature over a grid of \eqn{\boldsymbol{\theta}} values.
#' Higher values indicate harder items.
#'
#' @param item       An extracted mirt item object (from [mirt::extract.item()]).
#' @param max_score  Maximum possible item score. Defaults to `item@ncat - 1`.
#' @param theta_mean Numeric vector of latent trait means (length = number of
#'   dimensions). Defaults to the zero vector.
#' @param theta_cov  Numeric positive-definite covariance matrix for the latent
#'   traits. Defaults to the identity matrix.
#'
#' @return A single numeric value in \eqn{[0, 1]}.
#'
#' @examples
#' \dontrun{
#' library(mirt)
#' mod  <- mirt(Science, 1, itemtype = "graded", verbose = FALSE)
#' item <- extract.item(mod, 1)
#' population_difficulty(item)
#' }
#'
#' @seealso [population_discrimination()], [item_metrics()]
#' @export
population_difficulty <- function(item,
                                  max_score  = NULL,
                                  theta_mean = NULL,
                                  theta_cov  = NULL) {
  if (!isS4(item) || !methods::.hasSlot(item, "nfact") || !methods::.hasSlot(item, "ncat")) {
    stop("`item` must be an extracted mirt item object from `mirt::extract.item()`.")
  }

  D <- item@nfact

  if (!is.null(max_score)) {
    if (!is.numeric(max_score) || length(max_score) != 1L || is.na(max_score) || max_score <= 0) {
      stop("`max_score` must be a single positive number.")
    }
  }

  if (!is.null(theta_mean)) {
    if (!is.numeric(theta_mean) || anyNA(theta_mean) || length(theta_mean) != D) {
      stop(sprintf("`theta_mean` must be a numeric vector of length %d (number of latent dimensions).",
                   D))
    }
  }

  if (!is.null(theta_cov)) {
    if (!is.numeric(theta_cov) || !is.matrix(theta_cov) || anyNA(theta_cov) ||
        nrow(theta_cov) != D || ncol(theta_cov) != D) {
      stop(sprintf("`theta_cov` must be a numeric %d x %d covariance matrix.", D, D))
    }
  }

  if (is.null(max_score))  max_score  <- item@ncat - 1L
  if (is.null(theta_mean)) theta_mean <- rep(0, D)
  if (is.null(theta_cov))  theta_cov  <- diag(D)

  n_pts  <- if (D == 1L) 100L else 40L
  g      <- make_grid(D, n_pts)
  dens   <- grid_density(g$grid, theta_mean, theta_cov)
  scores <- mirt::expected.item(item, g$grid)
  sum((1 - scores / max_score) * dens) * (g$dx^D)
}


#' Population discrimination of an IRT item
#'
#' Computes the population discrimination by integrating the gradient of the
#' expected item score with respect to \eqn{\boldsymbol{\theta}} over the
#' population distribution.
#'
#' Two forms are available:
#'
#' * **`"vector"`** — the D-dimensional vector
#'   \deqn{\mathbf{A}_j = \int \frac{\nabla_{\boldsymbol{\theta}} E(X_j \mid \boldsymbol{\theta})}{\mathrm{max}_j} f(\boldsymbol{\theta})\, d\boldsymbol{\theta}}
#'
#' * **`"scalar"`** — the covariance-weighted norm
#'   \deqn{A_j = \int \frac{\sqrt{\nabla_{\boldsymbol{\theta}} E(X_j \mid \boldsymbol{\theta})^\top \boldsymbol{\Sigma}_\theta\, \nabla_{\boldsymbol{\theta}} E(X_j \mid \boldsymbol{\theta})}}{\mathrm{max}_j} f(\boldsymbol{\theta})\, d\boldsymbol{\theta}}
#'   which reduces to the absolute slope for unidimensional models and is
#'   rotation-invariant for multidimensional models.
#'
#' The gradient is computed via `mirt:::DerivTheta`, which provides the analytic
#' derivative of each category probability w.r.t. \eqn{\boldsymbol{\theta}}.
#'
#' @param item       An extracted mirt item object (from [mirt::extract.item()]).
#' @param max_score  Maximum possible item score. Defaults to `item@ncat - 1`.
#' @param type       Either `"scalar"` (default) or `"vector"`.
#' @param theta_mean Numeric vector of latent trait means. Defaults to zeros.
#' @param theta_cov  Numeric covariance matrix for the latent traits. Defaults
#'   to the identity matrix.
#'
#' @return For `type = "scalar"`: a single non-negative numeric value.
#'   For `type = "vector"`: a numeric vector of length D.
#'
#' @examples
#' \dontrun{
#' library(mirt)
#' mod  <- mirt(Science, 1, itemtype = "graded", verbose = FALSE)
#' item <- extract.item(mod, 1)
#' population_discrimination(item)
#' population_discrimination(item, type = "vector")
#' }
#'
#' @seealso [population_difficulty()], [item_metrics()]
#' @export
population_discrimination <- function(item,
                                      max_score  = NULL,
                                      type       = "scalar",
                                      theta_mean = NULL,
                                      theta_cov  = NULL) {
  type <- match.arg(type, c("scalar", "vector"))

  if (!isS4(item) || !methods::.hasSlot(item, "nfact") || !methods::.hasSlot(item, "ncat")) {
    stop("`item` must be an extracted mirt item object from `mirt::extract.item()`.")
  }

  D <- item@nfact

  if (!is.null(max_score)) {
    if (!is.numeric(max_score) || length(max_score) != 1L || is.na(max_score) || max_score <= 0) {
      stop("`max_score` must be a single positive number.")
    }
  }

  if (!is.null(theta_mean)) {
    if (!is.numeric(theta_mean) || anyNA(theta_mean) || length(theta_mean) != D) {
      stop(sprintf("`theta_mean` must be a numeric vector of length %d (number of latent dimensions).",
                   D))
    }
  }

  if (!is.null(theta_cov)) {
    if (!is.numeric(theta_cov) || !is.matrix(theta_cov) || anyNA(theta_cov) ||
        nrow(theta_cov) != D || ncol(theta_cov) != D) {
      stop(sprintf("`theta_cov` must be a numeric %d x %d covariance matrix.", D, D))
    }
  }

  if (is.null(max_score))  max_score  <- item@ncat - 1L
  if (is.null(theta_mean)) theta_mean <- rep(0, D)
  if (is.null(theta_cov))  theta_cov  <- diag(D)

  n_pts <- if (D == 1L) 100L else 40L
  g     <- make_grid(D, n_pts)
  dens  <- grid_density(g$grid, theta_mean, theta_cov)

  # Gradient of E[X|theta] w.r.t. theta via mirt's analytic DerivTheta.
  # grad[[k]] is an N x D matrix of dP(X = k-1 | theta) / d(theta).
  # dE/d(theta) = sum_k (k-1) * grad[[k]] (the min offset cancels: sum dP = 0).
  deriv <- mirt:::DerivTheta(item, g$grid)  # nolint: triple_colon_linter
  dE    <- matrix(0, nrow(g$grid), D)
  for (k in seq_len(item@ncat)) {
    dE <- dE + (k - 1L) * deriv$grad[[k]]
  }

  if (type == "scalar") {
    quad_forms <- rowSums((dE %*% theta_cov) * dE)
    val <- (sqrt(quad_forms) / max_score) * dens
    return(sum(val) * (g$dx^D))
  }

  # type == "vector"
  vapply(
    seq_len(D),
    function(k) sum((dE[, k] / max_score) * dens) * (g$dx^D),
    numeric(1L)
  )
}
