#' Population item metrics for a fitted mirt model
#'
#' A convenience wrapper that computes [population_difficulty()] and
#' [population_discrimination()] for every item (or a subset) in a fitted
#' [mirt::mirt()] model and returns the results as a `data.frame`.
#'
#' The latent distribution parameters (`theta_mean`, `theta_cov`) are taken
#' from the model's estimated group parameters unless overridden by the caller.
#'
#' @param mod        A fitted object of class `"SingleGroupClass"` or
#'   `"MultipleGroupClass"` as returned by [mirt::mirt()].
#' @param items      Integer vector of item indices to evaluate. Defaults to all
#'   items.
#' @param max_score  Named integer vector or single integer giving the maximum
#'   possible score for each item. If a single value is supplied it is recycled.
#'   Defaults to `item@ncat - 1` for each item.
#' @param theta_mean Numeric vector overriding the latent mean. Defaults to the
#'   model's estimated group means.
#' @param theta_cov  Numeric matrix overriding the latent covariance. Defaults
#'   to the model's estimated group covariance.
#' @param se         Logical. If `TRUE`, Delta-method standard errors are
#'   appended for each metric. SEs are only supported for 2PL items
#'   (`itemtype = "2PL"` or `"Rasch"`); other item types produce `NA` SE
#'   columns. Requires the model to have been fitted with `SE = TRUE`.
#'
#' @return A `data.frame` with columns:
#'   * `item` — item index (integer)
#'   * `item_name` — item name from the model
#'   * `B` — population difficulty
#'   * `A` — scalar population discrimination
#'   * `A_1, A_2, ...` — vector population discrimination components (one per
#'     latent dimension)
#'
#'   When `se = TRUE`, additional columns `B_se`, `A_se`, `A_1_se`, `A_2_se`,
#'   ... are appended.
#'
#' @examples
#' \dontrun{
#' library(mirt)
#' mod <- mirt(Science, 1, itemtype = "graded", verbose = FALSE)
#' item_metrics(mod)
#'
#' # With Delta-method SEs (requires SE = TRUE during fitting)
#' data <- expand.table(LSAT7)
#' mod_se <- mirt(data, 1, itemtype = "2PL", verbose = FALSE, SE = TRUE)
#' item_metrics(mod_se, se = TRUE)
#' }
#'
#' @seealso [population_difficulty()], [population_discrimination()],
#'   [se_pop_difficulty()], [se_pop_discrimination_scalar()],
#'   [se_pop_discrimination_vector()]
#' @importFrom methods .hasSlot
#' @importFrom stats4 vcov
#' @export
item_metrics <- function(mod,
                         items      = NULL,
                         max_score  = NULL,
                         theta_mean = NULL,
                         theta_cov  = NULL,
                         se         = FALSE) {

  if (!inherits(mod, c("SingleGroupClass", "MultipleGroupClass"))) {
    stop("`mod` must be a fitted mirt model of class 'SingleGroupClass' or ",
         "'MultipleGroupClass', as returned by `mirt::mirt()`.")
  }

  n_items <- mod@Data$nitems
  D       <- mod@Model$nfact

  if (is.null(items)) {
    items <- seq_len(n_items)
  } else {
    if (!is.numeric(items) || anyNA(items)) {
      stop("`items` must be a numeric vector of item indices with no NA values.")
    }
    if (any(items < 1) || any(items > n_items)) {
      stop(sprintf("`items` must contain indices between 1 and %d (the number of items in the model).",
                   n_items))
    }
  }
  items <- as.integer(items)

  if (!is.null(max_score)) {
    if (!is.numeric(max_score) || anyNA(max_score) || any(max_score <= 0)) {
      stop("`max_score` must be a positive numeric vector with no NA values.")
    }
    if (length(max_score) > 1L && length(max_score) != n_items) {
      stop(sprintf("`max_score` must be either length 1 or length %d (number of items).",
                   n_items))
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

  if (!is.logical(se) || length(se) != 1L || is.na(se)) {
    stop("`se` must be a single logical value (TRUE or FALSE).")
  }

  # Latent distribution from the model unless overridden
  if (is.null(theta_mean) || is.null(theta_cov)) {
    cc <- mirt::coef(mod, simplify = TRUE)
    if (is.null(theta_mean)) theta_mean <- as.numeric(cc$means)
    if (is.null(theta_cov))  theta_cov  <- as.matrix(cc$cov)
  }

  # Optionally retrieve full vcov for SE
  full_vcov     <- NULL
  vcov_parnums  <- NULL
  if (se) {
    if (!methods::.hasSlot(mod, "vcov") || nrow(mod@vcov) <= 1L) {
      warning("Model has no vcov matrix. Refit with SE = TRUE. Returning NA for all SEs.",
              immediate. = TRUE)
    } else {
      full_vcov    <- stats4::vcov(mod)
      vcov_parnums <- as.integer(sub(".*\\.", "", colnames(full_vcov)))
    }
  }

  rows <- vector("list", length(items))

  # Collect item names for SE warnings so we can emit a single message per cause.
  unsupported_type_items   <- character(0L)
  constrained_slope_items  <- character(0L)
  vcov_missing_items       <- character(0L)

  for (i_idx in seq_along(items)) {
    j    <- items[i_idx]
    item <- mirt::extract.item(mod, j)
    item_label <- colnames(mod@Data$data)[j]

    ms <- if (!is.null(max_score)) {
      if (length(max_score) == 1L) max_score else max_score[j]
    } else {
      item@ncat - 1L
    }

    B      <- population_difficulty(item, ms, theta_mean, theta_cov)
    A_sc   <- population_discrimination(item, ms, "scalar", theta_mean, theta_cov)
    A_vec  <- population_discrimination(item, ms, "vector", theta_mean, theta_cov)

    row <- c(list(item = j,
                  item_name = item_label),
             list(B = B, A = A_sc),
             stats::setNames(as.list(A_vec),
                             paste0("A_", seq_len(D))))

    if (se) {
      is_dich     <- inherits(item, "dich")
      has_guessing <- is_dich && item@par[D + 2L] > 0  # 3PL or similar
      all_slopes_free <- is_dich && !has_guessing &&
                           all(item@est[seq_len(D + 1L)])

      if (!is_dich || has_guessing) {
        unsupported_type_items <- c(unsupported_type_items,
                                    sprintf("%s (class: %s)",
                                            item_label, class(item)))
        se_row <- stats::setNames(
          as.list(rep(NA_real_, 2L + D)),
          c("B_se", "A_se", paste0("A_", seq_len(D), "_se"))
        )
      } else if (!all_slopes_free) {
        constrained_slope_items <- c(constrained_slope_items, item_label)
        se_row <- stats::setNames(
          as.list(rep(NA_real_, 2L + D)),
          c("B_se", "A_se", paste0("A_", seq_len(D), "_se"))
        )
      } else if (is.null(full_vcov)) {
        se_row <- stats::setNames(
          as.list(rep(NA_real_, 2L + D)),
          c("B_se", "A_se", paste0("A_", seq_len(D), "_se"))
        )
      } else {
        free_parnums <- item@parnum[item@est]
        vcov_pos     <- match(free_parnums, vcov_parnums)

        if (anyNA(vcov_pos)) {
          vcov_missing_items <- c(vcov_missing_items, item_label)
          se_row <- stats::setNames(
            as.list(rep(NA_real_, 2L + D)),
            c("B_se", "A_se", paste0("A_", seq_len(D), "_se"))
          )
        } else {
          vcov_j <- full_vcov[vcov_pos, vcov_pos, drop = FALSE]
          a_pars <- item@par[seq_len(D)]
          d_par  <- item@par[D + 1L]

          B_se   <- se_pop_difficulty(a_pars, d_par, vcov_j, theta_mean, theta_cov)
          A_se   <- se_pop_discrimination_scalar(a_pars, d_par, vcov_j, theta_mean, theta_cov)
          Av_se  <- se_pop_discrimination_vector(a_pars, d_par, vcov_j, theta_mean, theta_cov)

          se_row <- c(list(B_se = B_se, A_se = A_se),
                      stats::setNames(as.list(Av_se),
                                      paste0("A_", seq_len(D), "_se")))
        }
      }
      row <- c(row, se_row)
    }

    rows[[i_idx]] <- row
  }

  if (length(unsupported_type_items) > 0L) {
    warning(
      "Delta-method SEs are only implemented for 2PL/Rasch (dichotomous, no ",
      "guessing parameter) items. SEs set to NA for: ",
      paste(unsupported_type_items, collapse = ", "),
      ".",
      call. = FALSE, immediate. = TRUE
    )
  }

  if (length(constrained_slope_items) > 0L) {
    warning(
      "Delta-method SEs require all slope and intercept parameters to be ",
      "freely estimated, but one or more items have constrained (fixed) slope ",
      "parameters (e.g., due to identification constraints in an exploratory ",
      "model). SEs set to NA for: ",
      paste(constrained_slope_items, collapse = ", "),
      ".",
      call. = FALSE, immediate. = TRUE
    )
  }

  if (length(vcov_missing_items) > 0L) {
    warning(
      "Item parameters could not be matched in the model vcov matrix. ",
      "SEs set to NA for: ",
      paste(vcov_missing_items, collapse = ", "),
      ". Consider refitting the model with SE = TRUE.",
      call. = FALSE, immediate. = TRUE
    )
  }

  do.call(rbind, lapply(rows, as.data.frame))
}
