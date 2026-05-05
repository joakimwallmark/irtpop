# Tests for item_metrics() — the high-level wrapper.
# Basic smoke tests live in test-population_metrics.R; this file focuses on
# structure, value consistency, argument handling, and SE behaviour.

# ── helpers ───────────────────────────────────────────────────────────────────

# Graded (non-2PL) model used only in this file.
suppressPackageStartupMessages(library(mirt))
set.seed(42)
.dat_gr <- mirt::simdata(
  a     = matrix(rep(1.2, 4)),
  d     = matrix(cbind(1.5, 0.5, -0.5, -1.5), nrow = 4),
  N     = 500,
  itemtype = "graded"
)
mod_graded <- mirt::mirt(.dat_gr, 1, itemtype = "graded", verbose = FALSE, SE = TRUE)

# ── structure ─────────────────────────────────────────────────────────────────

test_that("item_metrics 2D has A_1 and A_2 columns", {
  res <- item_metrics(mod_2d)
  expect_true(all(c("A_1", "A_2") %in% names(res)))
  expect_false("A_3" %in% names(res))
})

test_that("item_metrics 2D returns correct number of rows", {
  res <- item_metrics(mod_2d)
  expect_equal(nrow(res), mod_2d@Data$nitems)
})

test_that("item_metrics item column equals items argument", {
  res <- item_metrics(mod_1d, items = c(2L, 4L, 5L))
  expect_equal(res$item, c(2L, 4L, 5L))
})

test_that("item_metrics item_name column is character", {
  res <- item_metrics(mod_1d)
  expect_type(res$item_name, "character")
  expect_equal(res$item_name, colnames(mod_1d@Data$data))
})

# ── value consistency ─────────────────────────────────────────────────────────

test_that("item_metrics A column matches population_discrimination scalar (1D)", {
  j   <- 3L
  res <- item_metrics(mod_1d, items = j)
  itm <- mirt::extract.item(mod_1d, j)
  cc  <- mirt::coef(mod_1d, simplify = TRUE)
  A_fn <- population_discrimination(itm, max_score = 1L,
                                    type = "scalar",
                                    theta_mean = as.numeric(cc$means),
                                    theta_cov  = as.matrix(cc$cov))
  expect_equal(res$A[1L], A_fn, tolerance = 1e-9)
})

test_that("item_metrics A_1, A_2 match population_discrimination vector (2D)", {
  j   <- 2L
  res <- item_metrics(mod_2d, items = j)
  itm <- mirt::extract.item(mod_2d, j)
  cc  <- mirt::coef(mod_2d, simplify = TRUE)
  A_vec <- population_discrimination(itm, max_score = 1L,
                                     type = "vector",
                                     theta_mean = as.numeric(cc$means),
                                     theta_cov  = as.matrix(cc$cov))
  expect_equal(res$A_1[1L], A_vec[1L], tolerance = 1e-9)
  expect_equal(res$A_2[1L], A_vec[2L], tolerance = 1e-9)
})

test_that("item_metrics B values are in [0, 1] for 1D model", {
  res <- item_metrics(mod_1d)
  expect_true(all(res$B >= 0 & res$B <= 1))
})

test_that("item_metrics B values are in [0, 1] for 2D model", {
  res <- item_metrics(mod_2d)
  expect_true(all(res$B >= 0 & res$B <= 1))
})

# ── argument handling ─────────────────────────────────────────────────────────

test_that("item_metrics theta_mean/theta_cov override changes results", {
  res_default  <- item_metrics(mod_1d, items = 1L)
  res_override <- item_metrics(mod_1d, items = 1L,
                               theta_mean = 0.5, theta_cov = matrix(2))
  expect_false(isTRUE(all.equal(res_default$B,   res_override$B)))
  expect_false(isTRUE(all.equal(res_default$A,   res_override$A)))
})

test_that("item_metrics scalar max_score override changes B", {
  res_default  <- item_metrics(mod_1d, items = 1L)
  # supply max_score = 2 (double the true maximum) — B must differ
  res_ms2 <- item_metrics(mod_1d, items = 1L, max_score = 2L)
  expect_false(isTRUE(all.equal(res_default$B, res_ms2$B)))
})

# ── SE columns ────────────────────────────────────────────────────────────────

test_that("item_metrics se=TRUE adds correct SE column names for 1D", {
  res <- item_metrics(mod_1d, se = TRUE)
  expect_true(all(c("B_se", "A_se", "A_1_se") %in% names(res)))
  # no A_2_se for 1D
  expect_false("A_2_se" %in% names(res))
})

test_that("item_metrics se=TRUE warns about constrained slopes for exploratory 2D model", {
  # Exploratory 2D models fix one slope per item for identification; the
  # Delta-method requires all slopes to be free, so a warning is expected.
  expect_warning(
    res <- item_metrics(mod_2d, se = TRUE),
    regexp = "constrained"
  )
  expect_true(all(c("B_se", "A_se", "A_1_se", "A_2_se") %in% names(res)))
})

test_that("item_metrics se=TRUE B_se and A_se match direct se_pop_* (1D, item 1)", {
  res  <- item_metrics(mod_1d, items = 1L, se = TRUE)
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)
  a_par  <- item_1d@par[1L]
  d_par  <- item_1d@par[2L]
  cc     <- mirt::coef(mod_1d, simplify = TRUE)
  mu  <- as.numeric(cc$means)
  sig <- as.matrix(cc$cov)

  expect_equal(res$B_se[1L],
               se_pop_difficulty(a_par, d_par, vcov_j, mu, sig),
               tolerance = 1e-9)
  expect_equal(res$A_se[1L],
               se_pop_discrimination_scalar(a_par, d_par, vcov_j, mu, sig),
               tolerance = 1e-9)
  expect_equal(res$A_1_se[1L],
               se_pop_discrimination_vector(a_par, d_par, vcov_j, mu, sig)[1L],
               tolerance = 1e-9)
})

test_that("item_metrics se=TRUE warns and produces NA SEs for graded items", {
  expect_warning(
    res <- item_metrics(mod_graded, se = TRUE),
    regexp = "only implemented for 2PL"
  )
  expect_true(all(is.na(res$B_se)))
  expect_true(all(is.na(res$A_se)))
  expect_true(all(is.na(res$A_1_se)))
})

test_that("item_metrics se=TRUE warns and returns NAs when model has no vcov", {
  mod_no_se <- mirt::mirt(mod_1d@Data$data, 1, itemtype = "2PL",
                           verbose = FALSE, SE = FALSE)
  expect_warning(
    res <- item_metrics(mod_no_se, se = TRUE),
    regexp = "SE = TRUE"
  )
  expect_true(all(is.na(res$B_se)))
})

# ── input validation ──────────────────────────────────────────────────────────

test_that("item_metrics errors when mod is not a mirt model", {
  expect_error(item_metrics(list()), regexp = "`mod` must be")
  expect_error(item_metrics("model"), regexp = "`mod` must be")
})

test_that("item_metrics errors when items index is out of range", {
  expect_error(item_metrics(mod_1d, items = 99L),
               regexp = "`items` must contain indices between")
  expect_error(item_metrics(mod_1d, items = 0L),
               regexp = "`items` must contain indices between")
})

test_that("item_metrics errors when items contains NA", {
  expect_error(item_metrics(mod_1d, items = c(1L, NA_integer_)),
               regexp = "`items` must be a numeric vector")
})

test_that("item_metrics errors when max_score is non-positive", {
  expect_error(item_metrics(mod_1d, max_score = 0),
               regexp = "`max_score` must be")
  expect_error(item_metrics(mod_1d, max_score = -1),
               regexp = "`max_score` must be")
})

test_that("item_metrics errors when max_score vector has wrong length", {
  n <- mod_1d@Data$nitems
  expect_error(item_metrics(mod_1d, max_score = rep(1, n + 1L)),
               regexp = "`max_score` must be either length 1")
})

test_that("item_metrics errors when theta_mean has wrong length", {
  expect_error(item_metrics(mod_1d, theta_mean = c(0, 0)),
               regexp = "`theta_mean` must be")
})

test_that("item_metrics errors when theta_cov has wrong dimensions", {
  expect_error(item_metrics(mod_1d, theta_cov = diag(2)),
               regexp = "`theta_cov` must be")
})

test_that("item_metrics errors when se is not a single logical", {
  expect_error(item_metrics(mod_1d, se = "yes"), regexp = "`se` must be")
  expect_error(item_metrics(mod_1d, se = NA),    regexp = "`se` must be")
})
