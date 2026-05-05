test_that("se_pop_difficulty returns positive scalar (1D)", {
  a_par <- item_1d@par[1L]
  d_par <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)

  se <- se_pop_difficulty(a_par, d_par, vcov_j)
  expect_true(is.numeric(se) && length(se) == 1L)
  expect_gte(se, 0)
})

test_that("se_pop_difficulty 1D matches reference value", {
  a_par <- item_1d@par[1L]
  d_par <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)

  se <- se_pop_difficulty(a_par, d_par, vcov_j)
  expect_equal(se, 0.004920, tolerance = 1e-4)
})

test_that("se_pop_discrimination_scalar 1D matches reference value", {
  a_par <- item_1d@par[1L]
  d_par <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)

  se <- se_pop_discrimination_scalar(a_par, d_par, vcov_j)
  expect_equal(se, 0.004964, tolerance = 1e-4)
})

test_that("se_pop_discrimination_vector 1D equals scalar SE (1D identity)", {
  a_par <- item_1d@par[1L]
  d_par <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)

  se_sc  <- se_pop_discrimination_scalar(a_par, d_par, vcov_j)
  se_vec <- se_pop_discrimination_vector(a_par, d_par, vcov_j)
  expect_length(se_vec, 1L)
  expect_equal(se_vec, se_sc, tolerance = 1e-6)
})

test_that("se_pop_discrimination_scalar returns NA when norm_a near zero", {
  a_zero <- c(0, 0)
  d_par  <- 0
  vcov_j <- diag(3)
  expect_true(is.na(se_pop_discrimination_scalar(a_zero, d_par, vcov_j)))
})

test_that("aux_integrals gamma is positive and kappa has opposite sign for typical params", {
  ints <- irtpop:::aux_integrals(a = 1.5, d = 0)
  expect_gt(ints$gamma, 0)
  expect_length(ints$mu_vec, 1L)
  expect_length(ints$lambda_vec, 1L)
})

# ── input validation ──────────────────────────────────────────────────────────

test_that("se_pop_difficulty errors on non-numeric a", {
  d      <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)
  expect_error(se_pop_difficulty("x", d, vcov_j), regexp = "`a` must be")
})

test_that("se_pop_difficulty errors on vector d", {
  a      <- item_1d@par[1L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)
  expect_error(se_pop_difficulty(a, c(0, 0), vcov_j), regexp = "`d` must be")
})

test_that("se_pop_difficulty errors on wrong-sized vcov_j", {
  a <- item_1d@par[1L]
  d <- item_1d@par[2L]
  expect_error(se_pop_difficulty(a, d, diag(3)), regexp = "`vcov_j` must be")
})

test_that("se_pop_difficulty errors on theta_mean with wrong length", {
  a      <- item_1d@par[1L]
  d      <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)
  expect_error(se_pop_difficulty(a, d, vcov_j, theta_mean = c(0, 0)),
               regexp = "`theta_mean` must be")
})

test_that("se_pop_difficulty errors on theta_cov with wrong dimensions", {
  a      <- item_1d@par[1L]
  d      <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)
  expect_error(se_pop_difficulty(a, d, vcov_j, theta_cov = diag(2)),
               regexp = "`theta_cov` must be")
})

test_that("se_pop_discrimination_scalar errors on non-numeric a", {
  d      <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)
  expect_error(se_pop_discrimination_scalar(TRUE, d, vcov_j),
               regexp = "`a` must be")
})

test_that("se_pop_discrimination_scalar errors on wrong-sized vcov_j", {
  a <- item_1d@par[1L]
  d <- item_1d@par[2L]
  expect_error(se_pop_discrimination_scalar(a, d, diag(3)),
               regexp = "`vcov_j` must be")
})

test_that("se_pop_discrimination_scalar errors on theta_cov with wrong dimensions", {
  a      <- item_1d@par[1L]
  d      <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)
  expect_error(se_pop_discrimination_scalar(a, d, vcov_j, theta_cov = diag(2)),
               regexp = "`theta_cov` must be")
})

test_that("se_pop_discrimination_vector errors on NA in a", {
  d      <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)
  expect_error(se_pop_discrimination_vector(NA_real_, d, vcov_j),
               regexp = "`a` must be")
})

test_that("se_pop_discrimination_vector errors on wrong-sized vcov_j", {
  a <- item_1d@par[1L]
  d <- item_1d@par[2L]
  expect_error(se_pop_discrimination_vector(a, d, diag(5)),
               regexp = "`vcov_j` must be")
})

test_that("se_pop_discrimination_vector errors on theta_mean with wrong length", {
  a      <- item_1d@par[1L]
  d      <- item_1d@par[2L]
  vcov_j <- .extract_item_vcov(mod_1d, item_1d)
  expect_error(se_pop_discrimination_vector(a, d, vcov_j, theta_mean = c(0, 0)),
               regexp = "`theta_mean` must be")
})

