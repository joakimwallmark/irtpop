test_that("population_difficulty returns value in [0,1] for 1D item", {
  B <- population_difficulty(item_1d, max_score = 1)
  expect_true(is.numeric(B) && length(B) == 1L)
  expect_gte(B, 0)
  expect_lte(B, 1)
})

test_that("population_difficulty 1D matches reference value", {
  B <- population_difficulty(item_1d, max_score = 1)
  expect_equal(B, 0.41192, tolerance = 1e-3)
})

test_that("population_difficulty 2D matches reference value", {
  B <- population_difficulty(item_2d, max_score = 1)
  expect_equal(B, 0.40338, tolerance = 1e-3)
})

test_that("population_discrimination scalar 1D matches reference value", {
  A <- population_discrimination(item_1d, max_score = 1, type = "scalar")
  expect_equal(A, 0.25396, tolerance = 1e-3)
})

test_that("population_discrimination vector == scalar for 1D (identity cov)", {
  A_sc  <- population_discrimination(item_1d, max_score = 1, type = "scalar")
  A_vec <- population_discrimination(item_1d, max_score = 1, type = "vector")
  expect_length(A_vec, 1L)
  expect_equal(abs(A_vec), A_sc, tolerance = 1e-4)
})

test_that("population_discrimination scalar 2D matches reference value", {
  A <- population_discrimination(item_2d, max_score = 1, type = "scalar")
  expect_equal(A, 0.25195, tolerance = 1e-3)
})

test_that("population_discrimination vector 2D has length 2 and correct signs", {
  A_vec <- population_discrimination(item_2d, max_score = 1, type = "vector")
  expect_length(A_vec, 2L)
  expect_equal(A_vec[1], -0.10645, tolerance = 1e-3)
  expect_equal(A_vec[2], -0.22836, tolerance = 1e-3)
})

test_that("population_difficulty defaults max_score from item", {
  B_explicit <- population_difficulty(item_1d, max_score = item_1d@ncat - 1L)
  B_default  <- population_difficulty(item_1d)
  expect_equal(B_explicit, B_default)
})

test_that("item_metrics returns a data.frame with correct structure (1D)", {
  res <- item_metrics(mod_1d)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), mod_1d@Data$nitems)
  expect_true(all(c("item", "item_name", "B", "A", "A_1") %in% names(res)))
})

test_that("item_metrics subset of items works", {
  res <- item_metrics(mod_1d, items = c(1L, 3L))
  expect_equal(nrow(res), 2L)
  expect_equal(res$item, c(1L, 3L))
})

test_that("item_metrics with se = TRUE adds SE columns", {
  res <- item_metrics(mod_1d, se = TRUE)
  expect_true(all(c("B_se", "A_se", "A_1_se") %in% names(res)))
  expect_true(all(res$B_se > 0, na.rm = TRUE))
  expect_true(all(res$A_se > 0, na.rm = TRUE))
})

test_that("item_metrics 1D B column matches population_difficulty", {
  res  <- item_metrics(mod_1d, items = 1L)
  B_fn <- population_difficulty(item_1d, max_score = 1)
  expect_equal(res$B[1L], B_fn, tolerance = 1e-6)
})

# ── input validation ──────────────────────────────────────────────────────────

test_that("population_difficulty errors when item is not an S4 mirt object", {
  expect_error(population_difficulty(list()), regexp = "`item` must be")
  expect_error(population_difficulty(42),     regexp = "`item` must be")
})

test_that("population_difficulty errors on non-positive max_score", {
  expect_error(population_difficulty(item_1d, max_score = 0),
               regexp = "`max_score` must be")
  expect_error(population_difficulty(item_1d, max_score = -2),
               regexp = "`max_score` must be")
})

test_that("population_difficulty errors on theta_mean with wrong length", {
  expect_error(population_difficulty(item_1d, theta_mean = c(0, 0)),
               regexp = "`theta_mean` must be")
})

test_that("population_difficulty errors on theta_cov with wrong dimensions", {
  expect_error(population_difficulty(item_1d, theta_cov = diag(2)),
               regexp = "`theta_cov` must be")
})

test_that("population_discrimination errors when item is not an S4 mirt object", {
  expect_error(population_discrimination(list()), regexp = "`item` must be")
})

test_that("population_discrimination errors on non-positive max_score", {
  expect_error(population_discrimination(item_1d, max_score = 0),
               regexp = "`max_score` must be")
})

test_that("population_discrimination errors on invalid type", {
  expect_error(population_discrimination(item_1d, type = "foo"),
               regexp = "'arg' should be one of")
})

test_that("population_discrimination errors on theta_mean with wrong length", {
  expect_error(population_discrimination(item_1d, theta_mean = c(0, 0)),
               regexp = "`theta_mean` must be")
})

test_that("population_discrimination errors on theta_cov with wrong dimensions", {
  expect_error(population_discrimination(item_1d, theta_cov = diag(2)),
               regexp = "`theta_cov` must be")
})
