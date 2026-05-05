# Shared fixtures — fitted once for the whole test suite.
# Models are cached in the session environment by testthat.

suppressPackageStartupMessages({
  library(mirt)
})

# ── 1D 2PL ────────────────────────────────────────────────────────────────────
set.seed(13)
.a1 <- matrix(rep(1.5, 5))
.d1 <- matrix(rnorm(5))
.dat1 <- mirt::simdata(.a1, .d1, N = 10000, itemtype = "2PL")

mod_1d <- mirt::mirt(.dat1, 1, itemtype = "2PL", verbose = FALSE, SE = TRUE)
item_1d <- mirt::extract.item(mod_1d, 1)

# ── 2D 2PL ────────────────────────────────────────────────────────────────────
set.seed(13)
.a2 <- cbind(c(rep(1.5, 3), rep(0.3, 3)),
             c(rep(0.3, 3), rep(1.5, 3)))
.d2 <- matrix(rnorm(6))
.dat2 <- mirt::simdata(.a2, .d2, N = 10000, itemtype = "2PL")

mod_2d <- mirt::mirt(.dat2, 2, itemtype = "2PL", verbose = FALSE, SE = TRUE)
item_2d <- mirt::extract.item(mod_2d, 1)

# ── Shared test helper ────────────────────────────────────────────────────────
.extract_item_vcov <- function(mod, item) {
  vc           <- stats4::vcov(mod)
  parnum_vc    <- as.integer(sub(".*\\.", "", colnames(vc)))
  free_parnums <- item@parnum[item@est]
  pos          <- match(free_parnums, parnum_vc)
  vc[pos, pos, drop = FALSE]
}
