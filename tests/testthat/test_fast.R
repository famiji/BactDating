context("fast internal-date update")

test_that("Optimized C++ internal-date update is bit-identical to the default loop", {
  # For each supported model, running with fast=TRUE must reproduce exactly the
  # same MCMC trajectory (record matrix) as fast=FALSE given the same seed.
  run_pair <- function(model, ntips, nbIts, seed, simmodel = model) {
    set.seed(seed)
    tree <- simcoaltree(dates = 1990:(1990 + ntips - 1))
    obsphy <- simobsphy(tree, mu = 10, sigma = 10, model = simmodel)
    date <- 1990:(1990 + ntips - 1)
    set.seed(seed + 999)
    slow <- bactdate(obsphy, date, nbIts = nbIts, model = model, fast = FALSE)
    set.seed(seed + 999)
    fast <- bactdate(obsphy, date, nbIts = nbIts, model = model, fast = TRUE)
    identical(slow$record, fast$record)
  }

  expect_true(run_pair("arc",          25, 2000, 1))
  expect_true(run_pair("carc",         25, 2000, 2))
  expect_true(run_pair("poisson",      25, 2000, 3))
  expect_true(run_pair("negbin",       25, 2000, 4))
  expect_true(run_pair("strictgamma",  25, 2000, 5))
  expect_true(run_pair("relaxedgamma", 25, 2000, 6))
  expect_true(run_pair("arc",          60, 1500, 7))  # larger tree
})
