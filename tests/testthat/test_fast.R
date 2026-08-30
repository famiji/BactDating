context("fast path correctness")

# The fast path replaces full likelihood/prior recomputations by incremental
# updates that are mathematically identical but not bit-identical in floating
# point: the accumulated l and p pick up O(1e-11) rounding differences, which
# occasionally flip a borderline MH decision and branch the trajectory (two
# valid draws from the SAME posterior, like two runs with different seeds).
# Correctness is therefore verified at two levels:
#   1. trajectory-level: records must agree to ~1e-6 absolute precision
#      (any logic bug shows up as macroscopic divergence, not as 1e-11 noise);
#   2. distribution-level: independent chains run with fast=TRUE and
#      fast=FALSE must estimate the same posterior within Monte Carlo error.

test_that("fast=TRUE trajectories match fast=FALSE to floating-point precision", {
  run_pair <- function(model, ntips, nbIts, seed, simmodel = model, ...) {
    set.seed(seed)
    tree <- simcoaltree(dates = 1990:(1990 + ntips - 1))
    obsphy <- simobsphy(tree, mu = 10, sigma = 10, model = simmodel)
    date <- 1990:(1990 + ntips - 1)
    set.seed(seed + 999)
    slow <- bactdate(obsphy, date, nbIts = nbIts, model = model, fast = FALSE, ...)
    set.seed(seed + 999)
    fast <- bactdate(obsphy, date, nbIts = nbIts, model = model, fast = TRUE, ...)
    expect_equal(slow$record, fast$record, tolerance = 1e-6, scale = 1)
    TRUE
  }

  # All likelihood models
  expect_true(run_pair("arc",          25, 2000, 1))
  expect_true(run_pair("carc",         25, 2000, 2))
  expect_true(run_pair("poisson",      25, 2000, 3))
  expect_true(run_pair("negbin",       25, 2000, 4))
  expect_true(run_pair("strictgamma",  25, 2000, 5))
  expect_true(run_pair("relaxedgamma", 25, 2000, 6))
  expect_true(run_pair("mixedgamma",   25, 2000, 11))
  expect_true(run_pair("mixedcarc",    25, 2000, 12))

  # Root-move and tuning variants (exercise the localized root moves)
  expect_true(run_pair("arc",          30, 2000, 7, updateRoot = FALSE))
  expect_true(run_pair("arc",          30, 2000, 8, updateRoot = "branch"))
  expect_true(run_pair("carc",         30, 2000, 9, tuning = FALSE))
  expect_true(run_pair("arc",          30, 2000, 10, useCoalPrior = FALSE))

  # Recombination-aware input (5-column tab)
  set.seed(13)
  tr <- simcoaltree(dates = 1990:2014)
  op <- simobsphy(tr, mu = 10, sigma = 10, model = "arc")
  op$unrec <- rep(0.85, length(op$edge.length))
  set.seed(1313)
  s1 <- bactdate(op, 1990:2014, nbIts = 1500, useRec = TRUE, fast = FALSE)
  set.seed(1313)
  s2 <- bactdate(op, 1990:2014, nbIts = 1500, useRec = TRUE, fast = TRUE)
  expect_equal(s1$record, s2$record, tolerance = 1e-6, scale = 1)

  # Larger tree, longer run: stresses the incremental prior accumulation
  expect_true(run_pair("arc", 100, 4000, 14))
})

test_that("fast=TRUE,exact=TRUE is bit-identical to the original fast=FALSE", {
  # exact=TRUE gives up the incremental prior and the localized root moves in
  # exchange for bit-identical trajectories: identical() must hold.
  run_pair <- function(model, ntips, nbIts, seed, simmodel = model, ...) {
    set.seed(seed)
    tree <- simcoaltree(dates = 1990:(1990 + ntips - 1))
    obsphy <- simobsphy(tree, mu = 10, sigma = 10, model = simmodel)
    date <- 1990:(1990 + ntips - 1)
    set.seed(seed + 999)
    slow <- bactdate(obsphy, date, nbIts = nbIts, model = model, fast = FALSE, ...)
    set.seed(seed + 999)
    fast <- bactdate(obsphy, date, nbIts = nbIts, model = model, fast = TRUE, exact = TRUE, ...)
    identical(slow$record, fast$record)
  }

  expect_true(run_pair("arc",          25, 2000, 21))
  expect_true(run_pair("carc",         25, 2000, 22))
  expect_true(run_pair("poisson",      25, 2000, 23))
  expect_true(run_pair("negbin",       25, 2000, 24))
  expect_true(run_pair("strictgamma",  25, 2000, 25))
  expect_true(run_pair("relaxedgamma", 25, 2000, 26))
  expect_true(run_pair("mixedgamma",   25, 2000, 27))
  expect_true(run_pair("mixedcarc",    25, 2000, 28))

  expect_true(run_pair("arc",          30, 2000, 29, updateRoot = FALSE))
  expect_true(run_pair("arc",          30, 2000, 30, updateRoot = "branch"))
  expect_true(run_pair("carc",         30, 2000, 31, tuning = FALSE))
  expect_true(run_pair("arc",          30, 2000, 32, useCoalPrior = FALSE))

  # Recombination-aware input (5-column tab)
  set.seed(33)
  tr <- simcoaltree(dates = 1990:2014)
  op <- simobsphy(tr, mu = 10, sigma = 10, model = "arc")
  op$unrec <- rep(0.85, length(op$edge.length))
  set.seed(3333)
  s1 <- bactdate(op, 1990:2014, nbIts = 1500, useRec = TRUE, fast = FALSE)
  set.seed(3333)
  s2 <- bactdate(op, 1990:2014, nbIts = 1500, useRec = TRUE, fast = TRUE, exact = TRUE)
  expect_true(identical(s1$record, s2$record))

  # Larger tree, longer run
  expect_true(run_pair("arc", 100, 4000, 34))
})

test_that("fast=TRUE and fast=FALSE chains estimate the same posterior", {
  # Six independent chains per path; the between-path difference of the
  # pooled posterior means must stay within ~3 standard errors of the
  # between-chain spread (a real bias would show up as many sds).
  ntips <- 40; nbIts <- 1.5e4
  set.seed(42)
  tree <- simcoaltree(dates = 1990:(1990 + ntips - 1))
  obsphy <- simobsphy(tree, mu = 10, sigma = 10, model = "arc")
  date <- 1990:(1990 + ntips - 1)

  parmeans <- function(fast) {
    vapply(991:996, function(sd) {
      set.seed(sd)
      res <- bactdate(obsphy, date, nbIts = nbIts, model = "arc", fast = fast)
      rec <- res$record
      keep <- floor(nrow(rec) / 2):nrow(rec)
      c(mu = mean(rec[keep, "mu"]), alpha = mean(rec[keep, "alpha"]),
        rootdate = mean(rec[keep, "root"]))
    }, numeric(3))
  }
  mfast <- parmeans(TRUE)
  mslow <- parmeans(FALSE)
  allm <- cbind(mfast, mslow)
  for (pname in rownames(allm)) {
    sd_mc <- sd(allm[pname, ])
    # std of the difference of two means of 6 chains ~ sd*sqrt(2/6) = 0.58*sd;
    # threshold 1.75*sd ~ 3 standard errors under the null hypothesis.
    expect_lt(abs(mean(mfast[pname, ]) - mean(mslow[pname, ])), 1.75 * sd_mc)
  }
})
