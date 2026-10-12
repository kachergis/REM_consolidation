# Run from rem_sleep/: Rscript -e 'testthat::test_file("tests/test_rem.R")' (testthat runs it from tests/)
library(testthat)
source("../R/rem.R")
source("../R/sleep.R")

par <- rem_par()

# Direct transcription of the REM likelihood ratio, one feature at a time.
ref_loglik <- function(probe, trace, par) {
  out <- 0
  for (f in which(trace > 0)) {
    out <- out + if (trace[f] == probe[f]) log_match_ratio(trace[f], par) else log(1 - par$c)
  }
  out
}

test_that("likelihood reproduces the Shiffrin & Steyvers (1997) worked example", {
  # g = .4, c = .7; traces (0,1,0,3) and (2,2,1,0), probe (2,3,4,3): mean LR = .92
  mem <- matrix(c(0, 1, 0, 3, 2, 2, 1, 0), 2, byrow = TRUE)
  p <- rem_par(g = .4, c = .7)
  expect_equal(mean(exp(loglik(c(2, 3, 4, 3), mem, p))), .919, tolerance = .001)
})

test_that("likelihood matches the per-feature definition on random traces", {
  set.seed(1)
  probe <- rgeom1(20, par$g)
  traces <- gen_items(50, 20, par$g) * (matrix(runif(1000), 50) < .5)
  traces[1:10, ] <- ifelse(traces[1:10, ] > 0, probe[col(traces[1:10, ])], 0)  # some true matches
  expect_equal(loglik(probe, traces, par), apply(traces, 1, ref_loglik, probe = probe, par = par))
})

test_that("repeated matches and pure mismatches are counted once each", {
  expect_equal(loglik(c(1, 1, 1, 1), matrix(c(1, 1, 1, 0), 1), par), 3 * log_match_ratio(1, par))
  expect_equal(loglik(c(1, 1, 1, 1), matrix(c(2, 2, 0, 0), 1), par), 2 * log(1 - par$c))
  expect_equal(loglik(c(1, 1), matrix(0, 1, 2), par), 0)
})

test_that("identification likelihood is the likelihood with the trace as the mask", {
  set.seed(2)
  lex <- gen_items(30, 20, par$g)
  tr <- encode(rep(0, 20), lex[7, ], par)
  expect_equal(identify_loglik(tr, lex, par), sapply(1:30, function(k) ref_loglik(lex[k, ], tr, par)))
  full <- lex[7, ]
  expect_equal(which.max(identify_loglik(full, lex, par)), 7)
})

test_that("encoding stores with probability u and copies correctly with probability c", {
  set.seed(3)
  x <- rgeom1(2e5, par$g)
  tr <- encode(rep(0, 2e5), x, par)
  expect_equal(mean(tr > 0), par$u, tolerance = .01)
  # stored value equals input if copied correctly, or by chance if not
  p_chance <- mean(par$g * (1 - par$g)^(x[tr > 0] - 1))
  expect_equal(mean(tr[tr > 0] == x[tr > 0]), par$c + (1 - par$c) * p_chance, tolerance = .01)
  # stored features are never overwritten
  expect_identical(encode(tr, x + 1, par)[tr > 0], tr[tr > 0])
})

test_that("sdt gives d' = 0 and AUC = .5 for identical distributions", {
  x <- c(-1, 0.5, 2, 3)
  expect_equal(unname(sdt(x, x)[c("dprime", "auc")]), c(0, .5))
  expect_equal(unname(sdt(c(2, 3), c(-1, -2))["auc"]), 1)
})

test_that("values start at 1", {
  expect_gt(min(gen_items(1000, 20, .45)), 0)
})

# Benchmark: recognition improves with study repetitions and worsens with list length.
recog <- function(n_study, reps, seed) {
  set.seed(seed)
  a <- new_agent(par, n_lex = 400)
  a <- study(a, 1:n_study, par, reps)
  sdt(log_odds(probes_of(a, 1:n_study), a$store, par),
      log_odds(probes_of(a, 201:(200 + n_study)), a$store, par))["dprime"]
}

test_that("d' increases with repetitions and decreases with list length", {
  d <- function(n, r) mean(sapply(1:20, function(s) recog(n, r, s)))
  expect_gt(d(40, 3), d(40, 1) + .3)
  expect_gt(d(20, 1), d(80, 1))
})

test_that("reactivation only fills empty features; re-encoding increments LS values", {
  set.seed(4)
  a <- new_agent(par, n_lex = 200)
  a <- study(a, 1:30, par)
  b <- sleep(a, par, n_replay = 30, mech = "reactivate")
  filled <- a$store > 0
  expect_identical(b$store[filled], a$store[filled])
  expect_gt(sum(b$store > 0), sum(filled))
  expect_identical(b$lex, a$lex)

  e <- sleep(a, par, n_replay = 30, mech = "reencode", f = 6)
  expect_equal(sum(e$lex - a$lex), 30 * 6)
  expect_true(all(e$lex >= a$lex))
  expect_true(all(e$store >= a$store))

  passive <- sleep(a, par)
  expect_identical(passive$store, a$store)
})

source("../R/recall.R")

test_that("cued recall retrieves targets from complete, correct traces", {
  set.seed(5)
  a <- new_agent(par, n_lex = 300, n_slots = 2)
  A <- 1:20; B <- 21:40
  a$store <- probes_of(a, cbind(A, B))
  a$src <- cbind(A, B)
  expect_identical(cued_recall(a, A, par), B)
})

test_that("recall rises with repetitions; A-C learning interferes and intrudes", {
  rec <- function(reps, ac, seed) {
    set.seed(seed)
    a <- new_agent(par, n_lex = 400, n_slots = 2)
    A <- 1:30; B <- 31:60; C <- 61:90
    a <- study(a, cbind(A, B), par, reps)
    if (ac) a <- study(a, cbind(A, C), par, reps)
    r <- cued_recall(a, A, par)
    c(score_recall(r, B), c_intr = mean(r %in% C))
  }
  m <- function(reps, ac) rowMeans(sapply(1:20, function(s) rec(reps, ac, s)))
  expect_gt(m(3, FALSE)["correct"], m(1, FALSE)["correct"] + .2)
  base <- m(2, FALSE); intf <- m(2, TRUE)
  expect_gt(base["correct"], intf["correct"] + .1)
  expect_gt(intf["c_intr"], .1)
})

test_that("a TMR cue on pair traces restores slot 1 from the cued item", {
  set.seed(6)
  a <- new_agent(par, n_lex = 200, n_slots = 2)
  a <- study(a, cbind(1:10, 11:20), par)
  b <- sleep(a, par, n_replay = 200, mech = "reactivate", cue = 3, p_cue = 1)
  # cued item 3's own trace gains slot-1 features copied from item 3
  expect_gt(sum(b$store[3, 1:20] > 0), sum(a$store[3, 1:20] > 0))
  expect_gt(mean(b$store[3, 1:20][b$store[3, 1:20] > 0] == a$lex[3, b$store[3, 1:20] > 0]), .6)
})
