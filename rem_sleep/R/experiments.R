# Simulated sleep-and-memory experiments. Each function simulates one participant and returns a
# data frame; run_subjects() runs many participants with per-subject seeds (reproducible
# regardless of the number of cores).

run_subjects <- function(fun, n_subj, seed, cores = max(1, parallel::detectCores() - 1), ...) {
  out <- parallel::mclapply(seq_len(n_subj), function(s) {
    set.seed(seed + s)
    cbind(subj = s, fun(...))
  }, mc.cores = cores)
  do.call(rbind, out)
}

# Retention-interval conditions shared by the experiments. Only wake stores new, unrelated
# events; "immediate" and passive "sleep" are therefore identical (see R/sleep.R).
retention <- function(a, cond, par, n_events, pool, n_replay, f, ...) {
  switch(cond,
         immediate = a,
         wake = wake(a, par, n_events, pool),
         sleep = a,
         `sleep+reactivate` = sleep(a, par, n_replay, "reactivate", f, ...),
         `sleep+reencode` = sleep(a, par, n_replay, "reencode", f, ...),
         `sleep+both` = sleep(a, par, n_replay, c("reactivate", "reencode"), f, ...))
}

CONDS <- c("immediate", "wake", "sleep", "sleep+reactivate", "sleep+reencode", "sleep+both")

# Exp. 1: sleep vs. wake retention x study repetitions x consolidation mechanism
# (single-item recognition). Generalizes the 2016 paper's simulation, adding a wake control.
exp_sleep_wake <- function(par, reps = c(1, 3, 5), n_study = 50, n_lex = 1000, n_events = 100,
                           n_replay = 50, f = 6, id_scope = "lexicon") {
  a0 <- new_agent(par, n_lex)
  targets <- 1:n_study
  foils <- n_study + 1:n_study
  pool <- setdiff(seq_len(n_lex), c(targets, foils))
  do.call(rbind, lapply(reps, function(r) {
    studied <- study(a0, targets, par, r)
    do.call(rbind, lapply(CONDS, function(cond) {
      a <- retention(studied, cond, par, n_events, pool, n_replay, f, id_scope = id_scope)
      res <- sdt(log_odds(probes_of(a, targets), a$store, par),
                 log_odds(probes_of(a, foils), a$store, par), par$crit)
      data.frame(id_scope = id_scope, reps = r, cond = cond, t(res))
    }))
  }))
}

# Exp. 2: interference challenge, after Ellenbogen et al. (2006): paired associates (A-B) are
# studied, followed by a sleep or wake interval, then (or not) an interfering A-C list, then
# associative recognition of intact A-B vs. rearranged A-B' pairs.
exp_interference <- function(par, n_pairs = 30, reps = 2, n_lex = 1000, n_events = 60,
                             n_replay = 30, f = 6, conds = c("wake", "sleep", "sleep+reactivate",
                                                            "sleep+reencode", "sleep+both"),
                             id_scope = "lexicon") {
  a0 <- new_agent(par, n_lex, n_slots = 2)
  A <- 1:n_pairs
  B <- n_pairs + A
  C <- 2 * n_pairs + A
  rearranged <- cbind(A, B[c(2:n_pairs, 1)])  # each A with another pair's B
  pool <- setdiff(seq_len(n_lex), c(A, B, C))
  studied <- study(a0, cbind(A, B), par, reps)
  do.call(rbind, lapply(conds, function(cond) {
    a <- retention(studied, cond, par, n_events, pool, n_replay, f, id_scope = id_scope)
    do.call(rbind, lapply(c("none", "A-C"), function(intf) {
      b <- if (intf == "A-C") study(a, cbind(A, C), par, reps) else a
      res <- sdt(log_odds(probes_of(b, cbind(A, B)), b$store, par),
                 log_odds(probes_of(b, rearranged), b$store, par), par$crit)
      data.frame(id_scope = id_scope, cond = cond, interference = intf, t(res))
    }))
  }))
}

# Exp. 3: targeted memory reactivation (e.g., Rudoy et al., 2009): half the studied items are
# cued during sleep. Cued and uncued d' share one false-alarm rate. The replay budget is fixed,
# so cued replays replace spontaneous ones.
exp_tmr <- function(par, n_study = 50, reps = 1, n_lex = 1000, n_events = 100, n_replay = 50,
                    f = 6, p_cue = .5, mechs = c("sleep+reactivate", "sleep+reencode", "sleep+both"),
                    id_scope = "lexicon") {
  a0 <- new_agent(par, n_lex)
  targets <- 1:n_study
  cued <- targets[seq(1, n_study, 2)]
  foils <- n_study + 1:n_study
  studied <- study(a0, targets, par, reps)
  do.call(rbind, lapply(mechs, function(cond) {
    do.call(rbind, lapply(c(FALSE, TRUE), function(tmr) {
      a <- retention(studied, cond, par, n_events, NULL, n_replay, f,
                     cue = if (tmr) cued, p_cue = if (tmr) p_cue else 0, id_scope = id_scope)
      old <- log_odds(probes_of(a, targets), a$store, par)
      new <- log_odds(probes_of(a, foils), a$store, par)
      is_cued <- targets %in% cued
      rbind(data.frame(id_scope = id_scope, cond = cond, tmr = tmr, items = "cued",
                       t(sdt(old[is_cued], new, par$crit))),
            data.frame(id_scope = id_scope, cond = cond, tmr = tmr, items = "uncued",
                       t(sdt(old[!is_cued], new, par$crit))))
    }))
  }))
}

# ---- Cued recall of paired associates (A -> B); see R/recall.R --------------------------------

# Exp. 4: sleep vs. wake x repetitions x mechanism, cued recall.
exp_recall_sleep_wake <- function(par, reps = 1:3, n_pairs = 40, n_lex = 1000, n_events = 80,
                                  n_replay = 40, f = 6, id_scope = "lexicon") {
  a0 <- new_agent(par, n_lex, n_slots = 2)
  A <- 1:n_pairs
  B <- n_pairs + A
  pool <- setdiff(seq_len(n_lex), c(A, B))
  do.call(rbind, lapply(reps, function(r) {
    studied <- study(a0, cbind(A, B), par, r)
    do.call(rbind, lapply(CONDS, function(cond) {
      a <- retention(studied, cond, par, n_events, pool, n_replay, f, id_scope = id_scope)
      data.frame(id_scope = id_scope, reps = r, cond = cond,
                 t(score_recall(cued_recall(a, A, par), B)))
    }))
  }))
}

# Exp. 5: interference challenge with cued recall, as in Ellenbogen et al. (2006): A-B, then a
# sleep or wake interval, then (or not) A-C learning, then recall of B given A.
exp_recall_interference <- function(par, n_pairs = 40, reps = 2, n_lex = 1000, n_events = 80,
                                    n_replay = 40, f = 6, id_scope = "lexicon",
                                    conds = c("wake", "sleep", "sleep+reactivate",
                                              "sleep+reencode", "sleep+both")) {
  a0 <- new_agent(par, n_lex, n_slots = 2)
  A <- 1:n_pairs
  B <- n_pairs + A
  C <- 2 * n_pairs + A
  pool <- setdiff(seq_len(n_lex), c(A, B, C))
  studied <- study(a0, cbind(A, B), par, reps)
  do.call(rbind, lapply(conds, function(cond) {
    a <- retention(studied, cond, par, n_events, pool, n_replay, f, id_scope = id_scope)
    do.call(rbind, lapply(c("none", "A-C"), function(intf) {
      b <- if (intf == "A-C") study(a, cbind(A, C), par, reps) else a
      r <- cued_recall(b, A, par)
      data.frame(id_scope = id_scope, cond = cond, interference = intf,
                 t(score_recall(r, B)), c_intrusion = mean(r %in% C))
    }))
  }))
}

# Exp. 6: TMR with cued recall: half the A cues are presented during sleep (fixed replay budget).
exp_recall_tmr <- function(par, n_pairs = 40, reps = 1, n_lex = 1000, n_replay = 40, f = 6,
                           p_cue = .5, id_scope = "lexicon",
                           mechs = c("sleep+reactivate", "sleep+reencode", "sleep+both")) {
  a0 <- new_agent(par, n_lex, n_slots = 2)
  A <- 1:n_pairs
  B <- n_pairs + A
  cued <- seq(1, n_pairs, 2)
  studied <- study(a0, cbind(A, B), par, reps)
  do.call(rbind, lapply(mechs, function(cond) {
    do.call(rbind, lapply(c(FALSE, TRUE), function(tmr) {
      a <- retention(studied, cond, par, 0, NULL, n_replay, f, id_scope = id_scope,
                     cue = if (tmr) A[cued], p_cue = if (tmr) p_cue else 0)
      r <- cued_recall(a, A, par)
      rbind(data.frame(id_scope = id_scope, cond = cond, tmr = tmr, items = "cued",
                       t(score_recall(r[cued], B[cued]))),
            data.frame(id_scope = id_scope, cond = cond, tmr = tmr, items = "uncued",
                       t(score_recall(r[-cued], B[-cued]))))
    }))
  }))
}

# Exp. 7: replay prioritization with a mixed-strength list (cf. Drosopoulos et al., 2007;
# Denis et al., 2021; Schapiro et al., 2018). Half the pairs are studied once (weak), half three
# times (strong); sleep replay is uniform (prio = 0) or favors weak traces (prio > 0).
# Cued recall is scored separately for weak and strong pairs.
exp_recall_strength <- function(par, n_pairs = 40, n_lex = 1000, n_events = 80, n_replay = 40,
                                f = 6, prios = c(0, 1, 4), id_scope = "lexicon") {
  a0 <- new_agent(par, n_lex, n_slots = 2)
  A <- 1:n_pairs
  B <- n_pairs + A
  strong <- seq(2, n_pairs, 2)
  pool <- setdiff(seq_len(n_lex), c(A, B))
  # strong pairs get two extra presentations; all pairs then get one (strong = 3, weak = 1)
  studied <- study(a0, cbind(A, B), par, 1)
  for (r in 1:2) for (i in strong) {
    studied$store[i, ] <- encode(studied$store[i, ], c(t(studied$lex[c(A[i], B[i]), ])), par)
  }
  conds <- rbind(data.frame(cond = c("immediate", "wake", "sleep"), prio = 0),
                 expand.grid(cond = c("sleep+reactivate", "sleep+both"), prio = prios,
                             stringsAsFactors = FALSE))
  is_strong <- A %in% strong
  do.call(rbind, lapply(seq_len(nrow(conds)), function(i) {
    cond <- conds$cond[i]
    prio <- conds$prio[i]
    a <- retention(studied, cond, par, n_events, pool, n_replay, f, id_scope = id_scope,
                   prio = prio)
    r <- cued_recall(a, A, par)
    rbind(data.frame(id_scope = id_scope, cond = cond, prio = prio, strength = "weak",
                     t(score_recall(r[!is_strong], B[!is_strong]))),
          data.frame(id_scope = id_scope, cond = cond, prio = prio, strength = "strong",
                     t(score_recall(r[is_strong], B[is_strong]))))
  }))
}
