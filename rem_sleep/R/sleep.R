# A simulated participant ("agent"): lexical-semantic (LS) traces and an episodic store. Each
# episodic trace holds n_slots content blocks of w features (1 = single item, 2 = paired
# associate).
#
# Context is not modeled. As in S&S (1997), test is assumed to activate the episodic traces
# from the relevant period (here, the whole day: study list plus any wake events); only content
# features enter the likelihood. So the only cost of wakefulness is interference from new traces,
# and passive sleep leaves memory exactly as it was at the end of study.

new_agent <- function(par, n_lex, n_slots = 1) {
  list(lex = gen_items(n_lex, par$w, par$g),
       store = matrix(0, 0, n_slots * par$w),
       src = matrix(0L, 0, n_slots),  # true lexical item(s) behind each trace; analysis only
       n_slots = n_slots)
}

slot_cols <- function(s, par) (s - 1) * par$w + seq_len(par$w)

# Probe for an event: the concatenated LS traces of its item(s); one row per event.
probes_of <- function(a, items) {
  items <- matrix(items, ncol = a$n_slots)
  t(apply(items, 1, function(i) c(t(a$lex[i, , drop = FALSE]))))
}

# Study each event (row of `items`: lexical indices, one column per slot) `reps` times, cycling
# through the list. Following S&S (1997), a repeated event adds to its existing trace.
study <- function(a, items, par, reps = 1) {
  items <- matrix(items, ncol = a$n_slots)
  rows <- nrow(a$store) + seq_len(nrow(items))
  a$store <- rbind(a$store, matrix(0, nrow(items), ncol(a$store)))
  a$src <- rbind(a$src, items)
  for (r in seq_len(reps)) for (i in seq_len(nrow(items))) {
    a$store[rows[i], ] <- encode(a$store[rows[i], ], c(t(a$lex[items[i, ], , drop = FALSE])), par)
  }
  a
}

# Wake retention interval: n_events new events drawn from `pool` (known items), each stored once.
wake <- function(a, par, n_events, pool) {
  study(a, matrix(sample(pool, n_events * a$n_slots), ncol = a$n_slots), par)
}

# Sleep retention interval: no new encoding. With mech = character(), sleep is purely passive.
# Otherwise n_replay traces are reactivated; for each content slot of a replayed trace the
# lexical item that generated it is inferred by sampling from the REM posterior, then:
#   "reactivate": empty features of the trace are re-stored from that LS trace, with the same
#                 storage (u) and copying (c) processes as study;
#   "reencode":   f randomly chosen features of that LS trace are incremented by 1 (rarer values
#                 are more diagnostic), and every trace slot whose best-matching item (at sleep
#                 onset) is that item and which stored the old value is updated to the new value.
# Targeted memory reactivation (TMR): with probability p_cue a replay is cued instead of
# spontaneous. A cued item's LS trace probes slot 1 of the store, a trace is retrieved in
# proportion to its likelihood ratio, and the cued item (rather than an inferred one) is the LS
# trace used for slot 1; other slots (e.g., the B of an A-B pair) are identified as usual.
# id_scope sets which LS traces are candidates when identifying a replayed trace: "lexicon" (all
# known items) or "today" (only items experienced since the last sleep, i.e., tagged at encoding).
# prio sets how strongly spontaneous replay favors weak traces: trace j is replayed with
# probability proportional to (fraction of its content features still empty)^prio, recomputed
# before each replay. prio = 0 is uniform replay.
sleep <- function(a, par, n_replay = 0, mech = character(), f = 6, cue = NULL, p_cue = 0,
                  id_scope = c("lexicon", "today"), prio = 0) {
  stopifnot(all(mech %in% c("reactivate", "reencode")))
  if (n_replay == 0 || !length(mech)) return(a)
  cand <- switch(match.arg(id_scope), lexicon = seq_len(nrow(a$lex)), today = sort(unique(c(a$src))))
  identify <- function(x) identify_loglik(x, a$lex[cand, , drop = FALSE], par)

  n <- nrow(a$store)
  if ("reencode" %in% mech) {
    map_id <- matrix(0L, n, a$n_slots)
    for (j in seq_len(n)) for (s in seq_len(a$n_slots)) {
      map_id[j, s] <- cand[which.max(identify(a$store[j, slot_cols(s, par)]))]
    }
  }
  for (r in seq_len(n_replay)) {
    k_cue <- NULL
    if (length(cue) && runif(1) < p_cue) {
      k_cue <- cue[sample.int(length(cue), 1)]
      ll <- loglik(a$lex[k_cue, ], a$store[, slot_cols(1, par), drop = FALSE], par)
      j <- sample.int(n, 1, prob = exp(ll - max(ll)))
    } else {
      j <- if (prio == 0) sample.int(n, 1) else {
        w <- rowMeans(a$store == 0)^prio
        if (sum(w) > 0) sample.int(n, 1, prob = w) else sample.int(n, 1)
      }
    }
    for (s in seq_len(a$n_slots)) {
      cols <- slot_cols(s, par)
      if (is.null(k_cue) || s > 1) {
        ll <- identify(a$store[j, cols])
        k <- cand[sample.int(length(cand), 1, prob = exp(ll - max(ll)))]
      } else {
        k <- k_cue
      }
      if ("reactivate" %in% mech) a$store[j, cols] <- encode(a$store[j, cols], a$lex[k, ], par)
      if ("reencode" %in% mech) {
        feats <- sample.int(par$w, f)
        old <- a$lex[k, feats]
        a$lex[k, feats] <- old + 1
        hits <- which(map_id == k, arr.ind = TRUE)  # any slot of any trace identified as k
        for (h in seq_len(nrow(hits))) {
          fc <- slot_cols(hits[h, 2], par)[feats]
          same <- a$store[hits[h, 1], fc] == old
          a$store[hits[h, 1], fc[same]] <- old[same] + 1
        }
      }
    }
  }
  a
}
