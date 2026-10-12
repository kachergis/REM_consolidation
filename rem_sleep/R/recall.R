# Cued recall, following the sampling-then-recovery scheme of ARC-REM (Diller, Nobel & Shiffrin,
# 2001), which substitutes REM likelihood ratios for SAM's activations.
# Sampling: the cue's LS trace is matched against the cue slot of every trace; up to k_max
# distinct traces are sampled in proportion to their likelihood ratios (each trace's posterior
# probability of being the cue's own trace, under a uniform prior). As in ARC-REM, a sampled
# trace is then accepted or rejected: it is used only if its cue log LR exceeds acc_crit.
# Recovery: the target slot of a sampled trace is identified against the whole lexicon by the
# same computation as in sleep replay (identify_loglik). If the best item's posterior is at least
# rec_crit, that item is output and search stops; otherwise the next sampled trace is tried.
# (This recovery rule is our choice; ARC-REM's own recovery equations were not consulted.)
# Returns the recalled lexical index for each cue, NA when nothing is recovered.
cued_recall <- function(a, cues, par, cue_slot = 1, target_slot = 2) {
  cue_cols <- slot_cols(cue_slot, par)
  tgt_cols <- slot_cols(target_slot, par)
  sapply(cues, function(k) {
    ll <- loglik(a$lex[k, ], a$store[, cue_cols, drop = FALSE], par)
    p <- exp(ll - max(ll))
    for (j in sample.int(length(p), min(par$k_max, sum(p > 0)), prob = p)) {
      if (ll[j] <= par$acc_crit) next
      id <- identify_loglik(a$store[j, tgt_cols], a$lex, par)
      post <- exp(id - max(id))
      best <- which.max(post)
      if (post[best] / sum(post) >= par$rec_crit) return(best)
    }
    NA_integer_
  })
}

# Proportion of cues whose target was recalled, and of cues yielding some other item.
score_recall <- function(recalled, targets) {
  hit <- !is.na(recalled) & recalled == targets
  c(correct = mean(hit), intrusion = mean(!is.na(recalled) & !hit))
}
