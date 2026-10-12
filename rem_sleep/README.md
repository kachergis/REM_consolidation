# rem_sleep: sleep-based consolidation in REM, re-implemented

A clean re-implementation of the two consolidation mechanisms proposed in Kachergis, de Kleijn &
Hommel (2016, CogSci; `../androids_dream_cogsci2016.tex`), built on REM (Shiffrin & Steyvers,
1997). The original code (`../modeling/`, also at github.com/kachergis/REM_consolidation) is
kept unchanged for reference.

```
Rscript -e 'testthat::test_file("tests/test_rem.R")'   # unit tests + REM benchmarks
Rscript run_experiments.R 500                           # ~8 min on 14 cores; writes results/, figures/
```

Files: `R/rem.R` (REM core), `R/sleep.R` (agent, study, wake, sleep and replay),
`R/recall.R` (cued recall), `R/experiments.R` (simulated experiments), `run_experiments.R`, `tests/test_rem.R`. See
`LIT_REVIEW.md` for the empirical targets.

## Why re-implement

The original simulation's results are artifacts of bugs:

1. `get_dprime()` computes the false-alarm rate from the hit rate (line 104), so FA is always
   1/(2N) and "d'" is a monotone function of the hit rate alone; foils never enter.
2. `likelihood_ratio()` assigns the whole product `(1-c)^n_q · ratio^n_im` to each matching
   feature and then multiplies over features: matches count as ratio^(n²), the mismatch penalty is
   applied once per match, and traces with only mismatches get λ = 1. It passes the S&S worked
   example only because every match there has a distinct value.
3. The "react+reencode" condition runs re-encoding twice and never reactivates.
4. `rgeom()` starts at 0, so ~45% of "non-zero" features are 0; reactivation's benefit came
   mostly from filling these artifactual gaps.
5. Traces were distorted copies (no storage probability u) while the likelihood assumed c = .7.

With bugs 1–3 fixed and the paper's settings, recognition is at ceiling (d' = 5.15) in every
condition.

## Model

Standard REM: items have w = 20 features with values ~ Geometric(g = .4) on {1, 2, ...}. Each
presentation stores each empty feature with probability u = 1 − (1 − .04)^10 = .335, copied
correctly with probability c = .7 (otherwise a random value). Repetitions add to the same trace.
Recognition: log odds = log mean λ over the episodic store, "old" if > 0. We report hit and
false-alarm rates, d' at that criterion, and AUC (criterion-free). All tests reproduce S&S's
worked example (mean λ = .919), per-feature definitions, storage/copy rates, and the
repetition and list-length benchmarks.

**Retention intervals.** Context is not modeled: the whole day's store is the search set.
*Wake* stores 100 new, unrelated events (known items). *Passive sleep* stores nothing, so it is
identical to an immediate test: in this model the only "passive" benefit of sleep is absence of
interference.

**Replay** (`sleep()`): n = 50 traces are sampled; for each, the lexical item that generated it
is inferred by sampling from the REM posterior over LS (lexical-semantic) traces. Then
- *reactivate*: empty features of the trace are re-stored from the inferred LS trace with the
  study process (u, c);
- *re-encode*: 6 features of the inferred LS trace are incremented by 1 (rarer values are more
  diagnostic), and traces whose best match (at sleep onset) is that item and which hold the old
  value are updated.

**Departures from the 2016 paper, each forced by a problem in the original:**
- Re-encoding selects LS items through replayed traces instead of sampling ∝ log familiarity
  (negative log familiarities were clamped to 0, and `sample()` failed once fewer than 50 items
  had positive familiarity).
- Episodic traces are updated only where they stored the old LS value, not "whatever is stored".
- Replayed traces are *identified* by sampling from the posterior. The original matched against
  a lexicon containing only the studied items, which made identification nearly perfect.
  Identification turns out to be the crux (below), so its scope is a factor: `id_scope =
  "lexicon"` (all 1000 known items) or `"today"` (only items experienced since the last sleep,
  i.e. tagged at encoding; cf. Stickgold & Walker, 2013).
- Cued (TMR) replays retrieve a trace ∝ λ given the cued item's LS trace, and use the cued item
  as the source.

## Results: recognition (500 simulated participants; AUC means, SE ≤ .003 throughout)

**Identification is the bottleneck.** Maximum a posteriori (MAP) identification of a
1-repetition trace is 60% correct against the full lexicon and 79% against today's items; 3
repetitions: 91% / 98%. Misidentified replays corrupt traces (reactivation) or change the wrong
LS trace, desynchronizing LS traces from their episodic traces (re-encoding).

### Exp. 1: sleep vs. wake × repetitions × mechanism (single-item recognition)

| AUC | 1 rep | 3 reps | 5 reps |
|---|---|---|---|
| immediate = passive sleep | .866 | .986 | .995 |
| wake (100 interfering events) | .815 | .964 | .982 |
| sleep+reactivate, lexicon / today | .806 / .888 | .980 / .988 | .995 / .996 |
| sleep+reencode, lexicon / today | .843 / .873 | .976 / .986 | .991 / .996 |
| sleep+both, lexicon / today | .801 / .877 | .972 / .987 | .990 / .996 |

- Passive sleep beats wake at every strength, purely by avoiding interference (the list-length
  effect). In d', the gap grows with repetitions (0.30, 0.65, 0.78), consistent with the
  meta-analytic moderator that multi-study designs show larger sleep benefits (Berres &
  Erdfelder, 2021). The size depends entirely on the arbitrary number of wake events.
- With full-lexicon identification, both mechanisms *hurt*. With "today" identification,
  reactivation helps weak memories only (+.022 AUC at 1 rep, ≈0 at 5), as the paper claimed.
  Re-encoding adds little in AUC. The mechanisms are *not* additive: both < reactivate alone.
- d' at the REM criterion tells a rosier story (today, 1 rep: +0.38 reactivate, +0.16 re-encode,
  +0.41 both) because reactivation lowers false alarms by storing more features overall
  (differentiation). Any empirical comparison needs hit and FA rates separately, not only d'.
- 3–5 repetitions are near ceiling; a harder regime (longer lists, lower u) is needed there.

### Exp. 2: interference challenge (A-B, sleep/wake, then optional A-C; associative recognition)

- Sleep and wake are equal: AUC .936 vs .935 without A-C, .918 vs .919 with A-C. Unrelated wake
  events don't interfere with associative recognition, and there is **no sleep × interference
  interaction**. This matches Bailes et al. (2020) and Pöhlchen et al. (2021), not Ellenbogen et
  al. (2006) or Abel et al. (2023).
- Both mechanisms *reduce* associative AUC (to .87–.93): they strengthen item information in each
  half of a pair trace, which also raises the odds for rearranged pairs. Concatenated item vectors
  are a poor model of associative memory. This experiment needs associative features or cued
  recall before it can say anything about Ellenbogen-type results.

### Exp. 3: targeted memory reactivation (half the items cued, fixed replay budget of 50)

| AUC, reactivate | no TMR | TMR cued | TMR uncued |
|---|---|---|---|
| lexicon | .804 | .900 | .817 |
| today | .888 | .928 | .870 |

- Cueing supplies the item's identity, so it bypasses the identification bottleneck: large
  cued benefits, largest when spontaneous identification is poor.
- Whether *uncued* items pay a cost depends on whether the cued replays displace useful
  spontaneous replays (today: −.017) or harmful ones (lexicon: +.013). This distinguishes the two
  scopes. TMR meta-analyses mostly report cued − uncued differences (Hu et al., 2020), which
  can't separate these; between-subjects no-cue controls can. Schechtman et al. (2021) found no
  set-size cost, which argues against a strictly fixed replay budget.
- Re-encoding gains little from cueing.

## Cued recall

Sleep benefits are reliable for recall but ≈0 for recognition (Newbury & Monaghan, 2019), so
the mechanisms were also tested with cued recall of A-B pairs (`R/recall.R`). The model follows
the sampling-then-recovery scheme of ARC-REM (Diller, Nobel & Shiffrin, 2001):
- **Sampling.** The cue A is matched to the A slot of every trace, and up to k_max = 3 distinct
  traces are sampled in proportion to their likelihood ratios.
- **Acceptance.** A sampled trace is used only if its cue log LR exceeds acc_crit = 2.
- **Recovery.** The B slot of an accepted trace is identified against the whole lexicon, the same
  computation sleep replay uses. The best item is output if its posterior is at least .5;
  otherwise the next trace is tried.

The recovery rule is ours; ARC-REM's own equations were not consulted. The acceptance threshold
was set to keep intrusions below omissions; without it (acc_crit = 0), 73% of intrusions are
other studied B's. Exp. 4 was rerun with acc_crit = 0 and the pattern is unchanged (correct
rates shift by ≤ .025). Benchmarks (tests): recall rises with repetitions (.28 / .66 / .83 at
1 / 2 / 3 reps), falls with list length, and A-C learning cuts B recall from .66 to .42 with
C intrusions of 44%.

Values below are paired differences in proportion correct (SE in parentheses), 500 simulated
participants, 40 pairs, 80 wake events, 40 replays.

### Exp. 4: sleep vs. wake × repetitions × mechanism

Change vs. immediate test, at 1 / 2 / 3 repetitions:

| Condition | 1 rep | 2 reps | 3 reps |
|---|---|---|---|
| wake | −.027 (.001) | −.025 (.001) | −.015 (.001) |
| passive sleep | 0 | 0 | 0 |
| sleep+reactivate, lexicon | −.054 (.002) | −.050 (.002) | −.025 (.002) |
| sleep+reactivate, today | **+.067 (.003)** | +.051 (.002) | +.030 (.002) |
| sleep+reencode, lexicon | −.006 (.002) | −.018 (.002) | −.014 (.001) |
| sleep+reencode, today | −.006 (.002) | +.021 (.003) | +.028 (.002) |
| sleep+both, today | +.033 (.003) | +.039 (.003) | +.035 (.002) |

- Passive sleep beats wake by only 1.5–2.7 points in recall, much less than in recognition:
  wake events mostly add intrusions (.30 vs .18 at 1 rep), not omissions.
- The recognition pattern holds. Reactivation needs good identification ("today"), and then
  helps weak pairs most. Re-encoding helps only stronger pairs. The two mechanisms are not
  additive at 1–2 reps.
- Sleep+reactivate vs. wake at 1 rep: +.094, about a third more pairs recalled.

### Exp. 5: interference challenge with cued recall (Ellenbogen et al., 2006 design)

Sleep benefit = condition − wake, without and with A-C learning before test (2 study reps).
Interaction = the wake group's A-C cost minus the sleep group's (positive = sleep protects).

| Condition | Benefit, no A-C | Benefit, A-C | Interaction |
|---|---|---|---|
| passive sleep | +.024 (.001) | +.009 (.003) | −.015 (.004) |
| sleep+reactivate, lexicon | −.023 (.002) | +.026 (.004) | +.049 (.003) |
| sleep+reactivate, today | +.080 (.002) | +.110 (.004) | +.030 (.004) |
| sleep+reencode, lexicon | +.005 (.002) | +.016 (.004) | +.012 (.003) |
| sleep+reencode, today | +.046 (.003) | +.041 (.004) | −.005 (.004) |
| sleep+both, lexicon | −.036 (.002) | +.033 (.003) | +.069 (.003) |
| sleep+both, today | +.067 (.003) | +.123 (.004) | +.055 (.004) |

- **Only active replay yields Ellenbogen's interaction**, and reactivation does so in both scopes.
  Passive sleep (less interference) predicts the opposite sign. So the interaction diagnoses
  active consolidation in this model, which raises the stakes of the replication dispute (Bailes
  et al., 2020 and Pöhlchen et al., 2021 found none; Abel et al., 2023 did).
- Mechanism: reactivation completes A-B traces, so they win the competition with A-C traces for
  cue A (C intrusions .42 after wake → .34 after sleep+reactivate, today).
- With full-lexicon identification, reactivation *crosses over*: it hurts without A-C and
  helps with it.

### Exp. 6: TMR with cued recall (half the A cues presented during sleep)

Change vs. the same items without TMR, sleep+reactivate:

| Scope | Cued, correct | Cued, intrusions | Uncued, correct | Uncued, intrusions |
|---|---|---|---|---|
| lexicon | +.035 (.004) | **+.156 (.006)** | +.009 (.004) | −.042 (.005) |
| today | +.049 (.005) | +.090 (.006) | −.030 (.004) | −.046 (.005) |

- Cueing helps cued pairs, but it also **raises their intrusions**. A cue sometimes retrieves the
  wrong pair trace and then writes its own features into that trace's A slot. This binds the cue
  to another pair's B, a falsifiable prediction: TMR of paired associates should raise cue-specific
  intrusions.
- As in recognition, uncued items pay a cost only when spontaneous replay is useful ("today").
- Re-encoding with TMR: cued +.042 (today), −.014 (lexicon).

## Limitations and next steps

- Recall is cued recall of pairs only; free recall needs a context model.
- No context model, so no temporal gradient of interference, and passive sleep = immediate test.
- The "today" scope assumes the day contained only the study list; realistic days add
  non-studied tagged items, which lowers identification accuracy.
- Replay is uniform over traces; evidence suggests weak memories are prioritized (Schapiro et
  al., 2018; Denis et al., 2021). This is a natural next manipulation.
