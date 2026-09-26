# rem_sleep: sleep-based consolidation in REM, re-implemented

A clean re-implementation of the two consolidation mechanisms proposed in Kachergis, de Kleijn &
Hommel (2016, CogSci; `../androids_dream_cogsci2016.tex`), built on REM (Shiffrin & Steyvers,
1997). The original code (`../modeling/`, also at github.com/kachergis/REM_consolidation) is
kept unchanged for reference.

```
Rscript -e 'testthat::test_file("tests/test_rem.R")'   # unit tests + REM benchmarks
Rscript run_experiments.R 500                           # ~1 min on 14 cores; writes results/, figures/
```

Files: `R/rem.R` (REM core), `R/sleep.R` (agent, study, wake, sleep and replay),
`R/experiments.R` (simulated experiments), `run_experiments.R`, `tests/test_rem.R`. See
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

## Results (500 simulated participants; AUC means, SE ≤ .003 throughout)

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

## Limitations and next steps

- Recognition only. Sleep benefits are reliable for recall but ≈0 for recognition
  (Newbury & Monaghan, 2019: veridical recognition g = 0.005), so the headline targets need
  REM-style cued/free recall (e.g., sampling + recovery as in SAM/REM variants).
- No context model, so no temporal gradient of interference, and passive sleep = immediate test.
- The "today" scope assumes the day contained only the study list; realistic days add
  non-studied tagged items, which lowers identification accuracy.
- Replay is uniform over traces; evidence suggests weak memories are prioritized (Schapiro et
  al., 2018; Denis et al., 2021). This is a natural next manipulation.
