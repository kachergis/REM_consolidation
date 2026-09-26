# Sleep and declarative memory: targets for simulation

Scope: behavioral findings a REM-based consolidation model could simulate, how robust each one
is, and competing computational accounts. Sources marked † were checked against the paper or
its abstract in this review (Sept 2026). The rest are standard citations recalled from memory:
verify before citing.

## 1. How big and how reliable is the sleep benefit?

- **Overall episodic benefit is moderate.** Berres & Erdfelder (2021, *Psych Bull*)†: 823
  effect sizes, 271 samples, sleep vs. equal wake g = 0.44. The benefit is larger when material
  is studied multiple times, and largest in night-sleep/day-wake designs, then naps, then
  deprivation designs.
- **Recall, not recognition.** Newbury & Monaghan (2019, *PB&R*)†, DRM meta-analysis:
  veridical recall g = 0.41 (CI includes 0), veridical recognition g = 0.005; false recall of
  critical lures g = 0.61, false recognition g = −0.15. Within recognition, sleep may help
  recollection but not familiarity (Drosopoulos, Wagner & Born, 2005, *Learn Mem*)†.
  *This matters most for us: the current model does recognition only.*
- **Smaller and less robust than early reports.** Cordi & Rasch (2021, *Curr Opin Neurobiol*)†
  argue that effects depend on boundary conditions such as prior knowledge. Several flagship
  results have failed to replicate (below).
- **Passive vs. active.** The contextual binding theory (Yonelinas et al., 2019, *Nat Rev
  Neurosci*)† holds that sleep helps only by reducing interference and context change. The
  active systems account (Rasch & Born, 2013, *Physiol Rev*; Klinzing, Niethard & Born, 2019,
  *Nat Neurosci*) holds that replay strengthens and transforms memories. Mednick et al. (2011,
  *TINS*) propose an "opportunistic" middle ground. The model's passive-sleep vs. replay
  conditions map directly onto this debate.

## 2. Target phenomena

| # | Phenomenon | Key studies | Robustness | What the model needs |
|---|---|---|---|---|
| A | Sleep > wake retention | Jenkins & Dallenbach 1924; Berres & Erdfelder 2021† | Robust for recall | Recall; a wake-interference model (have) |
| B | Interference challenge: sleep protects A-B from later A-C learning | Ellenbogen et al. 2006, *Curr Biol*†; Abel et al. 2023, *PB&R*† (N = 120, interaction replicated) vs. Bailes et al. 2020, *PLOS One*† and Pöhlchen et al. 2021, *J Sleep Res*† (failures) | **Contested** | Cued recall of paired associates; associative features |
| C | Encoding strength × sleep | Weak benefit more: Drosopoulos et al. 2007; Denis et al. 2021, *J Neurosci*† (spindles). Strong benefit more: Tucker & Fishbein 2008. Wernette et al. 2024, *J Sleep Res*†: depth matters, repetition doesn't. Berres & Erdfelder: repetition helps | **Mixed**, and "strength" is operationalized inconsistently | Have (repetitions); add depth ≈ u or c |
| D | Replay prioritizes weak items | Schapiro et al. 2018, *Nat Commun*† (fMRI replay during rest predicts gains after sleep) | Single study | Replay sampling ∝ weakness (easy) |
| E | TMR (targeted memory reactivation): cued > uncued | Rasch et al. 2007, *Science*; Rudoy et al. 2009, *Science*; meta: Hu et al. 2020, *Psych Bull*† (g = 0.29, NREM only; not REM sleep or wake) | Robust (meta) | Have |
| E2 | TMR helps low-accuracy items; needs a direct cue–memory link | Cairney et al. 2016, *Sleep*† | Moderate | Have (identification bottleneck predicts this) |
| E3 | No capacity cost: cueing groups of 1, 2 or 6 items is equally effective | Schechtman et al. 2021, *Commun Biol*† | Single study | Tests the fixed-replay-budget assumption |
| E4 | Competition: cueing impairs well-learned competing pairs | Antony et al. 2018, *Neurobiol Learn Mem*† | Single study | Competing traces sharing a cue |
| F | Sleep increases false recall of DRM lures | Payne et al. 2009, *NLM*†; registered report (Mak et al.)† found more lure and veridical recall but fewer intrusions; recognition: no effect (Newbury & Monaghan) | Recall: moderate | Semantic structure among LS traces; recall |
| G | Lexical integration: new words compete with neighbors only after sleep | Dumay & Gaskell 2007, *Psych Sci*†; Tamminen et al. 2010 (spindles; PDF in `../sleep consolidation papers/`) | Replicated, with boundary conditions | Re-encoding of LS traces is the natural mechanism; needs a lexical-decision / competition measure |
| H | Category structure is gained while exemplar details are preserved | Schapiro et al. 2017, *Sci Rep*† | Single study | Shared vs. unique features across LS traces |
| I | Future relevance gates the benefit | Wilhelm et al. 2011, *J Neurosci*†; a later study found future-relevant memories are *not* selectively strengthened (PMC8568116)† | **Contested** | Tagging (≈ our `id_scope = "today"`) |
| J | Schema/gist formation via overlapping replay | Lewis & Durrant 2011, *TiCS*† (iOtA model); Stickgold & Walker 2013, *Nat Neurosci*† (memory triage) | Theory | Overlapping LS features |

## 3. Computational models to position against

- **CLS (complementary learning systems):** McClelland, McNaughton & O'Reilly 1995 (interleaved
  replay prevents catastrophic interference); Norman, Newman & Perotte 2005 (oscillating
  inhibition during REM sleep); Káli & Dayan 2004 (replay maintains hippocampus-dependent
  memories; cited in the 2016 paper).
- **C-HORSE / NREM–REM alternation:** Singh, Norman & Schapiro 2022, *PNAS*†. The hippocampus
  teaches cortex during NREM; cortex replays its own knowledge during REM, which protects old
  knowledge. Alternation outperforms NREM alone.
- **iOtA:** Lewis & Durrant 2011†.
- **Generative-model accounts:** e.g., hippocampal replay training generative cortical models
  (PMC10963272)†.
- **Global matching / multitrace:** SARKAE (Nelson & Shiffrin 2013) is the closest REM relative
  that accumulates knowledge from events, and is the obvious home for LS updating during sleep.
  MINERVA 2 (Hintzman 1984). I found no REM-family model of sleep besides the 2016 paper. That
  is a niche, but the paper should cite CLS and C-HORSE and say what a REM-level model adds:
  closed-form likelihoods, and hit/FA predictions.

## 4. What the simulations so far imply, and priorities

1. **Implement recall** (REM/SAM-style sampling + recovery) before anything else. Recall carries
   the robust sleep effects (A, B, F), and recognition barely shows them.
2. **Wake vs. passive sleep** already reproduces a sleep advantage that grows with repetitions
   (cf. the Berres & Erdfelder moderator). Magnitude is set by the arbitrary amount of wake
   interference, so fit it to data (e.g., hours awake → events).
3. **Identification bottleneck → testable predictions:**
   - TMR benefit is largest for weakly encoded items (E2 is consistent).
   - Uncued items pay a cost under a fixed replay budget only if spontaneous replay is useful
     (E3 argues against a fixed budget).
   - Spontaneous replay of very weak traces can harm them. This is a counterintuitive,
     falsifiable prediction.
4. **Interference challenge (B)** is contested; simulate it with cued recall and report both the
   passive and active predictions. Current associative-recognition results: no sleep ×
   interference interaction.
5. **Lexical integration (G)** and **category structure (H)** are where LS re-encoding could earn
   its keep. It currently adds little to episodic recognition.
6. **Report hits and FAs separately.** Reactivation lowers FAs via differentiation, which
   inflates criterion-based d' relative to AUC.

## Sources checked in this review

- Ellenbogen et al. 2006: https://pubmed.ncbi.nlm.nih.gov/16824917/
- Bailes et al. 2020: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0220419
- Pöhlchen et al. 2021: https://onlinelibrary.wiley.com/doi/10.1111/jsr.13042
- Abel et al. 2023: https://pmc.ncbi.nlm.nih.gov/articles/PMC10728269/
- Hu et al. 2020: https://pubmed.ncbi.nlm.nih.gov/32027149/
- Berres & Erdfelder 2021: https://www.semanticscholar.org/paper/53a678d237a397ff6a46876285423bfff550718c
- Newbury & Monaghan 2019: https://pmc.ncbi.nlm.nih.gov/articles/PMC6488554/
- Drosopoulos et al. 2005: https://pubmed.ncbi.nlm.nih.gov/15687230/
- Cordi & Rasch 2021: https://www.sciencedirect.com/science/article/pii/S0959438820300994
- Yonelinas et al. 2019: https://pubmed.ncbi.nlm.nih.gov/30872808/
- Payne et al. 2009: https://pubmed.ncbi.nlm.nih.gov/19348959/
- Mak et al. (DRM registered report): https://eprints.whiterose.ac.uk/id/eprint/206094/1/Mak_FalseMemory.pdf
- Dumay & Gaskell 2007: https://journals.sagepub.com/doi/10.1111/j.1467-9280.2007.01845.x
- Denis et al. 2021: https://www.jneurosci.org/content/41/18/4088
- Wernette et al. 2024: https://pubmed.ncbi.nlm.nih.gov/38196130/
- Schapiro et al. 2017: https://www.nature.com/articles/s41598-017-12884-5
- Schapiro et al. 2018: https://www.nature.com/articles/s41467-018-06213-1
- Singh, Norman & Schapiro 2022: https://www.pnas.org/doi/10.1073/pnas.2123432119
- Wilhelm et al. 2011: https://pubmed.ncbi.nlm.nih.gov/21289163/ ; non-replication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8568116/
- Cairney et al. 2016: https://pubmed.ncbi.nlm.nih.gov/26856905/
- Schechtman et al. 2021: https://pubmed.ncbi.nlm.nih.gov/33398075/
- Antony et al. 2018: https://pubmed.ncbi.nlm.nih.gov/30092311/
- Lewis & Durrant 2011: https://pubmed.ncbi.nlm.nih.gov/21764357/
- Stickgold & Walker 2013: https://www.nature.com/articles/nn.3303
- Generative consolidation model: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10963272/
