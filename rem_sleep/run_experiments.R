# Run all simulated experiments; writes results/*.csv and figures/*.pdf.
# Usage (from rem_sleep/): Rscript run_experiments.R [n_subjects]
library(ggplot2)
for (f in c("R/rem.R", "R/sleep.R", "R/recall.R", "R/experiments.R")) source(f)

args <- commandArgs(TRUE)
n_subj <- if (length(args)) as.integer(args[1]) else 200
par <- rem_par()
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

both_scopes <- function(fun, seed, p = par) {
  do.call(rbind, lapply(c("lexicon", "today"), function(sc) {
    run_subjects(fun, n_subj, seed, par = p, id_scope = sc)
  }))
}

# Mean and SE over subjects of each measure, and paired differences in d' from a baseline
# condition. Figures show AUC, which (unlike d' at the fixed REM criterion) is unaffected by
# shifts of the whole log-odds distribution.
summarise <- function(d, by, base_col, base) {
  ag <- aggregate(cbind(hr, fa, dprime, auc) ~ ., d[c(by, "hr", "fa", "dprime", "auc")], mean)
  se <- function(x) sd(x) / sqrt(length(x))
  rhs <- paste(by, collapse = "+")
  ag$se <- aggregate(as.formula(paste("dprime ~", rhs)), d, se)$dprime
  ag$se_auc <- aggregate(as.formula(paste("auc ~", rhs)), d, se)$auc
  keys <- setdiff(c("subj", by), base_col)
  b <- d[d[[base_col]] == base, c(keys, "dprime")]
  names(b)[names(b) == "dprime"] <- "d_base"
  m <- merge(d, b, by = keys)
  m$diff <- m$dprime - m$d_base
  df <- aggregate(as.formula(paste("diff ~", paste(by, collapse = "+"))), m,
                  function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
  ag <- merge(ag, do.call(data.frame, df), by = by)
  names(ag)[names(ag) %in% c("diff.mean", "diff.se")] <- paste0(c("diff_", "diff_se_"), base)
  ag
}

dodge <- position_dodge(width = .3)

# Exp. 1: sleep vs. wake x repetitions x mechanism
e1 <- both_scopes(exp_sleep_wake, 1e4)
write.csv(e1, "results/exp1_sleep_wake.csv", row.names = FALSE)
s1 <- summarise(e1, c("id_scope", "reps", "cond"), "cond", "immediate")
s1$cond <- factor(s1$cond, levels = CONDS)
print(s1[order(s1$id_scope, s1$reps, s1$cond), ], digits = 3, row.names = FALSE)
ggplot(s1, aes(reps, auc, colour = cond, shape = cond)) +
  geom_point(position = dodge) + geom_line(position = dodge) +
  geom_errorbar(aes(ymin = auc - se_auc, ymax = auc + se_auc), width = .2, position = dodge) +
  facet_wrap(~ id_scope, labeller = label_both) + scale_x_continuous(breaks = c(1, 3, 5)) +
  labs(x = "Study repetitions per item", y = "AUC", colour = NULL, shape = NULL) + theme_bw()
ggsave("figures/exp1_sleep_wake.pdf", width = 8, height = 4)

# Exp. 2: interference challenge (associative recognition)
e2 <- both_scopes(exp_interference, 2e4)
write.csv(e2, "results/exp2_interference.csv", row.names = FALSE)
s2 <- summarise(e2, c("id_scope", "interference", "cond"), "cond", "wake")
print(s2[order(s2$id_scope, s2$interference, s2$cond), ], digits = 3, row.names = FALSE)
s2$interference <- factor(s2$interference, levels = c("none", "A-C"))
ggplot(s2, aes(interference, auc, colour = cond, shape = cond, group = cond)) +
  geom_point(position = dodge) + geom_line(position = dodge) +
  geom_errorbar(aes(ymin = auc - se_auc, ymax = auc + se_auc), width = .2, position = dodge) +
  facet_wrap(~ id_scope, labeller = label_both) +
  labs(x = "Interference list before test", y = "Associative recognition AUC", colour = NULL,
       shape = NULL) + theme_bw()
ggsave("figures/exp2_interference.pdf", width = 8, height = 4)

# Exp. 3: targeted memory reactivation
e3 <- both_scopes(exp_tmr, 3e4)
write.csv(e3, "results/exp3_tmr.csv", row.names = FALSE)
e3$cue <- paste(ifelse(e3$tmr, "TMR", "no TMR"), e3$items)
s3 <- summarise(e3, c("id_scope", "cond", "cue"), "cue", "no TMR cued")
print(s3[order(s3$id_scope, s3$cond, s3$cue), ], digits = 3, row.names = FALSE)
ggplot(s3, aes(cond, auc, colour = cue, shape = cue)) +
  geom_point(position = dodge) +
  geom_errorbar(aes(ymin = auc - se_auc, ymax = auc + se_auc), width = .2, position = dodge) +
  facet_wrap(~ id_scope, labeller = label_both) +
  labs(x = NULL, y = "AUC", colour = NULL, shape = NULL) + theme_bw()
ggsave("figures/exp3_tmr.pdf", width = 10, height = 4)

# ---- Cued recall -------------------------------------------------------------------------------

# Mean correct and intrusion rates with SEs, and paired differences in correct recall from a
# baseline condition.
summarise_recall <- function(d, by, base_col, base) {
  rhs <- paste(by, collapse = "+")
  se <- function(x) sd(x) / sqrt(length(x))
  ms <- intersect(c("correct", "intrusion", "c_intrusion"), names(d))
  ag <- aggregate(as.formula(paste0("cbind(", paste(ms, collapse = ","), ") ~ ", rhs)), d, mean)
  ag$se <- aggregate(as.formula(paste("correct ~", rhs)), d, se)$correct
  keys <- setdiff(c("subj", by), base_col)
  b <- d[d[[base_col]] == base, c(keys, "correct")]
  names(b)[names(b) == "correct"] <- "c_base"
  m <- merge(d, b, by = keys)
  m$diff <- m$correct - m$c_base
  df <- aggregate(as.formula(paste("diff ~", rhs)), m, function(x) c(mean = mean(x), se = se(x)))
  ag <- merge(ag, do.call(data.frame, df), by = by)
  names(ag)[names(ag) %in% c("diff.mean", "diff.se")] <- paste0(c("diff_", "diff_se_"), base)
  ag
}

# Exp. 4: sleep vs. wake x repetitions x mechanism
e4 <- both_scopes(exp_recall_sleep_wake, 4e4)
write.csv(e4, "results/exp4_recall_sleep_wake.csv", row.names = FALSE)
s4 <- summarise_recall(e4, c("id_scope", "reps", "cond"), "cond", "immediate")
s4$cond <- factor(s4$cond, levels = CONDS)
print(s4[order(s4$id_scope, s4$reps, s4$cond), ], digits = 3, row.names = FALSE)
ggplot(s4, aes(reps, correct, colour = cond, shape = cond)) +
  geom_point(position = dodge) + geom_line(position = dodge) +
  geom_errorbar(aes(ymin = correct - se, ymax = correct + se), width = .2, position = dodge) +
  facet_wrap(~ id_scope, labeller = label_both) + scale_x_continuous(breaks = 1:3) +
  labs(x = "Study repetitions per pair", y = "Cued recall (proportion correct)", colour = NULL,
       shape = NULL) + theme_bw()
ggsave("figures/exp4_recall_sleep_wake.pdf", width = 8, height = 4)

# Robustness: Exp. 4 with no acceptance threshold on sampled traces (acc_crit = 0)
e4r <- both_scopes(exp_recall_sleep_wake, 4e4, p = modifyList(par, list(acc_crit = 0)))
write.csv(e4r, "results/exp4_recall_sleep_wake_acc0.csv", row.names = FALSE)
s4r <- summarise_recall(e4r, c("id_scope", "reps", "cond"), "cond", "immediate")
print(s4r[order(s4r$id_scope, s4r$reps, s4r$cond), ], digits = 3, row.names = FALSE)

# Exp. 5: interference challenge, cued recall
e5 <- both_scopes(exp_recall_interference, 5e4)
write.csv(e5, "results/exp5_recall_interference.csv", row.names = FALSE)
s5 <- summarise_recall(e5, c("id_scope", "interference", "cond"), "cond", "wake")
print(s5[order(s5$id_scope, s5$interference, s5$cond), ], digits = 3, row.names = FALSE)
s5$interference <- factor(s5$interference, levels = c("none", "A-C"))
ggplot(s5, aes(interference, correct, colour = cond, shape = cond, group = cond)) +
  geom_point(position = dodge) + geom_line(position = dodge) +
  geom_errorbar(aes(ymin = correct - se, ymax = correct + se), width = .2, position = dodge) +
  facet_wrap(~ id_scope, labeller = label_both) +
  labs(x = "Interference list before test", y = "Cued recall of B (proportion correct)",
       colour = NULL, shape = NULL) + theme_bw()
ggsave("figures/exp5_recall_interference.pdf", width = 8, height = 4)

# Exp. 6: TMR, cued recall
e6 <- both_scopes(exp_recall_tmr, 6e4)
write.csv(e6, "results/exp6_recall_tmr.csv", row.names = FALSE)
e6$cue <- paste(ifelse(e6$tmr, "TMR", "no TMR"), e6$items)
s6 <- summarise_recall(e6, c("id_scope", "cond", "cue"), "cue", "no TMR cued")
print(s6[order(s6$id_scope, s6$cond, s6$cue), ], digits = 3, row.names = FALSE)
ggplot(s6, aes(cond, correct, colour = cue, shape = cue)) +
  geom_point(position = dodge) +
  geom_errorbar(aes(ymin = correct - se, ymax = correct + se), width = .2, position = dodge) +
  facet_wrap(~ id_scope, labeller = label_both) +
  labs(x = NULL, y = "Cued recall (proportion correct)", colour = NULL, shape = NULL) + theme_bw()
ggsave("figures/exp6_recall_tmr.pdf", width = 10, height = 4)

# Exp. 7: replay prioritization of weak memories (mixed-strength list), cued recall
e7 <- both_scopes(exp_recall_strength, 7e4)
write.csv(e7, "results/exp7_recall_strength.csv", row.names = FALSE)
e7$cond <- ifelse(e7$prio > 0, paste0(e7$cond, ", prio ", e7$prio), e7$cond)
s7 <- summarise_recall(e7, c("id_scope", "strength", "cond"), "cond", "wake")
s7$cond <- factor(s7$cond, levels = c("immediate", "wake", "sleep",
  sapply(c("sleep+reactivate", "sleep+both"), function(m) c(m, paste0(m, ", prio ", c(1, 4))))))
print(s7[order(s7$id_scope, s7$strength, s7$cond), ], digits = 3, row.names = FALSE)
ggplot(s7, aes(cond, correct, colour = strength, shape = strength)) +
  geom_point(position = dodge) +
  geom_errorbar(aes(ymin = correct - se, ymax = correct + se), width = .2, position = dodge) +
  facet_wrap(~ id_scope, labeller = label_both) +
  labs(x = NULL, y = "Cued recall (proportion correct)", colour = "Pair strength",
       shape = "Pair strength") + theme_bw() + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave("figures/exp7_recall_strength.pdf", width = 10, height = 5)

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
