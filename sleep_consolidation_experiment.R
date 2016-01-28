source("REM.R")

par = list("M"=100, "geom"=.45, "base_rate"=.4, "copy_prob"=.7)

# REM example (from Shiffrin & Steyvers, 1997)
mem = matrix(c(0,1,0,3, 
               2,2,1,0), ncol=4, byrow=T)
rat = likelihood_ratio(c(2,3,4,3), mem, par)
sum(exp(rat))/2 # .92

# want to see the difference between retrieval/recognition after sleep (no input) vs. a day of new memories
# and then compare to REM with the simple new proposed consolidation mechanisms

# experiment proposal 12/26/15
# simulate 24 hours of study, with or without a period of sleep consolidation in the middle
# - there should be some events (slight distortions of some ideal) that happen both days,
#   and some only on day 1 (for comparison: do we forget these?)
# - maybe some events are more frequent, too?

# ways to use reactivation to strengthen memories (random rows from M):
# 1) fill in some missing features in M (compare against lexicon, then use the strongest lexical match and fill in zeros from lexical trace)
# 2) 

reencode_lexicon <- function(lexicon, M, par, num_lex_reps_updated=5, num_updated_feats=5) {
  # go through lexicon, reactivate traces per lexical item and sample features (randomly or w prob prop to rarity (high values)) to increment
  lex_str = log(test_memory(M, lexicon, par)) # vector of familiarities for each lexical item in M
  # we're choosing lexical items to update that were experienced today
  lex_to_update = sample(1:nrow(lexicon$full_vecs), num_lex_reps_updated, prob=lex_str)
  # choose features to increment and increment them in the lexicon. but we also have to probe memory and reencode matching traces
  for(i in lex_to_update) {
    old_lex = lexicon$full_vecs[i,]
    feats_to_incr = sample(1:length(old_lex), num_updated_feats, prob=old_lex+1) # or could select features uniformly randomly
    lexicon$full_vecs[i,feats_to_incr] = old_lex[feats_to_incr] + 1 
    # lexicon updated. now probe memory with old_lex and update those traces
    match = match_memory_item(old_lex, M, par)
    for(m in match$item) {
      M[m,feats_to_incr] = M[m,feats_to_incr] + 1
    }
  }
  return(list(lexicon=lexicon, M=M))
}

# reactivate N memories for the day and try to fill their empty feature values - GK: DONE!
fill_missing_memory_features <- function(lexicon, M, par, N=5) {
  react = sample(1:nrow(M), N) # randomly sample memories
  for(mi in react) {
    match = match_lexical_item(M[mi,], lexicon, par) # find which lexical item the image is...
    # sometimes we get more than one match, so we let's just sample one (could do all..)
    if(length(match$item)>1) match$item = sample(match$item, 1)
    empty = which(M[mi,]==0)
    if(length(empty)==length(lexicon$full_vecs[match$item,empty])) {
      M[mi,empty] = lexicon$full_vecs[match$item,empty]
    } else {
      cat("mismatch!\n")
      cat("match item: ",match$item,"\n")
      cat("empty: ",empty,"\n")
      cat("mem: ",M[mi,],"\n")
      cat("lex: ",lexicon$full_vecs[match$item,],"\n")
    }
  }
  return(M)
}


simulate_day <- function(lexicon, study_reps, distort_prop, base_rate) {
  memory = matrix(0, nrow=nrow(lexicon$full_vecs)*study_reps, ncol=ncol(lexicon$full_vecs))
  # we take in a raw lexicon, select some items to appear day 1 and 2, and some to appear only day 2, and
  # for each study rep, create distortion of each item and store it in appropriate 
  ind = 1
  for(r in 1:study_reps) {
    for(it in 1:nrow(lexicon$full_vecs)) {
      memory[ind,] = distort_vector(lexicon$full_vecs[it,], distort_prop, base_rate)
      ind = ind + 1
    }
  }
  return(memory)
}

simulate_days <- function(num_days, lexicon, study_reps, distort_prop, base_rate) {
  days = list() 
  for(d in 1:num_days) {
    days[[d]] = simulate_day(lexicon, study_reps, distort_prop, base_rate)
    # consolidate(days[[d]]) # test the consolidation mechanisms we propose
  }
  return(days)
}

get_dprime <- function(fam_liks, unfam_liks, crit=1) {
  # looking at histogram of fam and unfam liks, 1 is not the best criterion placement:
  # median unfam=.03, med fam=.53
  hr = log(fam_liks) / length(fam_liks)
  hr = length(which(hr>crit)) / length(fam_liks)
  fa = log(unfam_liks) / length(unfam_liks)
  fa = length(which(hr>crit)) / length(fam_liks)
  eps = 1/(2*length(fam_liks))
  if(hr==0) hr=eps
  if(hr==1) hr=1-eps
  if(fa==1) fa=1-eps
  if(fa==0) fa=eps
  return(qnorm(hr) - qnorm(fa))
}

sleep_experiment <- function(par, fname, distort=.1, reps=10) {
  #conds = c("day","day+sleep")
  #dat = data.frame(list(cond=rep(conds,reps), dprime=rep(NA,length(conds)*reps)))
  dat = data.frame(list(sim=NA, cond=NA, strength=NA, dprime=NA))
  strengths = c(1,5)
  sim = 1
  for(rep in 1:reps) {
    lexicon = generate_vectors(100, par$M, 30, par$geom) # 100 items, 80 nonzero features
    foils = generate_vectors(100, par$M, 30, par$geom)
    
    for(str in strengths) {
      day = simulate_day(lexicon, str, distort, par$base_rate) # str reps/item
      day_consol = fill_missing_memory_features(lexicon, day, par, N=50) # reactivation of random memories from the day
      reencode_result = reencode_lexicon(lexicon, day, par, num_lex_reps_updated=50, num_updated_feats=6)
      new_lex = reencode_result$lexicon
      day_reencode = reencode_result$M
      
      fam_liks = test_memory(day, lexicon, par)
      unfam_liks = test_memory(day, foils, par)
      dprime = get_dprime(fam_liks, unfam_liks) # .25 
      dat = rbind(dat, c(sim, "none", str, dprime))
      
      # 1 day + sleep:
      fam_liks = test_memory(day_consol, lexicon, par)
      unfam_liks = test_memory(day_consol, foils, par)
      dprime = get_dprime(fam_liks, unfam_liks)
      dat = rbind(dat, c(sim, "reactivate", str, dprime))
      
      fam_liks = test_memory(day_reencode, new_lex, par) # could compare to old lexicon, too...(should obv be worse)
      unfam_liks = test_memory(day_reencode, foils, par)
      dprime = get_dprime(fam_liks, unfam_liks)
      dat = rbind(dat, c(sim, "reencode", str, dprime))
      sim = sim + 1
    }
  }
  #cat(mean(fam_liks), mean(unfam_liks)) # 2.991252e+50 24.97003
  #mean(fam_liks) / mean(unfam_liks)
  dat = na.omit(dat)
  dat$dprime = as.numeric(dat$dprime)
  save(dat, file=paste(fname,".RData",sep=''))
  
  print(t.test(subset(dat, cond=="none")$dprime - subset(dat, cond=="reactivate")$dprime))
  print(t.test(subset(dat, cond=="none")$dprime - subset(dat, cond=="reencode")$dprime))
  print(t.test(subset(dat, cond=="reactivate")$dprime - subset(dat, cond=="reencode")$dprime))
  
  ag = aggregate(dprime ~ cond + strength, mean, data=dat)
  ag$sd = aggregate(dprime ~ cond + strength, sd, data=dat)$dprime
  ag$SE = ag$sd / sqrt(reps-1)
  print(ag)
  limits <- with(ag, aes(ymax=dprime+SE, ymin=dprime-SE, width=.1))
  ggplot(ag, aes(y=dprime, x=cond, group=strength, colour=strength, shape=strength)) + geom_point() +
    geom_line() + geom_errorbar(limits) + theme_bw() + ylab("Mean Sensitivity (d')") + xlab("Consolidation Condition")
  ggsave(paste(fname,".pdf",sep=''), width=4,height=4)
  
  return(dat)
  # day5x = simulate_day(lexicon, 5, .1, par$base_rate) # 5 reps/item
  # day5x_consol = fill_missing_memory_features(lexicon, day5x, par, N=50) 
  #m2days = simulate_days(2, lexicon, 2, .1, .4)
  
  # less benefit of replay-based consolidation of more oft-repeated things (makes sense)
}


par = list("M"=40, "geom"=.45, "base_rate"=.4, "copy_prob"=.7)


sleep_experiment(par, "sleep_sim_reactivateM_distort.1_reps60", distort=.1, reps=60)
sleep_experiment(par, "sleep_sim_reactivateM_distort.15_reps60", distort=.15, reps=60)
sleep_experiment(par, "sleep_sim_reactivateM_distort.2_reps50", distort=.2, reps=50)
sleep_experiment(par, "sleep_sim_reactivateM_distort.25_reps50", distort=.2, reps=50)

run_model <- function(lexicon, par=list()) {
  if(length(par)==0) { # default parameters -- from paper
    par = list()
    par$M = 20 # number of non-zero feature values per trace (expected or guaranteed? rgeom generates 0s..)
    par$geom = .45 # g for high frequency (HF) words > g for LF words
    par$base_rate = .4 # long-term base rate (of features)
    par$copy_prob = .7 # probability of correctly copying a feature value
  }
  orig_probe = lexicon$full_vecs[2,]
  distort_props = c(.1,.2,.3,.4,.5,.6,.8,1)
  fams = c()
  for(d in distort_props) {
    probe = distort_vector(orig_probe, d, par$base_rate)
    fam = sum(exp(likelihood_ratio(probe, lexicon$full_vecs, par)))
    print(likelihood_ratio(probe, lexicon$full_vecs, par))
    fams = c(fams, fam)
  }
  return(fams)
}

