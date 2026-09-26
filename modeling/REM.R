# REM (Retrieving Effectively from Memory; Shiffrin & Steyvers, 1997)

generate_vectors <- function(num_items, num_features, num_nonzero, geom_param) {
  # geom_param can be a vector of different parameter values (for each item)
  #if(length(geom_param)==1) rep(geom_param, num_items) 
  compact_vecs <- matrix(0, nrow=num_items, ncol=num_nonzero)
  nonzero_indices <- matrix(0, nrow=num_items, ncol=num_nonzero)
  full_vecs <- matrix(0, nrow=num_items, ncol=num_features)
  for(i in 1:num_items) { 
    compact_vecs[i,] <- rgeom(num_nonzero, geom_param) # problem: generates many 0s (should we +1)?
    nonzero_indices[i,] <- sample(1:num_features, num_nonzero, replace=F)
    full_vecs[i,nonzero_indices[i,]] <- compact_vecs[i,]
  }
  return(list(compact_vecs=compact_vecs, nonzero_indices=nonzero_indices, full_vecs=full_vecs))
}


likelihood_ratio <- function(probe, vectors, par) {
  # calculate lambda_j: the prob that D_j would have been observed if image I_j was a same-image
  # divided by the prob that D_j would have been observed if I_j was a different-image
  c = par$copy_prob
  g = par$base_rate
  lik = matrix(1, nrow=nrow(vectors), ncol=ncol(vectors))
  for(j in 1:nrow(vectors)) {
    mismatch_inds = which(probe!=vectors[j,])
    n_q = length(which(vectors[j,mismatch_inds]>0)) # number of all nonzero mismatching features in image j
    #print(vectors[j,])
    #print(paste("n_jq =",n_q))
    for(i in 1:max(vectors[j,])) {
      inds = which(vectors[j,]==probe & vectors[j,]==i)
      n_im = length(inds) # number of nonzero features in image j that match the probe and have value i
      #print(paste("n_ijm =",n_im))
      lik[j,inds] =  (1-c)^n_q * ( (c+(1-c)*g*(1-g)^(i-1)) / (g*(1-g)^(i-1)) )^n_im # * lik[j,inds] 
    }
  }
  log_lik = rowSums(log(lik))
  #print(lik)
  return(log_lik)
}

# test each lexical item against memory
test_memory <- function(M, lex, par, verbose=F) {
  fams = c()
  for(i in 1:nrow(lex$full_vecs)) {
    #probe = distort_vector(lexicon$full_vecs[i,], .1, par$base_rate) # should we distort before probing?
    probe = lex$full_vecs[i,]
    fam = sum(exp(likelihood_ratio(probe, M, par))) 
    #fam = sum(likelihood_ratio(probe, lex$full_vecs, par))
    if(verbose) {
      print(likelihood_ratio(probe, M, par))
    }
    fams = c(fams, fam)
  }
  return(fams)
}

# probes the lexicon and returns best match (even if it's not great..)
match_lexical_item <- function(probe, lex, par) {
  liks = likelihood_ratio(probe, lex$full_vecs, par)
  return(list(item=which(liks==max(liks)), likelihood=max(liks)))
}

match_memory_item <- function(probe, M, par) {
  liks = likelihood_ratio(probe, M, par)
  return(list(item=which(liks==max(liks)), likelihood=max(liks)))
}

distort_vector <- function(vec, prop_distort, g) {
  num_to_change = round(prop_distort*length(vec))
  change_inds = sample(1:length(vec), num_to_change)
  new_vec = vec
  new_vec[change_inds] = rgeom(num_to_change, g)
  return(new_vec)
}

run_once <- function(par=list()) {
  if(length(par)==0) { # default parameters -- from paper
    par = list()
    par$M = 20 # number of non-zero feature values per trace (expected or guaranteed? rgeom generates 0s..)
    par$geom = .45 # g for high frequency (HF) words > g for LF words
    par$base_rate = .4 # long-term base rate (of features)
    par$copy_prob = .7 # probability of correctly copying a feature value
  }
  lexicon = generate_vectors(10, 20, 10, par$geom)
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


